/* Native file operations shared by the Fortran checkpoint stores.
 * Temporary files are claimed exclusively, including with Intel runtimes.
 * Replacement stays on one filesystem and never deletes the old file first.
 */
#include <stdio.h>

#ifdef _WIN32
#include <windows.h>
#include <stdint.h>
#include <stdlib.h>

/* The Fortran runtimes open narrow names in the ANSI code page, so the claimed and
 * replaced names are converted the same way before the UTF-16 file operations. */
static wchar_t *wide_path(const char *path)
{
   int length = MultiByteToWideChar(CP_ACP, MB_ERR_INVALID_CHARS, path, -1, NULL, 0);
   wchar_t *wide;
   if (length <= 0 || (size_t)length > SIZE_MAX / sizeof(wchar_t)) return NULL;
   wide = malloc((size_t)length * sizeof(wchar_t));
   if (wide == NULL) return NULL;
   if (MultiByteToWideChar(CP_ACP, MB_ERR_INVALID_CHARS, path, -1, wide, length) == 0) {
      free(wide);
      return NULL;
   }
   return wide;
}
#else
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#endif

/* 0: claimed, 1: already exists, -1: another failure. */
int moist_file_create_exclusive(const char *path)
{
#ifdef _WIN32
   wchar_t *wide = wide_path(path);
   HANDLE file;
   DWORD status;
   if (wide == NULL) return -1;
   file = CreateFileW(wide, GENERIC_WRITE, 0, NULL, CREATE_NEW, FILE_ATTRIBUTE_NORMAL, NULL);
   if (file == INVALID_HANDLE_VALUE) {
      status = GetLastError();
      free(wide);
      return status == ERROR_FILE_EXISTS || status == ERROR_ALREADY_EXISTS ? 1 : -1;
   }
   if (!CloseHandle(file)) {
      DeleteFileW(wide);
      free(wide);
      return -1;
   }
   free(wide);
   return 0;
#else
   int file = open(path, O_CREAT | O_EXCL | O_WRONLY, 0666);
   if (file < 0) return errno == EEXIST ? 1 : -1;
   if (close(file) != 0) {
      unlink(path);
      return -1;
   }
   return 0;
#endif
}

/* 0: replaced, -1: failed; no copy/delete fallback is permitted. */
int moist_file_replace(const char *source, const char *target)
{
#ifdef _WIN32
   wchar_t *wide_source = wide_path(source);
   wchar_t *wide_target = wide_path(target);
   BOOL replaced;
   if (wide_source == NULL || wide_target == NULL) {
      free(wide_source);
      free(wide_target);
      return -1;
   }
   replaced = MoveFileExW(wide_source, wide_target, MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH);
   free(wide_source);
   free(wide_target);
   return replaced ? 0 : -1;
#else
   return rename(source, target) == 0 ? 0 : -1;
#endif
}
