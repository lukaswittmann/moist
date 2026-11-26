/* moist replacement for the CMake GenerateExportHeader output of libmsym
 *
 * libmsym is compiled into libmoist and called only through moist's own shim
 * (msym_shim.c), so no libmsym symbol is exported: MSYM_EXPORT is empty
 */
#ifndef MSYM_EXPORT_H
#define MSYM_EXPORT_H

#define MSYM_EXPORT

#endif /* MSYM_EXPORT_H */
