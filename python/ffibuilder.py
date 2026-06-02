"""
FFI builder module for moist for usage from meson and from setup.py.

Meson handles preprocessing the header file and passes the resulting C
definitions to this script. Outside Meson we preprocess the header ourselves.
"""

import os

import cffi


library = "moist"
include_header = '#include "moist.h"'
prefix_var = "MOIST_PREFIX"
if prefix_var not in os.environ:
    prefix_var = "CONDA_PREFIX"

if __name__ == "__main__":
    import sys

    kwargs = dict(libraries=[library])

    header_file = sys.argv[1]
    module_name = sys.argv[2]

    with open(header_file) as f:
        cdefs = f.read()
else:
    import subprocess

    try:
        import pkgconfig

        if not pkgconfig.exists(library):
            raise ModuleNotFoundError("Unable to find pkg-config package 'moist'")
        if pkgconfig.installed(library, "< 0.5"):
            raise Exception("Installed 'moist' version is too old, 0.5 or newer is required")

        kwargs = pkgconfig.parse(library)
        cflags = pkgconfig.cflags(library).split()

    except ModuleNotFoundError:
        kwargs = dict(libraries=[library])
        cflags = []
        if prefix_var in os.environ:
            prefix = os.environ[prefix_var]
            kwargs.update(
                include_dirs=[os.path.join(prefix, "include")],
                library_dirs=[os.path.join(prefix, "lib")],
                runtime_library_dirs=[os.path.join(prefix, "lib")],
            )
            cflags.append("-I" + os.path.join(prefix, "include"))

    cc = os.environ["CC"] if "CC" in os.environ else "cc"
    module_name = "moist._libmoist"

    process = subprocess.Popen(
        [cc, *cflags, "-E", "-"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    out, err = process.communicate(include_header.encode())
    if process.returncode != 0:
        raise RuntimeError(err.decode())

    cdefs = out.decode()

ffibuilder = cffi.FFI()
ffibuilder.set_source(module_name, include_header, **kwargs)
ffibuilder.cdef(cdefs)

if __name__ == "__main__":
    ffibuilder.distutils_extension(".")
