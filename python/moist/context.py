"""Shared run settings and native context lifetime."""

from . import library


class Context:
    """Run context shared by cavities, models and their components.

    ``nthreads`` is fixed for the context's lifetime: a positive count is
    used as given, ``0`` (the default) takes the calling thread's OpenMP
    setting at construction and a negative count raises. The count sizes
    MOIST's own OpenMP regions; the host's OpenMP runtime and BLAS/LAPACK
    threading are left alone. Every model requires a context; a cavity or
    component takes one optionally and otherwise runs on its model's. Each
    retains its context and takes its verbosity and debug settings from it.
    """

    def __init__(self, nthreads: int = 0, verbosity: int = 0, debug: bool = False) -> None:
        self._handle = library.new_context(nthreads, verbosity, debug)

    @property
    def nthreads(self) -> int:
        """Thread count fixed at construction."""
        return library.get_context_num_threads(self._handle)

    def _as_handle(self) -> library.ContextHandle:
        return self._handle


def _resolve_context(context, *, optional=False):
    if context is None and optional:
        return None
    if not isinstance(context, Context):
        raise TypeError("context must be a Context")
    return context


def _context_handle(context):
    """Native handle of an optional context."""
    return None if context is None else context._as_handle()
