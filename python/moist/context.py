"""Shared run settings and native context lifetime."""

from contextlib import contextmanager, suppress
from typing import Iterator

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

    ``verbosity`` runs from ``0`` (silent) to ``6`` (debug), ``2`` by default, see
    the API reference for the levels; ``debug=True`` raises it to ``6``.
    """

    def __init__(self, nthreads: int = 0, verbosity: int = 2, debug: bool = False) -> None:
        self._handle = library.new_context(nthreads, verbosity, debug)

    @property
    def nthreads(self) -> int:
        """Thread count fixed at construction."""
        return library.get_context_num_threads(self._handle)

    def print_timer(self, simple: bool = False) -> None:
        """Print the recorded timings to the native output at this context's verbosity.

        Nothing prints at verbosity 0, the top-level phases at 1 and the
        whole tree from 4. ``simple=True`` prints a compact list of seconds
        per step, at most two levels deep. The output is written by the native
        library to the process standard output, not to ``sys.stdout``.
        """
        if not isinstance(simple, bool):
            raise TypeError("simple must be a bool")
        library.print_context_timer(self._handle, simple)

    @contextmanager
    def timer(self, name: str) -> Iterator[None]:
        """Time the enclosed host work as ``name`` in this context's timer report.

        Every host timer lands under a top-level ``Host`` node, apart from
        moist's own timers, and nests only inside other host timers.  Time
        only host work: a moist call inside the block is counted both here
        and in moist's own timers.  A repeated name under the same parent
        accumulates into one node.  The timer stops even when the block
        raises, and that exception stays the one raised: a failing stop does
        not replace it.

        A name has 1 to 32 bytes of UTF-8, no control character and no leading
        or trailing space. A non-``str`` raises
        :class:`TypeError`, a NUL or an unencodable character
        :class:`ValueError`, and a name breaking the other rules
        :class:`RuntimeError`.
        """
        self.start_timer(name)
        try:
            yield
        except BaseException:
            with suppress(RuntimeError):
                self.stop_timer(name)
            raise
        self.stop_timer(name)

    def start_timer(self, name: str) -> None:
        """Start the host timer ``name`` (rules as in :meth:`timer`); prefer :meth:`timer`."""
        library.start_context_timer(self._handle, name)

    def stop_timer(self, name: str) -> None:
        """Stop the innermost open host timer, which must be ``name``; anything else raises."""
        library.stop_context_timer(self._handle, name)

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
