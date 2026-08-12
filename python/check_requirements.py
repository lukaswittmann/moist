"""Fail the testsuite when a Python test dependency is missing

The python test targets are only defined when the extension module is built, so
at that point pytest and pyscf are requirements rather than optional extras;
this is done to make sure the CI and tests actualy test these instead of silently 
skipping these.
"""

import importlib
import sys


def main(specs):
    missing = []
    for spec in specs:
        package, _, subs = spec.partition(":")
        names = [package] + [f"{package}.{sub}" for sub in subs.split(",") if sub]
        for name in names:
            try:
                importlib.import_module(name)
            except ImportError as exc:
                missing.append(f"  {name}: {exc}")
                if name == package:
                    break

    if missing:
        sys.exit(
            "The Python testsuite requires these modules, which are unusable:\n"
            + "\n".join(missing)
            + "\nInstall them, or configure with -Dpython=false to drop the"
            + " python test targets."
        )

    print("all Python test requirements are importable")


if __name__ == "__main__":
    main(sys.argv[1:])
