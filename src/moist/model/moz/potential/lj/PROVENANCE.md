# Permissive sources for Lennard-Jones parameter tables

Citations identify the models; the linked source-code licenses supply the redistribution terms. Upstream notices are retained in the parameter modules and summarized in `THIRD_PARTY_NOTICES.md`, which is included in installed distributions.

## UFF

Source: Peter Boyd's [`lammps_interface/uff.py`](https://github.com/peteboyd/lammps_interface/blob/255f027cb76142d39c050a6810404debc6a06562/lammps_interface/uff.py), commit `255f027cb76142d39c050a6810404debc6a06562` (2023-05-24).
The repository's [`LICENSE`](https://github.com/peteboyd/lammps_interface/blob/255f027cb76142d39c050a6810404debc6a06562/LICENSE) is MIT, `Copyright (c) 2017 Peter Boyd`.

Citation: A. K. Rappe, C. J. Casewit, K. S. Colwell, W. A. Goddard III, and W. M. Skiff, *UFF, a full periodic table force field for molecular mechanics and molecular dynamics simulations*, **J. Am. Chem. Soc. 114** (1992), 10024-10035, [doi:10.1021/ja00051a040](https://doi.org/10.1021/ja00051a040).

## DREIDING

Source: Peter Boyd's [`lammps_interface/dreiding.py`](https://github.com/peteboyd/lammps_interface/blob/255f027cb76142d39c050a6810404debc6a06562/lammps_interface/dreiding.py), at the same revision and under the same MIT license as UFF, `Copyright (c) 2017 Peter Boyd`.

Citation: S. L. Mayo, B. D. Olafson, and W. A. Goddard III, *DREIDING: a generic force field for molecular simulations*, **J. Phys. Chem. 94** (1990), 8897-8909, [doi:10.1021/j100389a010](https://doi.org/10.1021/j100389a010), Table II and equations 36a/36c.

## GAFF2

Permissive reference: [`openmmforcefields/gaff-2.2.20.dat`](https://github.com/openmm/openmmforcefields/blob/f0dfefff34e8af3cd5742686b6a9015661d4d5b6/openmmforcefields/ffxml/amber/gaff/dat/gaff-2.2.20.dat), commit `f0dfefff34e8af3cd5742686b6a9015661d4d5b6` (2026-10-06). Its [license](https://github.com/openmm/openmmforcefields/blob/f0dfefff34e8af3cd5742686b6a9015661d4d5b6/LICENSE) is MIT, `Copyright (c) 2016-2019` Chodera lab / Memorial Sloan Kettering Cancer Center and Pande group / Stanford University.

GAFF citation: J. Wang, R. M. Wolf, J. W. Caldwell, P. A. Kollman, and D. A. Case, *Development and testing of a general amber force field*, **J. Comput. Chem. 25** (2004), 1157-1174, [doi:10.1002/jcc.20035](https://doi.org/10.1002/jcc.20035).

## OPLS-AA

Sources: [`OPLSAA-DB`](https://github.com/leelasd/OPLSAA-DB/tree/b09d51bf1bf5abc698a5b32e6c1715a625190a0b), commit `b09d51bf1bf5abc698a5b32e6c1715a625190a0b` (2018-08-02), supplies 526 native `db:` keys from OpenMM XML exports. Its [license](https://github.com/leelasd/OPLSAA-DB/blob/b09d51bf1bf5abc698a5b32e6c1715a625190a0b/LICENSE) is Apache-2.0; attribution: Leela S. Dodda, Jorgensen Lab, Yale University. [`TUK-FFDat`, version 1.0 (2023)](https://doi.org/10.5281/zenodo.8116422) supplies 125 `tuk:` chemical tags from `TUK-FFDat_OPLS-AA.xlsx`, intermolecular sheet, rows 2-126, under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Authors: G. Kanagalingam, S. Schmitt, F. Fleckenstein, and S. Stephan; the workbook SHA256 matches the module's source notice.

Force-field citation: W. L. Jorgensen, D. S. Maxwell, and J. Tirado-Rives, *Development and Testing of the OPLS All-Atom Force Field on Conformational Energetics and Properties of Organic Liquids*, **J. Am. Chem. Soc. 118** (1996), 11225-11236, [doi:10.1021/ja9621760](https://doi.org/10.1021/ja9621760). Data-scheme citation: G. Kanagalingam, S. Schmitt, F. Fleckenstein, and S. Stephan, *Data scheme and data format for transferable force fields for molecular simulation*, **Scientific Data 10** (2023), 495, [doi:10.1038/s41597-023-02369-8](https://doi.org/10.1038/s41597-023-02369-8). Individual parameter-reference DOIs remain beside the TUK rows.

## Transition metals

Source:
[`VeloxChem/src/pymodule/tmparameters.py`](https://github.com/VeloxChem/VeloxChem/blob/d4290334938fb5cab7e2e9954bd6dea5865d7b1c/src/pymodule/tmparameters.py), commit `d4290334938fb5cab7e2e9954bd6dea5865d7b1c` (2026-10-05). The file carries the complete BSD-3-Clause notice, `Copyright 2018-2025 VeloxChem developers`.

Citation: F. Sebesta, V. Slama, J. Melcr, Z. Futera, and J. V. Burda, *Estimation of Transition-Metal Empirical Parameters for Molecular Mechanical Force Fields*, **J. Chem. Theory Comput. 12** (2016), 3681-3688, [doi:10.1021/acs.jctc.6b00416](https://doi.org/10.1021/acs.jctc.6b00416).

## SPC/E water

Permissive reference: the `cspce` entry in [`VeloxChem/src/pymodule/waterparameters.py`](https://github.com/VeloxChem/VeloxChem/blob/d4290334938fb5cab7e2e9954bd6dea5865d7b1c/src/pymodule/waterparameters.py), at the same revision and with the same BSD-3-Clause notice as the TM source. The nonzero hydrogen LJ terms identify a modified SPC/E model.

Base model citation: H. J. C. Berendsen, J. R. Grigera, and T. P. Straatsma, *The Missing Term in Effective Pair Potentials*, **J. Phys. Chem. 91** (1987), 6269-6271, [doi:10.1021/j100308a038](https://doi.org/10.1021/j100308a038). Upstream's `cspce` entry cites T. Luchko, S. Gusarov, D. R. Roe, C. Simmerling, D. A. Case, J. Tuszynski, and A. Kovalenko, *Three-Dimensional Molecular Theory of Solvation Coupled with Molecular Dynamics in Amber*, **J. Chem. Theory Comput. 6** (2010), 607-624, [doi:10.1021/ct900460m](https://doi.org/10.1021/ct900460m).
