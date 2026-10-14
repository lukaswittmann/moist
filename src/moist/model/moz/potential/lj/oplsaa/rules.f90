! SPDX-License-Identifier: MPL-2.0 AND MIT
!> Foyer OPLS-AA environment definitions and independent LJ lookup keys
!>
!> Only names, SMARTS definitions and precedence are imported from Foyer
!> Numerical parameters continue to come from OPLSAA-DB and TUK-FFDat
!> Keys with different type IDs select equal element/sigma/epsilon values
!> They are LJ aliases, not equivalent charges or bonded force-field types
!> Missing parameter matches retain a blank key and fail automatic typing
module moist_model_moz_potential_lj_oplsaa_rules
   implicit none(type, external)
   private

   public :: oplsaa_rule_type, oplsaa_rules, oplsaa_type_len

   !> Length of the original Foyer OPLS-AA type label
   integer, parameter :: oplsaa_type_len = 10

   !> Atom environment definition with precedence and an LJ-only mapping
   type :: oplsaa_rule_type
      !> Original rule name
      character(len=oplsaa_type_len) :: name
      !> Modified SMARTS atom environment
      character(len=96) :: smarts
      !> Comma-separated names superseded by this rule
      character(len=96) :: overrides
      !> Existing LJ key; blank when its parameter pair is unavailable
      character(len=32) :: key
   end type oplsaa_rule_type

   ! Foyer: https://github.com/mosdef-hub/foyer
   ! Source: foyer/forcefields/xml/oplsaa.xml, retrieved 2026-10-08
   ! Source SHA256: c78ccb763cda33e3456a3f8b4ed8a8f0361b10c92c33aa02f2a9abf618c163e4
   ! Rules have their own Foyer SMARTS definitions; the XML parameter sections
   ! identify GROMACS as their upstream source and are not copied here
   ! Modifications: extracted definitions, converted to Fortran, added LJ keys
   ! Keys were verified against the independently sourced parameter table
   ! Zero-epsilon sites are equivalent regardless of the inert sigma value
   ! Matching criteria: element and sigma/epsilon within 1e-7 nm and kJ/mol
   ! Prefer the original type ID, then the lowest numerical compatible ID
   ! Preserve the original name separately to avoid implying full-FF identity
   !
   ! The MIT License (MIT)
   !
   ! Copyright (c) 2015 Vanderbilt University
   !
   ! Permission is hereby granted, free of charge, to any person obtaining a copy of
   ! this software and associated documentation files (the "Software"), to deal in
   ! the Software without restriction, including without limitation the rights to
   ! use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
   ! the Software, and to permit persons to whom the Software is furnished to do so,
   ! subject to the following conditions:
   ! The above copyright notice and this permission notice shall be included in all
   ! copies or substantial portions of the Software.
   ! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   ! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
   ! FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
   ! COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
   ! IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   ! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

   !> Environment rules 1-60
   type(oplsaa_rule_type), parameter :: rules_1(60) = [ &
      & oplsaa_rule_type("opls_135", &
      & "[C;X4](C)(H)(H)H", &
      & "", "db:opls_135"), &
      & oplsaa_rule_type("opls_136", &
      & "[C;X4](C)(C)(H)H", &
      & "", "db:opls_136"), &
      & oplsaa_rule_type("opls_137", &
      & "[C;X4](C)(C)(C)H", &
      & "", "db:opls_137"), &
      & oplsaa_rule_type("opls_138", &
      & "[C;X4](H)(H)(H)H", &
      & "", "db:opls_58"), &
      & oplsaa_rule_type("opls_139", &
      & "[C;X4](C)(C)(C)C", &
      & "", "db:opls_139"), &
      & oplsaa_rule_type("opls_140", &
      & "H[C;X4]", &
      & "", "db:opls_140"), &
      & oplsaa_rule_type("opls_141", &
      & "[C;X3;!r6](C)(C)C", &
      & "", "db:opls_141"), &
      & oplsaa_rule_type("opls_142", &
      & "[C;X3](C)(C)H", &
      & "", "db:opls_142"), &
      & oplsaa_rule_type("opls_143", &
      & "[C;X3](C)(H)H", &
      & "", "db:opls_143"), &
      & oplsaa_rule_type("opls_144", &
      & "[H][C;X3]", &
      & "", "db:opls_144"), &
      & oplsaa_rule_type("opls_145", &
      & "[C;X3;r6]1[C;X3;r6][C;X3;r6][C;X3;r6][C;X3;r6][C;X3;r6]1", &
      & "opls_141,opls_142", "db:opls_145"), &
      & oplsaa_rule_type("opls_146", &
      & "[H][C;%opls_145]", &
      & "opls_144", "db:opls_146"), &
      & oplsaa_rule_type("opls_147", &
      & "[C;R2;X3]([C,%opls_263])([C;R2;X3])", &
      & "opls_145", "db:opls_147"), &
      & oplsaa_rule_type("opls_148", &
      & "[C;X4]([C;%opls_145])(H)(H)H", &
      & "opls_149,opls_135", "db:opls_148"), &
      & oplsaa_rule_type("opls_149", &
      & "[C;X4]([C;%opls_145])(H)(H)*", &
      & "opls_136", "db:opls_149"), &
      & oplsaa_rule_type("opls_151", &
      & "[Cl][C;X4](C)(H)(H)", &
      & "opls_401", "db:opls_151"), &
      & oplsaa_rule_type("opls_151a", &
      & "[Cl][C;X4]([Cl])([C,H])(H)", &
      & "opls_401", "db:opls_151"), &
      & oplsaa_rule_type("opls_151b", &
      & "[Cl][C;X4]([Cl])([Cl])([C,H])", &
      & "opls_401", "db:opls_151"), &
      & oplsaa_rule_type("opls_151c", &
      & "[Cl][C;X4]([Cl])([F])([C,H])", &
      & "opls_401", "db:opls_151"), &
      & oplsaa_rule_type("opls_152", &
      & "[C;X4]([Cl])([C,H])(H)(H)", &
      & "", "db:opls_152"), &
      & oplsaa_rule_type("opls_152a", &
      & "[C;X4]([I])(C)([C,H])H", &
      & "", "db:opls_58"), &
      & oplsaa_rule_type("opls_152b", &
      & "[C;X4]([Cl])([Cl])([C,H])(H)", &
      & "", "db:opls_58"), &
      & oplsaa_rule_type("opls_152c", &
      & "[C;X4]([Cl])([Cl])([Cl])([C,H])", &
      & "", "db:opls_58"), &
      & oplsaa_rule_type("opls_152d", &
      & "[C;X4]([Br])([C,H])(H)(H)", &
      & "", "db:opls_58"), &
      & oplsaa_rule_type("opls_152e", &
      & "[C;X4]([Br])([Br])([C,H])(H)", &
      & "", "db:opls_58"), &
      & oplsaa_rule_type("opls_153", &
      & "[H][C;%opls_152,%opls_152d]", &
      & "opls_140", "db:opls_153"), &
      & oplsaa_rule_type("opls_153a", &
      & "[H]([C;X4;%opls_152a])", &
      & "opls_140", "db:opls_140"), &
      & oplsaa_rule_type("opls_153b", &
      & "[H]([C;X4;%opls_152b,%opls_152e,%opls_959])", &
      & "opls_140", "db:opls_140"), &
      & oplsaa_rule_type("opls_153c", &
      & "[H]([C;X4;%opls_152c])", &
      & "opls_140", "db:opls_140"), &
      & oplsaa_rule_type("opls_154", &
      & "[O;X2](H)([C])", &
      & "", "db:opls_154"), &
      & oplsaa_rule_type("opls_155", &
      & "H[O;%opls_154]", &
      & "", "db:opls_155"), &
      & oplsaa_rule_type("opls_156", &
      & "HC(H)(H)OH", &
      & "opls_140", "db:opls_140"), &
      & oplsaa_rule_type("opls_157", &
      & "[C;X4]([H])([H])([*])[O;%opls_154]", &
      & "opls_136,opls_159,opls_158", "db:opls_157"), &
      & oplsaa_rule_type("opls_158", &
      & "[C;X4]([H])[O;%opls_154]", &
      & "opls_136,opls_159", "db:opls_158"), &
      & oplsaa_rule_type("opls_159", &
      & "[C;X4][O;%opls_154]", &
      & "opls_136", "db:opls_159"), &
      & oplsaa_rule_type("opls_166", &
      & "[C;X3;r6]1(OH)[C;X3;r6][C;X3;r6][C;X3;r6][C;X3;r6][C;X3;r6]1", &
      & "opls_145", "db:opls_166"), &
      & oplsaa_rule_type("opls_167", &
      & "O(H)[C;%opls_166]", &
      & "opls_154", "db:opls_167"), &
      & oplsaa_rule_type("opls_168", &
      & "H[O;%opls_167]", &
      & "opls_155", "db:opls_168"), &
      & oplsaa_rule_type("opls_171", &
      & "[O;X2]([%opls_173,%opls_174])(H)", &
      & "opls_154", "db:opls_162"), &
      & oplsaa_rule_type("opls_172", &
      & "[H;X1]([O;%opls_171])", &
      & "opls_155", "db:opls_155"), &
      & oplsaa_rule_type("opls_173", &
      & "[C;X4]([H])([H])([O;X2]([H]))([C;X4]([O;X2]([H]))([C;X4]([O;X2]([H]))))", &
      & "opls_157,opls_158", "db:opls_173"), &
      & oplsaa_rule_type("opls_174", &
      & "[C;X4]([H])([O;X2]([H]))([C;X4]([O;X2]))([C;X4]([O;X2]))", &
      & "opls_158", "db:opls_174"), &
      & oplsaa_rule_type("opls_176", &
      & "[H][%opls_173,%opls_174]", &
      & "opls_140", "db:opls_176"), &
      & oplsaa_rule_type("opls_179", &
      & "[O;X2]([C;X3;r6])([C;X4])", &
      & "opls_180", "db:opls_179"), &
      & oplsaa_rule_type("opls_180", &
      & "[O;X2](C)C", &
      & "", "db:opls_180"), &
      & oplsaa_rule_type("opls_181", &
      & "[C;X4]([O;%opls_180])(H)(H)H", &
      & "opls_182", "db:opls_181"), &
      & oplsaa_rule_type("opls_182", &
      & "[C;X4]([O;%opls_180])(H)(H)C", &
      & "", "db:opls_182"), &
      & oplsaa_rule_type("opls_183", &
      & "[C;X4]([O;%opls_180])(C)(C)H", &
      & "", "db:opls_183"), &
      & oplsaa_rule_type("opls_184", &
      & "[C;X4]([O;%opls_180])(C)(C)(C)", &
      & "", "db:opls_184"), &
      & oplsaa_rule_type("opls_185", &
      & "HC[O;%opls_180]", &
      & "opls_144,opls_176,opls_140", "db:opls_185"), &
      & oplsaa_rule_type("opls_186", &
      & "[O;X2][C;%opls_193](O)(C)H", &
      & "opls_180", "db:opls_186"), &
      & oplsaa_rule_type("opls_189", &
      & "[C;X4]([O;%opls_180])([O;%opls_180])(H)(H)", &
      & "", "db:opls_189"), &
      & oplsaa_rule_type("opls_190", &
      & "H([%opls_189])", &
      & "opls_185", "db:opls_190"), &
      & oplsaa_rule_type("opls_193", &
      & "[C;X4](O)(O)(C)H", &
      & "", "db:opls_193"), &
      & oplsaa_rule_type("opls_194", &
      & "H[C;%opls_193][O;%opls_186]", &
      & "opls_185", "db:opls_194"), &
      & oplsaa_rule_type("opls_199", &
      & "[C;X3;r6]([O;X2]([C;X4](H)(H)(H)))([C;X3;r6])([C;X3;r6])", &
      & "opls_145", "db:opls_199"), &
      & oplsaa_rule_type("opls_200", &
      & "[S;X2]H", &
      & "opls_202", "db:opls_200"), &
      & oplsaa_rule_type("opls_201", &
      & "[S;X2](H)(H)", &
      & "opls_200,opls_202", ""), &
      & oplsaa_rule_type("opls_202", &
      & "[S;X2]", &
      & "", "db:opls_202"), &
      & oplsaa_rule_type("opls_203", &
      & "[S;X2](S)", &
      & "opls_200,opls_202", "db:opls_474")]

   !> Environment rules 61-120
   type(oplsaa_rule_type), parameter :: rules_2(60) = [ &
      & oplsaa_rule_type("opls_204", &
      & "[H;X1]([S;%opls_200])", &
      & "", "db:opls_204"), &
      & oplsaa_rule_type("opls_205", &
      & "[H;X1]([S;%opls_201])", &
      & "opls_204", "db:opls_155"), &
      & oplsaa_rule_type("opls_206", &
      & "[C;X4]([S;%opls_200])(H)(H)", &
      & "opls_207,opls_208,opls_210", "db:opls_206"), &
      & oplsaa_rule_type("opls_207", &
      & "[C;X4]([S;%opls_200])(H)([!H])([!H])", &
      & "opls_208,opls_211", "db:opls_58"), &
      & oplsaa_rule_type("opls_208", &
      & "[C;X4]([S;%opls_200])", &
      & "opls_212", "db:opls_58"), &
      & oplsaa_rule_type("opls_209", &
      & "[C;X4]([S;%opls_202])(H)(H)(H)", &
      & "opls_210,opls_211,opls_212", "db:opls_209"), &
      & oplsaa_rule_type("opls_210", &
      & "[C;X4]([S;%opls_202])(H)(H)", &
      & "opls_211,opls_212", "db:opls_210"), &
      & oplsaa_rule_type("opls_211", &
      & "[C;X4]([S;%opls_202])(H)", &
      & "opls_212", "db:opls_58"), &
      & oplsaa_rule_type("opls_212", &
      & "[C;X4]([S;%opls_202])", &
      & "", "db:opls_58"), &
      & oplsaa_rule_type("opls_213", &
      & "[C;X4]([S;%opls_203])(H)(H)(H)", &
      & "opls_214,opls_215,opls_216,opls_209", "db:opls_213"), &
      & oplsaa_rule_type("opls_214", &
      & "[C;X4]([S;%opls_203])(H)(H)", &
      & "opls_215,opls_216,opls_210", "db:opls_214"), &
      & oplsaa_rule_type("opls_215", &
      & "[C;X4]([S;%opls_203])(H)", &
      & "opls_216,opls_211", "db:opls_58"), &
      & oplsaa_rule_type("opls_216", &
      & "[C;X4][S;%opls_203]", &
      & "opls_212", "db:opls_58"), &
      & oplsaa_rule_type("opls_217", &
      & "[C;X4]([S;%opls_200])(H)(H)(H)", &
      & "opls_211,opls_207,opls_209,opls_206", "db:opls_58"), &
      & oplsaa_rule_type("opls_218", &
      & "[C;X4]([C;%opls_221])(O)(H)(H)", &
      & "opls_149,opls_157,opls_158", "db:opls_218"), &
      & oplsaa_rule_type("opls_221", &
      & "[C;X3;r6]([C;X3;r6])([C;X3;r6])([C;X4](O)(H)(H))", &
      & "opls_145", "db:opls_221"), &
      & oplsaa_rule_type("opls_226", &
      & "[Cl;X1]([C;%opls_227])", &
      & "opls_401", "db:opls_226"), &
      & oplsaa_rule_type("opls_227", &
      & "[C;X3;!r6;!r5]([C;X3;!r6;!r5](H)(H))(Cl)(Cl)", &
      & "", "db:opls_227"), &
      & oplsaa_rule_type("opls_232", &
      & "[C;X3]([O;X1])([C;X3;r6])[H]", &
      & "opls_277", "db:opls_232"), &
      & oplsaa_rule_type("opls_233", &
      & "[C;X3]([O;X1])([C;X3;r6])(C)", &
      & "opls_280", "db:opls_233"), &
      & oplsaa_rule_type("opls_235", &
      & "[C;X3]([O;X1])[N;X3]", &
      & "opls_277", "db:opls_235"), &
      & oplsaa_rule_type("opls_236", &
      & "O[C;%opls_235]", &
      & "opls_278", "db:opls_236"), &
      & oplsaa_rule_type("opls_237", &
      & "[N;X3](H)(H)[C;%opls_235]", &
      & "opls_900", "db:opls_237"), &
      & oplsaa_rule_type("opls_238", &
      & "[N;X3](H)(C)[C;%opls_235]", &
      & "", "db:opls_238"), &
      & oplsaa_rule_type("opls_239", &
      & "[N;X3](C)(C)[C;%opls_235]", &
      & "", "db:opls_239"), &
      & oplsaa_rule_type("opls_240", &
      & "[H;X1][N;%opls_237]", &
      & "opls_909", "db:opls_240"), &
      & oplsaa_rule_type("opls_241", &
      & "[H;X1][N;%opls_238]", &
      & "", "db:opls_241"), &
      & oplsaa_rule_type("opls_242", &
      & "[C;X4](H)(H)(H)[N;%opls_238]", &
      & "", "db:opls_242"), &
      & oplsaa_rule_type("opls_243", &
      & "[C;X4](H)(H)(H)[N;%opls_239]", &
      & "", "db:opls_243"), &
      & oplsaa_rule_type("opls_245", &
      & "[C;X4](C)(H)(H)[N;%opls_239]", &
      & "", "db:opls_245"), &
      & oplsaa_rule_type("opls_260", &
      & "[C;X3;r6]([C;X2;%opls_261])", &
      & "opls_145", "db:opls_260"), &
      & oplsaa_rule_type("opls_261", &
      & "[C;X2]([C;X3;r6])", &
      & "opls_754", "db:opls_261"), &
      & oplsaa_rule_type("opls_262", &
      & "[N;X1][C;%opls_261]", &
      & "opls_753", "db:opls_262"), &
      & oplsaa_rule_type("opls_263", &
      & "[C;X3;r6]1(Cl)[C;X3;r6][C;X3;r6][C;X3;r6][C;X3;r6][C;X3;r6]1", &
      & "opls_145", "db:opls_263"), &
      & oplsaa_rule_type("opls_264", &
      & "[Cl;X1]([C;%opls_263])", &
      & "opls_401", "db:opls_264"), &
      & oplsaa_rule_type("opls_267", &
      & "[C;X3]([O;X1])OH", &
      & "opls_277,opls_271,opls_465", "db:opls_267"), &
      & oplsaa_rule_type("opls_268", &
      & "[O;X2]([C;%opls_267])H", &
      & "opls_154", "db:opls_268"), &
      & oplsaa_rule_type("opls_269", &
      & "[O;X1]([C;%opls_267](OH))", &
      & "opls_278,opls_272", "db:opls_269"), &
      & oplsaa_rule_type("opls_270", &
      & "H([O;%opls_268])", &
      & "opls_155", "db:opls_270"), &
      & oplsaa_rule_type("opls_271", &
      & "[C;X3]([O;X1])([O;X1])", &
      & "", "db:opls_231"), &
      & oplsaa_rule_type("opls_272", &
      & "[O;X1]([C;%opls_271])", &
      & "", "db:opls_236"), &
      & oplsaa_rule_type("opls_277", &
      & "[C;X3]([O;X1])H", &
      & "", "db:opls_277"), &
      & oplsaa_rule_type("opls_278", &
      & "[O;X1][C;%opls_277]", &
      & "", "db:opls_278"), &
      & oplsaa_rule_type("opls_279", &
      & "H([C;X3]([O;X1])([N,C,O,H]))", &
      & "opls_185,opls_144,opls_485", "db:opls_279"), &
      & oplsaa_rule_type("opls_280", &
      & "[C;X3]([O;X1])(C)C", &
      & "", "db:opls_280"), &
      & oplsaa_rule_type("opls_281", &
      & "[O;X1]([C;%opls_280])", &
      & "", "db:opls_281"), &
      & oplsaa_rule_type("opls_282", &
      & "HC[C;%opls_277,%opls_280,%opls_465;!%opls_267;!%opls_233]", &
      & "opls_140,opls_144", "db:opls_282"), &
      & oplsaa_rule_type("opls_401", &
      & "[Cl]", &
      & "", ""), &
      & oplsaa_rule_type("opls_406", &
      & "Li", &
      & "", ""), &
      & oplsaa_rule_type("opls_440", &
      & "[P;X4]([O;X1])([O;X2])([O;X2])", &
      & "", "db:opls_781"), &
      & oplsaa_rule_type("opls_441", &
      & "[O;X1][P;X4]", &
      & "", ""), &
      & oplsaa_rule_type("opls_442", &
      & "[O;X2]([P;X4]([O;X1])([O;X2])([O;X2]))[#6]", &
      & "", "db:opls_177"), &
      & oplsaa_rule_type("opls_443", &
      & "[C;X4][O;X2][P;X4]([O;X2])([O;X2])(O)", &
      & "", "db:opls_58"), &
      & oplsaa_rule_type("opls_444", &
      & "[H][C;%opls_443;X4]", &
      & "opls_140,opls_176", "db:opls_140"), &
      & oplsaa_rule_type("opls_465", &
      & "[C;X3]([O;X1])([O;X2])", &
      & "opls_277", "db:opls_465"), &
      & oplsaa_rule_type("opls_466", &
      & "[O;X1]([C;%opls_465]([O;%opls_467]))", &
      & "opls_278", "db:opls_466"), &
      & oplsaa_rule_type("opls_467", &
      & "[O;X2]([C;%opls_465])([!H])", &
      & "opls_180", "db:opls_467"), &
      & oplsaa_rule_type("opls_468", &
      & "[C;X4]([O;%opls_467])(H)(H)H", &
      & "opls_181", "db:opls_468"), &
      & oplsaa_rule_type("opls_469", &
      & "H[C;%opls_468,%opls_490]", &
      & "opls_185,opls_140", "db:opls_469"), &
      & oplsaa_rule_type("opls_471", &
      & "[C;X3]([O;X2]([C;X4](H)(H)(H)))([O;X1])([C;X3;r6])", &
      & "opls_465", "db:opls_471")]

   !> Environment rules 121-180
   type(oplsaa_rule_type), parameter :: rules_3(60) = [ &
      & oplsaa_rule_type("opls_472", &
      & "[C;X3;r6]([O;%opls_473])", &
      & "opls_145,opls_147", "db:opls_472"), &
      & oplsaa_rule_type("opls_473", &
      & "[O;X2]([C;X3;r6])([C;X3;r6])", &
      & "opls_180", "db:opls_473"), &
      & oplsaa_rule_type("opls_484", &
      & "[C;X4]([S;X4](O)(O))(C)(H)(H)", &
      & "", "db:opls_58"), &
      & oplsaa_rule_type("opls_485", &
      & "H[!%opls_136;%opls_484]", &
      & "opls_140,opls_146", "db:opls_140"), &
      & oplsaa_rule_type("opls_490", &
      & "[C;X4]([O;%opls_467])(H)(H)C", &
      & "opls_182", "db:opls_490"), &
      & oplsaa_rule_type("opls_493", &
      & "[S;X4]([#6])([#6])([O;X1])([O;X1])", &
      & "", "db:opls_493"), &
      & oplsaa_rule_type("opls_494", &
      & "[O;X1][S;X4]([O;X1])([#6])([#6])", &
      & "", "db:opls_494"), &
      & oplsaa_rule_type("opls_496", &
      & "[S;X3]([O;%opls_497])(C)(C)", &
      & "", "db:opls_496"), &
      & oplsaa_rule_type("opls_497", &
      & "[O;X1][S;X3]", &
      & "", "db:opls_497"), &
      & oplsaa_rule_type("opls_498", &
      & "[C;X4]([S;X3])(H)(H)(H)", &
      & "", "db:opls_498"), &
      & oplsaa_rule_type("opls_518", &
      & "[C;X3;!r6]([O;X2])([C;X3;!r6])([H])", &
      & "", "db:opls_141"), &
      & oplsaa_rule_type("opls_520", &
      & "[N;X2;r6]1[C;X3;r6][C;X3;r6][C;X3;r6][C;X3;r6][C;X3;r6]1", &
      & "", "db:opls_520"), &
      & oplsaa_rule_type("opls_521", &
      & "[C;X3;r6][N;%opls_520]", &
      & "", "db:opls_521"), &
      & oplsaa_rule_type("opls_522", &
      & "[C;X3;r6][C;%opls_521]", &
      & "opls_142", "db:opls_522"), &
      & oplsaa_rule_type("opls_523", &
      & "[C;X3;r6]([C;%opls_522])[C;%opls_522]", &
      & "opls_142,opls_141", "db:opls_523"), &
      & oplsaa_rule_type("opls_524", &
      & "H[C;%opls_521]", &
      & "opls_144", "db:opls_524"), &
      & oplsaa_rule_type("opls_525", &
      & "H[C;%opls_522]", &
      & "opls_144", "db:opls_525"), &
      & oplsaa_rule_type("opls_526", &
      & "H[C;%opls_523]", &
      & "opls_144", "db:opls_526"), &
      & oplsaa_rule_type("opls_530", &
      & "[N;X2;r6]1[C;X3;r6][C;X3;r6][C;X3;r6][N;X2;r6][C;X3;r6]1", &
      & "", "db:opls_530"), &
      & oplsaa_rule_type("opls_531", &
      & "[C;X3;r6]([N;%opls_530])[N;%opls_530]", &
      & "", "db:opls_531"), &
      & oplsaa_rule_type("opls_532", &
      & "[C;X3;r6]([N;%opls_530])[C;X3;r6]", &
      & "", "db:opls_532"), &
      & oplsaa_rule_type("opls_533", &
      & "[C;X3;r6]([C;%opls_532])[C;%opls_532]", &
      & "opls_142", "db:opls_533"), &
      & oplsaa_rule_type("opls_534", &
      & "H[C;%opls_531]", &
      & "opls_144", "db:opls_534"), &
      & oplsaa_rule_type("opls_535", &
      & "H[C;%opls_532]", &
      & "opls_144", "db:opls_535"), &
      & oplsaa_rule_type("opls_536", &
      & "H[C;%opls_533]", &
      & "opls_144", "db:opls_536"), &
      & oplsaa_rule_type("opls_542", &
      & "[N;X3;r5]1[C;X3;r5][C;X3;r5][C;X3;r5][C;X3;r5]1H", &
      & "", "db:opls_542"), &
      & oplsaa_rule_type("opls_543", &
      & "[C;X3;r5]([N;%opls_542])", &
      & "", "db:opls_543"), &
      & oplsaa_rule_type("opls_544", &
      & "[C;X3;r5]([C;%opls_543])", &
      & "opls_142,opls_141", "db:opls_544"), &
      & oplsaa_rule_type("opls_545", &
      & "H[N;%opls_542]", &
      & "", "db:opls_545"), &
      & oplsaa_rule_type("opls_546", &
      & "H[C;%opls_543]", &
      & "opls_144", "db:opls_546"), &
      & oplsaa_rule_type("opls_547", &
      & "H[C;%opls_544]", &
      & "opls_144", "db:opls_547"), &
      & oplsaa_rule_type("opls_566", &
      & "[O;X2;r5]([C;X3;r5])([C;X3;r5])", &
      & "opls_180", "db:opls_566"), &
      & oplsaa_rule_type("opls_567", &
      & "[C;X3;r5]([%opls_633,%opls_566])", &
      & "opls_518", "db:opls_567"), &
      & oplsaa_rule_type("opls_568", &
      & "[C;X3;r5]([C;%opls_567])([C;X3;r5])", &
      & "opls_142", "db:opls_568"), &
      & oplsaa_rule_type("opls_569", &
      & "[H]([%opls_567])", &
      & "opls_185,opls_144", "db:opls_569"), &
      & oplsaa_rule_type("opls_570", &
      & "[H]([%opls_568])", &
      & "opls_144", "db:opls_570"), &
      & oplsaa_rule_type("opls_633", &
      & "[S;r5;X2]([C;r5;X3])([C;r5;X3])", &
      & "opls_202", "db:opls_474"), &
      & oplsaa_rule_type("opls_670", &
      & "[C;X4]([C;%opls_521])(H)(H)(H)", &
      & "opls_135", "db:opls_58"), &
      & oplsaa_rule_type("opls_672", &
      & "[C;X4]([C;%opls_522])(H)(H)(H)", &
      & "opls_135", "db:opls_58"), &
      & oplsaa_rule_type("opls_674", &
      & "[C;X4]([C;X3;r6]([C;X3;r6]([C;X3;r6](N))))(H)(H)(H)", &
      & "opls_135,opls_148", "db:opls_58"), &
      & oplsaa_rule_type("opls_678", &
      & "[C;X4](H)(H)(H)[C;%opls_543]", &
      & "opls_679", "db:opls_58"), &
      & oplsaa_rule_type("opls_679", &
      & "[C;X4](H)(H)[C;%opls_543]", &
      & "opls_136", "db:opls_679"), &
      & oplsaa_rule_type("opls_711", &
      & "[C;X4;r3]1(H)(H)[C;X4;r3][C;X4;r3]1", &
      & "opls_136,opls_712", "db:opls_58"), &
      & oplsaa_rule_type("opls_712", &
      & "[C;X4;r3]1(H)[C;X4;r3][C;X4;r3]1", &
      & "opls_137,opls_713", "db:opls_58"), &
      & oplsaa_rule_type("opls_713", &
      & "[C;X4;r3]1[C;X4;r3][C;X4;r3]1", &
      & "", "db:opls_58"), &
      & oplsaa_rule_type("opls_718", &
      & "[C;X3;r6](F)([C;%opls_145])([C;%opls_145])", &
      & "opls_145", "db:opls_718"), &
      & oplsaa_rule_type("opls_719", &
      & "F[C;%opls_718]", &
      & "opls_965", "db:opls_719"), &
      & oplsaa_rule_type("opls_720", &
      & "[C;X3;r6]([F])([C;r6](F))([C;r6](F))", &
      & "opls_145,opls_718,opls_727", "db:opls_145"), &
      & oplsaa_rule_type("opls_721", &
      & "[F;X1]([C;%opls_720])", &
      & "opls_719,opls_728", "db:opls_719"), &
      & oplsaa_rule_type("opls_722", &
      & "[Br;X1]([C;%opls_152d])", &
      & "", "db:opls_730"), &
      & oplsaa_rule_type("opls_722a", &
      & "[Br;X1]([C;%opls_152e])", &
      & "", "db:opls_730"), &
      & oplsaa_rule_type("opls_724", &
      & "[C;X3;r6]([C;%opls_725])([C;X3;r6])([C;X3;r6])", &
      & "opls_145", "db:opls_724"), &
      & oplsaa_rule_type("opls_725", &
      & "[C;X4](F)(F)(F)([C;X3;r6])", &
      & "opls_672,opls_961", "db:opls_725"), &
      & oplsaa_rule_type("opls_726", &
      & "[F;X1]([C;%opls_725])", &
      & "opls_965", "db:opls_726"), &
      & oplsaa_rule_type("opls_727", &
      & "[C;X3;r6](F)([C;X3;r6](F))([%opls_145])", &
      & "opls_145,opls_718", "db:opls_727"), &
      & oplsaa_rule_type("opls_728", &
      & "[F;X1]([C;%opls_727])", &
      & "opls_719", "db:opls_728"), &
      & oplsaa_rule_type("opls_732", &
      & "[I;X1][C;X4](C)(C)(H)", &
      & "", ""), &
      & oplsaa_rule_type("opls_734", &
      & "[S;X2]([C;X3;r6])(H)", &
      & "opls_200", "db:opls_474"), &
      & oplsaa_rule_type("opls_735", &
      & "[C;X3;r6]([S]([H]))([C;X3;r6])([C;X3;r6])", &
      & "opls_145", "db:opls_735"), &
      & oplsaa_rule_type("opls_753", &
      & "[N;X1]C", &
      & "", "db:opls_753")]

   !> Environment rules 181-228
   type(oplsaa_rule_type), parameter :: rules_4(48) = [ &
      & oplsaa_rule_type("opls_754", &
      & "[C;X2][N;%opls_753]", &
      & "", "db:opls_754"), &
      & oplsaa_rule_type("opls_755", &
      & "[C;X4](H)(H)(H)[C;%opls_754]", &
      & "opls_135,opls_756", "db:opls_754"), &
      & oplsaa_rule_type("opls_756", &
      & "[C;X4](H)(H)[C;%opls_754]", &
      & "opls_136,opls_757", "db:opls_756"), &
      & oplsaa_rule_type("opls_757", &
      & "[C;X4](H)[C;%opls_754]", &
      & "opls_137,opls_758", "db:opls_757"), &
      & oplsaa_rule_type("opls_758", &
      & "[C;X4][C;%opls_754]", &
      & "", "db:opls_754"), &
      & oplsaa_rule_type("opls_759", &
      & "HC[C;%opls_754]", &
      & "opls_140", "db:opls_759"), &
      & oplsaa_rule_type("opls_760", &
      & "[N;X3]([O;X1])([O;X1])", &
      & "opls_239", "db:opls_760"), &
      & oplsaa_rule_type("opls_761", &
      & "[O;X1][N][O;X1]", &
      & "", "db:opls_761"), &
      & oplsaa_rule_type("opls_762", &
      & "[C;X4](H)(H)(H)[N;%opls_760]", &
      & "", "db:opls_762"), &
      & oplsaa_rule_type("opls_763", &
      & "H([C;X4][N;%opls_760])", &
      & "opls_140", "db:opls_763"), &
      & oplsaa_rule_type("opls_764", &
      & "[C;X4](C)(H)(H)[N;%opls_760]", &
      & "opls_136", "db:opls_764"), &
      & oplsaa_rule_type("opls_765", &
      & "[C;X4]([N;X3](O)(O))([C;X4])([C;X4])", &
      & "", "db:opls_58"), &
      & oplsaa_rule_type("opls_767", &
      & "[N;X3]([O;X1])([O;X1])[C;X3;r6]", &
      & "opls_760", "db:opls_767"), &
      & oplsaa_rule_type("opls_768", &
      & "[C;X3;r6][N;%opls_767]", &
      & "opls_145", "db:opls_768"), &
      & oplsaa_rule_type("opls_771", &
      & "[O;X1][C;X3;%opls_772]", &
      & "opls_466", "db:opls_771"), &
      & oplsaa_rule_type("opls_772", &
      & "[C;X3;r5]1[O;X2;r5][C;X4;r5][C;X4;r5][O;X2;r5]1", &
      & "opls_465", "db:opls_772"), &
      & oplsaa_rule_type("opls_773", &
      & "[O;X2;r5]1[C;X4;r5][C;X4;r5][O;X2;r5][C;X3;r5]1", &
      & "opls_180,opls_467", "db:opls_773"), &
      & oplsaa_rule_type("opls_774", &
      & "[C;X4;r5]1(H)(H)[C;X4;r5][O;X2;r5][C;X3;r5][O;X2;r5]1", &
      & "opls_182,opls_490", "db:opls_774"), &
      & oplsaa_rule_type("opls_775", &
      & "[C;X4;r5]1(H)(C)[C;X4;r5][O;X2;r5][C;X3;r5][O;X2;r5]1", &
      & "opls_182,opls_183", "db:opls_775"), &
      & oplsaa_rule_type("opls_776", &
      & "[C;X4](H)(H)(H)[C;%opls_775]", &
      & "opls_135", "db:opls_776"), &
      & oplsaa_rule_type("opls_777", &
      & "H[C;%opls_774]", &
      & "opls_185,opls_140,opls_469", "db:opls_777"), &
      & oplsaa_rule_type("opls_778", &
      & "H[C;%opls_775]", &
      & "opls_185,opls_140", "db:opls_778"), &
      & oplsaa_rule_type("opls_779", &
      & "H[C;%opls_776]", &
      & "opls_185,opls_140", "db:opls_779"), &
      & oplsaa_rule_type("opls_900", &
      & "[N;X3](H)(H)C", &
      & "", "db:opls_900"), &
      & oplsaa_rule_type("opls_901", &
      & "[N;X3](H)([C;!%opls_235;!%opls_543])([C;!%opls_235;!%opls_543])", &
      & "", "db:opls_901"), &
      & oplsaa_rule_type("opls_902", &
      & "[N;X3]([C;X4])([C;X4])([C;X4])", &
      & "", "db:opls_902"), &
      & oplsaa_rule_type("opls_903", &
      & "[C;X4](H)(H)(H)([N;%opls_900])", &
      & "opls_906", "db:opls_903"), &
      & oplsaa_rule_type("opls_904", &
      & "[C;X4](H)(H)(H)([N;%opls_901])", &
      & "opls_906", "db:opls_904"), &
      & oplsaa_rule_type("opls_906", &
      & "[C;X4]([N;%opls_900])(H)(H)", &
      & "opls_136", "db:opls_906"), &
      & oplsaa_rule_type("opls_907", &
      & "[C;X4]([N;%opls_901])(C)(H)(H)", &
      & "", "db:opls_907"), &
      & oplsaa_rule_type("opls_908", &
      & "[C;X4]([H])([H])([C;X4])([N;X3]([C;X4])([C;X4]))", &
      & "", "db:opls_908"), &
      & oplsaa_rule_type("opls_909", &
      & "H[N;%opls_900]", &
      & "", "db:opls_909"), &
      & oplsaa_rule_type("opls_910", &
      & "H[N;%opls_901]", &
      & "", "db:opls_910"), &
      & oplsaa_rule_type("opls_911", &
      & "H([C;%opls_903,%opls_904,%opls_905,%opls_906,%opls_908,%opls_907,%opls_912,%opls_914])", &
      & "opls_140", "db:opls_911"), &
      & oplsaa_rule_type("opls_912", &
      & "[C;X4](H)([C;X4])([C;X4])([N;%opls_900])", &
      & "", "db:opls_58"), &
      & oplsaa_rule_type("opls_913", &
      & "[C;X4]([N;%opls_900](H)(H))(C)(C)(C)", &
      & "", "db:opls_58"), &
      & oplsaa_rule_type("opls_914", &
      & "[C;X4]([N;%opls_901])([C;%opls_135])([C;%opls_135])(H)", &
      & "", "db:opls_58"), &
      & oplsaa_rule_type("opls_916", &
      & "[C;X3;r6]([N;X3](H)(H))([C;X3;r6]([Cl]))([C;X3;r6](H))", &
      & "opls_145,opls_147", "db:opls_916"), &
      & oplsaa_rule_type("opls_917", &
      & "[C;X3;r6]([N;%opls_901])", &
      & "opls_145", "db:opls_917"), &
      & oplsaa_rule_type("opls_925", &
      & "[C;X2](C)H", &
      & "", "db:opls_925"), &
      & oplsaa_rule_type("opls_926", &
      & "H[C;%opls_925]", &
      & "", "db:opls_926"), &
      & oplsaa_rule_type("opls_927", &
      & "[C;X2](C(H)H)C", &
      & "opls_928", "db:opls_927"), &
      & oplsaa_rule_type("opls_930", &
      & "H([*][C][C;%opls_925])", &
      & "opls_140,opls_144", "db:opls_930"), &
      & oplsaa_rule_type("opls_959", &
      & "[C;X4](F)(H)([Cl,Br])", &
      & "", "db:opls_959"), &
      & oplsaa_rule_type("opls_960", &
      & "[C;X4](F)([!H])([Cl,Br])([Cl,Br])", &
      & "", "db:opls_58"), &
      & oplsaa_rule_type("opls_961", &
      & "[C;X4](F)(F)(F)*", &
      & "opls_962", "db:opls_961"), &
      & oplsaa_rule_type("opls_962", &
      & "[C;X4](F)(F)(*)*", &
      & "", "db:opls_58"), &
      & oplsaa_rule_type("opls_965", &
      & "[F;X1]([C;X4])", &
      & "", "db:opls_965")]

   !> Bundled rules; 223 of 228 have a compatible LJ key
   type(oplsaa_rule_type), parameter :: oplsaa_rules(228) = [rules_1, rules_2, rules_3, rules_4]

end module moist_model_moz_potential_lj_oplsaa_rules
