## Fission product yields: ENDF-6 → G4NDL

Converts ENDF-6 fission product yield data (MF=8) into the file format Geant4
expects under `Fission/FF/`, so that the neutron High Precision (HP) package can
generate real fission fragments for isotopes that G4NDL does not ship.

The immediate motivation was ²⁴¹Am. G4NDL 4.7.1 contains only 11 files in
`Fission/FF/` (²²⁷·²²⁹·²³²Th, ²³¹Pa, ²³²⁻²³⁸U); every other fissile isotope
silently produces no fragments. The yields here come from [UKFY-4.1][UKFY],
distributed as part of the JEFF decay data sub-library.

### Contents

| Path | Description |
| --- | --- |
| `endf2g4ndl_ff.py` | The converter. No dependencies beyond the standard library. |
| `neutron/am241n.out` | UKFY-4.1 yields for ²⁴¹Am(n,f), ENDF-6 format, MAT 3605 |
| `neutron/am242mn.out` | UKFY-4.1 yields for ²⁴²ᵐAm(n,f) |
| `neutron/am243n.out` | UKFY-4.1 yields for ²⁴³Am(n,f) |
| `Fission/FF/` | Generated output, laid out ready to drop into a G4NDL tree |
| `endf2g4ndl_fission.py` | Converter for every other `Fission/` directory — see [below](#the-other-fission-directories) |
| `prepro_fission.sh` | Runs the PREPRO chain that `endf2g4ndl_fission.py` expects its input to have been through |
| `g4ndl_compare.py` | Token-by-token (and cross-section-as-function) comparison of G4NDL files or trees |
| `jeff33/` | JEFF-3.3 neutron tapes, Z ≥ 88, from the [IAEA-NDS mirror][JEFF33-IAEA] |
| `Fission/{CrossSection,FS,FC,SC,TC,LC}/` | Output of `endf2g4ndl_fission.py` for `jeff33/` |
| `92_235_uranium` | Reference file, the decompressed `92_235_Uranium.z` from G4NDL 4.7.1 |
| `Source.asc` | Where the UKFY files were downloaded from |

### Usage

```
    python3 endf2g4ndl_ff.py neutron/am241n.out Fission/FF/95_241_Americium --compress
```

The script prints a summary to stderr — one line per incident-energy group with
the product count and ∑Y, which should be 2 to within rounding:

```
    ZA=95241 (Z=95, A=241)  AWR=238.99000
    37 incident-energy points -> Fission/FF/95_241_Americium
      E =  1.00000e-05 eV    1545 products   sum(Y) = 1.999999
      E =  2.50000e-02 eV    1545 products   sum(Y) = 1.999999
      ...
```

Options:

- `--mt 459` reads cumulative rather than independent yields. Leave it at the
  default 454; `G4WendtFissionFragmentGenerator` asks for independent yields and
  ignores any section with a different MT.
- `--library NAME` sets the label written on line 2. Defaults to `UKFY4.1`.
- `--energies E1,E2,...` keeps only the listed incident energies (in eV, nearest
  available point wins). See [Trimming the energy grid](#trimming-the-energy-grid).
- `--compress` additionally writes `<output>.z`, the raw zlib stream that G4NDL
  ships and that Geant4 prefers when both files are present.

### Installing into G4NDL

The filename is not free-form. `G4ParticleHPNames::GetName()` builds the path as
`<G4NEUTRONHPDATA>/Fission/FF/<Z>_<A>_<Element>`, with the element name spelled
exactly as in the table at `G4ParticleHPNames.cc:71` — so ²⁴¹Am must be
`95_241_Americium`, capital A, no extension. A metastable *target* would be
`<Z>_<A>m<M>_<Element>`, e.g. `95_242m1_Americium`.

G4NDL installs are usually read-only, so copy the tree first:

```
    cp -r /opt/Geant4/Libraries/G4NDL4.7.1 ~/G4NDL4.7.1-local
    cp Fission/FF/95_241_Americium.z ~/G4NDL4.7.1-local/Fission/FF/
    export G4NEUTRONHPDATA=~/G4NDL4.7.1-local
```

Then enable the model that actually reads these files:

```
    /process/had/particle_hp/use_Wendt_fission_model true
```

:warning: This is **not** the same as `/process/had/particle_hp/produce_fission_fragment`.
The two are mutually exclusive — `G4ParticleHPChannel.cc:58` switches the latter
off when the Wendt model is enabled — and only the Wendt path reads
`Fission/FF/`. Without this command the generated file is never opened, and
because `G4WendtFissionFragmentGenerator::InitializeANucleus()` returns quietly
on a stream error, a wrong filename or an unset `G4NEUTRONHPDATA` fails silently
too.

### File format

The `FF/` section of *G4NDL Final State Guide* (Ford, 2016) describes the layout;
the authoritative version is `G4ENDFTapeRead::ReadInData()`, at
`source/processes/hadronic/models/particle_hp/src/G4ENDFTapeRead.cc:409-456` in
geant4-v11.4.2. The two agree.

```
    G4NDL                     <- optional, consumed by GetDataStream()
    UKFY4.1                   <- library label, only read if line 1 is "G4NDL"
           454                <- MT: 454 = independent yields, 459 = cumulative
             8                <- MF: read into a dummy, never used
      2.389900e+02            <- AWR: also read into a dummy, never used
            37                <- number of incident-energy blocks
      1.000000e-05            <- E_inc [eV]
          1545                <- number of products at this energy
             2                <- interpolation scheme (parsed, currently ignored)
         21060  0  0.000000e+00  ...   <- ZZZAAA, metastate, yield; x1545
      ...                                <- and so on, x37
```

Every value is read with `operator>>`, so whitespace and line breaks carry no
meaning. The layout written here (14-column fields, two product triplets per
line) simply matches the files G4NDL ships.

A few consequences of that parser worth knowing:

- The `G4NDL` / library-name preamble is optional. `G4ParticleHPManager::GetDataStream()`
  peeks at the first token; if it is not the literal string `G4NDL` it rewinds
  and the file is parsed from byte zero.
- Products are keyed on `ZZZAAA*10 + metastate`, and yields are renormalised at
  load time (`G4FissionProductYieldDist.cc:1101`). Raw ENDF yields summing to 2
  are therefore fine as-is — do not pre-normalise to 1.
- Multiple energy blocks are genuinely used: `G4FissionProductYieldDist.cc:661-824`
  brackets the incident energy and interpolates the yields linearly between the
  two neighbouring groups.

### Trimming the energy grid

UKFY-4.1 tabulates 37 incident energies reaching 150 MeV. Geant4's HP package
stops at 20 MeV, so the 13 groups above that are unreachable; they cost a little
memory in `MakeTrees()` and nothing else. To keep only the useful range:

```
    python3 endf2g4ndl_ff.py neutron/am241n.out Fission/FF/95_241_Americium --compress \
        --energies 2.5e-2,1e5,5e5,1e6,2e6,5e6,1e7,1.4e7,2e7
```

For comparison, the G4NDL files converted from ENDF/B-VII.1 carry just two
groups, thermal and 500 keV.

:warning: UKFY places its thermal point at **2.5×10⁻² eV**, not the 2.53×10⁻² eV
used by ENDF/B. The converter passes the value through unchanged rather than
relabelling it.

### Verification

Three checks were run on the ²⁴¹Am output:

1. A standalone C++ replica of `GetDataStream()` + `ReadInData()` — 37 energy
   groups, 1545 distinct products, 57 165 items, clean EOF.
2. An independent re-parse of the generated file compared triplet-by-triplet
   against the ENDF source: 0 mismatches across all 57 165 values, with all 301
   metastable entries preserved.
3. ∑Y = 2.000000 ± 5×10⁻⁶ in every group.

The reference file `92_235_uranium` is byte-identical to the decompressed
`Fission/FF/92_235_Uranium.z` from G4NDL 4.7.1, which is what pinned down the
expected layout.

### Relation to PREPRO

None, directly. PREPRO (`../source`) linearises and reconstructs cross sections
in MF=2/MF=3, which is the step described in [INDC(NDS)-0758][INDC] §2 for
producing a full G4NDL library. Fission product yields in MF=8 need no such
pre-processing — they are already tabulated per incident energy — so this
converter reads the ENDF-6 tape directly. The other `Fission/` directories do
need it, for their cross sections; see below.

## The other Fission directories

`endf2g4ndl_fission.py` writes the rest of `Fission/` from an ENDF-6 neutron
evaluation. As with `FF/`, PREPRO only prepares the ENDF-6 data; the
translation into the token streams Geant4 reads is this script.

| Directory | ENDF-6 source | Read by |
| --- | --- | --- |
| `CrossSection/` | MF=3 MT=18 | `G4ParticleHPIsoData::Init` |
| `FS/` | MF=1 MT=452/455/456/458, MF=4/5 MT=18, MF=5 MT=455, MF=12/14/15 MT=18 | `G4ParticleHPFSFissionFS::Init` |
| `FC/` `SC/` `TC/` `LC/` | MF=3/4/5 for MT=19/20/21/38 (first, second, third, last chance) | `G4ParticleHPFissionBaseFS::Init` |

### Usage

```
    ./prepro_fission.sh jeff33/n_9228_92-U-235.dat /tmp/u235.endf
    python3 endf2g4ndl_fission.py /tmp/u235.endf Fission --compress
```

Several tapes can be given at once; each file lands at
`Fission/<dir>/<Z>_<A>[m<M>]_<Element>`, with `M` taken from LISO. A tape that
has no fission data for a directory simply produces no file there.
`--only FS,FC` restricts the directories; `--digits N` changes the number of
significant digits written (default 7, as the IAEA files — see below).

To regenerate the whole JEFF-3.3 set (about 5 minutes on 16 cores):

```
    mkdir -p jeff33/prepro
    ls jeff33/n_*.dat | xargs -P $(nproc) -I{} sh -c './prepro_fission.sh {} jeff33/prepro/$(basename {})'
    python3 endf2g4ndl_fission.py jeff33/prepro/n_*.dat Fission --compress
```

### What PREPRO has to do first

`prepro_fission.sh` runs, in order:

| Step | Why |
| --- | --- |
| LINEAR 0.1 % | Geant4 reads cross sections as bare (E, σ) pairs, i.e. lin-lin only |
| RECENT 0.1 % | adds the resolved/unresolved resonance contribution (MF=2) to MF=3 |
| LINEAR 0.1 % | thins the reconstructed MF=3 |
| FIXUP | tabulates polynomial ν̄ (LNU=1 → 3 points) and inserts a 0.0253 eV point in every table |
| DICTIN | rebuilds the MF=1/MT=451 index |
| SIXPAK | only when fission is given in MF=6 (JEFF-3.3: ²³²Th, ²³¹Pa, ²³³Pa): writes MF=4/5/12/14/15, appended to the tape |

INDC(NDS)-0758 gives 1 % for the second LINEAR pass. The JEFF-3.3 files the IAEA
distributes are only reproduced with 0.1 % — at 1 % the ²⁴⁰Cm resonance
integral in 1–10 eV comes out 16 % low. Tolerances can be overridden with
`LIN1_ERR`, `RECENT_ERR`, `LIN2_ERR`.

The converter refuses a tape whose fission cross sections are not lin-lin, and
warns when MT=451 still has LRP=1 (resonances not yet reconstructed).

### File layout

Everything is again read with `operator>>`. Each ENDF TAB1 becomes
`NP NR (NBT INT)… (x y)…` — note NP *before* NR — which is what
`G4ParticleHPVector::Init(istream)` reads. Distributions keep their ENDF
interpolation laws; only cross sections have to be lin-lin, since they are
read as pairs with no interpolation table.

- **CrossSection**: `0 0 NP` then the pairs.
- **FC/SC/TC/LC**: `QM LR NP` (QM in eV, as an integer), the pairs, and — only if
  the channel has MF=5 — `1 4` + angular block (isotropic if MF=4 is absent),
  `3 5` + energy block. Without final-state data Geant4 stops after the cross
  section.
- **FS**: a sequence of `(infoType, dataType)`-tagged blocks, in tape order:

  | Tag | ENDF | Geant4 |
  | --- | --- | --- |
  | `2 1` | MF=1 MT=452 | `G4ParticleHPParticleYield::InitMean` |
  | `3 1` | MF=1 MT=455 | `…::InitDelayed` (decay constants, then ν̄_d) |
  | `4 1` | MF=1 MT=456 | `…::InitPrompt` |
  | `5 1` | MF=1 MT=458 | `G4ParticleHPFissionERelease::Init` (constant terms only) |
  | `1 4` | MF=4 MT=18 | `G4ParticleHPAngular::Init` |
  | `1 5` | MF=5 MT=18 | `G4ParticleHPEnergyDistribution::Init` |
  | `3 5` | MF=5 MT=455 | same, one partial per delayed group |
  | `1 12` | MF=12 MT=18 | `G4ParticleHPPhotonDist::InitMean` (total yield skipped) |
  | `1 14` | MF=14 MT=18 | `…::InitAngular` |
  | `1 15` | MF=15 MT=18 | `…::InitEnergies` — written only if some MF=12 photon has LF=1, because the reader returns without consuming anything otherwise |

Supported representations: MF=4 LTT 0/1/2; MF=5 LF 1/5/7/9/11/12; MF=12 LO=1;
MF=14 LI=1 or LTT 1/2; ν̄ LNU 1/2. Anything else (MF=4 LTT=3, MF=1/455 LDG=1,
MF=6 without SIXPAK) stops with an error naming the section.

### Verification against the IAEA JEFF-3.3 files

All 72 JEFF-3.3 actinide tapes (Z ≥ 88) were converted and compared with
`/opt/Geant4/Libraries/JEFF-3.3/Fission` using `g4ndl_compare.py --tree`:

- **FS/** — 67 of 70 files token-identical (ints exact, reals to 7 digits). The
  other three are the SIXPAK cases (²³²Th, ²³¹Pa, ²³³Pa): same structure and
  token count, reals differing by ≤ 5×10⁻⁵ — the PREPRO-2017 SIXPAK used by the
  IAEA printed 5 significant digits.
- **FC/ SC/ TC/ LC/** — every final-state block and every header token-identical.
- **Cross sections** (CrossSection/ and the heads of FC–LC, 217 tables) — 157
  with exactly the IAEA point count, 43 fully token-identical; the rest differ
  in the resonance-region grid (PREPRO-2023 against the IAEA's PREPRO-2017
  RECENT). Integrals over every energy decade agree to ≤ 4×10⁻⁴.
- Every file the IAEA ships is produced, and no extra ones, except
  `CrossSection/0_0_Zero` — a Z=0 placeholder, not converted from any tape.

A pointwise comparison of ²³⁸U flags ~37 % of points in 1–20 keV. That comes
from the format itself: at 7 significant digits the energies inside narrow
resonances collapse — 27 044 of the 94 867 IAEA points in that range repeat the
previous energy. `--digits 9` avoids this in new files.

### Coding aids
Written with the help of Claude Code

### References

- W. Ford, *G4NDL4.2 Final State Guide*, McMaster University (2016) — `G4NDLFinalStateDecryptionCERNv1.pdf`
- E. Mendoza and D. Cano-Ott, *Update of the Evaluated Neutron Cross Section
  Libraries for the Geant4 Code*, [INDC(NDS)-0758][INDC] (2018) — `indc-nds-0758.pdf`
- [G4NDL format description][G4NDL-format], Geant4 Physics Reference Manual
- R.W. Mills, *UKFY-4.1 fission yield evaluation*, NexiaSol (2007)

[UKFY]: https://www.oecd-nea.org/dbdata/jeff/decay/UKFY4.1/UKFY4.1files/
[INDC]: https://doi.org/10.61092/iaea.5knd-4xdd
[JEFF33-IAEA]: https://www-nds.iaea.org/public/download-endf/JEFF-3.3/n/
[G4NDL-format]: https://geant4.web.cern.ch/documentation/pipelines/master/prm_html/PhysicsReferenceManual/hadronic/LowEnergyNeutron/G4NDL_format.html
