# Trajectum OO flow study

Use Trajectum to generate oxygen–oxygen (OO) events, convert them to a **JCORRAN Tree**, then run **two independent analyses** on that same tree.

**Install and run:** [INSTALL.md](INSTALL.md). **JCORRAN postprocess** (`AnalysisResults.root` → graphs): [POSTPROCESS.md](POSTPROCESS.md). This README is the physics workflow and what differs from upstream JCORRAN / OxygenSPC.

## Goal

Produce anisotropic-flow results for OO collisions from Trajectum:

1. **JCORRAN** — multi-particle flow (`vn`, fluctuations, `vn`–`pT` correlations)
2. **SPC** — symmetric/asymmetric cumulants (OxygenSPC)

Do not mix posterior draws. For this study, use **only posterior sample `08`** and **only OO** (not NeNe or PbPb).

## Codes

| Role | Upstream | What this study uses |
| --- | --- | --- |
| Event generator | [Trajectum-Releases](https://codeberg.org/Trajectum/Trajectum-Releases) | local `Trajectum/` plus `parfiles/2509.04299` |
| JCORRAN flow | [MaximVirta/JCORRAN](https://github.com/MaximVirta/JCORRAN) branch [`DJ`](https://github.com/MaximVirta/JCORRAN/tree/DJ) | local `JCORRAN/` on `DJ`, plus a Tree driver not on GitHub |
| SPC | [nmallick19/OxygenSPC](https://github.com/nmallick19/OxygenSPC) branch [`main`](https://github.com/nmallick19/OxygenSPC/tree/main) | local `OxygenSPC/analysis_pythia` (the folder name is leftover from Pythia samples) |

Paper for the Trajectum settings: [arXiv:2509.04299](https://arxiv.org/abs/2509.04299).

### Local changes vs those GitHub copies

Neither analysis fork is a full rewrite. The libraries stay as upstream; the study-specific edits are the Tree conversion, the Tree-input driver, and how SPC picks a centrality class.

**JCORRAN** ([MaximVirta/JCORRAN](https://github.com/MaximVirta/JCORRAN), branch `DJ` at `141d3df`). Clone with `git clone -b DJ https://github.com/MaximVirta/JCORRAN.git`. GitHub `main` is ahead of `DJ` (SMASH PDG helpers, dependency cleanup) and still has no Tree driver.

| Change | File | Against `DJ` |
| --- | --- | --- |
| **New** Tree-input executable | `Example_JCorran/main_tree.cpp` | not in the GitHub repo; reads `jTree` / `JTrackList` / `JEventHeaderList` via `JTreeDataManager`, uses `hdr->GetCentrality()`, default `-b OO_cent`, track cuts \(0.2 < p_T < 3.0\) GeV/\(c\), \(\lvert\eta\rvert\le 0.8\) |
| Build that executable | `Example_JCorran/Makefile` | `all` also builds `main_tree` |

HDF5 `main_hdf5` is unchanged from `DJ`. The Oscar → Tree converter is not part of JCORRAN; it lives in this study as [`tools/oscar_to_jtree.cpp`](tools/oscar_to_jtree.cpp).

**OxygenSPC** ([nmallick19/OxygenSPC](https://github.com/nmallick19/OxygenSPC), branch `main` at `004dc35`). Clone with `git clone https://github.com/nmallick19/OxygenSPC.git`.

| Change | File | Against `main` |
| --- | --- | --- |
| Bin on stored centrality % | `analysis_pythia/JCorrSPCRun3.C` | upstream used `GetYVertexMC()` as \(N_{\mathrm{ch}}\) and `FindMultiplicityBin`; this study uses `FindCentralityBin(hdr->GetCentrality())` with edges `{0, 5, 10, 20, 30, 40, 50, 60, 70, 100}` so SPC and JCORRAN share the converter percentile. QA `hmult` still fills `GetXVertexMC()` (raw \(N_{\mathrm{ch}}\)) |
| Same binning (unused here) | `analysis_WS/JCorrSPCRun3.C` | same `GetCentrality()` change; do **not** run `analysis_WS` for this Trajectum sample |
| macOS build | `analysis_pythia/Makefile`, `analysis_WS/Makefile` | SDK from `xcrun --show-sdk-path`; O2 paths `?=` so `alienv` can override; default O2 tag updated |

Do not run `analysis_WS` for this study. Name SPC output after the Trajectum sample, not Pythia.

## Trajectum setting (OO, sample 08 only)

That folder contains 20 posterior draws (`00`–`19`) for several systems. **Run only these OO files:**

- [`collisionOONLEFT5360_08.par`](https://codeberg.org/Trajectum/Trajectum-Releases/src/branch/master/parfiles/2509.04299/collisionOONLEFT5360_08.par) — NLEFT nuclear structure
- [`collisionOOPGCM5360_08.par`](https://codeberg.org/Trajectum/Trajectum-Releases/src/branch/master/parfiles/2509.04299/collisionOOPGCM5360_08.par) — PGCM nuclear structure

Do **not** run `collisionNeNe*` or `collisionPbPb*` from this folder, and do **not** run other `_00`–`_19` indices.

Both files already set `numevents=1000`, `output=OO`, and `format=smash`. Keep `numevents=1000` for the test run.

Pick **one** nuclear-structure variant per test (NLEFT or PGCM). Do not merge them.

## Workflow

```
Trajectum collide (OO, 08)
        → SMASH afterburner
        → convert particle lists to a JCORRAN Tree (preferred shared format)
              ↳ JCORRAN flow  (jTree  or  HDF5)
              ↳ SPC           (OxygenSPC tree analysis → NLEFT-Param08)
```

JCORRAN accepts **either** a JCORRAN Tree or HDF5. Prefer a **JCORRAN Tree** as the shared sample so both analyses read the same file. Do not regenerate Trajectum events between the two analyses.

### 1. Generate Trajectum events

Install and run Trajectum as in its [readme](https://codeberg.org/Trajectum/Trajectum-Releases). From the Trajectum `src/` directory:

```bash
./collide ../parfiles/2509.04299/collisionOONLEFT5360_08.par
```

(or the PGCM `_08` file). Expected files (prefix `OO`):

- `OO.main` — freeze-out particles for SMASH
- `OO.bypass`, `OO.geometry`, `OO.extra`

### 2. Run the SMASH afterburner

Match `Nevents` in `parfiles/smash.yaml` to `numevents` in the collide file (`1000` for the test). Rename `OO.main` to the prefix SMASH expects (`filename0` with the default yaml), then:

```bash
smash -i ../parfiles/smash.yaml
```

Rename the Oscar output to `OO.afterburned` next to the other `OO.*` files.

### 3. Convert to a JCORRAN Tree (shared sample)

Convert SMASH/`OO.afterburned` particles into a **JCORRAN Tree** and keep that file as the canonical event sample.

Tree layout (`JTreeDataManager`):

| Item | Name |
| --- | --- |
| `TTree` | `jTree` |
| Tracks | branch `JTrackList` — `TClonesArray` of `AliJBaseTrack` (`px, py, pz, E`, ID, PID, charge) |
| Event header | branch `JEventHeaderList` — `TClonesArray` of `AliJBaseEventHeader` (event ID, centrality, vertex) |

`JTreeDataManager::ChainInputStream` reads a **text list of ROOT files**, one path per line.

HDF5 is the other valid JCORRAN input (`main_hdf5`): per-event datasets with `pT`, `eta`, `phi`, `ID`, `charge`, `sample`, plus attributes `dNch_deta`, `e2`, `e3`. Use it only if you are not writing a Tree. Do not maintain two divergent samples.

Keep charged hadrons, PID, and kinematics; do not apply analysis-level `pT`/`η` cuts at conversion.

### Centrality

Centrality is **not** taken from Trajectum impact parameter. It is charged multiplicity at mid-rapidity, converted to a percentile, and stored in the JCORRAN Tree header. Measure it on **each new sample**. Do not reuse the 1k table for a 10k or 100k run.

**Estimator.** Charged \(N_\mathrm{ch}\) with \(|\eta|\le 0.5\). High \(N_\mathrm{ch}\) is most central (0%). The converter writes that percentile into `AliJBaseEventHeader::SetCentrality()` and also stores the raw \(N_\mathrm{ch}\) in `SetXVertexMC` / `SetYVertexMC` so later macros can re-derive cuts.

**Write it (production).** [`tools/oscar_to_jtree.cpp`](tools/oscar_to_jtree.cpp) is the converter. With no CSV it calibrates percentiles `{0, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}` from this Oscar file and writes `OO.jcorran.root.cent.csv`:

```bash
tools/oscar_to_jtree OO.afterburned OO.jcorran.root
```

To apply a table already measured on the **same campaign** (same generator settings and event count class):

```bash
tools/oscar_to_jtree OO.afterburned OO.jcorran.root path/to/OO.jcorran.root.cent.csv
```

Optional fourth argument is the CSV row index (`paramID`, default `0`). The CSV is one row of \(N_\mathrm{ch}\) edges, highest multiplicity first.

Do **not** pass `JCORRAN/Example_JCorran/dependencies/OO5360_main_cents.csv` into `oscar_to_jtree`. That file is for HDF5 `main_hdf5 --system`, not for this Tree.

**Check / re-derive cuts on an existing tree.** [`OxygenSPC/centralityCalibration.C`](OxygenSPC/centralityCalibration.C) chains `jTree` files, fills `nChFT0M` from `GetYVertexMC()`, writes `centralityCalibrationHist.root`, and prints \(N_\mathrm{ch}\) cuts for bins `{0, 30, 40, 50, 60, 70, 80, 90, 95, 100}`. Point `infile_path` at the directory of ROOT trees (the default path is leftover from another sample), then:

```bash
source ~/softwares/root_install/bin/thisroot.sh
cd OxygenSPC
# build dictionaries once if src/AliJBase*.so are missing:
#   cd analysis_pythia && make
root -l centralityCalibration.C
```

This macro does **not** rewrite `GetCentrality()` in the tree. Use it to inspect the \(N_\mathrm{ch}\) distribution. Production centrality still comes from `oscar_to_jtree`.

**Who reads `GetCentrality()`.**

| Analysis | Macro | Binning of the stored percentile |
| --- | --- | --- |
| JCORRAN Tree | [`JCORRAN/Example_JCorran/main_tree.cpp`](JCORRAN/Example_JCorran/main_tree.cpp) | `-b OO_cent`: 0–1, 1–2, 2–5, 5–10, then 5% bins up to 50% (`CentBin_OO_central`) |
| SPC | [`OxygenSPC/analysis_pythia/JCorrSPCRun3.C`](OxygenSPC/analysis_pythia/JCorrSPCRun3.C) `FindCentralityBin` | 0–5, 5–10, 10–20, …, 60–70, 70–100 → `Centrality_0` … `Centrality_8` |
| JCORRAN HDF5 only | `Example_JCorran/main_hdf5` | maps `dNch_deta` with `--system dependencies/OO5360_main_cents.csv`; skip this if you use the Tree |

JCORRAN and OxygenSPC both read the header; they do not recompute \(N_\mathrm{ch}\) percentiles. After conversion, pass `-b OO_cent` to `main_tree` and do not point it at a cents CSV from a different run.

**Weights / oversampling.** These collide files use non-uniform `entropyacceptanceprobability` and `oversampling`. Trajectum’s own `analyze`/`collect` correct for that. JCORRAN and SPC do not. For a first test, either flatten those sections to a constant (weights `1`, oversampling `1`) or propagate the weights from `OO.extra`.

### 4. JCORRAN flow analysis

Clone [MaximVirta/JCORRAN](https://github.com/MaximVirta/JCORRAN) branch `DJ`, put `JCORRAN/bin` on `PATH`, and build the library (`createJCORRANlib`).

**Tree input (preferred):** chain the JCORRAN Tree from step 3 with `JTreeDataManager` (`jTree` / `JTrackList` / `JEventHeaderList`) and run the same flow classes as in `Example_JCorran` (`AliJFFlucAnalysisTProfile`, `AliAnalysisPtVn`, …). For OO use binning `OO_cent`. Track cuts: \(0.2 < p_T < 3.0\) GeV/\(c\), \(|\eta| \le 0.8\). Centrality is already in the tree from this sample; do not point `main_tree` at the JCORRAN `OO5360_main_cents.csv` from a different run.

**HDF5 input:** if the sample was written as HDF5 instead of a Tree:

```bash
cd Example_JCorran
make
./main_hdf5 \
  -o AnalysisResults.root \
  -b OO_cent \
  --system dependencies/OO5360_main_cents.csv \
  --param 8 \
  path/to/events.hdf
```

Optional: `--addObs spc` runs JCORRAN’s own SPC; the official SPC result for this study still comes from OxygenSPC.

### 5. SPC analysis

Run OxygenSPC’s generic tree analysis ([`analysis_pythia`](https://github.com/nmallick19/OxygenSPC/tree/main/analysis_pythia) — the directory name is leftover from Pythia samples; it is not a second Trajectum initial condition) on the **same JCORRAN Tree** as step 4. Do **not** run `analysis_WS`. Name the output after this Trajectum sample, not after Pythia:

```bash
cd OxygenSPC/analysis_pythia
./JCorrSPCRun3 path/to/jtrees spc-results-OO-NLEFT-Param08-3SPC 0 6
./JCorrSPCRun3 path/to/jtrees spc-results-OO-NLEFT-Param08-4SPC 1 8
```

`0 6` is three-harmonic SPC (min multiplicity 6); `1 8` is four-harmonic SPC (min multiplicity 8). Track cuts stay \(0.2 < p_T < 5.0\) GeV/\(c\), \(|\eta| \le 0.8\) (not the JCORRAN \(p_T < 3\) GeV/\(c\)). Do not regenerate Trajectum events for SPC.

## Test run

Use **1000 events** (`numevents=1000` is already in the `_08` collide files). Check:

- [ ] `collide` finishes and writes `OO.main` / `OO.extra`
- [ ] SMASH writes `OO.afterburned` for 1000 events
- [ ] Conversion produces a non-empty JCORRAN Tree (`jTree` with `JTrackList` / `JEventHeaderList`)
- [ ] JCORRAN flow writes `AnalysisResults.root` (from Tree or HDF5)
- [ ] OxygenSPC runs on the **same** Tree without reprocessing Trajectum

After the 1000-event test works, increase `numevents` (and SMASH `Nevents`) for production. Keep posterior index `08` and the OO-only files.
