# Install and run

This GitHub repo is the **study overlay** (scripts, converter, Tree driver, SPC centrality patches). It does **not** contain Trajectum, SMASH, the full JCORRAN library, or O2Physics. Clone those next to (or as) `Trajectum/`, `deps/smash/`, `JCORRAN/`, and `OxygenSPC/`, then copy the overlay files on top.

Physics choices (posterior `08`, OO only, centrality) are in [README.md](README.md). This file is how to build and run.

## Prerequisites

| Tool | Used for | Notes |
| --- | --- | --- |
| CERN ROOT | converter, JCORRAN, SPC | Source **this** install: `source ~/softwares/root_install/bin/thisroot.sh` |
| C++ compiler | all builds | `clang++` on macOS, `g++` on Linux |
| Trajectum | `collide` | [Trajectum-Releases](https://codeberg.org/Trajectum/Trajectum-Releases) |
| SMASH | afterburner | built with Trajectum; binary on `PATH` |
| JCORRAN | flow | clone branch `DJ`, then overlay this repo’s Tree files |
| OxygenSPC + O2Physics | SPC | clone `main`, overlay this repo’s macros; enter `O2Physics/latest-master-o2` |

Trajectum also needs its own stack (GSL, HDF5, LHAPDF, FastJet, Pythia, Gemini, …). Follow Trajectum’s installer, not this file. On this machine those live under `deps/` and `$HOME/lhapdf`, `$HOME/softwares/fastjet-install`.

Edit paths in `setup_env.sh` (the `STUDY=` line) to your clone. `run/n100000/run_chain.sh` has the same hard-coded `STUDY` and ALICE `alienv` paths.

## 1. Clone this study and the analysis codes

```bash
git clone git@github.com:dongjokim/Trajectum-study.git
cd Trajectum-study
STUDY="$PWD"

git clone -b DJ https://github.com/MaximVirta/JCORRAN.git
# Overlay the Tree driver from this repo (already in JCORRAN/Example_JCorran/ if you cloned into STUDY).
# If JCORRAN was cloned elsewhere, copy:
#   JCORRAN/Example_JCorran/main_tree.cpp
#   JCORRAN/Example_JCorran/Makefile

git clone https://github.com/nmallick19/OxygenSPC.git
# Overlay from this repo:
#   OxygenSPC/analysis_pythia/JCorrSPCRun3.C
#   OxygenSPC/analysis_pythia/Makefile
#   OxygenSPC/centralityCalibration.C
```

Install Trajectum into `Trajectum/` and SMASH so `collide` and `smash` are on `PATH` (see Trajectum’s [readme](https://codeberg.org/Trajectum/Trajectum-Releases)). OO collide files: `parfiles/2509.04299/collisionOONLEFT5360_08.par` (or PGCM). Copies with `numevents` already set live in `run/`.

## 2. Environment

```bash
source ~/softwares/root_install/bin/thisroot.sh
# then either:
source "$STUDY/setup_env.sh"
```

Check: `root-config --version`, `which collide`, `which smash`, `jcorran-config-path`.

## 3. Build

### JCORRAN library

```bash
export PATH="$STUDY/JCORRAN/bin:$PATH"
createJCORRANlib
```

### Tree converter

```bash
source ~/softwares/root_install/bin/thisroot.sh
export PATH="$STUDY/JCORRAN/bin:$PATH"
g++ -O2 -o "$STUDY/tools/oscar_to_jtree" "$STUDY/tools/oscar_to_jtree.cpp" \
  $(root-config --cflags --libs) \
  $(jcorran-config-inc) $(jcorran-config-lib)
```

At runtime the JCORRAN dictionary (`.so` + `jcorranDict_rdict.pcm`) must be found. `run_chain.sh` does `ln -sfn "$STUDY/JCORRAN/include/jcorranDict_rdict.pcm" jcorranDict_rdict.pcm` in the working directory.

### JCORRAN Tree analysis

```bash
cd "$STUDY/JCORRAN/Example_JCorran"
make main_tree
```

### OxygenSPC (`JCorrSPCRun3`)

Needs O2Physics. Do **not** mix Trajectum/Pythia ROOT with the O2 `DYLD_LIBRARY_PATH`.

```bash
alias o2="alienv enter O2Physics/latest-master-o2"
o2
cd "$STUDY/OxygenSPC/analysis_pythia"
make
```

On macOS the Makefile takes the SDK from `xcrun --show-sdk-path` and O2 prefixes from `O2_ROOT` / `alienv` (`?=`). Override if your aliBuild tag differs.

Do not build or run `analysis_WS` for this Trajectum sample.

## 4. Run (1000-event test)

Match SMASH `Nevents` to Trajectum `numevents` (`1000` in `run/smash.yaml` and `run/collisionOONLEFT5360_08_test.par`).

```bash
source "$STUDY/setup_env.sh"
cd "$STUDY/run"          # or any empty work directory
cp "$STUDY/run/collisionOONLEFT5360_08_test.par" .
cp "$STUDY/run/smash.yaml" .

collide collisionOONLEFT5360_08_test.par
ln -f OO.main filename0
smash -i smash.yaml
ln -sfn data/0/particle_lists.oscar OO.afterburned

ln -sfn "$STUDY/JCORRAN/include/jcorranDict_rdict.pcm" jcorranDict_rdict.pcm
"$STUDY/tools/oscar_to_jtree" OO.afterburned OO.jcorran.root
# writes OO.jcorran.root.cent.csv (self-calibrated Nch |eta|<=0.5 percentiles)

mkdir -p jtrees
ln -sfn ../OO.jcorran.root jtrees/OO.jcorran.root

cd "$STUDY/JCORRAN/Example_JCorran"
ln -sfn "$STUDY/JCORRAN/include/jcorranDict_rdict.pcm" jcorranDict_rdict.pcm
./main_tree -o /path/to/AnalysisResults.root \
  -b OO_cent --pTmin 0.2 --pTmax 3.0 \
  /path/to/OO.jcorran.root

# SPC: fresh O2 env (see run_chain.sh). First argument is a directory of ROOT trees.
cd "$STUDY/OxygenSPC/analysis_pythia"
./JCorrSPCRun3 /path/to/jtrees spc-results-OO-NLEFT-Param08-3SPC 0 6
./JCorrSPCRun3 /path/to/jtrees spc-results-OO-NLEFT-Param08-4SPC 1 8
```

`0 6` = 3-harmonic SPC (min multiplicity 6). `1 8` = 4-harmonic (min multiplicity 8).

## 5. Production (100k)

Edit `STUDY` in `run/n100000/run_chain.sh`, then:

```bash
mkdir -p "$STUDY/run/n100000"
cp "$STUDY/run/n100000/smash.yaml" "$STUDY/run/n100000/"   # Nevents: 100000
# collide file: run/collisionOONLEFT5360_08_n100000.par
zsh "$STUDY/run/n100000/run_chain.sh"
```

That script: collide → SMASH → `oscar_to_jtree` (no CSV) → `main_tree` → SPC 3+4. Centrality is measured on **this** sample; do not reuse a 1k `.cent.csv`.

Turn `AnalysisResults.root` into `vn` / SC / ρ graphs with [POSTPROCESS.md](POSTPROCESS.md).

## Checks

- `collide` wrote `OO.main` / `OO.extra`
- SMASH wrote `data/0/particle_lists.oscar` (linked as `OO.afterburned`)
- `OO.jcorran.root` has `jTree` with `JTrackList` / `JEventHeaderList`
- `main_tree` wrote `AnalysisResults.root`
- `JCorrSPCRun3` wrote `spc-results-OO-NLEFT-Param08-3SPC.root` (and 4SPC)

If `oscar_to_jtree` or `main_tree` fails on `AliJBaseTrack` / dictionary, the `jcorranDict_rdict.pcm` symlink or `DYLD_LIBRARY_PATH` / `LD_LIBRARY_PATH` to `JCORRAN/lib` is missing.
