# Postprocess JCORRAN `AnalysisResults.root`

After `main_tree` (or `main_hdf5`) writes `AnalysisResults.root`, turn the raw TProfiles into graphs. Macros live in [`tools/flowPP/`](tools/flowPP/). Official SPC for this study is still OxygenSPC (`JCorrSPCRun3`); the optional JCORRAN SPC list in the same file is not the result to quote.

Build and generation: [INSTALL.md](INSTALL.md). Physics: [README.md](README.md).

## What is in `AnalysisResults.root`

`main_tree` writes under `hydro/`:

| Path | Source | Postprocess |
| --- | --- | --- |
| `hydro/jfluc/` | `AliJFFlucAnalysisTProfile` | `jfflucresults_TProfile.cpp` |
| `hydro/vnptCorr/OutputPtVnCorr` | `AliAnalysisPtVn` | `PP_RhoCorr.cxx` → `calculate` |
| `hydro/SPC/OutputListSPC` | optional `--addObs spc` | ignore for this study |

JFluc OO bins in the file are the native `OO_cent` edges: 0–1, 1–2, 2–5, 5–10, then 5% steps to 50% (`CentBin00` … `CentBin11`). The JFluc postprocessor rebins some of those for plots.

ROOT must be the CERN install:

```bash
source ~/softwares/root_install/bin/thisroot.sh
```

## 1. JFluc: `vn`, SC, ρ, χ

[`tools/flowPP/jfflucresults_TProfile.cpp`](tools/flowPP/jfflucresults_TProfile.cpp) is a ROOT macro (`#define OO true`). Run it from `tools/flowPP` so `observables.h` is found.

```bash
source ~/softwares/root_install/bin/thisroot.sh
STUDY="/Users/djkim/Documents/GitHub/Trajectum-study"   # edit
IN="$STUDY/run/n100000/AnalysisResults.root"
OUTDIR="$STUDY/run/n100000/PostProcessed/Jfluc"
mkdir -p "$OUTDIR"

cd "$STUDY/tools/flowPP"
root -l -b -q 'jfflucresults_TProfile.cpp("'"$IN"'","hydro/jfluc","'"$OUTDIR"'/FR_JflucPP_OO_NLEFT-Param08.root")'
```

Arguments:

1. input `AnalysisResults.root`
2. container path — **`hydro/jfluc`** for this Tree driver (not the ALICE default in the source)
3. output ROOT file (no `.root` extra suffix; the macro also writes `outfile.txt`)

Example already in the tree: `run/n100000/PostProcessed/Jfluc/FR_JflucPP_OO_NLEFT-Param08.root.txt`.

The output ROOT file holds `TGraphErrors`:

- `gr_v2`, `gr_v3`, … and `_QC`, `_2`, `_QC4`, `_nat`, `_wde` variants
- `gr_rho422`, `gr_chi422`, `gr_vnm422`, … (see `observables.h`)
- `gr_sc32_QC`, `gr_sc32N_QC`, …
- `gr_pTmean_*`, `gr_mult_*`, `gr_ecc2`, `gr_ecc3`

If histograms are missing, the container path is wrong (`hydro/jfluc` vs `jfluc`) or the file was not produced with `-b OO_cent`.

## 2. `vn`–`pT`: ρ, cumulants of pT

[`tools/flowPP/PP_RhoCorr.cxx`](tools/flowPP/PP_RhoCorr.cxx) reads `hydro/vnptCorr/OutputPtVnCorr`.

```bash
source ~/softwares/root_install/bin/thisroot.sh
cd "$STUDY/tools/flowPP"
bash compile_PP_RhoCorr.sh          # builds ./calculate

mkdir -p "$STUDY/run/n100000/PostProcessed/vnpt"
cd "$STUDY/run/n100000"             # so the output name stays short
"$STUDY/tools/flowPP/calculate" AnalysisResults.root PostProcessed/vnpt
```

That writes `PostProcessed/vnpt/output_AnalysisResults.root` (the code prefixes `output_` to the **input file name**, so pass a basename, not an absolute path).

Graphs inside: `c2sub`, `k2sub`, `stdskewsub`, `intskewsub`, `exkursub`, `cv22pt`, `ck`, `varv22`, `rhov22pt`, `cv24pt`.

## Checks

- `hydro/jfluc/hvna/hvnaK01CentBin00` exists → JFluc filled
- `hydro/vnptCorr/OutputPtVnCorr` is a `TList` → pT–vn filled
- JFluc output ROOT has `gr_v2` / `gr_v2_QC`
- pT–vn output has `rhov22pt`

Do not feed OxygenSPC `spc-results-*.root` into these macros. Those go to `OxygenSPC/analysis_WS/CalculateSPC.C` if you postprocess SPC separately.
