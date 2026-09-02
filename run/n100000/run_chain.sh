#!/bin/zsh
# Full 100k chain: collide → SMASH → self-calibrated tree → JCORRAN (0.2–3 GeV) → SPC 3+4.
# Centrality is measured on this sample (oscar_to_jtree with no CSV).
set -euo pipefail
STUDY="/Users/djkim/Documents/GitHub/Trajectum-study"
OUT="$STUDY/run/n100000"
PAR="$STUDY/run/collisionOONLEFT5360_08_n100000.par"

source "$STUDY/setup_env.sh"
unset PYTHIA8DATA
cd "$OUT"

echo "==== $(date) collide 100000 events ===="
collide "$PAR"
echo "==== $(date) collide done ===="
ls -lh OO.main OO.extra OO.geometry OO.bypass

ln -f OO.main filename0
echo "==== $(date) SMASH ===="
smash -i smash.yaml
ln -sfn data/0/particle_lists.oscar OO.afterburned
echo "==== $(date) SMASH done ===="
ls -lh OO.afterburned data/0/particle_lists.oscar

echo "==== $(date) convert (self-calibrate centrality) ===="
ln -sfn "$STUDY/JCORRAN/include/jcorranDict_rdict.pcm" jcorranDict_rdict.pcm
"$STUDY/tools/oscar_to_jtree" OO.afterburned OO.jcorran.root
mkdir -p jtrees
ln -sfn ../OO.jcorran.root jtrees/OO.jcorran.root
echo "==== $(date) convert done ===="
cat OO.jcorran.root.cent.csv

echo "==== $(date) JCORRAN pT 0.2-3.0 ===="
cd "$STUDY/JCORRAN/Example_JCorran"
ln -sfn "$STUDY/JCORRAN/include/jcorranDict_rdict.pcm" jcorranDict_rdict.pcm
./main_tree -o "$OUT/AnalysisResults.root" \
  -b OO_cent --addObs spc --pTmin 0.2 --pTmax 3.0 \
  "$OUT/OO.jcorran.root"
echo "==== $(date) JCORRAN done ===="

echo "==== $(date) OxygenSPC 3SPC + 4SPC ===="
TREEDIR="$OUT/jtrees"
run_spc() {
  local outbase="$1" which="$2" multmin="$3"
  env -i HOME="$HOME" USER="$USER" LOGNAME="${LOGNAME:-$USER}" TMPDIR="${TMPDIR:-/tmp}" \
    PATH="/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin:/sbin" \
    ALIBUILD_WORK_DIR=/Users/djkim/alice/sw \
    STUDY="$STUDY" TREEDIR="$TREEDIR" OUTBASE="$outbase" WHICH="$which" MULTMIN="$multmin" \
    /bin/bash --noprofile --norc -c '
      cd /Users/djkim/alice
      export STUDY TREEDIR OUTBASE WHICH MULTMIN
      alienv setenv O2Physics/latest-master-o2 -c /bin/bash --noprofile --norc -c "
        export DYLD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:/opt/homebrew/opt/libidn2/lib:/opt/homebrew/lib:\${DYLD_LIBRARY_PATH}
        cd \$STUDY/OxygenSPC/analysis_pythia
        ./JCorrSPCRun3 \$TREEDIR \$OUTBASE \$WHICH \$MULTMIN
      "
    '
}
run_spc "$OUT/spc-results-OO-NLEFT-Param08-3SPC" 0 6
run_spc "$OUT/spc-results-OO-NLEFT-Param08-4SPC" 1 8
echo "==== $(date) all done ===="
ls -lh "$OUT"/OO.jcorran.root "$OUT"/AnalysisResults.root \
  "$OUT"/spc-results-OO-NLEFT-Param08-*.root "$OUT"/OO.jcorran.root.cent.csv
