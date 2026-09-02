#!/bin/zsh
# Environment for Trajectum + JCORRAN. Source this from zsh or bash.
STUDY="/Users/djkim/Documents/GitHub/Trajectum-study"

export PATH="/opt/homebrew/opt/gsl/bin:/opt/homebrew/opt/hdf5/bin:/opt/homebrew/bin:$PATH"
source "$HOME/softwares/root_install/bin/thisroot.sh"

export PATH="$STUDY/JCORRAN/bin:$STUDY/Trajectum/src:$STUDY/deps/smash/build:$HOME/softwares/root_install/bin:$PATH"
export PKG_CONFIG_PATH="/opt/homebrew/lib/pkgconfig:/opt/homebrew/opt/eigen/share/pkgconfig:/opt/homebrew/opt/hdf5/lib/pkgconfig:${PKG_CONFIG_PATH:-}"

export LHAPDF_DATA_PATH="$HOME/lhapdf/share/LHAPDF"
export GINPUT="$STUDY/deps/gemini_2026.4.21/install/share/gemini/"

PYTHIA_LIB="$STUDY/deps/pythia8317/lib"
FASTJET_LIB="$HOME/softwares/fastjet-install/lib"
LHAPDF_LIB="$HOME/lhapdf/lib"
GEMINI_LIB="$STUDY/deps/gemini_2026.4.21/install/lib"
HDF5_LIB="/opt/homebrew/opt/hdf5/lib"
GSL_LIB="/opt/homebrew/opt/gsl/lib"

export DYLD_LIBRARY_PATH="$PYTHIA_LIB:$FASTJET_LIB:$LHAPDF_LIB:$GEMINI_LIB:$HDF5_LIB:$GSL_LIB:${DYLD_LIBRARY_PATH:-}"

echo "Trajectum-study env: ROOT=$(root-config --version)  JCORRAN=$(jcorran-config-path)"
