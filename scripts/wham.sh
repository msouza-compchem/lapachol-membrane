#!/usr/bin/env bash
# WHAM profile from the umbrella windows.
# NOT RUN YET: it will be tested when stage 1 (5 ns per window) finishes.
# -b 500 drops the first 500 ps of each window; -zprof0 3.5 puts the zero of the profile in bulk water.
set -euo pipefail
ls umbrella/umb_*.tpr        > tpr-files.dat
ls umbrella/umb_*_pullf.xvg  > pullf-files.dat
gmx wham -it tpr-files.dat -if pullf-files.dat \
  -o profile.xvg -hist histo.xvg -unit kJ \
  -temp 303.15 -b 500 -zprof0 3.5 \
  -nBootstrap 100 -bs-method b-hist \
  -bsres bsResult.xvg -bsprof bsProfs.xvg
