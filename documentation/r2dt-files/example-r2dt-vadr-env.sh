# Example site configuration for 'v-annotate.pl --draw_r2dt'.
#
# If the file
#
#   $R2DT_DIR/r2dt-vadr-env.sh
#
# exists, v-annotate.pl --draw_r2dt sources it with POSIX '.' immediately
# before invoking r2dt.py. It exists so that a site can supply whatever PATH
# and environment its R2DT installation needs without VADR hardcoding any
# install-specific paths.
#
# If that file does NOT exist, VADR proceeds on the assumption that
# 'python $R2DT_DIR/r2dt.py' already works in the environment v-annotate.pl
# was launched from.
#
# To enable it for an installation:
#
#   cp <vadr>/documentation/r2dt-files/example-r2dt-vadr-env.sh $R2DT_DIR/r2dt-vadr-env.sh
#   # then edit the paths below to match the local installation
#
# It is sourced with '.' from /bin/sh, so keep it POSIX compatible: use
# 'export VAR=...', not bash-only constructs, and do not rely on 'source'.
#
# The only environment VADR sets on its own is thread-count pinning
# (OPENBLAS_NUM_THREADS, OMP_NUM_THREADS and MKL_NUM_THREADS are set to 1 on
# the r2dt.py command line). Everything else is up to this file.
#
# See documentation/r2dt-drawing.md for more.
#
# ---------------------------------------------------------------------------
# The values below are PLACEHOLDERS. Replace every <...> with a real path, or
# delete the lines you do not need.
# ---------------------------------------------------------------------------

# Put the python that has R2DT's dependencies installed first on PATH, so that
# 'python' resolves to it. If R2DT's requirements were installed into a
# virtualenv, prepend that virtualenv's bin directory rather than activating
# it, which keeps this file POSIX-safe: a virtualenv python works by virtue of
# where it sits on PATH.
#
# Also put the external programs R2DT calls on PATH:
#   - Infernal          (cmalign, cmbuild, cmemit, esl-* miniapps)
#   - Bio-Easel scripts
#   - the jiffy Infernal/HMMER scripts
#   - Traveler          (the program that renders the diagrams)

export PATH=<path to python environment>/bin:<path to infernal>/bin:<path to Bio-Easel>/scripts:<path to jiffy-infernal-hmmer-scripts>:<path to traveler>/bin:$PATH

# Add any further environment your R2DT installation requires here.
