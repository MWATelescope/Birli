#!/bin/bash
set -ex
shopt -s globstar
shopt -s extglob

if [[ -z "$BIRLI_TEST_DIR" ]]; then
    echo "env var \$BIRLI_TEST_DIR not set."
    exit 1
fi

export FLAGS_OUT_DIR="$BIRLI_TEST_DIR/1247842824_vis/flags_cotter"
mkdir -p "$FLAGS_OUT_DIR"

for i in $(seq 10); do
    if [ -n "$(ls "$BIRLI_TEST_DIR/*")" ]; then
        rm "$BIRLI_TEST_DIR/*"
    fi
    time cotter \
        -m "$BIRLI_TEST_DIR/1247842824_vis/1247842824.metafits" \
        -o "$FLAGS_OUT_DIR/Flagfile_Cotter_%%.mwaf" \
        -nostats \
        -nogeom \
        -noantennapruning \
        -sbpassband tests/data/subband-passband-128ch-unitary.txt \
        -noflagautos \
        -noflagdcchannels \
        -nocablelength \
        -edgewidth 0 \
        -initflag 0 \
        -endflag 0 \
        -flag-strategy "/usr/local/share/aoflagger/strategies/mwa-default.lua" \
        $BIRLI_TEST_DIR/1247842824_vis/1247842824*gpubox*_00.fits |
        tee "output-$i.log"
done
