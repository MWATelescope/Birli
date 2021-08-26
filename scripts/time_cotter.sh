#!/bin/bash
set -ex
shopt -s globstar
shopt -s extglob

if [[ -z "$BIRLI_TEST_DIR" ]]; then
    echo "env var \$BIRLI_TEST_DIR not set."
    exit 1
fi
# export OBS_ID="1247842824" # Flags
export OBS_ID="1254670392" # Smol
# export OBS_ID="1196175296" # Big 

# export OUT_PATH="Flagfile.cotter.%%.mwaf"

# NO CORRECTIONS:
# export OUT_PATH="${OBSID}.cotter.none.uvfits"
# export GEOM_FLAG="-nogeom"
# export CABLE_FLAG="-nocablelength"

# BOTH CORRECTIONS:
export OUT_PATH="${OBSID}.cotter.both.uvfits"
export GEOM_FLAG=""
export CABLE_FLAG=""

cd "${BIRLI_TEST_DIR}/${OBS_ID}_vis/"
for i in $(seq 10); do
    [ -f "${OUT_PATH}" ] && rm "${OUT_PATH}"
    time cotter \
        -m "${OBS_ID}.metafits" \
        -o "${OUT_PATH}" \
        -nostats \
        ${GEOM_FLAG} \
        -noantennapruning \
        -noflagautos \
        -noflagdcchannels \
        ${CABLE_FLAG} \
        -edgewidth 0 \
        -initflag 0 \
        -endflag 0 \
        -flag-strategy "/usr/local/share/aoflagger/strategies/mwa-default.lua" \
        ${OBS_ID}*gpubox*_00.fits |
        tee "output-$i.log"
done
