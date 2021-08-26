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

# NO CORRECTIONS:
# export OUT_PATH="${OBSID}.birli.none.uvfits"
# export GEOM_FLAG="--no-geometric-delay"
# export CABLE_FLAG="--no-cable-delay"

# BOTH CORRECTIONS:
export OUT_PATH="${OBSID}.birli.both.uvfits"
export GEOM_FLAG=""
export CABLE_FLAG=""


# export RUST_LOG="trace"

cd "${BIRLI_TEST_DIR}/${OBS_ID}_vis/"

for i in $(seq 10); do
    time birli aoflagger \
        -m "${OBS_ID}.metafits" \
        -u "${OUT_PATH}" \
        ${GEOM_FLAG} ${CABLE_FLAG} \
        ${OBS_ID}*gpubox*_00.fits \
        | tee "output-birli-uvfits-$i.log"
done
