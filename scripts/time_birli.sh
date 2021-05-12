#!/bin/bash
set -ex
shopt -s globstar
shopt -s extglob

if [[ -z "$BIRLI_TEST_DIR" ]]; then
    echo "env var \$BIRLI_TEST_DIR not set."
    exit 1
fi

export RUST_LOG="trace"
export FLAGS_OUT_DIR="$BIRLI_TEST_DIR/1247842824_vis/flags_birli"
mkdir -p "$FLAGS_OUT_DIR"

for i in $(seq 10); do
    if [ -n "$(ls "$BIRLI_TEST_DIR/*")" ]; then
        rm "$BIRLI_TEST_DIR/*"
    fi
    time birli aoflagger \
        -m "$BIRLI_TEST_DIR/1247842824_vis/1247842824.metafits" \
        -f "$FLAGS_OUT_DIR/Flagfile_Birli_%%.mwaf" \
        $BIRLI_TEST_DIR/1247842824_vis/1247842824*gpubox*_00.fits \
        | tee "output-birli-mwax-$i.log"
    time birli aoflagger \
        -m "$BIRLI_TEST_DIR/1196175296_vis/1196175296.metafits" \
        -f "$FLAGS_OUT_DIR/Flagfile_Birli_%%.mwaf" \
        $BIRLI_TEST_DIR/1196175296_vis/1196175296*gpubox*_00.fits \
        | tee "output-birli-ord-$i.log"
done
