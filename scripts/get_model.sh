#!/bin/bash
# Download the default Hcor bicycle-gene GLM model from Zenodo.
#
# Usage:
#   scripts/get_model.sh                  # → $HOME/.bicycle-classifier/models/
#   scripts/get_model.sh /custom/dir      # → /custom/dir/
#
# Once the lab has minted a Zenodo deposit for the v5.5.6 model file, set
# MODEL_URL below to the direct download link and (optionally) pin
# EXPECTED_SHA256 to the checksum reported on the Zenodo record page.

set -euo pipefail

MODEL_NAME="Hcor.glm.full_v5.5.6"
MODEL_VERSION="v5.5.6"

# TODO: Fill in once Zenodo deposit is created. Example shape:
#   https://zenodo.org/records/XXXXXXX/files/Hcor.glm.full_v5.5.6
MODEL_URL=""
EXPECTED_SHA256=""   # leave empty to skip checksum verification

DEST_DIR="${1:-${HOME}/.bicycle-classifier/models}"
DEST_PATH="${DEST_DIR}/${MODEL_NAME}"

if [[ -z "${MODEL_URL}" ]]; then
    cat >&2 <<EOF
ERROR: Model download URL is not yet configured.

The default Hcor model (${MODEL_NAME}, ${MODEL_VERSION}) has not yet been
uploaded to Zenodo. Once it is, edit scripts/get_model.sh and set the
MODEL_URL variable to the direct file URL (plus EXPECTED_SHA256 if you
want checksum verification).

In the meantime you can:
  - Ask the lab for the model file directly, then point the classifier at it:
      bicycle_classifier -m /path/to/${MODEL_NAME} -g your_genes.gff3
  - Or set:
      export BICYCLE_MODEL=/path/to/${MODEL_NAME}
      bicycle_classifier -g your_genes.gff3
EOF
    exit 1
fi

mkdir -p "${DEST_DIR}"

if [[ -f "${DEST_PATH}" ]]; then
    echo "Model already present at: ${DEST_PATH}"
    echo "Delete and re-run if you want to re-download."
    exit 0
fi

echo "Downloading bicycle-gene-classifier model ${MODEL_VERSION}..."
echo "  source: ${MODEL_URL}"
echo "  dest:   ${DEST_PATH}"

curl -L --fail -o "${DEST_PATH}.partial" "${MODEL_URL}"
mv "${DEST_PATH}.partial" "${DEST_PATH}"

if [[ -n "${EXPECTED_SHA256}" ]]; then
    echo "Verifying SHA256..."
    ACTUAL_SHA256=$(sha256sum "${DEST_PATH}" | awk '{print $1}')
    if [[ "${ACTUAL_SHA256}" != "${EXPECTED_SHA256}" ]]; then
        echo "ERROR: SHA256 mismatch!" >&2
        echo "  expected: ${EXPECTED_SHA256}" >&2
        echo "  actual:   ${ACTUAL_SHA256}" >&2
        rm -f "${DEST_PATH}"
        exit 1
    fi
    echo "SHA256 OK."
fi

echo ""
echo "Done. To use this model, either:"
echo "  - Let bicycle_classifier find it automatically (it checks \$HOME/.bicycle-classifier/models/ by default)"
echo "  - Or set: export BICYCLE_MODEL=${DEST_PATH}"
echo "  - Or pass: bicycle_classifier -m ${DEST_PATH} ..."
