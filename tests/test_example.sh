#!/bin/bash
# tests/test_example.sh — end-to-end smoke test.
#
# Runs bicycle_classifier on data/example.gff3 and confirms the three expected
# output files appear and are non-empty. Does NOT validate biological accuracy
# (the example GFF3 is synthetic).
#
# Requirements:
#   - R + dependencies installed (see install.R)
#   - A bicycle GLM model available via $BICYCLE_MODEL env var OR
#     ~/.bicycle-classifier/models/Hcor.glm.full_v5.5.6
#
# Usage:
#   tests/test_example.sh
#   BICYCLE_MODEL=/path/to/model tests/test_example.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

WORK_DIR="$(mktemp -d)"
trap "rm -rf ${WORK_DIR}" EXIT

echo "=== bicycle_classifier smoke test ==="
echo "Work dir: ${WORK_DIR}"
echo ""

# If no model resolvable, fail with an actionable message rather than going through R.
if [[ -z "${BICYCLE_MODEL:-}" ]] \
        && [[ ! -f "${HOME}/.bicycle-classifier/models/Hcor.glm.full_v5.5.6" ]] \
        && [[ ! -f "${REPO_ROOT}/models/Hcor.glm.full_v5.5.6" ]]; then
    cat >&2 <<EOF
SKIP: No model file available.
This test needs the trained GLM model. Either:
  - Set BICYCLE_MODEL=/path/to/Hcor.glm.full_v5.5.6
  - Or run: bin/bicycle_classifier --download-model
  - Or drop the model file at ${REPO_ROOT}/models/
EOF
    exit 77   # POSIX "test skipped" convention
fi

cd "${WORK_DIR}"
"${REPO_ROOT}/bin/bicycle_classifier" \
    -g "${REPO_ROOT}/data/example.gff3" \
    -o smoketest \
    -d "${WORK_DIR}/out"

echo ""
echo "=== Assertions ==="
EXPECTED=(
    "${WORK_DIR}/out/smoketest_classifier_all_transcripts_response.txt"
    "${WORK_DIR}/out/smoketest_classifier_bicycle_gene_names.txt"
    "${WORK_DIR}/out/smoketest_classifier_response_histogram.pdf"
)
fail=0
for f in "${EXPECTED[@]}"; do
    if [[ ! -f "${f}" ]]; then
        echo "  FAIL: missing ${f}"
        fail=1
    elif [[ ! -s "${f}" ]]; then
        # Histogram PDF must be non-empty; the bicycle-gene-names file may be
        # empty if no genes pass cutoff on synthetic data — only flag the
        # response table + PDF as required-non-empty.
        case "${f}" in
            *_classifier_response_histogram.pdf|*_classifier_all_transcripts_response.txt)
                echo "  FAIL: empty ${f}"; fail=1 ;;
            *)
                echo "  WARN: empty ${f} (acceptable for this file)" ;;
        esac
    else
        echo "  OK:   ${f}  ($(wc -c < "${f}") bytes)"
    fi
done

if [[ ${fail} -ne 0 ]]; then
    echo ""
    echo "=== TEST FAILED ==="
    exit 1
fi

echo ""
echo "=== TEST PASSED ==="
