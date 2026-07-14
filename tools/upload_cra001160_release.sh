#!/usr/bin/env bash
# Upload the FULL CRA001160 vignette bundle as a GitHub Release asset and
# print the URL to store in repository secret PHENOMAPR_CRA001160_RDS_URL.
#
# Prerequisites:
#   - gh auth login
#   - local file (default): /tmp/phenomapr_data/PAAD_CRA001160_full.rds
#     build with: Rscript tools/make_cra001160_full_rds.R
#
set -euo pipefail
REPO="${REPO:-brooksbenard/PhenoMapR}"
TAG="${TAG:-vignette-data-cra001160}"
FILE="${1:-/tmp/phenomapr_data/PAAD_CRA001160_full.rds}"

if [[ ! -f "$FILE" ]]; then
  echo "Missing $FILE — run: Rscript tools/make_cra001160_full_rds.R" >&2
  exit 1
fi

if ! gh auth status >/dev/null 2>&1; then
  echo "Run: gh auth login" >&2
  exit 1
fi

if gh release view "$TAG" -R "$REPO" >/dev/null 2>&1; then
  echo "Release $TAG exists; uploading/replacing asset..."
  gh release upload "$TAG" "$FILE" -R "$REPO" --clobber
else
  gh release create "$TAG" "$FILE" -R "$REPO" \
    --title "CRA001160 full vignette data" \
    --notes "Full CRA001160 expression+metadata RDS for pkgdown (21,066 genes × 57,443 cells)."
fi

ASSET="$(basename "$FILE")"
URL="https://github.com/${REPO}/releases/download/${TAG}/${ASSET}"
echo
echo "Set repository secret PHENOMAPR_CRA001160_RDS_URL to:"
echo "  $URL"
echo
echo "gh secret set PHENOMAPR_CRA001160_RDS_URL -R $REPO -b '$URL'"
