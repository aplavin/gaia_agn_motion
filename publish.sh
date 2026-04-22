#!/usr/bin/env bash
set -euo pipefail

REPO="aplavin/gaia_agn_motion"
MSG="Code for \"Gaia Sees Blazars Move: Locating Optical Flares Using Astrometry\", A. Plavin, 2026"

# --- Phase 1: Push clean orphan commit ---
COMMIT=$(git commit-tree HEAD^{tree} -m "$MSG")
git push origin "$COMMIT":refs/heads/master --force
echo "Published as $COMMIT"

# --- Opt-in: Create GitHub release for Zenodo ---
if [[ "${1:-}" == "--release" ]]; then
    VERSION="${2:?Usage: $0 --release <version>}"

    git tag "$VERSION" "$COMMIT"
    git push origin "$VERSION"

    gh release create "$VERSION" \
        --repo "$REPO" \
        --target "$COMMIT" \
        --title "$VERSION: $MSG"

    echo "Release: https://github.com/$REPO/releases/tag/$VERSION"
fi
