#!/usr/bin/env bash
# One-shot publication to GitHub + GitHub Pages.
#
# Run this from the repository root on a machine where you are logged in to
# GitHub (either the `gh` CLI, or git with a credential helper / PAT).
#
#   bash publish.sh                 # uses the defaults below
#   OWNER=someone-else bash publish.sh
#
# What it does:
#   1. initialises the repo and makes the first commit if that has not happened
#   2. creates the public repository on GitHub (gh) or asks you to create it
#   3. pushes `main`
#   4. turns on GitHub Pages serving from /docs on main
#
# The data are NOT committed (4.2 GB). A clean clone reproduces them with
#   bash scripts/download_data.sh && make all

set -euo pipefail

OWNER="${OWNER:-ashokdwivedi}"
REPO="${REPO:-xenium-breast-tme}"
DESC="${DESC:-Reproducible QC, reference-based cell typing and spatial neighbourhood analysis of a public 10x Xenium breast-cancer dataset}"

if [ ! -d .git ]; then
  git init -b main
  git add -A
  git commit -m "Reproducible Xenium breast-tumour TME analysis"
fi
git branch -M main

if command -v gh >/dev/null 2>&1; then
  if ! gh repo view "$OWNER/$REPO" >/dev/null 2>&1; then
    gh repo create "$OWNER/$REPO" --public --description "$DESC" \
      --homepage "https://$OWNER.github.io/$REPO/"
  fi
  git remote get-url origin >/dev/null 2>&1 || \
    git remote add origin "https://github.com/$OWNER/$REPO.git"
  git push -u origin main
  # Pages from the /docs folder on main. Create, or update if it already exists.
  gh api -X POST "repos/$OWNER/$REPO/pages" \
      -f "source[branch]=main" -f "source[path]=/docs" >/dev/null 2>&1 || \
  gh api -X PUT  "repos/$OWNER/$REPO/pages" \
      -f "source[branch]=main" -f "source[path]=/docs" >/dev/null
  echo
  echo "Repository: https://github.com/$OWNER/$REPO"
  echo "Website:    https://$OWNER.github.io/$REPO/   (live in ~1 minute)"
else
  cat <<EOF

The GitHub CLI (gh) is not installed, so finish the last two steps by hand:

  1. Create an empty PUBLIC repository named "$REPO" at
     https://github.com/new  (no README, no .gitignore, no licence).

  2. Then run:

       git remote add origin https://github.com/$OWNER/$REPO.git
       git push -u origin main

  3. In the repository: Settings -> Pages -> Build and deployment ->
     Source: "Deploy from a branch", Branch: "main", Folder: "/docs" -> Save.

     The site appears at https://$OWNER.github.io/$REPO/ within a minute.

EOF
fi
