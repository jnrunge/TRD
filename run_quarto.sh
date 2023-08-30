export XDG_RUNTIME_DIR=$(mktemp -d)
quarto convert Draft.ipynb
quarto render Draft.qmd
