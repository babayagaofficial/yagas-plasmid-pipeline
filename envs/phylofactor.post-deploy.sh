#!/usr/bin/env bash
set -euo pipefail
Rscript -e 'remotes::install_github("reptalex/phylofactor@b87f652ca7a6422e5c6da37563f69ea45e69f60e", dependencies = FALSE, upgrade = "never")'
