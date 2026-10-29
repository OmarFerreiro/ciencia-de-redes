set -e
source ~/.cargo/env
cd "$(dirname "$0")"
CARGO_TARGET_DIR=$HOME/rt/es cargo build --release 2>&1 | tail -2
time $HOME/rt/es/release/explosive_sync starscan 0.01 resultados/starscan.csv
column -s, -t resultados/starscan.csv
