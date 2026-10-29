# Control de convergencia con el generador del modelo de 2006 (8 redes por caso, dl = 0.02):
# base (transitorio 40, dt estándar), transitorio x5 y dt/2, para α = 1 (panel a) y α = 0 (BA, panel d).
source ~/.cargo/env
cd "$(dirname "$0")"
CARGO_TARGET_DIR=$HOME/rt/es cargo build --release 2>&1 | tail -1
B=$HOME/rt/es/release/explosive_sync
R=resultados/conv_final
mkdir -p $R
t() { /usr/bin/time -f "  -> %es" "$@"; }
echo "a base";   t $B sweep 1 1000 6 8 0.4 0.9 0.02 $R/a_base.csv
echo "a trans";  TTRANS=200 TMEAS=20 t $B sweep 1 1000 6 8 0.4 0.9 0.02 $R/a_trans200.csv
echo "a dt/2";   DTF=0.5 t $B sweep 1 1000 6 8 0.4 0.9 0.02 $R/a_dt05.csv
echo "d base";   t $B sweep 0 1000 6 8 1.1 1.9 0.02 $R/d_base.csv
echo "d dt/2";   DTF=0.5 t $B sweep 0 1000 6 8 1.1 1.9 0.02 $R/d_dt05.csv
echo "d trans";  TTRANS=200 TMEAS=20 t $B sweep 0 1000 6 8 1.1 1.9 0.02 $R/d_trans200.csv
echo TERMINADO
