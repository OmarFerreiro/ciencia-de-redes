# Comparación con la Fig. 1 del artículo (ventanas de λ de cada panel), con el generador del modelo
# BA–ER de 2006 (PRE 73, 056124). 4 redes por caso, paso 0.04, N=1000, <k>=6.
source ~/.cargo/env
cd /mnt/c/Users/jorge/proyectos/ciencia-de-redes/codigo/explosive_sync
CARGO_TARGET_DIR=$HOME/rt/es cargo build --release 2>&1 | tail -1
B=$HOME/rt/es/release/explosive_sync
R=resultados/fig1
mkdir -p $R
t() { /usr/bin/time -f "  -> %es" "$@"; }
echo "(a) alpha=1";   t $B sweep 1   1000 6 4 0.2 1.2 0.04 $R/a_alpha1.csv
echo "(b) alpha=0.6"; t $B sweep 0.6 1000 6 4 0.4 1.4 0.04 $R/b_alpha06.csv
echo "(c) alpha=0.2"; t $B sweep 0.2 1000 6 4 0.8 1.8 0.04 $R/c_alpha02.csv
echo "(d) alpha=0";   t $B sweep 0   1000 6 4 0.8 1.8 0.04 $R/d_alpha0.csv
echo TERMINADO
