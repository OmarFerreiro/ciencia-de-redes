# Control de convergencia: ventana de la transición, 8 redes por caso.
source ~/.cargo/env
cd /mnt/c/Users/jorge/proyectos/ciencia-de-redes/codigo/explosive_sync
B=$HOME/rt/es/release/explosive_sync
R=resultados/conv
mkdir -p $R
t() { /usr/bin/time -f "  -> %es" "$@"; }
echo "ER base";        t $B sweep 1 1000 6 8 1.2 1.6 0.02 $R/er_base.csv
echo "ER transitorio"; TTRANS=200 TMEAS=20 t $B sweep 1 1000 6 8 1.2 1.6 0.02 $R/er_trans200.csv
echo "ER dt/2";        DTF=0.5 t $B sweep 1 1000 6 8 1.2 1.6 0.02 $R/er_dt05.csv
echo "ER N=4000";      t $B sweep 1 4000 6 8 1.2 1.6 0.02 $R/er_n4000.csv
echo "BA transitorio"; TTRANS=200 TMEAS=20 t $B sweep 0 1000 6 8 1.1 1.9 0.02 $R/ba_trans200.csv
echo TERMINADO
