set -e
cd /mnt/c/Users/jorge/proyectos/ciencia-de-redes/codigo/explosive_sync
B=$HOME/rt/es/release/explosive_sync
R=resultados
$B star 10 0.02 $R/star_K10.csv
$B plot $R/star_K10.csv $R/star_K10.svg "Estrella K=10 (teoría: λc=0.818, rc=0.909)"
$B sweep 1 1000 6 8 1.0 2.2 0.02 $R/er_alpha1_fine.csv
$B plot $R/er_alpha1_fine.csv $R/er_alpha1_fine.svg "ER (α=1), N=1000, ⟨k⟩=6, 8 redes"
$B sweep 0 1000 6 8 1.0 2.2 0.02 $R/ba_alpha0_fine.csv
$B plot $R/ba_alpha0_fine.csv $R/ba_alpha0_fine.svg "BA (α=0), N=1000, ⟨k⟩=6, 8 redes"
echo TERMINADO