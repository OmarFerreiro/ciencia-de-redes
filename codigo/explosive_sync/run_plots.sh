# Gráficas SVG a partir de los CSV del generador final.
source ~/.cargo/env
cd "$(dirname "$0")"
B=$HOME/rt/es/release/explosive_sync
R=resultados
$B plot $R/conv_final/d_base.csv $R/conv_final/d_base.svg "BA (α=0), N=1000, ⟨k⟩=6, media de 8 redes"
$B plot $R/conv_final/a_base.csv $R/conv_final/a_base.svg "Red aleatoria (α=1), N=1000, ⟨k⟩=6, media de 8 redes"
