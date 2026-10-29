# Reproduje en Rust la sincronización explosiva: el generador de la red es parte del resultado

*Serie "Ciencia de Redes": lectura y reproducción de trabajos del grupo BIFI/COSNET (Universidad de Zaragoza).*

> **Borrador.** Antes de publicar: (1) hacer público el repositorio y sustituir el marcador del enlace, (2) poner una etiqueta (`v0.1`) en el commit que se cite, (3) revisar cifras y citas, (4) decidir si se adjunta la gráfica `codigo/explosive_sync/resultados/conv_final/d_base.svg` (BA, 8 redes; es de mis simulaciones, las figuras del artículo original no se reproducen). Decisión tuya: la introducción menciona que se usaron herramientas de IA para leer, programar y revisar; es cierto y encaja con el tema, pero puedes quitarlo o ampliarlo.

---

Últimamente parece que las nuevas tecnologías sirven sobre todo para producir contenido: más textos, más imágenes, más vídeos, casi siempre iguales. Sin embargo, leo pocas referencias a su utilidad para la comprensión: de un problema científico, de una técnica, de cómo funciona algo por dentro.

Así que me puse un ejercicio: reproducir desde cero un artículo científico serio, con un lenguaje moderno (Rust) y con ayuda de herramientas de IA para leer, programar y revisar, comprobando cada resultado contra la teoría y contra lo publicado. Si algo no coincide, se cuenta.

El artículo es "Explosive Synchronization Transitions in Scale-free Networks", de Jesús Gómez-Gardeñes, Sergio Gómez, Alex Arenas y Yamir Moreno (Physical Review Letters, 2011). Les agradezco un trabajo que se deja reproducir. Lo que sigue es una reproducción, no un resultado nuevo.

**Qué dice.** Osciladores acoplados sobre una red (modelo de Kuramoto): si la frecuencia de cada nodo es su número de conexiones y la red es libre de escala, la sincronización aparece de golpe y con histéresis; en redes más homogéneas, de forma gradual. Un grafo estrella (un hub y K hojas) tiene solución analítica.

**Qué hice.** Lo programé en Rust y lo comparé con la teoría y con la Fig. 1.

- **Estrella:** para K de 5 a 70, el acoplamiento crítico simulado coincide con (K−1)/(K+1) dentro de la resolución del barrido (0.01).
- **Fig. 1** (N = 1000, ⟨k⟩ = 6): las transiciones de los paneles (a), (b) y (c) caen donde las publicadas, y la rama de vuelta de la red libre de escala también. Los valores del artículo los leí a ojo, con un error de ±0.03. El salto de ida en esa red varía de una red a otra (λ ≈ 1.44–1.64 en mis 8 redes); el publicado (1.42, una sola red) queda en el borde inferior.

**Lo que el artículo no cuenta.** Mi primera versión daba, para la red aleatoria, una transición abrupta y al doble de λ. La causa: yo hacía crecer la red eligiendo destinos solo entre los nodos ya presentes, y eso crea hubs que una red aleatoria no tiene. El modelo que el artículo cita (Jesús Gómez-Gardeñes y Yamir Moreno, "From scale-free to Erdős–Rényi networks", Phys. Rev. E 73, 056124, 2006) fija el número de nodos desde el principio y elige los destinos entre todos. Con ese generador, todo cae donde debe. En redes complejas, el generador es parte del resultado. (Con dt a la mitad y con un transitorio cinco veces mayor, las curvas de la red aleatoria no cambian.)

**Rendimiento.** Rust y C++ rinden igual; Python con NumPy, unas 2,4 veces más lento; Python puro, unas 20 (una máquina, N = 1000, un hilo). Detalles y limitaciones en el repositorio.

**Código y datos:** <ENLACE AL REPOSITORIO: pendiente de hacerlo público> · etiqueta `v0.1`

Siguiente paso: la Fig. 3. Si alguien del grupo ve algo mal en mi lectura, me encantará que me lo diga.

#CienciaDeRedes #RedesComplejas #Rust #Sincronización #FísicaDeSistemasComplejos
