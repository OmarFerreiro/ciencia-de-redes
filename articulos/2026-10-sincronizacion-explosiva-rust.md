# Reproduje en Rust la sincronización explosiva: el generador de la red es parte del resultado

*Serie "Ciencia de Redes": lectura y reproducción de trabajos del grupo BIFI/COSNET (Universidad de Zaragoza).*

> **Borrador.** Antes de publicar: (1) hacer público el repositorio y sustituir el marcador del enlace, (2) poner una etiqueta (`v0.1`) en el commit que se cite, (3) revisar cifras y citas, (4) decidir si se adjunta la gráfica `ba_alpha0_fine.svg` (es de mis simulaciones; las figuras del artículo original no se reproducen). Palabras: ≈ 620 (el texto de LinkedIn rinde mejor con menos; candidatos a recortar: el párrafo de rendimiento y el de "qué aporta"). Decisión tuya: la introducción menciona que se usaron herramientas de IA para leer, programar y revisar; es cierto, y encaja con el tema, pero puedes quitarlo o ampliarlo.

---

Últimamente parece que las nuevas tecnologías sirven sobre todo para producir contenido: más textos, más imágenes, más vídeos, a toda velocidad y casi siempre iguales. A mí me interesa lo contrario: usar estas herramientas en algo que se pueda comprobar y que sirva de verdad. Para entender un problema científico, escribir código que funcione y dejarlo abierto para que otra persona lo repita y lo corrija.

Así que me puse un ejercicio concreto: coger un artículo científico serio y reproducirlo desde cero, con un lenguaje moderno (Rust) y con ayuda de herramientas de IA para leer, programar y revisar, comprobando cada resultado contra la teoría y contra lo publicado. Sin atajos: si algo no coincide, se cuenta.

Y quiero agradecer a Gómez-Gardeñes, Gómez, Arenas y Moreno un trabajo que se deja reproducir. Lo que sigue es una reproducción, no un resultado nuevo.

**El trabajo.** En un modelo de Kuramoto sobre una red, cada nodo oscila con una frecuencia propia y se acopla a sus vecinos con intensidad λ. Si la frecuencia de cada nodo es igual a su grado y la red es libre de escala, la sincronización aparece de golpe y con histéresis (primer orden). En redes más homogéneas, la transición es continua. Además, el grafo estrella (un hub y K hojas) tiene solución analítica.

> J. Gómez-Gardeñes, S. Gómez, A. Arenas, Y. Moreno, "Explosive Synchronization Transitions in Scale-free Networks", Phys. Rev. Lett. 106, 128701 (2011). DOI: 10.1103/PhysRevLett.106.128701 · arXiv:1102.4823

**Lo que hice.** Lo escribí en Rust (Runge–Kutta de 4.º orden, semillas fijas) y lo comparé con la teoría y con la Fig. 1.

- **Estrella:** para K entre 5 y 70, el acoplamiento crítico simulado coincide con λ_c = (K−1)/(K+1) dentro de la resolución del barrido (0.01), y la rama enganchada reproduce la curva analítica con error cuadrático medio menor que 5·10⁻⁶.
- **Fig. 1:** con N = 1000 y ⟨k⟩ = 6, la posición de las transiciones de los paneles (a), (b) y (c) coincide con la publicada, y la rama hacia atrás de la red Barabási–Albert también. (Los valores del artículo los leí a ojo de la figura, con un error de ±0.03.) El salto hacia adelante en BA varía de una red a otra (entre λ ≈ 1.44 y 1.64 en mis 8 redes), y el publicado (1.42, una sola red) queda en el borde inferior de ese intervalo.

**Lo que no estaba en el artículo.** Mi primera versión daba, para la red "aleatoria", una transición abrupta con histéresis y al doble de λ que el artículo. La causa: yo hacía crecer la red eligiendo destinos al azar solo entre los nodos ya presentes. Eso concentra enlaces en los nodos antiguos y crea hubs que una red aleatoria real no tiene. El modelo que cita el artículo (Gómez-Gardeñes y Moreno, Phys. Rev. E 73, 056124, 2006) fija el número de nodos desde el principio y elige los destinos entre todos. Con ese generador, las transiciones caen donde deben.

Moraleja: en redes complejas, el generador es parte del resultado. El artículo tampoco fija el integrador, el paso ni el transitorio; con dt a la mitad y con un transitorio cinco veces mayor, las curvas de la red aleatoria no cambian.

**Rendimiento, sin exagerar.** Con el mismo problema y un solo hilo, Rust y C++ rinden igual; NumPy con matriz dispersa queda unas 2,4 veces por detrás y Python puro unas 20 (una máquina, una red, N = 1000). Los detalles y las limitaciones están en el repositorio. Mi impresión es que lo que Rust me aporta aquí es seguridad de memoria más que velocidad, pero eso no lo he medido.

**Qué aporta y qué no.** Aporta una reproducción independiente, abierta y con semillas, y una lista de detalles del protocolo que condicionan el resultado. No aporta ningún resultado científico nuevo.

**Código y datos:** <ENLACE AL REPOSITORIO: pendiente de hacerlo público> · etiqueta `v0.1`

Si alguien del grupo ve algo mal en mi lectura, me encantará que me lo diga. Siguiente paso: la Fig. 3 (distintos exponentes γ).

#CienciaDeRedes #RedesComplejas #Rust #Sincronización #FísicaDeSistemasComplejos
