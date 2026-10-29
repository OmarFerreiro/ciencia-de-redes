# Reproduje la sincronización explosiva en Rust: lo que el artículo no dice y cambia el resultado

*Serie "Ciencia de Redes": lectura y reproducción de trabajos del grupo BIFI/COSNET (Universidad de Zaragoza).*

> **Borrador.** Antes de publicar: (1) hacer público el repositorio y poner el enlace real, (2) revisar con calma las cifras y las citas, (3) decidir si se adjunta la gráfica `ba_alpha0_fine.svg` (es de mis simulaciones; las figuras del artículo original no se reproducen).

---

## En una frase

Reproduje en Rust los resultados centrales de *Explosive Synchronization Transitions in Scale-free Networks* (Gómez-Gardeñes, Gómez, Arenas y Moreno, PRL 2011) y aprendí que un detalle que el artículo no explica, cómo se genera la red, cambia por completo el resultado.

## El trabajo original

En un modelo de Kuramoto sobre una red, cada nodo oscila con una frecuencia propia y se acopla a sus vecinos con intensidad λ. El artículo muestra que, si la frecuencia de cada nodo es igual a su grado (ω = k) y la red es libre de escala (Barabási–Albert), la sincronización aparece de golpe, con histéresis: la transición es de primer orden. En redes más homogéneas (Erdős–Rényi), la transición es continua. Además, resuelven analíticamente un grafo estrella: un hub de grado K y K hojas.

> Referencia: J. Gómez-Gardeñes, S. Gómez, A. Arenas, Y. Moreno, "Explosive Synchronization Transitions in Scale-free Networks", Phys. Rev. Lett. 106, 128701 (2011). DOI: 10.1103/PhysRevLett.106.128701 · arXiv:1102.4823

## Qué hice

Escribí el modelo en Rust (integrador Runge–Kutta de 4.º orden, redes en memoria contigua, semillas fijas) y lo comparé con la teoría y con la Fig. 1 del artículo.

**1. El grafo estrella coincide con la teoría.** Para K entre 5 y 70, el acoplamiento crítico simulado coincide con λ_c = (K−1)/(K+1) dentro de la resolución del barrido (diferencias de 0.002 a 0.009 con paso 0.01). En la rama enganchada, la simulación reproduce la curva analítica exacta con un error cuadrático medio menor que 5·10⁻⁶.

Un matiz de definición: el artículo da r_c = K/(K+1), calculado como promedio de cosenos. El módulo exacto del parámetro de orden es √(K²+1)/(K+1), un 0.5 % mayor para K = 10, y la simulación sigue el módulo. No es un error del artículo, solo una cuestión de cómo se define r.

**2. Redes BA, ER y mezclas.** Con redes de N = 1000 y ⟨k⟩ = 6, mis curvas reproducen la posición de las transiciones de los paneles (a), (b) y (c) de la Fig. 1, y la rama hacia atrás de BA (caída en λ ≈ 1.3). El salto hacia adelante de BA me sale algo más tarde que el publicado y con dispersión entre redes. Esa diferencia sigue abierta.

(Los valores del artículo los leí a ojo de las figuras, con un error de ±0.03.)

## Lo que no estaba en el artículo

Mi primera versión daba un Erdős–Rényi con transición **abrupta** y con histéresis, en λ ≈ 1.3, el doble que en el artículo. La causa: generaba la red haciendo crecer los nodos y eligiendo destinos al azar **solo entre los nodos ya presentes**. Eso concentra enlaces en los nodos antiguos y crea hubs que un ER real no tiene.

El artículo remite al modelo de Gómez-Gardeñes y Moreno (2006), en el que el número de nodos es fijo desde el principio y los destinos uniformes se eligen **entre todos los nodos**. Con ese generador, las transiciones caen donde dice el artículo.

> Referencia: J. Gómez-Gardeñes, Y. Moreno, "From scale-free to Erdős–Rényi networks", Phys. Rev. E 73, 056124 (2006). arXiv:cond-mat/0604100

Moraleja: en redes complejas, el generador es parte del resultado. Y el artículo tampoco especifica el integrador, el paso ni el tiempo de transitorio, así que esos también hay que anotarlos y comprobar que no cambian nada. Con dt a la mitad los resultados no cambiaron. El transitorio sí influye cerca de la transición.

## Una nota sobre rendimiento (sin exagerar)

Con el mismo problema (red BA con N = 1000, mismo algoritmo, un solo hilo) en mi portátil:

| Versión | µs por paso RK4 |
|---|---|
| Rust | ≈ 170 |
| C++ (`-O3`, con `sincos`) | ≈ 177 |
| Python con NumPy + SciPy (matriz dispersa) | ≈ 403 |
| Python puro | ≈ 3260 |

Rust y C++ rinden lo mismo (una primera diferencia del 20 % venía de cómo llamaba a seno y coseno). NumPy queda unas 2,4 veces por detrás con N = 1000, y Python puro unas 20 veces. Es una máquina, una red y un hilo: no extrapolo a otros casos. Lo que Rust me aporta aquí no es velocidad frente a C++, sino seguridad de memoria y un paralelismo más fácil de escribir sin errores.

## Qué aporta esto (y qué no)

- **Aporta:** una reproducción independiente, abierta y con semillas; un control numérico explícito; y una lista de los detalles del protocolo que el artículo no fija.
- **No aporta:** un resultado científico nuevo. Es una reproducción de un trabajo de 2011. No afirmo ningún resultado inédito.

## Código

Repositorio: **<ENLACE AL REPOSITORIO: pendiente de hacerlo público>** · Commit: `00652d3` · Semillas 1 a 8 en los barridos finos.

## Siguiente paso

Cerrar la diferencia del salto hacia adelante en BA con un transitorio mayor y un barrido más fino, reproducir la Fig. 3 (distintos exponentes γ) y pasar a otro trabajo del grupo.

Agradezco el trabajo original a sus autores. Si alguien del grupo ve algo mal en mi lectura, me encantaría que me lo dijera.

#CienciaDeRedes #RedesComplejas #Rust #Sincronización #FísicaDeSistemasComplejos
