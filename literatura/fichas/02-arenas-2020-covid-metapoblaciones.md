# Ficha 02: Modeling the Spatiotemporal Epidemic Spreading of COVID-19 and the Impact of Mobility and Social Distancing Interventions

- **Referencia:** A. Arenas, W. Cota, J. Gómez-Gardeñes, S. Gómez, C. Granell, J. T. Matamalas, D. Soriano-Paños, B. Steinegger, Phys. Rev. X 10, 041055 (2020). DOI 10.1103/PhysRevX.10.041055 (publicado el 18 de diciembre de 2020; acceso abierto, CC BY 4.0).
- **Estado:** leído (PDF de la revista, texto extraído con pypdf; las fórmulas pueden tener errores de extracción y habría que contrastarlas con el PDF antes de implementarlas) · sin reproducir
- **Fecha:** 2026-10-29
- **Nota de autoría:** Yamir Moreno **no** es autor. Según las afiliaciones del artículo, Cota, Gómez-Gardeñes, Granell y Soriano-Paños pertenecen al GOTHAM Lab-BIFI (Universidad de Zaragoza); Arenas, Gómez y Steinegger, a la Universitat Rovira i Virgili.

## Resumen

Modelo de metapoblaciones por municipios, con tres grupos de edad y diez compartimentos clínicos, escrito con el método MMCA (Microscopic Markov Chain Approach): ecuaciones deterministas en tiempo discreto (1 paso = 1 día) sobre probabilidades, no sobre fracciones continuas. Se aplica a España con movilidad real (INE). Contiene tres aportaciones: el modelo y su validación con muertes y casos diarios (Fig. 2), una expresión analítica para el número reproductivo efectivo R(t) (Sec. III) y un confinamiento crítico κ₀ᶜ que separa "aplanar" de "doblar" la curva (Sec. IV).

## Modelo

- **Compartimentos** por municipio i y edad g: S, E, A (asintomático infeccioso), I, PH (pre-UCI), PD (pre-fallecimiento), HR y HD (en UCI, recuperará o fallecerá), R, D. Edades: Y (hasta 25), M (26–65), O (más de 65).
- **Dinámica (Ecs. 2–11):** transiciones lineales entre compartimentos con tasas η, α, μ, ζ, λ, ψ, χ y probabilidades θ, γ, ω. El único término no lineal es Π, la probabilidad de contagio.
- **Contagio (Ecs. 12–19):** Π = (1 − p)·P_i + p·Σ_j R_ij·P_j (primero en la residencia, luego en el destino de movilidad). P depende de β_A, β_I, de los contactos por persona ⟨k⟩, de una función de densidad f(x) = 1 + (1 − e^(−ξx)), de la matriz de contactos entre edades C y de los infecciosos presentes en el municipio (tras el desplazamiento). Solo los adultos se mueven entre municipios: jóvenes y mayores no salen del suyo (R = δ_ij).
- **Confinamiento (Ecs. 20–27):** fracción confinada κ₀(t), distanciamiento δ en los no confinados, y permeabilidad de los hogares φ. Se añade un compartimento CH (confinados en el hogar).
- **Calibración:** ABC (approximate Bayesian computation) con muertes diarias hasta el 8 de abril de 2020. Se calibran β_I, η, α, μ, φ y δ; el resto sale de la literatura. κ₀(t) se toma de la reducción de movilidad laboral del INE (base: 19 de noviembre de 2019).
- **Parámetros calibrados (Tablas II y IV):** β_I = 0.075, β_A = 0.5·β_I (supuesto), η = 2.444 d⁻¹, α⁻¹ = 2.756 d (M, O) y 1/5.671 d (Y), μ⁻¹ = 3.915 d (M, O) y 1 d (Y), φ = 0.174, δ = 0.207. Contactos medios ⟨k⟩ = (11.8, 13.3, 6.8); ξ = 0.01; tamaño medio del hogar σ = 2.5.

## Resultados analíticos

- **R(t) (Ecs. 33–43):** expresión explícita a partir de las probabilidades del modelo. Versión estática: R = (β_A/α + β_I/μ)·Q, con Q el número de contactos susceptibles. Se valida contra EpiEstim (Apéndice E); la diferencia en el día en que R cae bajo 1 es de unos 2 días.
- **Aproximación de un parche (Ec. 44–51):** R depende cuadráticamente de κ₀ y se obtiene κ₀ᶜ resolviendo A·κ₀² − B·κ₀ + C = 0. La transición R = 1 separa "aplanar" de "doblar" la curva y desaparece para φ ≳ 0.83.
- **España:** κ₀ entre 0.6 y 0.7 mantiene R < 1 en mayo; R(t) crece al relajar el confinamiento. Los adultos pueden tener R_M > 1 con R global < 1.

## Datos y código

- **Código oficial:** `MMCAcovid19.jl` (Julia, licencia GNU AGPL-3), de J. T. Matamalas y S. Gómez; documentación en jtmatamalas.github.io/MMCAcovid19.jl. Ya existe, así que cualquier reproducción en Rust sería una verificación independiente, no una primicia.
- **Datos:** población, superficie y movilidad (encuesta de desplazamientos de 2011) del INE a nivel de municipio; matriz de contactos de un estudio social internacional (ref. [42] del artículo). **Por localizar y verificar** el acceso a los datos concretos y su formato.

## Cómo reproducirlo en Rust (borrador de plan, sin hacer)

1. **Núcleo con datos sintéticos:** implementar las Ecs. 2–19 con unos pocos parches y 3 edades, comprobando que cada fila de probabilidades suma 1 (invariante dado por el artículo).
2. **Verificar R(t)** (Ec. 42) contra una estimación independiente a partir de la incidencia simulada, como en el Apéndice E.
3. **Confinamiento crítico:** comparar la fórmula cerrada κ₀ᶜ (Ec. 51) con el barrido numérico de κ₀ (Fig. 5), de forma parecida a lo hecho con la estrella en la ficha 01.
4. **Cotejo con el código en Julia**, si se consigue ejecutar, en un caso pequeño.
5. **Solo después:** datos reales de España a nivel de municipio, sujeto a que los datos sean de acceso libre.

## Posibles aportes (a valorar, sin afirmar novedad)

- Implementación independiente y abierta en Rust, con comprobación numérica de las identidades del modelo.
- Estudio de rendimiento frente a Julia en el barrido de κ₀ y φ (sin asumir de antemano quién gana).
- Control de convergencia respecto a la discretización diaria y a las condiciones iniciales (el artículo siembra la epidemia en un patrón concreto; por verificar en el material suplementario).

## Pendiente por verificar

- Fórmulas exactas (Ecs. 13, 14, 23) contra el PDF original, no contra el texto extraído.
- Cómo se siembran los primeros casos y la condición inicial, y la forma exacta de la función de densidad f.
- Material suplementario (escenarios de contención).
- Disponibilidad y formato de los datos del INE y de las matrices de contactos.
