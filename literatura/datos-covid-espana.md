# Datos necesarios para reproducir Arenas et al. (PRX 10, 041055) con España

Estado a 2026-10-30. "Verificado" = visto en una página de la fuente. "Por contrastar" = solo en un resumen de búsqueda.

| Dato | Qué dice el artículo | Qué se ha encontrado | Estado |
|---|---|---|---|
| Población, pirámide y tamaño medio del hogar por municipio | INE (ref. [50], solo la web del INE, sin tabla concreta) | No se ha localizado la tabla concreta (probablemente Padrón/Censo) | **por localizar** |
| Superficie por municipio (s_i) | INE | Sin localizar | **por localizar** |
| Matriz de movilidad de adultos entre municipios R_ij | Encuesta de desplazamientos del INE de 2011 (Censo 2011), municipal | Según un resumen de búsqueda de un artículo de la Revista Índice (no abierto: la URL devolvió 404): matriz de unos 8.116 × 8.116 municipios, no publicada completa; el INE ofrecería una explotación específica para investigación. El Censo 2011 preguntó el lugar de trabajo o estudio a una muestra de ~9 % | **por contrastar; probablemente requiere solicitud al INE** |
| Matriz de contactos por edad C_gh, y ⟨k⟩ | Ref. [42]: Prem, Cook, Jit, "Projecting social contact matrices in 152 countries using contact surveys and demographic data", PLoS Comput. Biol. 13, e1005697 (2017) | Verificada la referencia (DOI 10.1371/journal.pcbi.1005697). España no es país POLYMOD, así que su matriz es una **proyección**. Repositorio de las matrices: github.com/kieshaprem/synthetic-contact-matrices (licencia GPL-3.0 según la página; la carpeta con los datos de España no se ha abierto) | referencia verificada; **fichero de España por localizar**; el paso a los 3 estratos (Tabla I) por reproducir |
| Reducción de movilidad κ₀(t) | "Mobility reduction" del INE (ref. [44], covid_movilidad_en.htm); el texto dice que la serie se estima como reducción de la movilidad laboral frente al 19 de noviembre de 2019 | La página del INE ofrece flujos origen-destino entre "áreas de movilidad" (telefonía móvil, desde el 16 de marzo de 2020, cada dos días y luego diarios; solo flujos de más de 15 personas), con CSV, shapefiles y tablas por CCAA y provincias, con licencia CC BY-SA 4.0 según su pie (el enlace apunta a CC BY 4.0). **La página no define un indicador de "reducción de movilidad" ni habla de movilidad laboral** | datos de flujos verificados; **la serie κ₀(t) exacta de los autores no está publicada que yo haya visto**, habría que reconstruirla o pedirla |
| Casos y muertes diarios por comunidad (Fig. 2) | Ref. [37], informe nº 31 de la Red Nacional de Vigilancia Epidemiológica (ISCIII, 14 de mayo de 2020) | Referencia leída en el artículo; sin localizar el archivo | **por localizar** |
| Seroprevalencia (89.4 % sin diagnosticar) | Ref. [38]: Pollán et al., ENE-COVID, Lancet 396, 535 (2020) | Referencia leída en el artículo; no abierta | por verificar |
| Fecha inicial y siembra | Modelo desde el 9 de febrero de 2020; casos notificados hasta el 3 de marzo retrotraídos al 9 de febrero, con A₀ calibrado | Texto del artículo (Apéndice C) | leído |

## Consecuencias para el plan

1. **Reproducir a escala de municipios de España depende de la matriz de movilidad de 2011**, que probablemente haya que solicitar al INE. Mientras no se tenga, la reproducción numérica de las Figs. 2 y 3 queda bloqueada.
2. **Alternativa sin esa matriz:** una versión a nivel de provincia o de las "áreas de movilidad" del INE (datos abiertos de telefonía móvil de 2020). Sería otro modelo (otra escala y otra fuente de movilidad), así que **no sería una reproducción del artículo**, sino una variante, y habría que decirlo.
3. **Código oficial (Julia):** la documentación de MMCAcovid19.jl menciona una sección "Model data"; no se ha podido leer desde la página inicial. Podría traer los datos de España ya preparados. Sin comprobar.
4. Lo que se puede hacer ya, sin datos: más pruebas con datos sintéticos (convergencia en el paso temporal, comparación con el código en Julia en un caso pequeño).
