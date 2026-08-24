# ZIP‑SSM

Modelo Zero‑Inflated Poisson State‑Space Model (ZIP‑SSM) implementado en R. Este paquete permite estimar series de tiempo con inflación de ceros bajo una estructura de espacio de estados, combinando:

- Inicialización mediante filtro local‑level

- Aproximación gaussiana

- Verosimilitud con parametrización estable (log(sigma))

- Estimación del parámetro de inflación de ceros mediante EM

- Muestreo por importancia para obtener la verosimilitud exacta

- Funciones de simulación para pruebas y validación

El objetivo del paquete es ofrecer una implementación modular, transparente y reproducible del ZIP‑SSM, adecuada para investigación, docencia y análisis aplicado.

 ## Instalación

Clona el repositorio y carga el paquete manualmente:
````
# Clonar repositorio
# git clone https://github.com/usuario/ZIP-SSM.git

# Cargar funciones
source("R/zip_initialization.R")
source("R/zip_densities.R")
source("R/zip_likelihood.R")
source("R/zip_EM.R")
source("R/zip_importance_sampling.R")
source("R/zip_simulation.R")
````
Estructura del paquete
````
ZIP-SSM/
├── R/
│   ├── zip_initialization.R
│   ├── zip_densities.R
│   ├── zip_likelihood.R
│   ├── zip_EM.R
│   ├── zip_importance_sampling.R
│   └── zip_simulation.R
├── tests/
│   └── prueba_minima.R
├── DESCRIPTION
├── NAMESPACE
└── README.md
````

Autor

Francisco Ariel Vázquez Chávez (@fcoavc)
Economista y estadístico.

Licencia

MIT
