# ============================================================
# Archivo: zip_likelihood.R
# Autor: Francisco Vázquez Chávez
# Descripción:
#   Funciones de verosimilitud para el modelo ZIP-SSM:
#   - loglg(): verosimilitud normal aproximada (filtro local-level)
#   - logv(): verosimilitud exacta del modelo SSM
# ============================================================


# ------------------------------------------------------------
# 1. Verosimilitud normal aproximada (loglg)
# ------------------------------------------------------------
# Esta verosimilitud se usa para optimizar psi = (sigma2_eta, beta)
# mediante un filtro local-level con covariable x.
#
# Modelo:
#   y*_t = theta_t + error_t
#   theta_t = theta_{t-1} + beta * x_t + eta_t
#
# Entrada:
#   psi   -> vector de parámetros (sigma2_eta, beta)
#   y     -> serie transformada (yt_)
#   x     -> covariable
#   Ht    -> varianza aproximada del error de observación
#
# Salida:
#   log-verosimilitud (normal)
# ------------------------------------------------------------

loglg <- function(psi, y, x, Ht) {

    n <- length(y)

    sigma2_eta <- psi[1]   # varianza del estado
    beta       <- psi[2]   # coeficiente de regresión

    # Inicialización del filtro
    at <- y[1]             # estado inicial
    Pt <- sigma2_eta + Ht[1]

    suma <- 0

    for (t in 2:n) {

        # Predicción del estado
        theta_t <- at + beta * x[t]

        # Innovación
        vt <- y[t] - theta_t

        # Varianza de la innovación
        Ft <- Pt + Ht[t]

        # Ganancia de Kalman
        Kt <- Pt / Ft

        # Acumular verosimilitud
        suma <- suma + log(Ft) + vt^2 / Ft

        # Actualización del estado
        at <- at + Kt * vt

        # Actualización de la varianza del estado
        Pt <- Pt * (1 - Kt) + sigma2_eta
    }

    return(0.5 * suma)
}


# ------------------------------------------------------------
# 2. Verosimilitud exacta del modelo SSM (logv)
# ------------------------------------------------------------
# Esta verosimilitud se usa para:
#   - comparar modelos
#   - corregir la verosimilitud normal con muestreo por importancia
#
# Modelo:
#   y*_t = theta_t + error_t
#   theta_t = theta_{t-1} + beta * x_t + eta_t
#
# Entrada:
#   psi   -> vector de parámetros (sigma2_eta, beta)
#   y     -> serie transformada (yt_)
#   x     -> covariable
#   Ht    -> varianza aproximada del error de observación
#
# Salida:
#   log-verosimilitud exacta
# ------------------------------------------------------------

logv <- function(psi, y, x, Ht) {

    n <- length(y)

    sigma2_eta <- psi[1]
    beta       <- psi[2]

    at <- y[1]
    Pt <- sigma2_eta + Ht[1]

    suma <- 0

    for (t in 2:n) {

        theta_t <- at + beta * x[t]
        vt <- y[t] - theta_t

        Ft <- Pt + Ht[t]
        Kt <- Pt / Ft

        suma <- suma + log(Ft) + vt^2 / Ft

        at <- at + Kt * vt
        Pt <- Pt * (1 - Kt) + sigma2_eta
    }

    return(-0.5 * suma)
}
