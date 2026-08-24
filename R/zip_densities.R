# ============================================================
# Archivo: zip_densities.R
# Autor: Francisco Vázquez Chávez
# Descripción:
#   Funciones de densidad para el modelo ZIP-SSM:
#   - f1(): densidad ZIP para un solo y_t
#   - f2(): densidad normal para aproximación gaussiana
#   - gaussian_density(): aproximación gaussiana iterada
# ============================================================


# ------------------------------------------------------------
# 1. Densidad ZIP para un solo y_t
# ------------------------------------------------------------
# Fórmulas:
#   Si y_t = 0:
#       f1 = p + (1 - p) * exp(-exp(theta))
#
#   Si y_t > 0:
#       f1 = (1 - p) * exp(y_t * theta - exp(theta) - log(y_t!))
#
# Entrada:
#   y      -> observación
#   theta  -> estado latente (log lambda)
#   p      -> probabilidad de cero estructural
#
# Salida:
#   f1(y_t | theta_t, p)
# ------------------------------------------------------------

f1 <- function(y, theta, p) {

    if (y == 0) {
        dens <- p + (1 - p) * exp(-exp(theta))
    } else {
        dens <- (1 - p) * exp(
            y * theta - exp(theta) - lfactorial(y)
        )
    }

    # Evitar ceros numéricos
    return(max(dens, 1e-12))
}


# ------------------------------------------------------------
# 2. Densidad normal para aproximación gaussiana
# ------------------------------------------------------------
# Entrada:
#   y      -> observación transformada
#   theta  -> estado latente
#   sd     -> desviación estándar
#
# Salida:
#   dnorm(y | theta, sd)
# ------------------------------------------------------------

f2 <- function(y, theta, sd) {
    return(dnorm(y, mean = theta, sd = sd))
}


# ------------------------------------------------------------
# 3. Aproximación gaussiana iterada
# ------------------------------------------------------------
# Esta función implementa la transformación:
#
#   y*_t = theta_t + exp(-theta_t) * (y_t - exp(theta_t))
#
# y calcula:
#   - yt_          -> serie transformada
#   - theta_suav   -> estado suavizado vía statespacer
#   - ee_theta     -> error estándar del estado suavizado
#   - mu_t         -> media aproximada
#   - H_t          -> varianza aproximada
#
# Entrada:
#   yt    -> serie observada ZIP
#   theta -> estado inicial (log lambda)
#   num   -> número de iteraciones de refinamiento
#
# Salida:
#   lista con:
#     yt_
#     theta_suav
#     ee_theta_suav
#     mu_t
#     H_t
# ------------------------------------------------------------

gaussian_density <- function(yt, theta_ini, num = 1) {
    n <- length(yt)

    # Transformación inicial
    yt_ <- log(yt + 1)

    # Suavizado gaussiano usando tu propio filtro
    s2e <- 1
    s2w <- 1

    theta_pred <- numeric(n)
    theta_suav <- numeric(n)
    P_pred <- numeric(n)
    P_suav <- numeric(n)

    theta_suav[1] <- theta_ini[1]
    P_suav[1] <- 10

    for (t in 2:n) {
        # Predicción
        theta_pred[t] <- theta_suav[t - 1]
        P_pred[t] <- P_suav[t - 1] + s2w

        # Actualización
        Kt <- P_pred[t] / (P_pred[t] + s2e)
        theta_suav[t] <- theta_pred[t] + Kt * (yt_[t] - theta_pred[t])
        P_suav[t] <- (1 - Kt) * P_pred[t]
    }

    # Varianza del estado suavizado
    ee_theta_suav <- sqrt(P_suav)

    # H_t para la verosimilitud normal
    H_t <- P_suav

    return(list(
        yt_ = yt_,
        theta_suav = theta_suav,
        ee_theta_suav = ee_theta_suav,
        H_t = H_t
    ))
}

