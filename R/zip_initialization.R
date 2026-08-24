# ============================================================
# Archivo: zip_initialization.R
# Autor: Francisco Vázquez Chávez
# Descripción:
#   Funciones de inicialización para el modelo ZIP-SSM:
#   - vals_ini(): estado inicial vía suavizado local-level
#   - ini_p(): estimación inicial de p
#   - E_wt(): esperanza del indicador de cero estructural
#   - params_ini(): parámetros iniciales para optimización
# ============================================================


# ------------------------------------------------------------
# 1. Inicialización del estado latente con statespacer
# ------------------------------------------------------------
# Obtiene:
#   - theta_ini: nivel suavizado (estado inicial)
#   - ee_theta_ini: error estándar del estado
#
# Entrada:
#   y  -> serie observada (ts o vector)
#   x  -> lista de covariables (opcional)
#
# Salida:
#   lista con:
#     theta_ini
#     ee_theta_ini
# ------------------------------------------------------------

vals_ini <- function(yt) {
    n <- length(yt)

    # Modelo local-level:
    # theta_{t+1} = theta_t + eta_t
    # y_t = theta_t + eps_t

    # Varianzas iniciales (puedes ajustarlas si quieres)
    s2e <- 1 # varianza del error de observación
    s2w <- 1 # varianza del error del estado

    # Inicialización
    theta_pred <- numeric(n)
    theta_filt <- numeric(n)
    P_pred <- numeric(n)
    P_filt <- numeric(n)

    # Valores iniciales
    theta_filt[1] <- yt[1]
    P_filt[1] <- 10 # varianza inicial grande

    # Filtro de Kalman
    for (t in 2:n) {
        # Predicción
        theta_pred[t] <- theta_filt[t - 1]
        P_pred[t] <- P_filt[t - 1] + s2w

        # Actualización
        Kt <- P_pred[t] / (P_pred[t] + s2e)
        theta_filt[t] <- theta_pred[t] + Kt * (yt[t] - theta_pred[t])
        P_filt[t] <- (1 - Kt) * P_pred[t]
    }

    return(list(
        theta_ini = theta_filt,
        ee_theta_ini = sqrt(P_filt)
    ))
}


# ------------------------------------------------------------
# 2. Estimación inicial de p
# ------------------------------------------------------------
# p = proporción de ceros observados
#
# Entrada:
#   yt -> serie observada
#
# Salida:
#   p inicial
# ------------------------------------------------------------

ini_p <- function(yt) {
    n0 <- sum(yt == 0)
    p_ini <- n0 / length(yt)
    return(p_ini)
}


# ------------------------------------------------------------
# 3. Esperanza del indicador de cero estructural
# ------------------------------------------------------------
# E[w_t | y_t, theta_t, p]
#
# Fórmula:
#   Si y_t = 0:
#       E[w_t] = p / (p + (1 - p) * exp(-exp(theta_t)))
#   Si y_t > 0:
#       E[w_t] = 0
#
# Entrada:
#   p      -> probabilidad de cero estructural
#   theta  -> estado latente (log lambda)
#   y      -> serie observada
#
# Salida:
#   vector E[w_t]
# ------------------------------------------------------------

E_wt <- function(p, theta, y) {
    n <- length(y)
    wt <- numeric(n)

    for (t in 1:n) {
        if (y[t] == 0) {
            wt[t] <- p / (p + (1 - p) * exp(-exp(theta[t])))
        } else {
            wt[t] <- 0
        }
    }

    return(wt)
}


# ------------------------------------------------------------
# 4. Parámetros iniciales para optimización
# ------------------------------------------------------------
# psi = (sigma2_eta, beta)
#
# Entrada:
#   y      -> serie transformada (yt_)
#   x      -> covariable
#   theta  -> estado suavizado
#
# Salida:
#   vector psi_ini = (var(diff(theta)), beta_ini)
# ------------------------------------------------------------

params_ini <- function(y, x, theta) {

    # Varianza inicial del ruido del estado
    sigma2_eta_ini <- var(diff(theta))

    # Estimación inicial de beta vía regresión lineal
    beta_ini <- as.numeric(solve(t(x) %*% x) %*% t(x) %*% y)

    return(c(sigma2_eta_ini, beta_ini))
}
