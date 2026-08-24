# ============================================================
# Archivo: zip_forecast.R
# Autor: Francisco Vázquez Chávez
# Descripción:
#   Funciones de pronóstico para el modelo ZIP-SSM.
# ============================================================


# ------------------------------------------------------------
# 1. Pronóstico del estado latente (theta)
# ------------------------------------------------------------
# Modelo:
#   theta_{t+1} = theta_t + beta * x_{t+1} + eta_{t+1}
#   eta ~ N(0, sigma2_eta)
#
# Entrada:
#   theta   -> estado suavizado (vector)
#   params  -> c(sigma2_eta, beta)
#   x       -> covariable futura (longitud h)
#   h       -> horizonte de pronóstico
#
# Salida:
#   vector con theta pronosticado
# ------------------------------------------------------------

forecast_theta <- function(theta, params, x, h = 12) {

    sigma2_eta <- params[1]
    beta       <- params[2]

    eta <- rnorm(h, mean = 0, sd = sqrt(sigma2_eta))

    n <- length(theta)
    theta_fore <- c(theta, numeric(h))

    for (j in 1:h) {
        theta_fore[n + j] <- theta_fore[n + j - 1] + beta * x[j] + eta[j]
    }

    return(theta_fore)
}


# ------------------------------------------------------------
# 2. Pronóstico del nivel (lambda = exp(theta))
# ------------------------------------------------------------
forecast_level <- function(theta_fore) {
    return(exp(theta_fore))
}


# ------------------------------------------------------------
# 3. Pronóstico ZIP (conteos)
# ------------------------------------------------------------
# Simula valores futuros ZIP usando:
#   y ~ ZIP(lambda, p)
#
# Entrada:
#   lambda_fore -> exp(theta_fore)
#   p_est       -> probabilidad estimada de cero estructural
#
# Salida:
#   serie pronosticada ZIP
# ------------------------------------------------------------

forecast_zip <- function(lambda_fore, p_est) {

    n <- length(lambda_fore)
    yt_fore <- numeric(n)

    for (t in 1:n) {
        w <- rbinom(1, 1, p_est)
        if (w == 1) {
            yt_fore[t] <- 0
        } else {
            yt_fore[t] <- rpois(1, lambda_fore[t])
        }
    }

    return(yt_fore)
}
