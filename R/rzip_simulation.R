# ============================================================
# Archivo: rzip_simulation.R
# Autor: Francisco Vázquez Chávez
# Descripción:
#   Funciones de simulación para el modelo ZIP-SSM:
#   - rzip(): simula datos ZIP con cero inflado
#   - sim_rw(): simula un proceso de caminata aleatoria
# ============================================================


# ------------------------------------------------------------
# 1. Simulación ZIP con cero inflado
# ------------------------------------------------------------
# yt ~ ZIP(lambda_t, p)
# p = probabilidad de cero estructural
# lambda_t = exp(theta_t)
#
# Entrada:
#   n      -> longitud de la serie
#   theta  -> vector de estados latentes (log lambda)
#   p      -> probabilidad de cero estructural
#
# Salida:
#   lista con:
#     yt -> serie simulada
#     wt -> indicador de cero estructural (1 = cero estructural)
# ------------------------------------------------------------

rzip <- function(n, theta, p) {
    yt <- numeric(n)
    wt <- numeric(n)

    for (t in 1:n) {
        wt[t] <- rbinom(1, 1, p)       # cero estructural
        if (wt[t] == 0) {
            yt[t] <- rpois(1, lambda = exp(theta[t]))
        } else {
            yt[t] <- 0
        }
    }

    return(list(yt = yt, wt = wt))
}


# ------------------------------------------------------------
# 2. Simulación ZIP con lambda explícito
# ------------------------------------------------------------
# Igual que rzip(), pero recibe lambda_t directamente
# ------------------------------------------------------------

rzip_lambda <- function(n, lambda, p) {
    yt <- numeric(n)
    wt <- numeric(n)

    for (t in 1:n) {
        wt[t] <- rbinom(1, 1, p)
        if (wt[t] == 0) {
            yt[t] <- rpois(1, lambda = lambda[t])
        } else {
            yt[t] <- 0
        }
    }

    return(list(yt = yt, wt = wt))
}


# ------------------------------------------------------------
# 3. Simulación del estado latente (caminata aleatoria)
# ------------------------------------------------------------
# theta_t = theta_{t-1} + eta_t
# eta_t ~ N(0, s2e)
#
# Entrada:
#   n      -> longitud de la serie
#   s2e    -> varianza del ruido del estado
#   alfa0  -> valor inicial del estado
#
# Salida:
#   vector theta_t
# ------------------------------------------------------------

sim_rw <- function(n, s2e, alfa0 = 0) {
    eta <- rnorm(n, mean = 0, sd = sqrt(s2e))
    theta <- numeric(n)

    theta[1] <- alfa0 + eta[1]

    for (t in 2:n) {
        theta[t] <- theta[t - 1] + eta[t]
    }

    return(theta)
}
