# ============================================================
# Archivo: zip_EM.R
# Autor: Francisco Vázquez Chávez
# Descripción:
#   Estimación de p (probabilidad de cero estructural)
#   mediante el algoritmo EM para el modelo ZIP-SSM.
#
#   Funciones:
#     - EM(): estima p usando E[w_t] y actualización iterativa
# ============================================================


# ------------------------------------------------------------
# 1. Algoritmo EM para estimar p
# ------------------------------------------------------------
# Modelo ZIP:
#   y_t = 0 con probabilidad p + (1 - p) exp(-exp(theta_t))
#   y_t > 0 con probabilidad (1 - p) Poisson(exp(theta_t))
#
# EM:
#   Paso E:
#       w_t = E[w_t | y_t, theta_t, p]
#
#   Paso M:
#       p_new = sum(w_t) / n
#
# Entrada:
#   y      -> serie observada
#   theta  -> estado latente (log lambda)
#   beta   -> coeficiente de regresión
#   x      -> covariable
#   tol    -> tolerancia de convergencia
#   nmax   -> máximo número de iteraciones
#   verbose-> imprimir progreso
#
# Salida:
#   p_est -> estimación final de p
# ------------------------------------------------------------

EM <- function(y, theta, beta = 0, x, tol = 0.01, nmax = 20, verbose = FALSE) {

    n <- length(y)

    # Nivel del estado
    nivel <- theta + beta * x

    # Inicialización
    p_old <- ini_p(y)
    dif <- Inf
    iter <- 0

    # Paso E inicial
    w_t <- E_wt(p = p_old, theta = nivel, y = y)

    # Iteración EM
    while (dif > tol && iter < nmax) {

        # Paso M: actualizar p
        p_new <- sum(w_t) / n

        # Diferencia para convergencia
        dif <- abs(p_new - p_old)

        # Paso E: actualizar w_t
        w_t <- E_wt(p = p_new, theta = nivel, y = y)

        # Actualizar iteración
        iter <- iter + 1
        p_old <- p_new
    }

    if (verbose) {
        cat("Iteraciones:", iter, "\n")
        cat("Diferencia final:", dif, "\n")
        cat("p estimado:", p_new, "\n")
    }

    return(p_new)
}
