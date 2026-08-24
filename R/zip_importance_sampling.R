# ============================================================
# Archivo: zip_importance_sampling.R
# Autor: Francisco Vázquez Chávez
# Descripción:
#   Muestreo por importancia para el modelo ZIP-SSM
# ============================================================


# ------------------------------------------------------------
# 1. Una iteración del muestreo por importancia
# ------------------------------------------------------------
importance_sampling_step <- function(yt, yt_, theta, ee_theta, p_est, beta, x) {

    n <- length(yt)

    # Simular estado
    theta_sim <- rnorm(n, mean = theta, sd = ee_theta)

    # Nivel simulado
    nivel_sim <- theta_sim + beta * x

    # Densidad ZIP
    py <- sapply(1:n, function(t) f1(yt[t], nivel_sim[t], p_est))

    # Densidad normal aproximada
    gy <- f2(yt_, theta_sim, ee_theta)

    # Peso de importancia
    wi <- sum(log(py)) / sum(log(gy))

    return(list(
        wi  = wi,
        xi  = nivel_sim * wi,
        xi2 = nivel_sim^2 * wi,
        xxi = (1 - p_est) * exp(nivel_sim) * wi,
        xxi2 = ((1 - p_est) * exp(nivel_sim))^2 * wi
    ))
}

# ------------------------------------------------------------
# 2. Muestreo por importancia completo
# ------------------------------------------------------------
zip_importance_sampling <- function(yt, yt_, theta, ee_theta, p_est, beta, x, m = 1000) {

    n <- length(yt)

    W   <- numeric(m)
    X   <- matrix(0, n, m)
    XX  <- matrix(0, n, m)
    X2  <- matrix(0, n, m)
    XX2 <- matrix(0, n, m)

    for (i in 1:m) {

        step <- importance_sampling_step(
            yt = yt,
            yt_ = yt_,
            theta = theta,
            ee_theta = ee_theta,
            p_est = p_est,
            beta = beta,
            x = x
        )

        W[i]    <- step$wi
        X[, i]  <- step$xi
        XX[, i] <- step$xi2
        X2[, i] <- step$xxi
        XX2[, i] <- step$xxi2
    }

    # Señal
    num <- sum(W)
    est_signal <- rowSums(X) / num
    var_signal <- rowSums(XX) / num - est_signal^2
    ee_signal <- sqrt(var_signal)

    # Nivel
    est_level <- rowSums(X2) / num
    var_level <- rowSums(XX2) / num - est_level^2
    ee_level <- sqrt(var_level)

    # Verosimilitud exacta corregida
    s2w <- var(W)
    lv_exact <- log(mean(W)) + s2w / (2 * m * mean(W)^2)

    return(list(
        signal_est = est_signal,
        level_est = est_level,
        ee_signal = ee_signal,
        ee_level = ee_level,
        lv_exact = lv_exact
    ))
}
