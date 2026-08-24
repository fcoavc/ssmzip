# 1. Cargar módulos
source("R/rzip_simulation.R")
source("R/zip_initialization.R")
source("R/zip_densities.R")
source("R/zip_likelihood.R")
source("R/zip_EM.R")
source("R/zip_importance_sampling.R")

# 2. Simular datos ZIP simples
set.seed(123)
n <- 50

theta_true <- sim_rw(n, s2e = 0.2)
datos <- rzip(n, theta_true, p = 0.3)
yt <- datos$yt

# 3. Inicialización
ini <- vals_ini(yt)
theta_ini <- ini$theta_ini
ee_theta_ini <- ini$ee_theta_ini

# 4. Aproximación gaussiana
gd <- gaussian_density(yt, theta_ini, num = 1)
yt_ <- gd$yt_
theta_suav <- gd$theta_suav
ee_theta_suav <- gd$ee_theta_suav
Ht <- gd$H_t

# 5. Estimación de parámetros
x <- rep(1, n)
psi_ini <- params_ini(yt_, x, theta_suav)

opt <- optim(
    par     = psi_ini,
    fn      = loglg,
    y       = yt_,
    x       = x,
    Ht      = Ht,
    method  = "BFGS"
)

Psi_est <- opt$par
beta_est <- Psi_est[2]

# 6. Estimación de p
p_est <- EM(yt, theta_suav, beta_est, x)

# 7. Muestreo por importancia
imp <- zip_importance_sampling(
    yt = yt,
    yt_ = yt_,
    theta = theta_suav,
    ee_theta = ee_theta_suav,
    p_est = p_est,
    beta = beta_est,
    x = x,
    m = 200
)

imp

