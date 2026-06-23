// ============================================================
// blca_stan.stan  —  BLCA in Stan (robustness check)
// ============================================================
// 
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

file <- file.path(cmdstan_path(), "examples", "bernoulli", "bernoulli.stan")
mod  <- cmdstan_model(file)

data_list <- list(N = 10, y = c(1,0,1,1,0,1,0,1,1,0))
fit <- mod$sample(data = data_list, chains = 4, parallel_chains = 4)

install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE)
print(fit$summary())
install_cmdstan()
cmdstanr::cmdstan_version()
## Carpetas
dir.create("data")
dir.create("output")
dir.create("scripts")


## Rtools & compilador

Sys.setenv(PATH = paste("C:/rtools45/usr/bin",
                        "C:/rtools45/x86_64-w64-mingw32.static.posix/bin",
                        Sys.getenv("PATH"), sep=";"))

Sys.which("g++")
Sys.which("make")


## Paquetes 

if (!require("rstan"))     install.packages("rstan", repos="https://cloud.r-project.org/")
if (!require("readxl"))    install.packages("readxl")
if (!require("tidyverse")) install.packages("tidyverse")

library(rstan)
library(readxl)
library(tidyverse)

rstan_options(auto_write = TRUE)
options(mc.cores = 3)

cat("Stan version:", as.character(stan_version()), "\n")
cat("rstan version:", as.character(packageVersion("rstan")), "\n")



df <- read_excel("data/DataBVF.xlsx", na = c("", "NA", "N/A", " "))

cat("Filas:", nrow(df), "\n")
cat("Columnas:", paste(names(df), collapse=", "), "\n")


## Matrices para datos Stan

prep <- function(Y_mat) {
  obs     <- ifelse(!is.na(Y_mat), 1, 0)
  Y_clean <- Y_mat
  Y_clean[is.na(Y_clean)] <- 0
  list(N   = nrow(Y_mat),
       K   = ncol(Y_mat),
       Y   = Y_clean,
       obs = obs)
}

Y_a <- df %>% select(MAT, elisaa, lifea, sda, lepa) %>% as.matrix()
Y_c <- df %>% select(MAT, elisac, lifec, sdc, lepc) %>% as.matrix()

datos_aguda <- prep(Y_a)
datos_conv  <- prep(Y_c)

cat("Datos aguda — N:", datos_aguda$N, "K:", datos_aguda$K, "\n")
cat("NA por columna aguda:\n")
print(colSums(is.na(Y_a)))

cat("\nDatos conv — N:", datos_conv$N, "K:", datos_conv$K, "\n")
cat("NA por columna conv:\n")
print(colSums(is.na(Y_c)))



# Verificar la formula conjunta para 3 tests dependientes
# usando valores de prueba
Se <- c(0.9, 0.92, 0.95)   # ELISA, LifeAssay, LeptoCheck
covS <- c(0.02, 0.01, 0.015)  # cov(2,3), cov(2,5), cov(3,5)

# P(111|D=1) bajo independencia
p111_indep <- Se[1]*Se[2]*Se[3]
cat("P(1,1,1|D=1) independencia:", p111_indep, "\n")


modelo_keddie_real <- "
data {
  int<lower=1> N;
  int<lower=0,upper=1> t1[N]; // MAT
  int<lower=0,upper=1> t2[N]; // ELISA
  int<lower=0,upper=1> t3[N]; // LifeAssay
  int<lower=0,upper=1> t4[N]; // SD Bioline
  int<lower=0,upper=1> t5[N]; // LeptoCheck
}
parameters {
  real<lower=0,upper=1> a11;
  real<lower=1-inv_logit(logit(a11)*2),upper=1> a12;
  real<lower=0,upper=1> a21;
  real<lower=1-inv_logit(logit(a21)*2),upper=1> a22;
  real<lower=0,upper=1> a31;
  real<lower=1-inv_logit(logit(a31)*2),upper=1> a32;
  real<lower=0,upper=1> a41;
  real<lower=1-inv_logit(logit(a41)*2),upper=1> a42;
  real<lower=0,upper=1> a51;
  real<lower=1-inv_logit(logit(a51)*2),upper=1> a52;
  real<lower=0,upper=1> prev;
}
transformed parameters {
  simplex[2] theta;
  vector[N] prob[5,2];
  theta[1] = 1-prev;
  theta[2] = prev;
  prob[1,1] = rep_vector(1-a11, N); prob[1,2] = rep_vector(a12, N);
  prob[2,1] = rep_vector(1-a21, N); prob[2,2] = rep_vector(a22, N);
  prob[3,1] = rep_vector(1-a31, N); prob[3,2] = rep_vector(a32, N);
  prob[4,1] = rep_vector(1-a41, N); prob[4,2] = rep_vector(a42, N);
  prob[5,1] = rep_vector(1-a51, N); prob[5,2] = rep_vector(a52, N);
}
model {
  real ps[2];
  a11~beta(47,1);  a12~beta(4,6);
  a21~beta(12,3);  a22~beta(8,4);
  a31~beta(12,3);  a32~beta(8,4);
  a41~beta(12,3);  a42~beta(4,8);
  a51~beta(12,3);  a52~beta(8,4);
  prev~beta(2,6);

  for(n in 1:N){
    for(k in 1:2){
      ps[k] = log(theta[k]) + binomial_lpmf(t1[n]|1,prob[1,k,n])
            + binomial_lpmf(t2[n]|1,prob[2,k,n])
            + binomial_lpmf(t3[n]|1,prob[3,k,n])
            + binomial_lpmf(t4[n]|1,prob[4,k,n])
            + binomial_lpmf(t5[n]|1,prob[5,k,n]);
    }
    target += log_sum_exp(ps);
  }
}
generated quantities {
  real Se_mean[5];
  real Sp_mean[5];
  for(m in 1:5){
    Se_mean[m] = mean(prob[m,2,]);
    Sp_mean[m] = mean(1-prob[m,1,]);
  }
}
"

#Compilando
m_keddie <- stan_model(model_code = modelo_keddie_real,
                       model_name = "BLCA_Keddie_real")
cat("Compilado OK\n")
cat("Modelo Keddie adaptado, listo para compilar\n")

##Prueba

df <- read_excel("data/DataBVF.xlsx", na = c("", "NA", "N/A", " "))

datos_stan <- list(
  N  = nrow(df),
  t1 = ifelse(is.na(df$MAT),    0L, as.integer(df$MAT)),
  t2 = ifelse(is.na(df$elisaa), 0L, as.integer(df$elisaa)),
  t3 = ifelse(is.na(df$lifea),  0L, as.integer(df$lifea)),
  t4 = ifelse(is.na(df$sda),    0L, as.integer(df$sda)),
  t5 = ifelse(is.na(df$lepa),   0L, as.integer(df$lepa))
)

## 1000

fit_prueba <- sampling(
  m_keddie, data = datos_stan,
  chains = 1, iter = 1000, warmup = 500,
  seed = 2024
)
print(summary(fit_prueba, pars=c("a11","a12","a21","a22","prev"))$summary[,"Rhat"])


## 5000

fit_aguda_v2 <- sampling(
  m_keddie, data = datos_stan,
  chains = 3, iter = 5000, warmup = 1500,
  seed = 2024,
  control = list(adapt_delta = 0.95)
)
print(round(summary(fit_aguda_v2)$summary[,"Rhat"], 4))


# Resolver faltantes

cat("NA en cada prueba (fase aguda):\n")
cat("MAT:", sum(is.na(df$MAT)), "\n")
cat("ELISA:", sum(is.na(df$elisaa)), "\n")
cat("LifeAssay:", sum(is.na(df$lifea)), "\n")
cat("SD Bioline:", sum(is.na(df$sda)), "\n")
cat("LeptoCheck:", sum(is.na(df$lepa)), "\n")



# Indices de pacientes CON dato observado, por prueba
idx_elisa <- which(!is.na(df$elisaa))
idx_life  <- which(!is.na(df$lifea))
idx_sd    <- which(!is.na(df$sda))
idx_lepto <- which(!is.na(df$lepa))
# MAT no tiene NA, así que todos los indices 1:N

datos_stan_v2 <- list(
  N = nrow(df),
  t1 = as.integer(df$MAT),
  
  N2 = length(idx_elisa),
  idx2 = idx_elisa,
  t2 = as.integer(df$elisaa[idx_elisa]),
  
  N3 = length(idx_life),
  idx3 = idx_life,
  t3 = as.integer(df$lifea[idx_life]),
  
  N4 = length(idx_sd),
  idx4 = idx_sd,
  t4 = as.integer(df$sda[idx_sd]),
  
  N5 = length(idx_lepto),
  idx5 = idx_lepto,
  t5 = as.integer(df$lepa[idx_lepto])
)

cat("N total:", datos_stan_v2$N, "\n")
cat("N2 (ELISA observados):", datos_stan_v2$N2, "\n")
cat("N3 (LifeAssay observados):", datos_stan_v2$N3, "\n")
cat("N4 (SD Bioline observados):", datos_stan_v2$N4, "\n")
cat("N5 (LeptoCheck observados):", datos_stan_v2$N5, "\n")

modelo_keddie_na <- "
data {
  int<lower=1> N;
  int<lower=0,upper=1> t1[N];

  int<lower=1> N2;
  int<lower=1,upper=N> idx2[N2];
  int<lower=0,upper=1> t2[N2];

  int<lower=1> N3;
  int<lower=1,upper=N> idx3[N3];
  int<lower=0,upper=1> t3[N3];

  int<lower=1> N4;
  int<lower=1,upper=N> idx4[N4];
  int<lower=0,upper=1> t4[N4];

  int<lower=1> N5;
  int<lower=1,upper=N> idx5[N5];
  int<lower=0,upper=1> t5[N5];
}
parameters {
  real<lower=0,upper=1> a11;
  real<lower=1-inv_logit(logit(a11)*2),upper=1> a12;
  real<lower=0,upper=1> a21;
  real<lower=1-inv_logit(logit(a21)*2),upper=1> a22;
  real<lower=0,upper=1> a31;
  real<lower=1-inv_logit(logit(a31)*2),upper=1> a32;
  real<lower=0,upper=1> a41;
  real<lower=1-inv_logit(logit(a41)*2),upper=1> a42;
  real<lower=0,upper=1> a51;
  real<lower=1-inv_logit(logit(a51)*2),upper=1> a52;
  real<lower=0,upper=1> prev;
}
transformed parameters {
  simplex[2] theta;
  theta[1] = 1-prev;
  theta[2] = prev;
}
model {
  real ps_base[N,2];

  a11~beta(47,1);  a12~beta(4,6);
  a21~beta(12,3);  a22~beta(8,4);
  a31~beta(12,3);  a32~beta(8,4);
  a41~beta(12,3);  a42~beta(4,8);
  a51~beta(12,3);  a52~beta(8,4);
  prev~beta(2,6);

  // Empezamos con MAT (siempre observado) + log(theta)
  for(n in 1:N){
    ps_base[n,1] = log(theta[1]) + binomial_lpmf(t1[n]|1, 1-a11);
    ps_base[n,2] = log(theta[2]) + binomial_lpmf(t1[n]|1, a12);
  }

  // Sumamos ELISA solo donde hay dato
  for(i in 1:N2){
    ps_base[idx2[i],1] += binomial_lpmf(t2[i]|1, 1-a21);
    ps_base[idx2[i],2] += binomial_lpmf(t2[i]|1, a22);
  }
  // Sumamos LifeAssay (siempre observado, pero por consistencia usamos N3=N)
  for(i in 1:N3){
    ps_base[idx3[i],1] += binomial_lpmf(t3[i]|1, 1-a31);
    ps_base[idx3[i],2] += binomial_lpmf(t3[i]|1, a32);
  }
  // Sumamos SD Bioline solo donde hay dato
  for(i in 1:N4){
    ps_base[idx4[i],1] += binomial_lpmf(t4[i]|1, 1-a41);
    ps_base[idx4[i],2] += binomial_lpmf(t4[i]|1, a42);
  }
  // Sumamos LeptoCheck solo donde hay dato
  for(i in 1:N5){
    ps_base[idx5[i],1] += binomial_lpmf(t5[i]|1, 1-a51);
    ps_base[idx5[i],2] += binomial_lpmf(t5[i]|1, a52);
  }

  for(n in 1:N){
    target += log_sum_exp(ps_base[n,1], ps_base[n,2]);
  }
}
generated quantities {
  real Se_MAT = a12;       real Sp_MAT = a11;
  real Se_ELISA = a22;     real Sp_ELISA = a21;
  real Se_LifeAssay = a32; real Sp_LifeAssay = a31;
  real Se_SDBioline = a42; real Sp_SDBioline = a41;
  real Se_LeptoCheck = a52; real Sp_LeptoCheck = a51;
}
"
cat("Modelo con manejo de NA escrito\n")

## Revisión de sintaxis

cat("Compilando modelo con manejo de NA... espera 2-3 minutos\n")
m_keddie_na <- stan_model(model_code = modelo_keddie_na,
                          model_name = "BLCA_Keddie_NA")
cat("Compilado OK\n")

##Prueba


cat("Prueba corta con manejo de NA: 1 cadena, 1000 iteraciones\n")
fit_na_prueba <- sampling(
  m_keddie_na, data = datos_stan_v2,
  chains = 1, iter = 1000, warmup = 500,
  seed = 2024
)

cat("\nR-hat:\n")
print(round(summary(fit_na_prueba,
                    pars=c("Se_MAT","Sp_MAT","Se_ELISA","Sp_ELISA",
                           "Se_LifeAssay","Sp_LifeAssay",
                           "Se_SDBioline","Sp_SDBioline",
                           "Se_LeptoCheck","Sp_LeptoCheck","prev"))$summary[,"Rhat"], 4))

cat("\nEstimados (mediana):\n")
print(round(summary(fit_na_prueba,
                    pars=c("Se_MAT","Sp_MAT","Se_ELISA","Sp_ELISA",
                           "Se_LifeAssay","Sp_LifeAssay",
                           "Se_SDBioline","Sp_SDBioline",
                           "Se_LeptoCheck","Sp_LeptoCheck","prev"))$summary[,"50%"], 4))



## Aumentando cadenas

cat("Corrida FINAL - fase aguda: 3 cadenas, 20000 iteraciones\n")
cat("Esto puede tardar 10-20 minutos\n")

set.seed(2024)
fit_aguda_final <- sampling(
  m_keddie_na, data = datos_stan_v2,
  chains = 3, iter = 20000, warmup = 5000,
  seed = 2024,
  control = list(adapt_delta = 0.95)
)

cat("\n========== RESULTADOS FINALES - FASE AGUDA ==========\n")

pars_finales <- c("Se_MAT","Sp_MAT","Se_ELISA","Sp_ELISA",
                  "Se_LifeAssay","Sp_LifeAssay",
                  "Se_SDBioline","Sp_SDBioline",
                  "Se_LeptoCheck","Sp_LeptoCheck","prev")

resumen_final <- summary(fit_aguda_final, pars = pars_finales,
                         probs = c(0.025, 0.5, 0.975))$summary

cat("\nR-hat maximo:", round(max(resumen_final[,"Rhat"]), 4), "\n\n")

print(round(resumen_final[,c("2.5%","50%","97.5%","Rhat")], 4))

saveRDS(fit_aguda_final, "output/fit_stan_aguda_FINAL.rds")
cat("\n✓ Guardado en output/fit_stan_aguda_FINAL.rds\n")

## FASE CONVALECIENTE

# Indices de pacientes CON dato observado - fase convaleciente
idx_elisac <- which(!is.na(df$elisac))
idx_lifec  <- which(!is.na(df$lifec))
idx_sdc    <- which(!is.na(df$sdc))
idx_lepc   <- which(!is.na(df$lepc))

datos_stan_conv <- list(
  N = nrow(df),
  t1 = as.integer(df$MAT),
  
  N2 = length(idx_elisac),
  idx2 = idx_elisac,
  t2 = as.integer(df$elisac[idx_elisac]),
  
  N3 = length(idx_lifec),
  idx3 = idx_lifec,
  t3 = as.integer(df$lifec[idx_lifec]),
  
  N4 = length(idx_sdc),
  idx4 = idx_sdc,
  t4 = as.integer(df$sdc[idx_sdc]),
  
  N5 = length(idx_lepc),
  idx5 = idx_lepc,
  t5 = as.integer(df$lepc[idx_lepc])
)

cat("N total:", datos_stan_conv$N, "\n")
cat("N2 (ELISA conv observados):", datos_stan_conv$N2, "\n")
cat("N3 (LifeAssay conv observados):", datos_stan_conv$N3, "\n")
cat("N4 (SD Bioline conv observados):", datos_stan_conv$N4, "\n")
cat("N5 (LeptoCheck conv observados):", datos_stan_conv$N5, "\n")



cat("Corrida FINAL - fase convaleciente: 3 cadenas, 20000 iteraciones\n")
cat("Esto puede tardar 10-20 minutos\n")

set.seed(2024)
fit_conv_final <- sampling(
  m_keddie_na, data = datos_stan_conv,
  chains = 3, iter = 20000, warmup = 5000,
  seed = 2024,
  control = list(adapt_delta = 0.95)
)

cat("\n========== RESULTADOS FINALES - FASE CONVALECIENTE ==========\n")

resumen_conv <- summary(fit_conv_final, pars = pars_finales,
                        probs = c(0.025, 0.5, 0.975))$summary

cat("\nR-hat maximo:", round(max(resumen_conv[,"Rhat"]), 4), "\n\n")

print(round(resumen_conv[,c("2.5%","50%","97.5%","Rhat")], 4))

saveRDS(fit_conv_final, "output/fit_stan_conv_FINAL.rds")
cat("\n✓ Guardado en output/fit_stan_conv_FINAL.rds\n")




###Modelo de dependencia


modelo_dep_na <- "
data {
  int<lower=1> N;
  int<lower=0,upper=1> t1[N];

  int<lower=1> N2;
  int<lower=1,upper=N> idx2[N2];
  int<lower=0,upper=1> t2[N2];

  int<lower=1> N3;
  int<lower=1,upper=N> idx3[N3];
  int<lower=0,upper=1> t3[N3];

  int<lower=1> N4;
  int<lower=1,upper=N> idx4[N4];
  int<lower=0,upper=1> t4[N4];

  int<lower=1> N5;
  int<lower=1,upper=N> idx5[N5];
  int<lower=0,upper=1> t5[N5];
}
parameters {
  real<lower=0,upper=1> prev;
  vector[N] RE;
  real<lower=0,upper=5> bpos;

  real<lower=0,upper=1> a11;
  real<lower=1-inv_logit(logit(a11)*2),upper=1> a12;
  real<lower=0,upper=1> a21;
  real<lower=1-inv_logit(logit(a21)-bpos*2),upper=1> a22;
  real<lower=0,upper=1> a31;
  real<lower=1-inv_logit(logit(a31)-bpos*2),upper=1> a32;
  real<lower=0,upper=1> a41;
  real<lower=1-inv_logit(logit(a41)*2),upper=1> a42;
  real<lower=0,upper=1> a51;
  real<lower=1-inv_logit(logit(a51)-bpos*2),upper=1> a52;
}
transformed parameters {
  simplex[2] theta;
  theta[1] = 1-prev;
  theta[2] = prev;
}
model {
  real ps_base[N,2];

  a11~beta(47,1);  a12~beta(4,6);
  a21~beta(12,3);  a22~beta(8,4);
  a31~beta(12,3);  a32~beta(8,4);
  a41~beta(12,3);  a42~beta(4,8);
  a51~beta(12,3);  a52~beta(8,4);
  prev~beta(2,6);
  RE~normal(0,1);
  bpos~gamma(1,1);

  for(n in 1:N){
    ps_base[n,1] = log(theta[1]) + binomial_lpmf(t1[n]|1, 1-a11);
    ps_base[n,2] = log(theta[2]) + binomial_lpmf(t1[n]|1, a12);
  }
  for(i in 1:N2){
    int n = idx2[i];
    ps_base[n,1] += binomial_lpmf(t2[i]|1, 1-a21);
    ps_base[n,2] += binomial_lpmf(t2[i]|1, inv_logit(logit(a22)+bpos*RE[n]));
  }
  for(i in 1:N3){
    int n = idx3[i];
    ps_base[n,1] += binomial_lpmf(t3[i]|1, 1-a31);
    ps_base[n,2] += binomial_lpmf(t3[i]|1, inv_logit(logit(a32)+bpos*RE[n]));
  }
  for(i in 1:N4){
    int n = idx4[i];
    ps_base[n,1] += binomial_lpmf(t4[i]|1, 1-a41);
    ps_base[n,2] += binomial_lpmf(t4[i]|1, a42);
  }
  for(i in 1:N5){
    int n = idx5[i];
    ps_base[n,1] += binomial_lpmf(t5[i]|1, 1-a51);
    ps_base[n,2] += binomial_lpmf(t5[i]|1, inv_logit(logit(a52)+bpos*RE[n]));
  }

  for(n in 1:N){
    target += log_sum_exp(ps_base[n,1], ps_base[n,2]);
  }
}
generated quantities {
  real Se_MAT = a12;        real Sp_MAT = a11;
  real Se_ELISA = a22;      real Sp_ELISA = a21;
  real Se_LifeAssay = a32;  real Sp_LifeAssay = a31;
  real Se_SDBioline = a42;  real Sp_SDBioline = a41;
  real Se_LeptoCheck = a52; real Sp_LeptoCheck = a51;
  real bpos_out = bpos;
}
"
cat("Modelo de dependencia (CDP) escrito\n")



cat("Compilando modelo de dependencia... espera 2-4 minutos\n")
m_dep_na <- stan_model(model_code = modelo_dep_na,
                       model_name = "BLCA_Dependencia_CDP")
cat("Compilado OK\n")



##Prueba usando datos anteriores para independence

cat("Prueba corta - modelo dependencia: 1 cadena, 1000 iteraciones\n")
set.seed(2024)
fit_dep_prueba <- sampling(
  m_dep_na, data = datos_stan_v2,
  chains = 1, iter = 1000, warmup = 500,
  seed = 2024,
  control = list(adapt_delta = 0.95)
)

pars_dep <- c("Se_MAT","Sp_MAT","Se_ELISA","Sp_ELISA",
              "Se_LifeAssay","Sp_LifeAssay",
              "Se_SDBioline","Sp_SDBioline",
              "Se_LeptoCheck","Sp_LeptoCheck",
              "prev","bpos_out")

cat("\nR-hat:\n")
print(round(summary(fit_dep_prueba, pars=pars_dep)$summary[,"Rhat"], 4))

cat("\nEstimados (mediana):\n")
print(round(summary(fit_dep_prueba, pars=pars_dep)$summary[,"50%"], 4))


cat("Corrida FINAL - DEPENDENCIA - fase aguda: 3 cadenas, 20000 iteraciones\n")
cat("Esto puede tardar 15-25 minutos (modelo mas complejo)\n")

set.seed(2024)
fit_dep_aguda_final <- sampling(
  m_dep_na, data = datos_stan_v2,
  chains = 3, iter = 20000, warmup = 5000,
  seed = 2024,
  control = list(adapt_delta = 0.95)
)

cat("\n========== RESULTADOS FINALES - DEPENDENCIA - FASE AGUDA ==========\n")

resumen_dep_aguda <- summary(fit_dep_aguda_final, pars = pars_dep,
                             probs = c(0.025, 0.5, 0.975))$summary

cat("\nR-hat maximo:", round(max(resumen_dep_aguda[,"Rhat"]), 4), "\n\n")
print(round(resumen_dep_aguda[,c("2.5%","50%","97.5%","Rhat")], 4))

saveRDS(fit_dep_aguda_final, "output/fit_stan_dependencia_aguda_FINAL.rds")
cat("\n✓ Guardado\n")


###Final

cat("Corrida FINAL - DEPENDENCIA - fase aguda: 3 cadenas, 20000 iteraciones\n")
cat("Esto puede tardar 15-25 minutos (modelo mas complejo)\n")

set.seed(2024)
fit_dep_aguda_final <- sampling(
  m_dep_na, data = datos_stan_v2,
  chains = 3, iter = 20000, warmup = 5000,
  seed = 2024,
  control = list(adapt_delta = 0.95)
)

cat("\n========== RESULTADOS FINALES - DEPENDENCIA - FASE AGUDA ==========\n")

resumen_dep_aguda <- summary(fit_dep_aguda_final, pars = pars_dep,
                             probs = c(0.025, 0.5, 0.975))$summary

cat("\nR-hat maximo:", round(max(resumen_dep_aguda[,"Rhat"]), 4), "\n\n")
print(round(resumen_dep_aguda[,c("2.5%","50%","97.5%","Rhat")], 4))

saveRDS(fit_dep_aguda_final, "output/fit_stan_dependencia_aguda_FINAL.rds")
cat("\n✓ Guardado\n")



cat("Corrida FINAL - DEPENDENCIA - fase convaleciente: 3 cadenas, 20000 iteraciones\n")
cat("Esto puede tardar 15-25 minutos\n")

set.seed(2024)
fit_dep_conv_final <- sampling(
  m_dep_na, data = datos_stan_conv,
  chains = 3, iter = 20000, warmup = 5000,
  seed = 2024,
  control = list(adapt_delta = 0.95)
)

cat("\n========== RESULTADOS FINALES - DEPENDENCIA - FASE CONVALECIENTE ==========\n")

resumen_dep_conv <- summary(fit_dep_conv_final, pars = pars_dep,
                            probs = c(0.025, 0.5, 0.975))$summary

cat("\nR-hat maximo:", round(max(resumen_dep_conv[,"Rhat"]), 4), "\n\n")
print(round(resumen_dep_conv[,c("2.5%","50%","97.5%","Rhat")], 4))



