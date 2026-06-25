# =============================================================
# Time-dependent sensitivity model (Bayesian latent class)
# Leptospirosis rapid tests vs MAT reference, Colombia
# Conditional independence; sensitivity as a function of
# days since symptom onset (logit-linear).
# =============================================================

library(readxl)
library(R2jags)

# ---- 1. Data ----
datos <- read_excel("data/DataBVF.xlsx", na = "NA")

# Keep patients with a recorded onset-to-sample interval
datos <- datos[!is.na(datos$dias_muestra1), ]

# Centre days (aids convergence and interpretability)
dias_c <- datos$dias_muestra1 - mean(datos$dias_muestra1)

# Test matrix: 4 index tests (acute) + MAT (combined)
T <- as.matrix(datos[, c("elisaa", "lifea", "sda", "lepa", "MAT")])

jags_data <- list(T = T, dias = dias_c, N = nrow(T), J = ncol(T))

# ---- 2. Model ----
model <- "
model {
  for (i in 1:N) {
    D[i] ~ dbern(pi)
    for (j in 1:J) {
      p[i,j] <- D[i] * Se[i,j] + (1 - D[i]) * (1 - Sp[j])
      T[i,j] ~ dbern(p[i,j])
    }
    # Index tests (j = 1..4): sensitivity varies with days
    for (j in 1:4) {
      logit(Se[i,j]) <- a[j] + b[j] * dias[i]
    }
    Se[i,5] <- Se5          # MAT: fixed sensitivity
  }

  # Priors
  pi ~ dbeta(2, 6)
  for (j in 1:4) {
    a[j] ~ dnorm(0, 0.25)   # intercept (logit scale)
    b[j] ~ dnorm(0, 1)      # slope per day
  }
  Se5 ~ dbeta(4, 6)         # MAT sensitivity

  Sp[1] ~ dbeta(12, 3)      # ELISA
  Sp[2] ~ dbeta(12, 3)      # LifeAssay
  Sp[3] ~ dbeta(12, 3)      # SD Bioline
  Sp[4] ~ dbeta(12, 3)      # Leptocheck-WB
  Sp[5] ~ dbeta(47, 1)      # MAT
}
"
writeLines(model, "blca_time_model.txt")

# ---- 3. Fit ----
params <- c("a", "b", "Sp", "Se5", "pi")

set.seed(2025)
fit <- jags(
  data               = jags_data,
  parameters.to.save = params,
  model.file         = "blca_time_model.txt",
  n.chains           = 3,
  n.iter             = 20000,
  n.burnin           = 5000,
  n.thin             = 5
)

# ---- 4. Output ----
print(fit)
cat("Max R-hat:", max(fit$BUGSoutput$summary[, "Rhat"]), "\n")

# Predicted sensitivity (%) at selected days
inv_logit <- function(x) 1 / (1 + exp(-x))
a_med <- fit$BUGSoutput$summary[paste0("a[", 1:4, "]"), "50%"]
b_med <- fit$BUGSoutput$summary[paste0("b[", 1:4, "]"), "50%"]
days  <- c(1, 3, 5, 7, 9, 12, 15)
days_c <- days - mean(datos$dias_muestra1)

pred <- sapply(days_c, function(d) round(100 * inv_logit(a_med + b_med * d), 1))
dimnames(pred) <- list(c("ELISA", "LifeAssay", "SD Bioline", "Leptocheck"),
                       paste0("day", days))
print(pred)
