# -------------------------------------------------------------------
# SCRIPT PRINCIPAL — FATOR DINÂMICO COMUM (FFBS EXTERNO)
# -------------------------------------------------------------------

# --- PASSO 0: CARREGAR PACOTES ---
if (!require(nimble)) install.packages("nimble")
if (!require(coda)) install.packages("coda")
if (!require(ggplot2)) install.packages("ggplot2")

library(nimble)
library(coda)
library(ggplot2)

# -------------------------------------------------------------------
# IMPORTANTE: NÃO USAR attach()
# -------------------------------------------------------------------

cat("--- Carregando dados e iniciais ---\n")

# Dados reais
source("_dataCaseStudy.r")      # Deve preencher um objeto, ex: dataCaseStudy
# NÃO usar attach(data)

# Dados simulados + iniciais para o fator comum
source("ValoresIniciais_FatorComum.R")

# FFBS customizado
source("ffbs_sampler_FatorComum.R")

# -------------------------------------------------------------------
# PASSO 1 — DEFINIÇÃO DO MODELO NIMBLE
# -------------------------------------------------------------------

cat("\n--- PASSO 1: Construindo modelo NIMBLE ---\n")

code_ffbs_externo <- nimbleCode({
  
  # -------------------------------------------------
  # Priors para beta
  # -------------------------------------------------
  for (j in 1:p) {
    beta[j] ~ dnorm(mu_beta[j], sd = 5)
  }
  
  # -------------------------------------------------
  # Priors para gamma (stick-breaking dependente)
  # -------------------------------------------------
  gamma[1] ~ dunif(a_unif[1], b_unif[1])
  for (j in 2:K) {
    gamma[j] ~ dunif(
      min = a_unif[j] * (1 - sum(gamma[1:(j - 1)])),
      max = b_unif[j] * (1 - sum(gamma[1:(j - 1)]))
    )
  }
  
  # epsilon[i] = 1 - soma_h[i,·] * gamma[·]
  for (i in 1:n_regions) {
    epsilon[i] <- 1 - sum(h[i, 1:K] * gamma[1:K])
  }
  
  # -------------------------------------------------
  # Determinação de g_it (componente espacial estática)
  # -------------------------------------------------
  for (i in 1:n_regions) {
    for (t in 1:n_times) {
      log_g_it[i, t] <- inprod(beta[1:p], x[i, t, 1:p])
      g_it[i, t] <- E[i, t] * epsilon[i] * exp(log_g_it[i, t])
    }
  }
  
  # -------------------------------------------------
  # Prior fictício para lambda (substituído pelo FFBS)
  # -------------------------------------------------
  for (t in 1:n_times) {
    lambda[t] ~ dgamma(1, 1)
  }
  
  # -------------------------------------------------
  # Verossimilhança
  # -------------------------------------------------
  for (i in 1:n_regions) {
    for (t in 1:n_times) {
      Y[i, t] ~ dpois(lambda[t] * g_it[i, t])
    }
  }
  
})


# -------------------------------------------------------------------
# PASSO 2 — CONFIGURAÇÃO DO MCMC
# -------------------------------------------------------------------

cat("\n--- PASSO 2: Configurando MCMC ---\n")

model_ext <- nimbleModel(code_ffbs_externo,
                         constants = constants_nimble,
                         data      = data_nimble,
                         inits     = inits_list_nimble[[1]],
                         check     = FALSE)

Cmodel_ext <- compileNimble(model_ext)

conf_ext <- configureMCMC(model = Cmodel_ext,
                          monitors = c("lambda", "beta", "gamma"))

cat("Substituindo amostrador de lambda pelo FFBS...\n")

conf_ext$removeSamplers("lambda")
conf_ext$addSampler(target = "lambda",
                    type = "ffbs_sampler_FatorComum",
                    control = list(
                      n_regions = constants_nimble$n_regions,
                      n_times   = constants_nimble$n_times,
                      p         = constants_nimble$p,
                      w         = constants_nimble$w,
                      a0        = constants_nimble$a0,
                      b0        = constants_nimble$b0
                    ))

conf_ext$printSamplers()

Rmcmc_ext <- buildMCMC(conf_ext)
Cmcmc_ext <- compileNimble(Rmcmc_ext, project = Cmodel_ext)

# -------------------------------------------------------------------
# PASSO 3 — EXECUTAR MCMC
# -------------------------------------------------------------------

cat("\n--- Executando MCMC com 2 cadeias ---\n")

samples_ext <- runMCMC(Cmcmc_ext,
                       niter            = 50000,
                       nburnin          = 0,
                       nchains          = 2,
                       thin             = 1,
                       inits            = inits_list_nimble,
                       samplesAsCodaMCMC = TRUE)

# -------------------------------------------------------------------
# PASSO 4 — GERAR TRACEPLOTS
# -------------------------------------------------------------------

cat("\n--- Gerando traceplots ---\n")

samples_df_ext <- do.call(
  rbind,
  lapply(1:length(samples_ext), function(i) {
    df <- as.data.frame(samples_ext[[i]])
    df$chain <- as.factor(i)
    df$iter  <- 1:nrow(df)
    df
  })
)

save_traceplot <- function(df, param_name, true_value) {
  p <- ggplot(df, aes_string(x = "iter", y = paste0("`", param_name, "`"), color = "chain")) +
    geom_line(alpha = 0.7) +
    geom_hline(yintercept = true_value, color = "blue", linetype = "dashed", size = 1) +
    labs(
      title    = paste("Traceplot de", param_name),
      subtitle = paste("Valor verdadeiro =", round(true_value, 3)),
      x = "Iteração",
      y = "Valor"
    ) +
    theme_minimal() +
    scale_color_manual(values = c("black", "red"))
  
  ggsave(paste0("traceplot_", gsub("\\[|,|\\]", "_", param_name), ".png"),
         p, width = 8, height = 4)
}

# Betas
for(i in 1:p) {
  save_traceplot(samples_df_ext, paste0("beta[", i, "]"), beta_true[i])
}

# Gammas
for(i in 1:K) {
  save_traceplot(samples_df_ext, paste0("gamma[", i, "]"), gamma_true[i])
}

# Algumas lambdas
selected_lambdas <- c(1, 10, 23)
for(t in selected_lambdas) {
  save_traceplot(samples_df_ext, paste0("lambda[", t, "]"), lambda_true[t])
}

cat("\n--- Análise FFBS concluída. ---\n")
