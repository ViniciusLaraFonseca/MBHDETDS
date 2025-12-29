###############################################################################
# run_sensitivity_parallel_w07_optimized.R
# Versão otimizada (mesmos outputs) - limita threads, parLapplyLB, export explícito
###############################################################################

# ---------------------------
# 0. Pacotes (Mestre)
# ---------------------------
pkgs <- c("parallel","nimble","coda","dplyr","stringr","tibble","ggplot2","purrr","readr")
for(p in pkgs) if(!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(parallel); library(nimble); library(coda)
library(dplyr); library(stringr); library(tibble); library(ggplot2); library(purrr); library(readr)

# Limitar threads do BLAS/OMP (evita oversubscription quando cada worker roda C/C++)
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")
if(requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
}

# ---------------------------
# 1. Carregar Dados e Funções (Mestre)
# ---------------------------
cat("--- PREPARAÇÃO DO AMBIENTE MESTRE ---\n")
print(Sys.time())
print("Checkpoint: source _dataCaseStudy.r")
source("_dataCaseStudy.r")
if(!file.exists("ffbs_sampler_Lambda_it.R")) stop("Sampler FFBS não encontrado.")

# A. Carregar T=100
print("Checkpoint: carregar ValoresIniciais_T100")
source("ValoresIniciais_Lambda_it_T100_w07.R")
data_nimble_T100       <- data_nimble
constants_nimble_T100  <- constants_nimble
inits_list_nimble_T100 <- inits_list_nimble
lambda_true_T100       <- lambda_true
rm(data_nimble, constants_nimble, inits_list_nimble, lambda_true) # Limpar genéricos

# B. Carregar T=23
print("Checkpoint: carregar ValoresIniciais_T23")
source("ValoresIniciais_Lambda_it_T23_w07.R")
# (Este script cria objetos _T23 automaticamente)

cat("Datasets carregados: T=100 e T=23 prontos na memória.\n")

# Segurança: garantir que beta_true / gamma_true existam (senão NULL)
if(!exists("beta_true")) { warning("beta_true não encontrado — definindo NULL"); beta_true <- NULL }
if(!exists("gamma_true")) { warning("gamma_true não encontrado — definindo NULL"); gamma_true <- NULL }

# Source do sampler localmente no mestre (validação)
source("ffbs_sampler_Lambda_it.R")

# --- Nimble Code (mantido igual) ---
code_ffbs_externo_it <- nimbleCode({
  for (j in 1:p) {
    beta[j] ~ dnorm(mu_beta[j], sd = 5)
  }
  gamma[1] ~ dunif(a_unif[1], b_unif[1])
  for (j in 2:K) {
    gamma[j] ~ dunif(min = a_unif[j] * (1 - sum(gamma[1:(j - 1)])),
                     max = b_unif[j] * (1 - sum(gamma[1:(j - 1)])))
  }
  for (i in 1:n_regions) {
    epsilon[i] <- 1 - sum(h[i, 1:K] * gamma[1:K])
  }
  for (i in 1:n_regions) {
    for (t in 1:n_times) {
      log_g_it[i, t] <- inprod(beta[1:p], x[i, t, 1:p])
      g_it[i, t] <- E[i, t] * epsilon[i] * exp(log_g_it[i, t])
    }
  }
  for (i in 1:n_regions) {
    for (t in 1:n_times) {
      lambda[i, t] ~ dgamma(1, 1) # Será substituído pelo FFBS
    }
  }
  for (i in 1:n_regions) {
    for (t in 1:n_times) {
      mu[i, t] <- lambda[i, t] * g_it[i, t]
      Y[i, t] ~ dpois(mu[i, t])
    }
  }
})

# ---------------------------
# 2. Funções Utilitárias (idênticas ao seu pipeline)
# ---------------------------
extract_lambda_summary <- function(samples_coda, regions_vec) {
  varn <- colnames(as.matrix(samples_coda[[1]]))
  lambda_names <- grep("^lambda\\[", varn, value = TRUE)
  if(length(lambda_names) == 0) return(NULL)
  parse_it <- function(name) {
    m <- str_match(name, "lambda\\[\\s*([0-9]+)\\s*,\\s*([0-9]+)\\s*\\]")
    c(as.integer(m[2]), as.integer(m[3]))
  }
  parsed <- t(sapply(lambda_names, parse_it))
  col_regions <- parsed[,1]; col_times <- parsed[,2]
  sel_idx <- which(col_regions %in% regions_vec)
  if(length(sel_idx) == 0) return(NULL)
  lambda_sel_names <- lambda_names[sel_idx]
  mat_all <- do.call(rbind, lapply(samples_coda, as.matrix))[, lambda_sel_names, drop = FALSE]
  out <- tibble(
    param = lambda_sel_names,
    region = col_regions[sel_idx],
    time   = col_times[sel_idx],
    Mean   = colMeans(mat_all),
    low95  = apply(mat_all, 2, quantile, 0.025, na.rm = TRUE),
    upp95  = apply(mat_all, 2, quantile, 0.975, na.rm = TRUE)
  ) %>% arrange(region, time)
  return(out)
}

matrix_to_tidy_df <- function(mat_true) {
  if(is.null(mat_true)) return(NULL)
  nr <- nrow(mat_true)
  nc <- ncol(mat_true)
  df <- data.frame(
    region   = rep(1:nr, times = nc),
    time     = rep(1:nc, each = nr),
    true_val = as.vector(mat_true) 
  )
  df$region <- as.integer(df$region)
  df$time   <- as.integer(df$time)
  return(df)
}

compute_diagnostics <- function(samples_mcmc, beta_true = NULL, gamma_true = NULL,
                                df_true_lambda = NULL) {
  mlist <- samples_mcmc
  gel <- tryCatch(gelman.diag(mlist, autoburnin = FALSE, multivariate = FALSE),
                  error = function(e) NULL)
  rhats <- NULL
  if(!is.null(gel))
    rhats <- tibble(
      param = names(gel$psrf[,1]),
      Rhat  = gel$psrf[,1]
    )
  ess_vals <- effectiveSize(mlist)
  ess_df <- tibble(
    param = names(ess_vals),
    ESS   = as.numeric(ess_vals)
  )
  mean_df <- tibble(
    param = colnames(as.matrix(do.call(rbind, lapply(mlist, as.matrix)))),
    Mean  = colMeans(do.call(rbind, lapply(mlist, as.matrix)))
  )
  diag_bg <- mean_df %>% 
    filter(grepl("beta\\[|gamma\\[", param)) %>%
    left_join(ess_df,  by = "param") %>%
    left_join(rhats,    by = "param")
  diag_bg$true_val <- NA_real_
  if(!is.null(beta_true)) {
    idx_beta <- grep("^beta\\[", diag_bg$param)
    for(k in idx_beta) {
      j <- as.integer(str_match(diag_bg$param[k], "beta\\[([0-9]+)\\]")[,2])
      diag_bg$true_val[k] <- beta_true[j]
    }
  }
  if(!is.null(gamma_true)) {
    idx_gamma <- grep("^gamma\\[", diag_bg$param)
    for(k in idx_gamma) {
      j <- as.integer(str_match(diag_bg$param[k], "gamma\\[([0-9]+)\\]")[,2])
      diag_bg$true_val[k] <- gamma_true[j]
    }
  }
  diag_bg <- diag_bg %>%
    mutate(
      bias = Mean - true_val,
      rmse = sqrt((Mean - true_val)^2)
    )
  diag_lambda <- NULL
  if(!is.null(df_true_lambda)) {
    df_lambda_names <- tibble(
      param  = names(ess_vals),
      ESS    = as.numeric(ess_vals)
    ) %>% filter(grepl("^lambda\\[", param))
    df_lambda_means <- mean_df %>%
      filter(grepl("^lambda\\[", param)) %>%
      select(param, Mean)
    df_lambda <- df_lambda_means %>%
      left_join(df_lambda_names, by = "param") %>%
      mutate(
        region = as.integer(str_match(param, "lambda\\[([0-9]+),")[,2]),
        time   = as.integer(str_match(param, ",([0-9]+)\\]")[,2])
      ) %>%
      left_join(df_true_lambda, by = c("region","time")) %>%
      mutate(
        bias = Mean - true_val,
        rmse = (Mean - true_val)^2
      )
    diag_lambda <- df_lambda
  }
  list(
    diag_beta_gamma = diag_bg,
    diag_lambda     = diag_lambda
  )
}

compute_EQM_lambda_robust <- function(lambda_summary, df_true_tidy, regions_vec) {
  if(is.null(lambda_summary) || is.null(df_true_tidy)) return(NA_real_)
  df_calc <- lambda_summary %>% 
    filter(region %in% regions_vec) %>%
    left_join(df_true_tidy, by = c("region", "time"))
  mean((df_calc$Mean - df_calc$true_val)^2, na.rm = TRUE)
}

save_lambda_panel_robust <- function(lambda_df, regions_to_plot, outfile, df_true_tidy = NULL) {
  if(is.null(lambda_df)) return(NULL)
  df <- lambda_df %>% filter(region %in% regions_to_plot)
  if(!is.null(df_true_tidy)) {
    df <- df %>% left_join(df_true_tidy, by = c("region", "time"))
  } else {
    df$true_val <- NA_real_
  }
  p <- ggplot(df, aes(x = time)) +
    geom_ribbon(aes(ymin = low95, ymax = upp95), alpha = 0.25) +
    geom_line(aes(y = Mean), size = 0.9, color = "black") +
    geom_line(aes(y = true_val), color = "red", size = 0.9, na.rm = TRUE) +
    facet_wrap(~ region, ncol = 3, nrow = 4, scales = "free_y") +
    labs(x = "Tempo", y = expression(lambda[t]), title = basename(outfile)) + 
    theme_bw(base_size = 12)
  ggsave(outfile, p, width = 12, height = 14, dpi = 150)
}

save_traceplots_bg <- function(mcmc_list, out_dir, beta_true = NULL, gamma_true = NULL, thin_factor = 10) {
  
  # Convertendo para mcmc.list com thinning
  mcmc_thin <- lapply(mcmc_list, function(ch) {
    window(ch, thin = thin_factor)
  })
  mcmc_thin <- mcmc.list(mcmc_thin)
  
  varn <- colnames(as.matrix(mcmc_thin[[1]]))
  params <- c(grep("^beta\\[", varn, value = TRUE),
              grep("^gamma\\[", varn, value = TRUE))
  
  if(length(params) == 0) return(NULL)
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  for(par in params){
    png(file.path(out_dir, paste0("traceplot_", gsub("\\[|\\]|,| ", "_", par), ".png")),
        width = 1000, height = 450)
    
    coda::traceplot(mcmc_thin[, par], main = par)
    
    # Linha do valor verdadeiro
    if(startsWith(par, "beta[") && !is.null(beta_true)) {
      idx <- as.integer(str_match(par, "beta\\[([0-9]+)\\]")[,2])
      if(!is.na(idx)) abline(h = beta_true[idx], col = "red", lwd = 2, lty = 2)
    }
    if(startsWith(par, "gamma[") && !is.null(gamma_true)) {
      idx <- as.integer(str_match(par, "gamma\\[([0-9]+)\\]")[,2])
      if(!is.na(idx)) abline(h = gamma_true[idx], col = "red", lwd = 2, lty = 2)
    }
    
    dev.off()
  }
}

# ---------------------------
# 3. WORKER FUNCTION (aceita datasets NULL -> usa objetos exportados no worker)
# ---------------------------
run_scenario_worker <- function(idx, config_df, 
                                # Argumentos opcionais (se NULL: usar objetos já exportados no worker)
                                data_T100 = NULL, const_T100 = NULL, inits_T100 = NULL, true_T100 = NULL,
                                data_T23  = NULL, const_T23  = NULL, inits_T23  = NULL, true_T23  = NULL,
                                # Globais
                                code_nimble = NULL, regions_roi = NULL, beta_true = NULL, gamma_true = NULL, wd_master = NULL,
                                niter = 50000, nchains = 2, nburnin = 0) {
  
  # Log curto
  start_time <- Sys.time()
  cat(sprintf("worker[%d] start %s\n", idx, format(start_time, "%Y-%m-%d %H:%M:%S")))
  
  # carregar libs (rápido) e garantir wd
  library(nimble); library(coda); library(dplyr); library(stringr); library(ggplot2); library(readr)
  if(!is.null(wd_master)) tryCatch(setwd(wd_master), error = function(e) warning("setwd falhou: ", e$message))
  
  # Se argumentos NULL, buscar objetos no ambiente do worker (exportados pelo mestre)
  if(is.null(data_T100))       data_T100       <- if(exists("data_nimble_T100")) data_nimble_T100 else NULL
  if(is.null(const_T100))      const_T100      <- if(exists("constants_nimble_T100")) constants_nimble_T100 else NULL
  if(is.null(inits_T100))      inits_T100      <- if(exists("inits_list_nimble_T100")) inits_list_nimble_T100 else NULL
  if(is.null(true_T100))       true_T100       <- if(exists("lambda_true_T100")) lambda_true_T100 else NULL
  
  if(is.null(data_T23))        data_T23        <- if(exists("data_nimble_T23")) data_nimble_T23 else NULL
  if(is.null(const_T23))       const_T23       <- if(exists("constants_nimble_T23")) constants_nimble_T23 else NULL
  if(is.null(inits_T23))       inits_T23       <- if(exists("inits_list_nimble_T23")) inits_list_nimble_T23 else NULL
  if(is.null(true_T23))        true_T23        <- if(exists("lambda_true_T23")) lambda_true_T23 else NULL
  
  if(is.null(code_nimble))     code_nimble     <- code_ffbs_externo_it
  if(is.null(regions_roi))     regions_roi     <- get0("regions_of_interest", ifnotfound = NULL)
  
  # assegurar sampler no worker
  if(!exists("ffbs_sampler_Lambda_it")) {
    source("ffbs_sampler_Lambda_it.R")
  }
  
  cfg <- config_df[idx, ]
  nome <- cfg$cenario; nT <- cfg$n_times; wval <- cfg$w
  out_root <- file.path(wd_master, "resultados_sensibilidade", nome)
  dir.create(file.path(out_root, "lambdas"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_root, "traceplots"), recursive = TRUE, showWarnings = FALSE)
  
  res <- tryCatch({
    if(nT == 100) {
      data_use  <- data_T100
      const_use <- const_T100
      inits_use <- inits_T100
      true_use  <- true_T100
    } else {
      data_use  <- data_T23
      const_use <- const_T23
      inits_use <- inits_T23
      true_use  <- true_T23
    }
    
    if(is.null(data_use) || is.null(const_use) || is.null(inits_use)) stop("Dados/inits/constantes do cenário não disponíveis no worker.")
    
    const_use$w <- wval
    df_true_tidy <- matrix_to_tidy_df(true_use)
    
    # construir e compilar modelo (necessário por cenário)
    model <- nimbleModel(code_nimble, constants = const_use, data = data_use, inits = inits_use[[1]], check = FALSE)
    Cmodel <- compileNimble(model)
    conf <- configureMCMC(model, monitors = c("lambda", "beta", "gamma"), enableWAIC = TRUE)
    conf$removeSamplers("lambda")
    conf$addSampler(target = "lambda",
                    type = "ffbs_sampler_Lambda_it",
                    control = list(n_regions = const_use$n_regions,
                                   n_times   = nT,
                                   p         = const_use$p,
                                   w         = wval,
                                   a0        = const_use$a0,
                                   b0        = const_use$b0))
    Rmcmc <- buildMCMC(conf)
    Cmcmc <- compileNimble(Rmcmc, project = Cmodel)
    
    cat(sprintf("[%s] Rodando MCMC (%d it) - inicio %s\n", nome, niter, format(Sys.time(), "%H:%M:%S")))
    samples_res <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin, nchains = nchains,
                           samplesAsCodaMCMC = TRUE, summary = TRUE, WAIC = TRUE)
    samples_mcmc <- samples_res$samples
    
    # diagnosticos e salvamento (idêntico ao original)
    diagn <- compute_diagnostics(samples_mcmc,
                                 beta_true   = beta_true,
                                 gamma_true  = gamma_true,
                                 df_true_lambda = df_true_tidy)
    mean_rhat <- mean(diagn$diag_beta_gamma$Rhat, na.rm = TRUE)
    mean_ess  <- mean(diagn$diag_beta_gamma$ESS,  na.rm = TRUE)
    mean_bias <- mean(diagn$diag_beta_gamma$bias, na.rm = TRUE)
    waic_val <- if(!is.null(samples_res$WAIC)) samples_res$WAIC$WAIC else NA_real_
    
    saveRDS(samples_mcmc, file = file.path(out_root, paste0("samples_", nome, ".rds")))
    
    lambda_summ <- extract_lambda_summary(samples_mcmc, regions_roi)
    if(!is.null(lambda_summ)) {
      save_lambda_panel_robust(lambda_summ, regions_roi,
                               file.path(out_root, "lambdas", paste0("painel_", nome, ".png")),
                               df_true_tidy)
    }
    save_traceplots_bg(samples_mcmc, file.path(out_root, "traceplots"), beta_true, gamma_true)
    eqm_val <- compute_EQM_lambda_robust(lambda_summ, df_true_tidy, regions_roi)
    
    cat(sprintf("[%s] Sucesso - fim %s\n", nome, format(Sys.time(), "%H:%M:%S")))
    
    tibble(
      cenario = nome,
      T = nT,
      w = wval,
      WAIC = waic_val,
      EQM = eqm_val,
      Rhat_mean = mean_rhat,
      ESS_mean  = mean_ess,
      Bias_mean = mean_bias,
      status = "OK"
    )
  }, error = function(e) {
    cat(sprintf("[%s] ERRO FATAL: %s\n", nome, e$message))
    tibble(cenario = nome, T = nT, w = wval, WAIC = NA, EQM = NA, status = paste("ERRO:", e$message))
  })
  
  end_time <- Sys.time()
  cat(sprintf("worker[%d] end %s (dur %s)\n", idx, format(end_time, "%Y-%m-%d %H:%M:%S"),
              round(difftime(end_time, start_time, units="secs"),1)))
  return(res)
}

# ---------------------------
# 4. Orquestração (EXPORT EXPLÍCITO DOS GRANDES OBJETOS ANTES DO parLapplyLB)
# ---------------------------
estudos <- tibble::tibble(
  cenario = c("C1_T100_w09","C2_T100_w07","C3_T23_w09","C4_T23_w07"),
  n_times = c(100, 100, 23, 23),
  w       = c(0.9, 0.7, 0.9, 0.7)
)
regions_of_interest <- c(1,8,15,19,22,31,34,40,46,55,65,75)

# --- DEFINIÇÃO DE ITERAÇÕES ATUALIZADA ---
niter_run <- 50000   # troque para 50000 quando quiser rodar full
nchains_run <- 2
nburnin_run <- 0

wd_master <- getwd()
n_cores <- min(parallel::detectCores() - 1, nrow(estudos))
if(n_cores < 1) n_cores <- 1

message("\n=== Iniciando Cluster com ", n_cores, " núcleos ===")
print(Sys.time())

# 0) Criar cluster
cl <- makeCluster(n_cores)

# Carregar pacotes nos workers (reduz overhead)
clusterEvalQ(cl, {
  library(nimble); library(coda); library(dplyr); library(stringr); library(tibble); library(ggplot2); library(readr)
  Sys.setenv(OMP_NUM_THREADS = "1"); Sys.setenv(MKL_NUM_THREADS = "1")
  if(requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1)
  NULL
})

# Exportar wd_master ANTES do setwd nos workers
clusterExport(cl, c("wd_master"), envir = environment())

# Sourcing do sampler nos workers (idempotente)
clusterEvalQ(cl, {
  setwd(wd_master)
  if(!file.exists("ffbs_sampler_Lambda_it.R")) stop("ffbs_sampler_Lambda_it.R não encontrado no worker")
  source("ffbs_sampler_Lambda_it.R")
  NULL
})

# 1) Exportar explicitamente os objetos grandes ANTES do parLapplyLB
big_objs <- c(
  "data_nimble_T100", "constants_nimble_T100", "inits_list_nimble_T100", "lambda_true_T100",
  "data_nimble_T23",  "constants_nimble_T23",  "inits_list_nimble_T23",  "lambda_true_T23"
)
cat("Exportando objetos grandes para workers (uma vez)...\n")
t0 <- Sys.time()
clusterExport(cl, varlist = big_objs, envir = environment())
cat("Término exportação: ", as.numeric(difftime(Sys.time(), t0, units = "secs")), " segundos\n")

# 2) Exportar funções/objetos leves
clusterExport(cl, varlist = c(
  "code_ffbs_externo_it",
  "compute_diagnostics","compute_EQM_lambda_robust","extract_lambda_summary",
  "matrix_to_tidy_df","save_lambda_panel_robust","save_traceplots_bg",
  "beta_true","gamma_true","run_scenario_worker"
), envir = environment())

# DEBUG breve para confirmar workers vivos (mostra alguns nomes)
cat(">>> workers list (head):\n")
print(clusterEvalQ(cl, head(ls(), 20)))

# 3) Chamada parLapplyLB - NÃO passamos os objetos grandes como parâmetros
tempo_inicio <- Sys.time()
res_list <- parLapplyLB(cl, X = seq_len(nrow(estudos)), fun = run_scenario_worker,
                        config_df = estudos,
                        # NOTA: NÃO passamos data_T100 / data_T23 aqui (os workers já têm via clusterExport)
                        code_nimble = code_ffbs_externo_it, regions_roi = regions_of_interest,
                        beta_true = beta_true, gamma_true = gamma_true, wd_master = wd_master,
                        niter = niter_run, nchains = nchains_run, nburnin = nburnin_run)
stopCluster(cl)

resumo_final <- bind_rows(res_list)
print(resumo_final)
dir.create("resultados_sensibilidade", showWarnings = FALSE)
write_csv(resumo_final, file.path("resultados_sensibilidade", "resumo_metricas_FINAL.csv"))
message("Tempo Total: ", round(difftime(Sys.time(), tempo_inicio, units="mins"), 2), " min.")
