###############################################################################
# postprocess_metrics_waic_lpml_clusters_FINAL.R
# VERSÃO OTIMIZADA: Mantém apenas arquivos essenciais (boxplots + tabelas consolidadas)
# CORREÇÃO: Cálculo correto do epsilon estimado a partir dos gammas
###############################################################################

library(coda)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)
library(stringr)

###############################################################################
# 1. Leitura MCMC
###############################################################################
read_mcmc_samples <- function(file, thin = 10) {
  obj <- readRDS(file)
  stopifnot(inherits(obj, "mcmc.list"))
  M <- do.call(rbind, lapply(obj, as.matrix))
  if (thin > 1) M <- M[seq(1, nrow(M), by = thin), ]
  M
}

###############################################################################
# 2. Loader do cenário - SIMPLIFICADO
###############################################################################
load_scenario_data <- function(scen_dir, main_dir = ".", w_generation) {
  
  cat("    Carregando dados para:", basename(scen_dir), "\n")
  cat("    w de geração:", w_generation, "\n")
  
  # PASSO 1: Carrega _dataCaseStudy.R do diretório principal
  old_wd <- getwd()
  setwd(main_dir)
  
  if (!file.exists("_dataCaseStudy.R")) {
    setwd(old_wd)
    stop("Arquivo _dataCaseStudy.R não encontrado em ", main_dir)
  }
  
  # Carrega _dataCaseStudy.R em um ambiente separado
  data_env <- new.env()
  source("_dataCaseStudy.R", local = data_env)
  
  # Extrai hAI
  if (exists("data", envir = data_env)) {
    hAI <- data_env$data$hAI
  } else {
    hAI <- data_env$hAI
  }
  
  if (is.null(hAI)) {
    setwd(old_wd)
    stop("Variável hAI não encontrada em _dataCaseStudy.R")
  }
  
  # PASSO 2: Encontra o arquivo de valores iniciais correto
  dir_name <- basename(scen_dir)
  
  # Identificar T do nome do diretório
  t_value <- ifelse(grepl("T23|t23", dir_name), "T23", "T100")
  
  cat("      T:", t_value, "\n")
  
  # Buscar por arquivo de valores iniciais
  init_pattern <- paste0("ValoresIniciais_Lambda_it_", t_value, ".*\\.R")
  init_files <- list.files(pattern = init_pattern, full.names = TRUE)
  
  if (length(init_files) == 0) {
    setwd(old_wd)
    stop("Nenhum arquivo ValoresIniciais_Lambda_it_*.R encontrado para T = ", t_value)
  }
  
  # Selecionar o arquivo que corresponde ao w de geração
  init_file <- NULL
  for (f in init_files) {
    if (w_generation == "w07" && grepl("w07|w0\\.7", f, ignore.case = TRUE)) {
      init_file <- f
      break
    } else if (w_generation == "w09" && grepl("w09|w0\\.9", f, ignore.case = TRUE)) {
      init_file <- f
      break
    }
  }
  
  # Se não encontrou, pegar o primeiro
  if (is.null(init_file)) {
    init_file <- init_files[1]
    cat("      AVISO: Usando arquivo de valores iniciais padrão\n")
  }
  
  cat("      Arquivo de valores iniciais:", basename(init_file), "\n")
  
  # Carrega os valores iniciais
  init_env <- new.env()
  source(init_file, local = init_env)
  
  setwd(old_wd)
  
  # Verifica e extrai variáveis necessárias
  if (!exists("x_true", envir = init_env)) {
    if (exists("x", envir = init_env)) {
      x <- init_env$x
    } else {
      stop("Variável x (ou x_true) não encontrada")
    }
  } else {
    x <- init_env$x_true
  }
  
  if (!exists("E_true", envir = init_env)) {
    if (exists("E", envir = init_env)) {
      E <- init_env$E
    } else {
      stop("Variável E (ou E_true) não encontrada")
    }
  } else {
    E <- init_env$E_true
  }
  
  if (!exists("Y_ini", envir = init_env)) {
    if (exists("y", envir = init_env)) {
      y <- init_env$y
    } else {
      stop("Variável y (ou Y_ini) não encontrada")
    }
  } else {
    y <- init_env$Y_ini
  }
  
  required_vars <- c("epsilon_true", "lambda_true", "beta_true", "gamma_true")
  missing_vars <- setdiff(required_vars, ls(envir = init_env))
  if (length(missing_vars) > 0) {
    stop("Variáveis faltando: ", paste(missing_vars, collapse = ", "))
  }
  
  epsilon_true <- init_env$epsilon_true
  lambda_true <- init_env$lambda_true
  beta_true <- init_env$beta_true
  gamma_true <- init_env$gamma_true
  
  # Converte lambda_true para matriz se necessário
  if (!is.matrix(lambda_true)) {
    if (is.array(lambda_true) && length(dim(lambda_true)) == 3) {
      lambda_true_mat <- lambda_true[1,,]
    } else {
      lambda_true_mat <- matrix(lambda_true, nrow = nrow(x), ncol = ncol(x))
    }
  } else {
    lambda_true_mat <- lambda_true
  }
  
  # Calcula cluster_id
  cluster_id <- rowSums(hAI)
  
  # Identificar w de inferência do nome do diretório
  w_inference <- ifelse(grepl("w07|w0\\.7", dir_name, ignore.case = TRUE), "w07", "w09")
  
  # Retorna lista com todos os dados
  list(
    x = x,
    E = E,
    y = y,
    epsilon_true = epsilon_true,
    lambda_true = lambda_true_mat,
    beta_true = beta_true,
    gamma_true = gamma_true,
    hAI = hAI,
    cluster_id = cluster_id,
    w_inference = w_inference,
    w_generation = w_generation,
    t_value = t_value
  )
}

###############################################################################
# 3. Extração de parâmetros - CORRIGIDA PARA INCLUIR GAMMA
###############################################################################
extract_params <- function(M, R, Tt) {
  beta_idx   <- grep("^beta\\[", colnames(M))
  lambda_idx <- grep("^lambda", colnames(M))
  gamma_idx  <- grep("^gamma\\[", colnames(M))
  
  lambda_arr <- array(
    M[, lambda_idx],
    dim = c(nrow(M), R, Tt)
  )
  
  list(
    beta   = M[, beta_idx, drop = FALSE],
    lambda = lambda_arr,
    gamma  = M[, gamma_idx, drop = FALSE]
  )
}

###############################################################################
# 4. Métricas básicas
###############################################################################
bias_eqm <- function(samples, true) {
  c(
    bias = mean(samples) - true,
    eqm  = mean((samples - true)^2)
  )
}

###############################################################################
# 5. Cálculo do epsilon estimado - NOVA FUNÇÃO
###############################################################################
compute_epsilon_estimated <- function(gamma_mat, hAI) {
  # gamma_mat: matriz S x K (amostras x número de clusters)
  # hAI: matriz R x K (regiões x clusters)
  
  S <- nrow(gamma_mat)  # número de amostras
  R <- nrow(hAI)        # número de regiões
  K <- ncol(hAI)        # número de clusters
  
  # Inicializar matriz de epsilon estimado
  epsilon_est <- matrix(NA, nrow = S, ncol = R)
  
  # Calcular epsilon estimado para cada amostra e região
  for (s in 1:S) {
    # epsilon_i = 1 - sum_{j=1}^{K} h_{ij} * gamma_j
    epsilon_est[s, ] <- 1 - rowSums(hAI * matrix(gamma_mat[s, ], nrow = R, ncol = K, byrow = TRUE))
  }
  
  return(epsilon_est)
}

###############################################################################
# 6. Lambda médio por cluster - CORRIGIDO
###############################################################################
lambda_cluster_samples <- function(lambda_arr, cluster_id) {
  
  S  <- dim(lambda_arr)[1]
  R  <- dim(lambda_arr)[2]
  Tt <- dim(lambda_arr)[3]
  
  map_dfr(sort(unique(cluster_id)), function(cl) {
    
    idx <- which(cluster_id == cl)
    
    tibble(
      cluster = cl,
      sample = sapply(1:S, function(s)
        mean(lambda_arr[s, idx, , drop = FALSE]))
    )
  })
}

###############################################################################
# 7. WAIC e LPML - CORREÇÃO PARA EVITAR -Inf
###############################################################################
compute_waic_lpml <- function(lambda_arr, beta_mat, x, E, epsilon_vec, y) {
  
  S  <- dim(lambda_arr)[1]
  R  <- dim(lambda_arr)[2]
  Tt <- dim(lambda_arr)[3]
  
  epsilon_mat <- matrix(epsilon_vec, nrow = R, ncol = Tt, byrow = TRUE)
  loglik <- array(NA_real_, dim = c(S, R, Tt))
  
  for (s in 1:S) {
    # Calcula x*beta de forma vetorizada
    xb <- array(0, dim = c(R, Tt))
    for (p in 1:ncol(beta_mat)) {
      xb <- xb + beta_mat[s, p] * x[,,p]
    }
    
    mu <- E * epsilon_mat * lambda_arr[s,,] * exp(xb)
    loglik[s,,] <- dpois(y, mu, log = TRUE)
  }
  
  # CORREÇÃO: Usar log-sum-exp com estabilidade numérica aprimorada
  # Calcula WAIC
  lppd <- 0
  pwaic <- 0
  
  for (r in 1:R) {
    for (t in 1:Tt) {
      log_lik_rt <- loglik[, r, t]
      valid_idx <- is.finite(log_lik_rt)
      
      if (sum(valid_idx) > 0) {
        # Log-sum-exp trick para lppd
        max_log <- max(log_lik_rt[valid_idx])
        lppd <- lppd + (max_log + log(mean(exp(log_lik_rt[valid_idx] - max_log))))
        
        # Calcular variância para pwaic
        if (sum(valid_idx) > 1) {
          pwaic <- pwaic + var(log_lik_rt[valid_idx])
        }
      }
    }
  }
  
  waic  <- -2 * (lppd - pwaic)
  
  # CORREÇÃO: Cálculo do LPML usando log-sum-exp para evitar underflow
  lpml <- 0
  
  for (r in 1:R) {
    for (t in 1:Tt) {
      log_lik_rt <- loglik[, r, t]
      valid_idx <- is.finite(log_lik_rt)
      
      if (sum(valid_idx) > 0) {
        # Calcular CPO usando log-sum-exp
        # CPO = 1 / mean(exp(-loglik))
        # log(CPO) = -log(mean(exp(-loglik)))
        
        neg_log_lik <- -log_lik_rt[valid_idx]
        max_neg <- max(neg_log_lik)
        
        # log(mean(exp(neg_log_lik))) = max_neg + log(mean(exp(neg_log_lik - max_neg)))
        log_mean_exp <- max_neg + log(mean(exp(neg_log_lik - max_neg)))
        
        # log(CPO) = -log_mean_exp
        log_cpo <- -log_mean_exp
        
        # Adicionar ao LPML (soma dos log(CPO))
        if (is.finite(log_cpo)) {
          lpml <- lpml + log_cpo
        }
      }
    }
  }
  
  list(WAIC = waic, LPML = lpml)
}

###############################################################################
# 8. Viés lambda_t (média espacial por tempo)
###############################################################################
lambda_t_bias_df <- function(lambda_arr, lambda_true) {
  
  S  <- dim(lambda_arr)[1]
  R  <- dim(lambda_arr)[2]
  Tt <- dim(lambda_arr)[3]
  
  # Garante que lambda_true seja matriz R x T
  if (!is.matrix(lambda_true)) {
    lambda_true <- matrix(lambda_true, nrow = R, ncol = Tt)
  }
  
  # Média verdadeira por tempo
  true_mean_t <- colMeans(lambda_true, na.rm = TRUE)
  
  # Calcula para cada amostra
  out <- vector("list", S)
  
  for (s in 1:S) {
    est_mean_t <- colMeans(lambda_arr[s,,], na.rm = TRUE)
    
    out[[s]] <- tibble(
      time = seq_len(Tt),
      bias = est_mean_t - true_mean_t
    )
  }
  
  bind_rows(out)
}

###############################################################################
# 9. BUSCA DIRETA DOS ARQUIVOS SAMPLES - BUSCA EXPLÍCITA NAS PASTAS w07 E w09
###############################################################################
find_sample_files <- function(root_dir = ".") {
  cat("\nBuscando arquivos samples_*.rds nas pastas de w de geração...\n")
  
  # Definir as pastas de w de geração
  w_generation_folders <- c(
    "resultados_sensibilidade_w07",
    "resultados_sensibilidade_w09"
  )
  
  # Verificar se as pastas existem
  existing_folders <- w_generation_folders[dir.exists(w_generation_folders)]
  
  if (length(existing_folders) == 0) {
    cat("AVISO: Nenhuma pasta de w de geração encontrada. Buscando recursivamente...\n")
    # Fallback: busca recursiva
    all_files <- list.files(root_dir, pattern = "^samples_.*\\.rds$", 
                            recursive = TRUE, full.names = TRUE)
    all_files <- all_files[!grepl("\\.git|\\.Rproj\\.user", all_files)]
    
    if (length(all_files) == 0) {
      cat("Nenhum arquivo samples_*.rds encontrado.\n")
      return(list())
    }
    
    # Tentar inferir w de geração do caminho do arquivo
    result <- list()
    for (f in all_files) {
      w_gen <- ifelse(grepl("sensibilidade_w07", f), "w07", 
                      ifelse(grepl("sensibilidade_w09", f), "w09", NA))
      dir_path <- dirname(f)
      result[[length(result) + 1]] <- list(
        dir = dir_path,
        files = f,
        w_generation = w_gen
      )
    }
    return(result)
  }
  
  cat("Pastas de w de geração encontradas:", paste(existing_folders, collapse = ", "), "\n")
  
  # Procurar em cada pasta de w de geração
  result <- list()
  
  for (w_folder in existing_folders) {
    # Determinar w de geração a partir do nome da pasta
    w_gen <- ifelse(grepl("w07", w_folder), "w07", "w09")
    
    cat("\n  Processando pasta:", w_folder, "(w de geração =", w_gen, ")\n")
    
    # Listar subdiretórios (cenários)
    subdirs <- list.dirs(w_folder, recursive = FALSE, full.names = TRUE)
    
    if (length(subdirs) == 0) {
      cat("    Nenhum subdiretório encontrado.\n")
      next
    }
    
    for (subdir in subdirs) {
      # Encontrar arquivos samples no subdiretório
      files <- list.files(subdir, pattern = "^samples_.*\\.rds$", full.names = TRUE)
      
      if (length(files) > 0) {
        cat("    Cenário:", basename(subdir), "-", length(files), "arquivo(s) encontrado(s)\n")
        
        # Adicionar à lista de resultados
        for (f in files) {
          result[[length(result) + 1]] <- list(
            dir = subdir,
            files = f,
            w_generation = w_gen
          )
        }
      }
    }
  }
  
  if (length(result) == 0) {
    cat("Nenhum arquivo samples_*.rds encontrado nas pastas de w de geração.\n")
    return(list())
  }
  
  cat("\nTotal de arquivos samples encontrados:", length(result), "\n")
  return(result)
}

###############################################################################
# 10. LOOP PRINCIPAL - MODIFICADO PARA PROCESSAR EXPLICITAMENTE w DE GERAÇÃO
###############################################################################
main_dir <- getwd()
cat("Diretório principal:", main_dir, "\n")

scenario_files <- find_sample_files(main_dir)

if (length(scenario_files) == 0) {
  stop("Nenhum arquivo samples_*.rds encontrado.")
}

results <- list()
waic_lpml_list <- list()
lambda_t_all <- list()  # Lista para armazenar todos os dados de lambda_t por cenário
individual_results <- list()  # Lista para armazenar resultados individuais

# Contadores
processed_count <- 0
error_count <- 0

for (i in seq_along(scenario_files)) {
  item <- scenario_files[[i]]
  scen_dir <- item$dir
  sfile <- item$files
  w_generation <- item$w_generation
  
  cat("\n=======================================\n")
  cat("Processando item", i, "de", length(scenario_files), "\n")
  cat("Diretório:", basename(scen_dir), "\n")
  cat("Arquivo:", basename(sfile), "\n")
  cat("w de geração:", w_generation, "\n")
  cat("=======================================\n")
  
  tryCatch({
    data_list <- load_scenario_data(scen_dir, main_dir, w_generation)
  }, error = function(e) {
    cat("  ERRO ao carregar dados do cenário:", e$message, "\n")
    error_count <<- error_count + 1
    return()
  })
  
  x <- data_list$x
  E <- data_list$E
  y <- data_list$y
  epsilon_true <- data_list$epsilon_true
  lambda_true <- data_list$lambda_true
  beta_true <- data_list$beta_true
  gamma_true <- data_list$gamma_true
  hAI <- data_list$hAI
  cluster_id <- data_list$cluster_id
  w_inference <- data_list$w_inference
  w_generation <- data_list$w_generation
  t_value <- data_list$t_value
  Tt <- ncol(x)  # Número de períodos de tempo
  
  cat("  Dimensões: R =", nrow(x), ", T =", Tt, "\n")
  cat("  Número de clusters:", length(unique(cluster_id)), "\n")
  cat("  Configuração: w (inferência) =", w_inference, 
      "| w (geração) =", w_generation, "| T =", t_value, "\n")
  
  cat("  --- Processando:", basename(sfile), "\n")
  
  M <- tryCatch({
    read_mcmc_samples(sfile)
  }, error = function(e) {
    cat("    ERRO ao ler arquivo:", e$message, "\n")
    error_count <<- error_count + 1
    return(NULL)
  })
  
  if (is.null(M)) next
  
  pars <- extract_params(M, nrow(x), Tt)
  
  # --- beta
  beta_tab <- map_dfr(
    1:ncol(pars$beta),
    ~ tibble(
      param = paste0("beta", .x),
      bias = bias_eqm(pars$beta[, .x], beta_true[.x])["bias"],
      eqm = bias_eqm(pars$beta[, .x], beta_true[.x])["eqm"]
    )
  )
  
  # --- epsilon por cluster - CORREÇÃO: calcular epsilon estimado a partir dos gammas
  # Calcular epsilon estimado para cada amostra e região
  epsilon_est <- compute_epsilon_estimated(pars$gamma, hAI)
  
  # Agrupar por cluster e calcular métricas
  eps_tab <- map_dfr(
    sort(unique(cluster_id)),
    function(cl) {
      idx <- which(cluster_id == cl)
      if (length(idx) == 0) return(NULL)
      
      # Calcular epsilon verdadeiro médio no cluster
      epsilon_true_cluster <- mean(epsilon_true[idx])
      
      # Calcular epsilon estimado médio no cluster para cada amostra
      epsilon_est_cluster <- rowMeans(epsilon_est[, idx, drop = FALSE])
      
      # Calcular métricas
      metrics <- bias_eqm(epsilon_est_cluster, epsilon_true_cluster)
      
      tibble(
        cluster = cl,
        bias = metrics["bias"],
        eqm = metrics["eqm"]
      )
    }
  )
  
  # --- lambda por cluster
  lam_samp <- lambda_cluster_samples(pars$lambda, cluster_id)
  
  # Calcular a média verdadeira por cluster de forma segura
  # Primeiro, criar um data.frame com lambda_true e cluster_id
  lambda_true_df <- as.data.frame(lambda_true)
  colnames(lambda_true_df) <- paste0("t", 1:ncol(lambda_true))
  lambda_true_df$region <- 1:nrow(lambda_true)
  lambda_true_df$cluster_id <- cluster_id
  
  # Calcular a média verdadeira por cluster
  true_means <- lambda_true_df %>%
    pivot_longer(cols = starts_with("t"), names_to = "time", values_to = "lambda") %>%
    group_by(cluster_id) %>%
    summarise(true_mean = mean(lambda))
  
  # Agora, juntar com lam_samp e calcular as métricas
  lam_tab <- lam_samp %>%
    left_join(true_means, by = c("cluster" = "cluster_id")) %>%
    group_by(cluster) %>%
    summarise(
      est_mean = mean(sample),
      bias = est_mean - first(true_mean),
      eqm = mean((sample - first(true_mean))^2),
      .groups = "drop"
    ) %>%
    select(cluster, bias, eqm)
  
  # --- WAIC / LPML
  crit <- tryCatch({
    compute_waic_lpml(pars$lambda, pars$beta, x, E, epsilon_true, y)
  }, error = function(e) {
    cat("    ERRO no cálculo WAIC/LPML:", e$message, "\n")
    list(WAIC = NA, LPML = NA)
  })
  
  waic_lpml_list[[sfile]] <- tibble(
    scenario = basename(scen_dir),
    file = basename(sfile),
    w_inference = w_inference,
    w_generation = w_generation,
    T_config = t_value,
    WAIC = crit$WAIC,
    LPML = crit$LPML
  )
  
  # --- Lambda por tempo
  tmp_lambda_t <- tryCatch({
    lambda_t_bias_df(pars$lambda, lambda_true)
  }, error = function(e) {
    cat("    ERRO no cálculo lambda_t:", e$message, "\n")
    tibble(time = integer(), bias = numeric())
  })
  
  if (nrow(tmp_lambda_t) > 0 && "bias" %in% names(tmp_lambda_t)) {
    tmp_lambda_t <- tmp_lambda_t %>%
      mutate(
        scenario = basename(scen_dir),
        file = basename(sfile),
        w_inference = w_inference,
        w_generation = w_generation,
        T_config = t_value
      )
    
    # Armazenar dados para boxplots com chave única
    key <- paste(basename(scen_dir), basename(sfile), w_generation, sep = "_")
    lambda_t_all[[key]] <- tmp_lambda_t
    cat("    lambda_t: OK (", nrow(tmp_lambda_t), " linhas)\n", sep = "")
  } else {
    cat("    AVISO: Dados de lambda_t vazios ou sem coluna 'bias'\n")
  }
  
  # Armazenar resultados individuais para este cenário
  individual_results[[sfile]] <- list(
    scenario = basename(scen_dir),
    file = basename(sfile),
    w_inference = w_inference,
    w_generation = w_generation,
    T_config = t_value,
    beta = beta_tab,
    epsilon = eps_tab,
    lambda = lam_tab,
    WAIC = crit$WAIC,
    LPML = crit$LPML
  )
  
  results[[sfile]] <- list(
    beta = beta_tab,
    epsilon = eps_tab,
    lambda = lam_tab
  )
  
  processed_count <- processed_count + 1
}

###############################################################################
# 11. OUTPUTS - VERSÃO OTIMIZADA (APENAS ARQUIVOS ESSENCIAIS)
###############################################################################
cat("\n\n=======================================\n")
cat("SALVANDO RESULTADOS (VERSÃO OTIMIZADA)\n")
cat("=======================================\n")

# Estatísticas de processamento
cat("Estatísticas de processamento:\n")
cat("- Total de itens a processar:", length(scenario_files), "\n")
cat("- Itens processados com sucesso:", processed_count, "\n")
cat("- Erros durante o processamento:", error_count, "\n")

# Criar diretórios para boxplots
boxplot_dir <- "boxplot_bias_lambda_t"
if (!dir.exists(boxplot_dir)) {
  dir.create(boxplot_dir)
  cat("✓ Diretório criado:", boxplot_dir, "\n")
}

# Criar subdiretórios para w de geração
boxplot_w07_dir <- file.path(boxplot_dir, "w07_generation")
boxplot_w09_dir <- file.path(boxplot_dir, "w09_generation")
if (!dir.exists(boxplot_w07_dir)) dir.create(boxplot_w07_dir, recursive = TRUE)
if (!dir.exists(boxplot_w09_dir)) dir.create(boxplot_w09_dir, recursive = TRUE)

# Criar diretórios para tabelas de métricas
tables_dir <- "tabelas_metricas"
if (!dir.exists(tables_dir)) {
  dir.create(tables_dir)
  cat("✓ Diretório criado:", tables_dir, "\n")
}

# Criar subdiretórios para w de geração
tables_w07_dir <- file.path(tables_dir, "w07_generation")
tables_w09_dir <- file.path(tables_dir, "w09_generation")
if (!dir.exists(tables_w07_dir)) dir.create(tables_w07_dir, recursive = TRUE)
if (!dir.exists(tables_w09_dir)) dir.create(tables_w09_dir, recursive = TRUE)

# Salvar WAIC/LPML consolidado
if (length(waic_lpml_list) > 0) {
  waic_lpml_df <- bind_rows(waic_lpml_list)
  write.csv(waic_lpml_df, "WAIC_LPML_summary.csv", row.names = FALSE)
  cat("✓ Arquivo WAIC_LPML_summary.csv salvo (", nrow(waic_lpml_df), " linhas)\n", sep = "")
  
  # Estatísticas resumidas por w de geração
  cat("\nResumo WAIC/LPML por w de geração:\n")
  summary_table <- waic_lpml_df %>%
    group_by(w_generation, T_config) %>%
    summarise(
      n = n(),
      WAIC_mean = mean(WAIC, na.rm = TRUE),
      WAIC_sd = sd(WAIC, na.rm = TRUE),
      LPML_mean = mean(LPML, na.rm = TRUE),
      LPML_sd = sd(LPML, na.rm = TRUE),
      .groups = "drop"
    )
  print(summary_table)
  
  # Salvar tabela de resumo
  write.csv(summary_table, "WAIC_LPML_resumo.csv", row.names = FALSE)
  cat("✓ Arquivo WAIC_LPML_resumo.csv salvo\n")
} else {
  cat("✗ AVISO: Nenhum dado WAIC/LPML para salvar\n")
}

# Criar boxplots individuais (deve gerar 8 no total)
if (length(lambda_t_all) > 0) {
  cat("\nCriando boxplots individuais para cada cenário...\n")
  
  # Contadores para w de geração
  count_w07 <- 0
  count_w09 <- 0
  
  # Gerar um boxplot para cada cenário
  for (i in seq_along(lambda_t_all)) {
    df <- lambda_t_all[[i]]
    scenario_name <- unique(df$scenario)
    file_name <- unique(df$file)
    w_gen <- unique(df$w_generation)
    w_inf <- unique(df$w_inference)
    t_val <- unique(df$T_config)
    
    if (length(scenario_name) > 0 && nrow(df) > 0) {
      # Determinar limites do eixo X baseado no tempo máximo
      max_time <- max(df$time, na.rm = TRUE)
      
      # Calcular limites do eixo Y baseado nos percentis 1% e 99%
      y_limits <- quantile(df$bias, c(0.01, 0.99), na.rm = TRUE)
      
      # Criar boxplot
      p <- ggplot(df, aes(x = factor(time), y = bias)) +
        geom_boxplot(outlier.shape = NA) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
        theme_bw() +
        labs(
          x = "Tempo (t)",
          y = expression("Viés de" ~ bar(lambda)[t]),
          title = paste("Distribuição do viés de λ por tempo"),
          subtitle = paste("Cenário:", scenario_name, 
                           "| w (geração) =", w_gen,
                           "| w (inferência) =", w_inf,
                           "| T =", t_val)
        ) +
        theme(
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 10)
        ) +
        scale_x_discrete(limits = factor(1:max_time)) +
        coord_cartesian(ylim = y_limits)
      
      # Criar nome do arquivo
      safe_name <- gsub("[^a-zA-Z0-9_]", "_", paste(scenario_name, w_gen, sep = "_"))
      filename <- paste0("boxplot_", safe_name, ".png")
      
      # Determinar diretório baseado no w de geração
      if (w_gen == "w07") {
        filepath <- file.path(boxplot_w07_dir, filename)
        count_w07 <- count_w07 + 1
      } else if (w_gen == "w09") {
        filepath <- file.path(boxplot_w09_dir, filename)
        count_w09 <- count_w09 + 1
      } else {
        filepath <- file.path(boxplot_dir, filename)
      }
      
      ggsave(filepath, plot = p, width = 14, height = 8, dpi = 300)
      cat("  ✓ Boxplot salvo:", filepath, "\n")
      cat("    - Cenário:", scenario_name, "\n")
      cat("    - w (geração):", w_gen, "| w (inferência):", w_inf, "| T:", t_val, "\n")
    }
  }
  
  # Salvar todos os dados combinados
  all_lambda_t <- bind_rows(lambda_t_all)
  if (nrow(all_lambda_t) > 0) {
    write.csv(all_lambda_t, "lambda_t_bias_data.csv", row.names = FALSE)
    cat("✓ Arquivo lambda_t_bias_data.csv salvo\n")
  }
} else {
  cat("\n✗ AVISO: Nenhum dado de lambda_t para processar\n")
}

# Criar tabelas de métricas otimizadas (apenas 2 arquivos por cenário)
if (length(individual_results) > 0) {
  cat("\nCriando tabelas de métricas otimizadas...\n")
  cat("NOTA: Cada cenário gerará apenas 2 arquivos:\n")
  cat("  1. consolidated_*.rds (todos os dados em formato R)\n")
  cat("  2. all_metrics_*.csv (todas as métricas em formato CSV)\n\n")
  
  # Contadores para w de geração
  count_tables_w07 <- 0
  count_tables_w09 <- 0
  
  # Gerar arquivos essenciais para cada cenário
  for (i in seq_along(individual_results)) {
    result_item <- individual_results[[i]]
    
    if (!is.null(result_item)) {
      scenario_name <- result_item$scenario
      file_name <- result_item$file
      w_gen <- result_item$w_generation
      w_inf <- result_item$w_inference
      t_val <- result_item$T_config
      
      # Criar nome do arquivo
      safe_name <- gsub("[^a-zA-Z0-9_]", "_", paste(scenario_name, w_gen, sep = "_"))
      
      # Determinar diretório baseado no w de geração
      if (w_gen == "w07") {
        table_dir <- tables_w07_dir
        count_tables_w07 <- count_tables_w07 + 1
      } else if (w_gen == "w09") {
        table_dir <- tables_w09_dir
        count_tables_w09 <- count_tables_w09 + 1
      } else {
        table_dir <- tables_dir
      }
      
      # Criar tabela consolidada para este cenário
      beta_table <- result_item$beta
      epsilon_table <- result_item$epsilon
      lambda_table <- result_item$lambda
      
      # Tabela de WAIC/LPML
      criteria_table <- tibble(
        Metrica = c("WAIC", "LPML"),
        Valor = c(result_item$WAIC, result_item$LPML)
      )
      
      # Metadados
      metadata_table <- tibble(
        Campo = c("Cenário", "Arquivo", "w (geração)", "w (inferência)", "T", "Data processamento"),
        Valor = c(
          as.character(scenario_name), 
          as.character(file_name), 
          as.character(w_gen), 
          as.character(w_inf), 
          as.character(t_val), 
          as.character(Sys.time())
        )
      )
      
      # 1. Criar relatório consolidado em um único arquivo RDS (para uso em R)
      consolidated_data <- list(
        Metadados = metadata_table,
        Beta = beta_table,
        Epsilon = epsilon_table,
        Lambda = lambda_table,
        Criteria = criteria_table
      )
      
      # Salvar como RDS (formato binário do R)
      rds_file <- file.path(table_dir, paste0("consolidated_", safe_name, ".rds"))
      saveRDS(consolidated_data, rds_file)
      cat("  ✓ Arquivo RDS salvo:", basename(rds_file), "\n")
      
      # 2. Criar um único arquivo CSV com todas as informações (para uso geral)
      all_tables <- bind_rows(
        # Metadados - já está no formato correto e character
        metadata_table %>% 
          mutate(Tabela = "Metadados") %>%
          select(Tabela, Campo, Valor),
        
        # Beta - converter Valor para character explicitamente
        beta_table %>%
          pivot_longer(cols = c(bias, eqm), names_to = "Campo", values_to = "Valor_num") %>%
          mutate(
            Tabela = "Beta",
            Campo = paste(param, Campo, sep = "_"),
            Valor = as.character(Valor_num)
          ) %>%
          select(Tabela, Campo, Valor),
        
        # Epsilon - converter Valor para character explicitamente
        epsilon_table %>%
          pivot_longer(cols = c(bias, eqm), names_to = "Campo", values_to = "Valor_num") %>%
          mutate(
            Tabela = "Epsilon",
            Campo = paste("cluster", cluster, Campo, sep = "_"),
            Valor = as.character(Valor_num)
          ) %>%
          select(Tabela, Campo, Valor),
        
        # Lambda - converter Valor para character explicitamente
        lambda_table %>%
          pivot_longer(cols = c(bias, eqm), names_to = "Campo", values_to = "Valor_num") %>%
          mutate(
            Tabela = "Lambda",
            Campo = paste("cluster", cluster, Campo, sep = "_"),
            Valor = as.character(Valor_num)
          ) %>%
          select(Tabela, Campo, Valor),
        
        # Criteria - converter Valor para character explicitamente
        criteria_table %>%
          mutate(
            Tabela = "Criteria",
            Campo = Metrica,
            Valor = as.character(Valor)
          ) %>%
          select(Tabela, Campo, Valor)
      )
      
      # Garantir que todos os valores sejam character
      all_tables <- all_tables %>%
        mutate(Valor = as.character(Valor))
      
      csv_file <- file.path(table_dir, paste0("all_metrics_", safe_name, ".csv"))
      write.csv(all_tables, csv_file, row.names = FALSE)
      cat("  ✓ Arquivo CSV salvo:", basename(csv_file), "\n")
      cat("    - Cenário:", scenario_name, "(w geração =", w_gen, ")\n")
      cat("    - Diretório:", table_dir, "\n\n")
    }
  }
  
  # Resumo das tabelas geradas
  cat("\nRESUMO DAS TABELAS GERADAS (VERSÃO OTIMIZADA):\n")
  cat("- Total de cenários processados:", length(individual_results), "\n")
  cat("- Cenários para w (geração) = 0.7:", count_tables_w07, "\n")
  cat("- Cenários para w (geração) = 0.9:", count_tables_w09, "\n")
  cat("- Arquivos por cenário: 2 (consolidated_*.rds + all_metrics_*.csv)\n")
  cat("- Total de arquivos de métricas:", 2 * length(individual_results), "\n")
  cat("- Tabelas salvas em:", tables_dir, "\n")
  
  # Salvar um resumo consolidado de todas as métricas
  if (length(individual_results) > 0) {
    # Extrair WAIC/LPML de todos os cenários
    all_criteria <- map_dfr(individual_results, function(x) {
      tibble(
        Scenario = x$scenario,
        File = x$file,
        w_generation = x$w_generation,
        w_inference = x$w_inference,
        T_config = x$T_config,
        WAIC = x$WAIC,
        LPML = x$LPML
      )
    })
    
    write.csv(all_criteria, file.path(tables_dir, "all_scenarios_criteria.csv"), row.names = FALSE)
    cat("✓ Arquivo all_scenarios_criteria.csv salvo\n")
  }
} else {
  cat("\n✗ AVISO: Nenhum resultado individual para salvar\n")
}

# Salvar resultados consolidados
if (length(results) > 0) {
  saveRDS(results, "all_results_consolidated.rds")
  cat("\n✓ Resultados consolidados salvos em all_results_consolidated.rds\n")
}

cat("\n=======================================\n")
cat("PROCESSAMENTO CONCLUÍDO!\n")
cat("=======================================\n")

# Resumo final
cat("\nRESUMO FINAL DA EXECUÇÃO:\n")
cat("- Total de arquivos samples encontrados:", length(scenario_files), "\n")
cat("- Arquivos processados com sucesso:", processed_count, "\n")
cat("- Erros durante o processamento:", error_count, "\n")
cat("- Total de boxplots gerados:", length(lambda_t_all), "\n")
cat("- Total de tabelas de métricas geradas (otimizadas):", 
    ifelse(length(individual_results) > 0, 2 * length(individual_results), 0), "\n")
cat("\nORGANIZAÇÃO DOS ARQUIVOS:\n")
cat("1. Boxplots por w de geração:\n")
cat("   *", boxplot_w07_dir, " (", count_w07, " arquivos PNG)\n", sep = "")
cat("   *", boxplot_w09_dir, " (", count_w09, " arquivos PNG)\n", sep = "")
cat("2. Tabelas de métricas por w de geração:\n")
cat("   *", tables_w07_dir, " (", count_tables_w07 * 2, " arquivos: RDS + CSV por cenário)\n", sep = "")
cat("   *", tables_w09_dir, " (", count_tables_w09 * 2, " arquivos: RDS + CSV por cenário)\n", sep = "")
cat("3. Arquivos consolidados (gerais):\n")
cat("   * WAIC_LPML_summary.csv\n")
cat("   * WAIC_LPML_resumo.csv\n")
cat("   * lambda_t_bias_data.csv\n")
cat("   * all_results_consolidated.rds\n")
cat("   *", file.path(tables_dir, "all_scenarios_criteria.csv"), "\n")
cat("\nARQUITETURA OTIMIZADA:\n")
cat("Para cada cenário são gerados APENAS 2 arquivos:\n")
cat("  1. consolidated_[CENARIO]_[WGEN].rds - Todos os dados em formato R (para análise posterior)\n")
cat("  2. all_metrics_[CENARIO]_[WGEN].csv - Todas as métricas em CSV (para visualização/exportação)\n")
cat("  (Total: 2 arquivos por cenário, em vez de 7 anteriormente)\n")