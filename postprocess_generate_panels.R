###############################################################
# postprocess_generate_panels_FINAL.R
# Versão final — robusta ao nome das colunas de lambda
# COM HPD (Highest Posterior Density) em vez de intervalos de credibilidade
# Fórmulas:
#   theta_it = lambda_it * exp(beta^T x_it)
#   mu_it    = E_it * epsilon_i * theta_it
###############################################################
regions_of_interest <- c(1,8,15,19,22,31,34,40,46,55,65,75)
start_time <- Sys.time()
logfile <- "postprocess_log.txt"
writeLines(paste("Início:", start_time), con = logfile)

safe_log <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - ",
                paste(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = logfile, append = TRUE)
}

safe_log("Starting postprocess_generate_panels_FINAL.R (com HPD)")

# Ajuste aqui se necessário
base_dir <- "C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/MBHDETDS"
base_dir <- gsub("\\\\", "/", base_dir)
setwd(base_dir)
safe_log("Working dir:", getwd())

# Detect root dirs
root_candidates <- c("resultados_sensibilidade_w07", "resultados_sensibilidade_w09")
root_dirs <- root_candidates[file.exists(root_candidates)]
if(length(root_dirs) == 0) stop("Nenhuma pasta resultados_sensibilidade_* encontrada no base_dir")
safe_log("Root dirs:", paste(root_dirs, collapse = " | "))

# load packages
pkgs <- c("coda","ggplot2","dplyr","stringr","tibble","purrr","readr")
for(p in pkgs) if(!requireNamespace(p, quietly=TRUE)) install.packages(p)
library(coda); library(ggplot2); library(dplyr); library(stringr); library(tibble); library(readr)

# map scenario -> init file
map_scenario_to_init <- function(scen_name) {
  scen <- basename(scen_name)
  if (grepl("T100", scen) && grepl("w09", scen)) return("ValoresIniciais_Lambda_it_T100_w09.R")
  if (grepl("T100", scen) && grepl("w07", scen)) return("ValoresIniciais_Lambda_it_T100_w07.R")
  if (grepl("T23",  scen) && grepl("w09", scen)) return("ValoresIniciais_Lambda_it_T23_w09.R")
  if (grepl("T23",  scen) && grepl("w07", scen)) return("ValoresIniciais_Lambda_it_T23_w07.R")
  stop("Cenário não reconhecido: ", scen)
}

# load scenario: sources and standardizes x, E, epsilon_vec
load_scenario_data <- function(scen_dir) {
  scen <- basename(scen_dir)
  safe_log("Loading scenario:", scen_dir)
  
  # try sourcing _dataCaseStudy.R from a few likely places
  tried <- c(file.path(scen_dir, "_dataCaseStudy.R"), "_dataCaseStudy.R", file.path("..", "_dataCaseStudy.R"))
  found <- FALSE
  for(p in tried) if(file.exists(p)) { source(p, local = .GlobalEnv); safe_log("Sourced", p); found = TRUE; break }
  if(!found) safe_log("Warning: _dataCaseStudy.R not found in tried paths:", paste(tried, collapse=" | "))
  
  # initial values file
  init_file <- map_scenario_to_init(scen)
  tried2 <- c(file.path(scen_dir, init_file), init_file, file.path("..", init_file))
  found2 <- NULL
  for(p in tried2) if(file.exists(p)) { source(p, local = .GlobalEnv); found2 = p; break }
  if(is.null(found2)) {
    safe_log("ERROR: initial-values file not found for", scen, "tried:", paste(tried2, collapse=" | "))
    return(FALSE)
  }
  safe_log("Sourced initial-values:", found2)
  
  # standardize x
  if (exists("x_true")) assign("x", x_true, envir = .GlobalEnv)
  else if (exists("x")) assign("x", x, envir = .GlobalEnv)
  else if (exists("X")) assign("x", X, envir = .GlobalEnv)
  
  # standardize E
  if (exists("E_true")) assign("E", E_true, envir = .GlobalEnv)
  else if (exists("E_raw")) assign("E", E_raw, envir = .GlobalEnv)
  else if (exists("E")) assign("E", E, envir = .GlobalEnv)
  
  # epsilon
  if (exists("epsilon")) assign("epsilon_vec", as.numeric(epsilon), envir = .GlobalEnv)
  else if (exists("epsilon_true")) assign("epsilon_vec", as.numeric(epsilon_true), envir = .GlobalEnv)
  else if (exists("hAI") && exists("gamma_true")) {
    assign("epsilon_vec", as.numeric(1 - as.vector(hAI %*% gamma_true)), envir = .GlobalEnv)
    safe_log("epsilon_vec computed from hAI %*% gamma_true")
  } else {
    safe_log("ERROR: epsilon not available and cannot be computed.")
    return(FALSE)
  }
  
  # checks
  if(!exists("x") || !exists("E") || !exists("epsilon_vec")) {
    safe_log("ERROR: required objects x/E/epsilon_vec not present after sourcing.")
    return(FALSE)
  }
  
  if (length(dim(x)) != 3) stop("x must be array [region, time, p]. Current dims: ", paste(dim(x), collapse=","))
  if (!all(dim(E)[1:2] == dim(x)[1:2])) stop("E dims must match x dims. E:", paste(dim(E), collapse=","), " x:", paste(dim(x)[1:2], collapse=","))
  
  safe_log("Scenario data ready. dims x:", paste(dim(x), collapse="x"), "E:", paste(dim(E), collapse="x"))
  return(TRUE)
}

# read mcmc.list and combine chains (thinned)
read_mcmc_samples <- function(sfile, thin = 10) {
  obj <- readRDS(sfile)
  if(!inherits(obj, "mcmc.list")) stop("RDS is not mcmc.list:", sfile)
  mats <- lapply(obj, function(ch) as.matrix(ch))
  M <- do.call(rbind, mats)
  if(thin > 1) {
    idx <- seq(1, nrow(M), by = thin)
    M <- M[idx, , drop = FALSE]
  }
  return(M)
}

# compute stats for one (region i, time t) - COM HPD
compute_stats_for_it <- function(lambda_vec, beta_mat, x_it, E_it, eps_i) {
  xb <- as.numeric(beta_mat %*% x_it)   # S-length
  exp_xb <- exp(xb)
  theta_vec <- lambda_vec * exp_xb
  mu_vec <- as.numeric(E_it * eps_i * theta_vec)
  
  # Calcular HPD (95%) usando coda::HPDinterval
  theta_hpd <- HPDinterval(mcmc(theta_vec), prob = 0.95)
  mu_hpd <- HPDinterval(mcmc(mu_vec), prob = 0.95)
  
  list(
    theta_mean = mean(theta_vec, na.rm=TRUE),
    theta_low  = theta_hpd[1, "lower"],
    theta_up   = theta_hpd[1, "upper"],
    mu_mean    = mean(mu_vec, na.rm=TRUE),
    mu_low     = mu_hpd[1, "lower"],
    mu_up      = mu_hpd[1, "upper"]
  )
}

# build panels robust to lambda column naming and ordering
build_panels <- function(mcmc_mat, x, E, epsilon_vec) {
  varn <- colnames(mcmc_mat)
  beta_idx <- grep("^beta\\[", varn, value = FALSE)
  lambda_idx <- grep("^lambda", varn, value = FALSE)  # any lambda*
  if(length(beta_idx) == 0) stop("No beta found in MCMC cols")
  if(length(lambda_idx) == 0) stop("No lambda found in MCMC cols")
  
  # infer R and T from x (most reliable)
  R <- dim(x)[1]; Tt <- dim(x)[2]
  L <- length(lambda_idx)
  if (L != R * Tt) stop("Number of lambda columns (", L, ") != R*T (", R * Tt, "). Cannot map lambdas.")
  
  # Because your head(...) showed order: lambda[1,1], lambda[2,1],...,lambda[75,1], lambda[1,2],...
  # we fill matrix column-wise (by column -> time), so that lambda_mat[r,t] = column index in mcmc_mat
  lambda_mat <- matrix(lambda_idx, nrow = R, ncol = Tt, byrow = FALSE)
  
  safe_log("Building panels: S=", nrow(mcmc_mat), " R=", R, " T=", Tt, " p=", length(beta_idx))
  
  beta_mat <- mcmc_mat[, beta_idx, drop = FALSE]  # S x p
  
  rows <- R * Tt
  out_theta <- tibble(region = integer(rows), time = integer(rows), 
                      Mean = double(rows), low95 = double(rows), upp95 = double(rows))
  out_mu <- tibble(region = integer(rows), time = integer(rows), 
                   Mean = double(rows), low95 = double(rows), upp95 = double(rows))
  
  k <- 1
  for(r in 1:R) {
    for(t in 1:Tt) {
      lam_col_idx <- lambda_mat[r, t]
      lambda_vec <- mcmc_mat[, lam_col_idx]
      
      x_it <- as.numeric(x[r, t, ])
      E_it <- as.numeric(E[r, t])
      eps_i <- as.numeric(epsilon_vec[r])
      
      stats <- compute_stats_for_it(lambda_vec, beta_mat, x_it, E_it, eps_i)
      
      out_theta[k, ] <- list(region = r, time = t, 
                             Mean = stats$theta_mean, 
                             low95 = stats$theta_low, 
                             upp95 = stats$theta_up)
      out_mu[k, ]    <- list(region = r, time = t,
                             Mean = stats$mu_mean,
                             low95 = stats$mu_low,
                             upp95 = stats$mu_up)
      k <- k + 1
    }
  }
  
  list(theta = out_theta %>% arrange(region, time), 
       mu = out_mu %>% arrange(region, time))
}

save_panel_from_df <- function(df, outfile, ylab = "") {
  if(!dir.exists(dirname(outfile))) dir.create(dirname(outfile), recursive = TRUE)
  
  df <- df %>% filter(region %in% regions_of_interest)
  
  p <- ggplot(df, aes(x = time)) +
    geom_ribbon(aes(ymin = low95, ymax = upp95), alpha = 0.25) +
    geom_line(aes(y = Mean), size = 0.9, color = "black") +
    geom_line(aes(y = true_val), color = "red", size = 0.9, na.rm = TRUE) +
    facet_wrap(~ region, ncol = 3, nrow = 4, scales = "free_y") +
    labs(x = "Tempo", y = ylab, title = basename(outfile)) +
    theme_bw(base_size = 12)
  
  ggsave(outfile, p, width = 12, height = 10, dpi = 150)
  safe_log("Saved", outfile)
}


# main loop
scenario_dirs <- unlist(lapply(root_dirs, function(r) list.dirs(r, recursive = FALSE)))
safe_log("Scenario dirs:", paste(scenario_dirs, collapse = " | "))

panels_created <- 0
thin_factor <- 10

for(sc in scenario_dirs) {
  safe_log("Processing scenario:", sc)
  ok <- tryCatch(load_scenario_data(sc), error = function(e){ safe_log("Error loading scenario:", e$message); FALSE })
  if(!isTRUE(ok)) { safe_log("Skipping scenario due to load failure:", sc); next }
  
  sfiles <- list.files(sc, pattern = "^samples_.*\\.rds$", full.names = TRUE)
  safe_log("Found sample files:", paste(sfiles, collapse = " | "))
  if(length(sfiles) == 0) next
  
  for(sfile in sfiles) {
    safe_log("Reading:", sfile)
    mcmc_mat <- tryCatch(read_mcmc_samples(sfile, thin = thin_factor), 
                         error = function(e){ safe_log("Error reading mcmc:", e$message); NULL })
    if(is.null(mcmc_mat)) next
    
    # quick column sanity
    if(!any(grepl("^lambda", colnames(mcmc_mat)))) { safe_log("No lambda columns in", sfile); next }
    if(!any(grepl("^beta\\[", colnames(mcmc_mat)))) { safe_log("No beta columns in", sfile); next }
    
    panels <- tryCatch(build_panels(mcmc_mat, x, E, epsilon_vec), 
                       error = function(e){ safe_log("Error building panels:", e$message); NULL })
    if(is.null(panels)) next
    
    # attach true values if available
    if (exists("lambda_true") && exists("beta_true")) {
      R <- dim(x)[1]; Tt <- dim(x)[2]
      
      theta_true_mat <- matrix(NA_real_, nrow = R, ncol = Tt)
      mu_true_mat    <- matrix(NA_real_, nrow = R, ncol = Tt)
      
      for(i in 1:R) for(t in 1:Tt) {
        xb_true <- sum(beta_true * as.numeric(x[i,t,]))
        theta_true_mat[i,t] <- lambda_true[i,t] * exp(xb_true)
        mu_true_mat[i,t]    <- E[i,t] * epsilon_vec[i] * theta_true_mat[i,t]
      }
      
      # long format for merging
      long_theta <- tibble(
        region = rep(1:R, times = Tt),
        time   = rep(1:Tt, each = R),
        true_val = as.vector(theta_true_mat)
      )
      long_mu <- tibble(
        region = rep(1:R, times = Tt),
        time   = rep(1:Tt, each = R),
        true_val = as.vector(mu_true_mat)
      )
      
      panels$theta <- panels$theta %>% left_join(long_theta, by = c("region","time"))
      panels$mu    <- panels$mu    %>% left_join(long_mu,    by = c("region","time"))
    } else {
      panels$theta$true_val <- NA_real_
      panels$mu$true_val    <- NA_real_
    }
    
    
    out_dir <- file.path(sc, "lambdas")
    if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    base <- tools::file_path_sans_ext(basename(sfile))
    out_theta <- file.path(out_dir, paste0("painel_theta_", base, ".png"))
    out_mu    <- file.path(out_dir, paste0("painel_mu_", base, ".png"))
    
    save_panel_from_df(panels$theta, out_theta, ylab = "theta[t]")
    save_panel_from_df(panels$mu, out_mu, ylab = "mu[t]")
    
    panels_created <- panels_created + 2
  }
}

safe_log("PROCESS COMPLETE. Panels created:", panels_created)
safe_log("Total time (mins):", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")),2))
cat("Finished. See", logfile, "for details.\n")