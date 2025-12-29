# Nome: ValoresIniciais_Lambda_it_T23.R
# Objetivo: Simular dados específicos para T=23
# para os cenários C3 e C4.

set.seed(987) # Mesmo seed para consistência de parêmetros, mas a geração aleatória mudará pela dimensão
cat("--- Iniciando Simulação (lambda_it, T=23) ---\n")

# --- PASSO 1: DEFINIR PARÂMETROS ---
n_regions <- 75 
n_times   <- 23 # T REDUZIDO
p         <- 3   
K         <- 4   

# Parâmetros "Verdadeiros"
beta_true  <- c(-1.0, 1.0, 0.5)
gamma_true <- c(0.05, 0.1, 0.1, 0.15)
w_true     <- 0.7
a0_true    <- 1.0
b0_true    <- 1.0

# --- PASSO 2: SIMULAR DADOS ESTÁTICOS ---
cat("--- PASSO 2: Simulando x, E, epsilon (T=23) ---\n")

x_true <- array(rnorm(n_regions * n_times * p), dim = c(n_regions, n_times, p))
E_raw  <- matrix(runif(n_regions * n_times, 150, 250), nrow = n_regions)
E_true <- E_raw / mean(E_raw)

# Carrega hAI do dataCaseStudy se não existir
if (!exists("hAI")) {
  source("_dataCaseStudy.r")
}
epsilon_true <- matrix(NA, nrow = n_regions, ncol = 1)
for(i in 1:n_regions) {
  epsilon_true[i] <- 1 - sum(hAI[i, ] * gamma_true)
}

g_it_true <- array(NA, dim = c(n_regions, n_times))
for(i in 1:n_regions) {
  for(t in 1:n_times) {
    prod_val <- sum(x_true[i, t, ] * beta_true)
    g_it_true[i, t] <- E_true[i, t] * epsilon_true[i] * exp(prod_val)
  }
}

# --- PASSO 3: SIMULAR DADOS DINÂMICOS ---
cat("--- PASSO 3: Simulando trajetória lambda_it e Y_it (T=23) ---\n")

lambda_true <- matrix(NA, ncol = n_times, nrow = n_regions)
Y_ini       <- matrix(NA, nrow = n_regions, ncol = n_times)

at_true  <- matrix(NA, nrow = n_regions, ncol = n_times + 1)
bt_true  <- matrix(NA, nrow = n_regions, ncol = n_times + 1)
att_true <- matrix(NA, nrow = n_regions, ncol = n_times)
btt_true <- matrix(NA, nrow = n_regions, ncol = n_times)

for(i in 1:n_regions) {
  at_true[i, 1] <- a0_true
  bt_true[i, 1] <- b0_true
  
  for(t in 1:n_times) {
    att_true[i, t] <- w_true * at_true[i, t]
    btt_true[i, t] <- w_true * bt_true[i, t]
    lambda_true[i, t] <- rgamma(1, shape = att_true[i, t], rate = btt_true[i, t])
    mu_it <- lambda_true[i, t] * g_it_true[i, t]
    Y_ini[i, t] <- rpois(1, mu_it)
    at_true[i, t+1] <- att_true[i, t] + Y_ini[i, t]
    bt_true[i, t+1] <- btt_true[i, t] + g_it_true[i, t]
  }
}

cat("Simulação T=23 concluída.\n")

# --- PASSO 4: DEFINIR OBJETOS COM SUFIXO _T23 ---
constants_nimble_T23 <- list(
  n_regions = n_regions, 
  n_times   = n_times, 
  p         = p, 
  K         = K, 
  h         = hAI,
  mu_beta   = rep(0, p),
  a_unif    = c(0.00, rep(0, K - 1)),
  b_unif    = c(0.1, rep(0.99, K - 1)),
  a0        = a0_true, 
  b0        = b0_true, 
  w         = w_true
)

data_nimble_T23 <- list(
  Y = Y_ini, 
  E = E_true, 
  x = x_true
)

# Iniciais T23
inits_nimble_1 <- list(beta = beta_true, gamma = gamma_true, lambda = lambda_true)
inits_nimble_2 <- list(
  beta   = rnorm(p, 0, 0.5),
  gamma  = c(0.02, 0.05, 0.05, 0.05),
  lambda = matrix(rgamma(n_regions * n_times, 1, 1), nrow = n_regions)
)

inits_list_nimble_T23 <- list(inits_nimble_1, inits_nimble_2)
lambda_true_T23 <- lambda_true # Salvar referência explícita

cat("--- Objetos '_T23' prontos. ---\n")