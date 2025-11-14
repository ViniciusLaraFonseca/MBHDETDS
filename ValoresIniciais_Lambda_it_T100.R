# Nome: ValoresIniciais_Lambda_it_T100.R
# Objetivo: Simular dados (T=100) de acordo com o modelo original (lambda_it)
# para testar a hipótese de identificabilidade com séries longas.

set.seed(987)
cat("--- Iniciando Simulação (lambda_it, T=100) ---\n")

# --- PASSO 1: DEFINIR PARÂMETROS DA SIMULAÇÃO ---
n_regions <- 75 
n_times   <- 100 # Nova dimensão de tempo
p         <- 3   # N° de covariáveis beta
K         <- 4   # N° de clusters gamma

# Parâmetros "Verdadeiros" para Simulação
beta_true  <- c(-1.0, 1.0, 0.5) # Valores do teste original 'gamma_fixed'
gamma_true <- c(0.05, 0.1, 0.1, 0.15)
w_true     <- 0.9
a0_true    <- 1.0
b0_true    <- 1.0

# --- PASSO 2: SIMULAR DADOS ESTÁTICOS (x, E, epsilon) ---
cat("--- PASSO 2: Simulando x, E, epsilon ---\n")

x_true <- array(rnorm(n_regions * n_times * p), dim = c(n_regions, n_times, p))
E_raw  <- matrix(runif(n_regions * n_times, 150, 250), nrow = n_regions)
E_true <- E_raw / mean(E_raw) # Normalizado

if (!exists("hAI")) {
  stop("Erro: Objeto 'hAI' não encontrado. Carregue '_dataCaseStudy.r'.")
}
epsilon_true <- matrix(NA, nrow = n_regions, ncol = 1)
for(i in 1:n_regions) {
  epsilon_true[i] <- 1 - sum(hAI[i, ] * gamma_true)
}

# Pré-calcular o componente de regressão g_it = E_it * epsilon_i * exp(x_it' * beta)
g_it_true <- array(NA, dim = c(n_regions, n_times))
for(i in 1:n_regions) {
  for(t in 1:n_times) {
    prod_val <- sum(x_true[i, t, ] * beta_true)
    g_it_true[i, t] <- E_true[i, t] * epsilon_true[i] * exp(prod_val)
  }
}

# --- PASSO 3: SIMULAR DADOS DINÂMICOS (lambda_it e Y_it) ---
cat("--- PASSO 3: Simulando trajetória lambda_it e Y_it ---\n")

# Buffers de filtragem (agora são matrizes 2D)
lambda_true <- matrix(NA, ncol = n_times, nrow = n_regions)
Y_ini       <- matrix(NA, nrow = n_regions, ncol = n_times)

at_true  <- matrix(NA, nrow = n_regions, ncol = n_times + 1)
bt_true  <- matrix(NA, nrow = n_regions, ncol = n_times + 1)
att_true <- matrix(NA, nrow = n_regions, ncol = n_times)
btt_true <- matrix(NA, nrow = n_regions, ncol = n_times)

# Simulação Forward (agora por área, i)
for(i in 1:n_regions) {
  at_true[i, 1] <- a0_true
  bt_true[i, 1] <- b0_true
  
  for(t in 1:n_times) {
    # 1. Preditiva (para tempo t)
    att_true[i, t] <- w_true * at_true[i, t]
    btt_true[i, t] <- w_true * bt_true[i, t]
    
    # 2. Amostrar lambda_it "verdadeiro" da preditiva
    lambda_true[i, t] <- rgamma(1, shape = att_true[i, t], rate = btt_true[i, t])
    
    # 3. Gerar Y_it
    mu_it <- lambda_true[i, t] * g_it_true[i, t]
    Y_ini[i, t] <- rpois(1, mu_it)
    
    # 4. Atualizar (para tempo t+1)
    at_true[i, t+1] <- att_true[i, t] + Y_ini[i, t]
    bt_true[i, t+1] <- btt_true[i, t] + g_it_true[i, t]
  }
}

cat("Simulação concluída.\n")
cat("Verdadeiros (para traceplots):\n")
print(beta_true)
print(gamma_true)
cat("Lambda[1,1] Verdadeiro:", lambda_true[1,1], "\n")

# --- PASSO 4: DEFINIR OBJETOS PARA O NIMBLE ---
constants_nimble <- list(
  n_regions = n_regions, 
  n_times   = n_times, 
  p         = p, 
  K         = K, 
  h         = hAI,
  mu_beta   = rep(0, p),
  a_unif    = c(0.00, rep(0, K - 1)),
  b_unif    = c(0.1, rep(0.99, K - 1)), # Prior informativa
  a0        = a0_true, 
  b0        = b0_true, 
  w         = w_true
)

data_nimble <- list(
  Y = Y_ini, 
  E = E_true, 
  x = x_true
)

# Iniciais (Cadeia 1 = valores verdadeiros, Cadeia 2 = perturbado)
inits_nimble_1 <- list(
  beta   = beta_true,
  gamma  = gamma_true,
  lambda = lambda_true 
)

inits_nimble_2 <- list(
  beta   = rnorm(p, 0, 0.5),
  gamma  = c(0.02, 0.05, 0.05, 0.05),
  lambda = matrix(rgamma(n_regions * n_times, 1, 1), nrow = n_regions)
)

inits_list_nimble <- list(inits_nimble_1, inits_nimble_2)

cat("--- Objetos 'constants_nimble', 'data_nimble' e 'inits_list_nimble' prontos. ---\n")