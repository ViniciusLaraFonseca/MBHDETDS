# Nome: ValoresIniciais_debug_FatorComum.R
# Objetivo: Simular dados de acordo com o modelo de FATOR DINÂMICO COMUM (lambda_t univariado)
# para testes de identificabilidade.

set.seed(42)
cat("--- Iniciando Simulação (Fator Comum) ---\n")

# --- PASSO 1: DEFINIR PARÂMETROS DA SIMULAÇÃO ---
n_regions <- 75 
n_times   <- 23
p         <- 3 # N° de covariáveis beta
K         <- 4 # N° de clusters gamma

# Parâmetros "Verdadeiros" para Simulação
beta_true  <- c(-1, 1, 0.5) # Valores do seu script v8
gamma_true <- c(0.05, 0.1, 0.1, 0.15)
w_true     <- 0.9
a0_true    <- 1.0
b0_true    <- 1.0

# --- PASSO 2: SIMULAR DADOS ESTÁTICOS (x, E, epsilon) ---
cat("--- PASSO 2: Simulando x, E, epsilon ---\n")

# Simular Covariáveis 'x' (n_regions, n_times, p)
x_true <- array(rnorm(n_regions * n_times * p), dim = c(n_regions, n_times, p))

# Simular e Normalizar Offset 'E'
E_raw <- matrix(runif(n_regions * n_times, 150, 250), nrow = n_regions)
E_true <- E_raw / mean(E_raw) # Normalizado

# Calcular Epsilon Fixo (baseado em gamma_true)
# (Assumindo que 'hAI' está carregado do _dataCaseStudy.r)
if (!exists("hAI")) {
  stop("Erro: Objeto 'hAI' não encontrado. Carregue '_dataCaseStudy.r' e use attach(data).")
}
epsilon_true <- matrix(NA, nrow = n_regions, ncol = 1)
for(i in 1:n_regions) {
  epsilon_true[i] <- 1 - sum(hAI[i, ] * gamma_true)
}

# Pré-calcular o componente espacial/regressão g_it
# g_it = E_it * epsilon_i * exp(x_it' * beta)
g_it_true <- array(NA, dim = c(n_regions, n_times))
for(i in 1:n_regions) {
  for(t in 1:n_times) {
    prod_val <- sum(x_true[i, t, ] * beta_true)
    g_it_true[i, t] <- E_true[i, t] * epsilon_true[i] * exp(prod_val)
  }
}

# --- PASSO 3: SIMULAR DADOS DINÂMICOS (lambda_t e Y_it) ---
cat("--- PASSO 3: Simulando trajetória lambda_t e Y_it ---\n")

# Inicialização de vetores (NÃO mais matrizes) para o filtro
lambda_true <- rep(NA, n_times)
Y_ini       <- matrix(NA, nrow = n_regions, ncol = n_times) # Matriz de observações

at_true  <- rep(NA, n_times + 1)
bt_true  <- rep(NA, n_times + 1)
att_true <- rep(NA, n_times)
btt_true <- rep(NA, n_times)

# Simulação Forward (Agora univariada)
at_true[1] <- a0_true
bt_true[1] <- b0_true

for(t in 1:n_times) {
  # 1. Passos Preditivos (para o tempo t)
  att_true[t] <- w_true * at_true[t]
  btt_true[t] <- w_true * bt_true[t]
  
  # 2. Amostrar o lambda_t "verdadeiro" da preditiva
  lambda_true[t] <- rgamma(1, shape = att_true[t], rate = btt_true[t])
  
  # 3. Gerar TODAS as observações Y_it para o tempo t
  for(i in 1:n_regions) {
    mu_it <- lambda_true[t] * g_it_true[i, t]
    Y_ini[i, t] <- rpois(1, mu_it)
  }
  
  # 4. Atualizar os parâmetros do filtro (para o tempo t+1)
  sum_Y_t <- sum(Y_ini[, t])
  sum_g_t <- sum(g_it_true[, t])
  
  at_true[t+1] <- att_true[t] + sum_Y_t
  bt_true[t+1] <- btt_true[t] + sum_g_t
}

cat("Simulação concluída.\n")
cat("Verdadeiros (para traceplots):\n")
print(beta_true)
print(gamma_true)
print(lambda_true[c(1, 10, 23)]) # Valores de exemplo

# --- PASSO 4: DEFINIR OBJETOS PARA O NIMBLE ---
# Estes são os objetos que você passará para o script do modelo
constants_nimble <- list(
  n_regions = n_regions, 
  n_times   = n_times, 
  p         = p, 
  K         = K, 
  h         = hAI,
  mu_beta   = rep(0, p), # Priori para beta
  a_unif    = c(0.00, rep(0, K - 1)),
  b_unif    = c(0.1, rep(0.99, K - 1)), # Priori informativa para gamma[1]
  a0        = a0_true, 
  b0        = b0_true, 
  w         = w_true
)

data_nimble <- list(
  Y = Y_ini, 
  E = E_true, 
  x = x_true
)

inits_nimble <- list(
  beta   = beta_true,     # Iniciar no valor real
  gamma  = gamma_true,    # Iniciar no valor real
  lambda = lambda_true    # Iniciar no valor real (agora um VETOR)
)

inits_nimble_cadeia2 <- list(
  beta   = rnorm(p, 0, 0.5),
  gamma  = c(0.02, 0.05, 0.05, 0.05), # Perturbado
  lambda = rgamma(n_times, 1, 1)      # Aleatório
)

inits_list_nimble <- list(inits_nimble, inits_nimble_cadeia2)

cat("--- Objetos 'constants_nimble', 'data_nimble' e 'inits_list_nimble' prontos. ---\n")