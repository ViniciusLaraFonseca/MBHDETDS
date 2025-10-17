# ------------------------------------------------------------------
# SCRIPT DE GERAÇÃO DE VALORES INICIAIS
# ------------------------------------------------------------------
# Objetivo: Gerar um conjunto de valores iniciais plausíveis para os parâmetros
# e variáveis latentes do modelo dinâmico. Um bom ponto de partida é
# crucial para acelerar a convergência do MCMC e evitar que as cadeias
# fiquem presas em regiões de baixa probabilidade.
# ------------------------------------------------------------------


# --- PASSO 1: NORMALIZAÇÃO DO OFFSET (E) ---

# O offset 'E' (população exposta ou número esperado de eventos) pode ter valores
# muito grandes, o que causa instabilidade numérica no MCMC, especialmente
# em termos que envolvem a função exp().
# A melhor prática é escalonar 'E' para que a sua média global seja 1.
# Isto estabiliza os cálculos e torna o risco relativo (theta) mais interpretável.

cat("A normalizar a matriz de offset 'E'...\n")
media_global_E <- mean(E)
E_escalonado <- E / media_global_E
cat("Normalização de 'E' concluída. Média de E agora é:", mean(E_escalonado), "\n\n")


# --- PASSO 2: OBTENÇÃO DE VALORES INICIAIS PARA BETA, a0 e b0 VIA GLM ---

# Para obter um ponto de partida razoável para os coeficientes de regressão (beta),
# ajustamos um Modelo Linear Generalizado (GLM) de Poisson simples a um corte 
# transversal dos dados (Y2009_2011). Este modelo ignora a estrutura dinâmica
# mas fornece uma estimativa rápida e informada.

cat("A ajustar GLM para obter valores iniciais de beta...\n")

# Define as dimensões do problema a partir dos dados já carregados
K <- 4 # Número de clusters de qualidade de dados
n_regions <- dim(E)[1]
n_times <- dim(E)[2]

# Ajusta o GLM. Usamos o offset log(E) para ser consistente com a teoria de modelos de contagem.
# Nota: Usamos a matriz 'E' original aqui, pois o GLM lida bem com a escala.
glm_fit <- glm(Y1999_2001 ~ x[, 1, 1] + x[, 1, 2] + x[, 1, 3] + offset(log(E[,1])), 
               family = poisson(link = "log"))

print(summary(glm_fit))

# Extrai o intercepto (nível base do log-risco) e os coeficientes de regressão
lambda0_log <- glm_fit$coefficients[1] 
beta_ini <- glm_fit$coefficients[2:4]
cat("\nValores iniciais para beta obtidos do GLM:\n"); print(beta_ini); cat("\n")

# Usa o intercepto para definir os hiperparâmetros (a0, b0) da distribuição Gamma
# para o estado inicial de lambda, usando o método dos momentos.
# Média = a / b  |  Variância = a / b^2
# Assumimos uma variância arbitrária, mas razoável, para lambda0.
mu_lambda0 <- exp(lambda0_log) # Converte o intercepto log para a escala original
var_lambda0 <- 2.0 # Variância inicial assumida para lambda

# CORREÇÃO: Fórmulas corretas do método dos momentos para Gamma(shape, rate)
b0 <- mu_lambda0 / var_lambda0  # rate
a0 <- mu_lambda0 * b0          # shape
cat(paste("Hiperparâmetros iniciais para a Gamma: a0 =", round(a0, 2), ", b0 =", round(b0, 2), "\n\n"))


# --- PASSO 3: DEFINIÇÃO DE PARÂMETROS FIXOS E SIMULAÇÃO DA TRAJETÓRIA ---

# Define valores iniciais para gamma e o fator de desconto w
gamma_ini <- c(0.05, 0.10, 0.10, 0.15)
w <- 0.9

# Calcula epsilon inicial a partir de gamma_ini
epsilon_ini <- 1 - hAI %*% gamma_ini

# Este bloco executa uma simulação forward-backward simplificada para gerar
# uma trajetória inicial coerente para lambda e Y, dados os parâmetros iniciais.
# Isto é muito mais eficaz do que simplesmente usar valores aleatórios.

cat("A simular trajetória inicial para lambda e Y via FFBS simplificado...\n")

# Inicialização de matrizes para o processo de simulação
lambda_at_fwd <- matrix(NA, ncol = n_times, nrow = n_regions)
y_at <- matrix(NA, nrow = n_regions, ncol = n_times + 1)
y_at[, 1] <- Y1999_2001 # Inicia com a média observada por região

att_ini <- matrix(NA, nrow = n_regions, ncol = n_times)
btt_ini <- matrix(NA, nrow = n_regions, ncol = n_times)
at_ini  <- matrix(NA, nrow = n_regions, ncol = n_times + 1)
bt_ini  <- matrix(NA, nrow = n_regions, ncol = n_times + 1)

# Simulação Forward
for(i in 1:n_regions) {
  at_ini[i, 1] <- a0
  bt_ini[i, 1] <- b0
  for(t in 2:(n_times + 1)) {
    # Passos de filtragem (forward)
    att_ini[i, t-1] <- w * at_ini[i, t-1]
    btt_ini[i, t-1] <- w * bt_ini[i, t-1]
    at_ini[i, t]    <- att_ini[i, t-1] + y_at[i, t-1]
    
    prod_val <- sum(x[i, t-1, ] * beta_ini)
    bt_ini[i, t]    <- btt_ini[i, t-1] + E_escalonado[i, t-1] * epsilon_ini[i] * exp(prod_val)
    
    # Amostra um valor de lambda a partir da distribuição filtrada (forward)
    lambda_at_fwd[i, t-1] <- rgamma(1, shape = at_ini[i, t], rate = bt_ini[i, t])
    
    # Simula um novo valor de Y para alimentar o próximo passo do filtro
    mu_at <- lambda_at_fwd[i, t-1] * exp(prod_val) * epsilon_ini[i] * E_escalonado[i, t-1]
    y_at[i, t] <- rpois(1, mu_at)
  }
}

# Simulação Backward (Suavização)
lambda_smooth <- matrix(NA, ncol = n_times, nrow = n_regions)
y_smooth <- matrix(NA, ncol = n_times, nrow = n_regions)

lambda_smooth[, n_times] <- lambda_at_fwd[, n_times] # O último ponto é o mesmo
y_smooth[, n_times] <- y_at[, n_times]

for(t in n_times:2) {
  # Amostra a "inovação" nu e calcula o lambda suavizado
  nu <- rgamma(n_regions, shape = (1-w) * at_ini[, t], rate = bt_ini[, t])
  lambda_smooth[, t-1] <- nu + w * lambda_smooth[, t]
  
  # Gera contagens Y consistentes com a trajetória suavizada de lambda
  prod_val_vec <- x[, t-1, ] %*% beta_ini
  mu_smooth <- lambda_smooth[, t-1] * exp(prod_val_vec) * epsilon_ini * E_escalonado[, t-1]
  y_smooth[, t-1] <- rpois(n_regions, mu_smooth)
}

# --- PASSO 4: ATRIBUIÇÃO FINAL DOS VALORES INICIAIS ---

# Os objetos finais 'y_ini' e 'lambda_ini' contêm as trajetórias completas
# e temporalmente consistentes que serão usadas para iniciar o MCMC.
y_ini <- y_at[,2:24]
lambda_ini <- lambda_at_fwd

# Valores iniciais para os priors uniformes de gamma
a_unif_ini <- c(0.00, rep(0, K - 1))
b_unif_ini <- c(0.1, rep(0.99, K - 1))

# Substitui a matriz E original pela versão escalonada para uso no modelo NIMBLE
E <- E_escalonado

cat("Geração de valores iniciais concluída com sucesso.\n")

# Limpeza de objetos intermédios (opcional)
rm(E_escalonado, glm_fit, lambda0_log, mu_lambda0, var_lambda0, a0, b0,
   epsilon_ini, lambda_at_fwd, y_at, att_ini, btt_ini, at_ini, bt_ini,
   lambda_smooth, y_smooth, prod_val_vec, mu_smooth, media_global_E)