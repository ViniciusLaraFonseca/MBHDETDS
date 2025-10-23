# --- PASSO 0: Carregar pacotes e limpar ambiente (se necessário) ---
if (!require(nimble)) { install.packages("nimble"); library(nimble) }
if (!require(coda)) { install.packages("coda"); library(coda) }

# Opcional: Limpar compilações anteriores do NIMBLE se houver problemas
# try(nimble::clearCompiled(), silent = TRUE)
# gc()

set.seed(123) # Para reprodutibilidade

# --- PASSO 1: Definir o código do modelo híbrido ---
print("Definindo o código do modelo NIMBLE (code_hibrido)...")
code_hibrido <- nimbleCode({
  # --- Priors (beta, gamma) ---
  for (j in 1:p) {
    beta[j] ~ dnorm(mu_beta[j], sd = 5) # Usando mu_beta de constants
  }
  # Usando a_unif e b_unif de constants
  gamma[1] ~ dunif(a_unif[1], b_unif[1])
  for (j in 2:K) {
    gamma[j] ~ dunif(
      min = a_unif[j] * (1 - sum(gamma[1:(j - 1)])),
      max = b_unif[j] * (1 - sum(gamma[1:(j - 1)]))
    )
  }
  for (i in 1:n_regions) {
    epsilon[i] <- 1 - sum(h[i, 1:K] * gamma[1:K])
  }
  
  # --- Lógica do Filtro Forward (declarada no nimbleCode) ---
  for (i in 1:n_regions) {
    # Inicialização t=1
    at[i, 1] <- a0
    bt[i, 1] <- b0
    
    # Evolução t >= 2
    for (t in 2:(n_times + 1)) {
      # Passos preditivos
      att[i, t-1] <- w * at[i, t-1]
      btt[i, t-1] <- w * bt[i, t-1] # Taxa preditiva (assumindo bt é taxa)
      
      # Atualização com observação Y[i, t-1]
      at[i, t] <- att[i, t-1] + Y[i, t-1] # Y são os dados
      prod_val_update <- inprod(beta[1:p], x[i, t-1, 1:p])
      # Atualiza taxa/scale bt
      bt[i, t] <- btt[i, t-1] + E[i, t-1] * epsilon[i] * exp(prod_val_update)
    }
    
    # --- Lambda Filtrado (Dummy prior) ---
    for (t in 1:n_times) {
      lambda_filtrado[i, t] ~ dgamma(shape = at[i, t+1], rate = bt[i, t+1]) # Dummy
    }
    
    # --- Lambda Suavizado (Nó alvo do sampler customizado) ---
    for (t in 1:n_times) {
      lambda_smooth[i, t] ~ dnorm(0, sd = 100) # Dummy
    }
    
    # --- Verossimilhança (depende de lambda_smooth) ---
    for(t in 1:n_times){
      # Cálculo intermediário explícito (opcional, pode ajudar na depuração)
      lambda_smooth_safe[i,t] <- max(lambda_smooth[i,t], 1e-10) # Evitar log(0)
      log_lambda_smooth[i,t] <- log(lambda_smooth_safe[i,t])
      inprod_beta_x[i,t] <- inprod(beta[1:p], x[i, t, 1:p])
      
      log_theta[i, t] <- log_lambda_smooth[i,t] + inprod_beta_x[i,t]
      theta[i, t] <- exp(log_theta[i, t])
      mu[i,t] <- E[i, t] * epsilon[i] * theta[i, t]
      # Adicionar verificação para mu > 0
      mu_safe[i,t] <- max(mu[i,t], 1e-10)
      Y[i, t] ~ dpois(mu_safe[i,t])
    }
  }
})
print("Código do modelo definido.")

# --- PASSO 2: Definir o sampler customizado para Backward Sampling ---
print("Definindo o sampler customizado (backward_sampler)...")
backward_sampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    n_regions <- control$n_regions
    n_times   <- control$n_times
    w         <- control$w
  },
  
  run = function() {
    for (i in 1:n_regions) {
      # Amostrar último estado
      a_T_plus_1 <- model$at[i, n_times + 1]
      b_T_plus_1 <- model$bt[i, n_times + 1]
      # Precaução numérica
      shape_T <- max(a_T_plus_1, 1e-10)
      rate_T  <- max(b_T_plus_1, 1e-10)
      model$lambda_smooth[i, n_times] <<- rgamma(1, shape = shape_T, rate = rate_T)
    }
    
    # Loop Backward (t = n_times-1 ... 1)
    if(n_times > 1) { # Evitar loop se n_times = 1
      for (t in (n_times - 1):1) {
        for (i in 1:n_regions) {
          at_t_plus_1 <- model$at[i, t + 1]
          bt_t_plus_1 <- model$bt[i, t + 1]
          # Parâmetros para nu
          shape_nu <- (1 - w) * at_t_plus_1
          rate_nu  <- bt_t_plus_1
          
          # Precaução numérica
          shape_nu <- max(shape_nu, 1e-10)
          rate_nu  <- max(rate_nu, 1e-10)
          
          nu <- rgamma(1, shape = shape_nu, rate = rate_nu)
          
          lambda_futuro <- model$lambda_smooth[i, t + 1]
          model$lambda_smooth[i, t] <<- nu + w * lambda_futuro
        }
      }
    }
    
    # Calcular dependências e copiar
    model$calculate(calcNodes) # Calcula Y, mu, theta, etc. que dependem de lambda_smooth
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)
print("Sampler customizado definido.")

# --- PASSO 3: Criar constantes, dados e inicializações (Adaptado de debug_script_v2.R) ---
print("Definindo constantes, dados e inicializações...")
# Parâmetros do cenário de depuração
n_regions <- 1
n_times   <- 20
p         <- 3 # Número de covariáveis beta
K         <- 4 # Número de clusters gamma

# Valores 'verdadeiros' ou iniciais do script de debug
beta_ini  <- c(0.5, 1, -0.3)
gamma_ini <- c(0.05, 0.10, 0.10, 0.15)
w         <- 0.9
a0        <- 1.0
b0        <- 1.0
a_unif_ini = c(0.00, rep(0, K - 1)) # Exemplo de ValoresIniciais.R
b_unif_ini = c(0.1, rep(0.99, K - 1))# Exemplo de ValoresIniciais.R

# Simular dados E, x, Y (como no debug_script_v2.R)
set.seed(123)
E_debug_at <- matrix(runif(n_regions * (n_times + 1), 100, 200), nrow = n_regions)
# IMPORTANTE: E e x devem ter dimensão n_times, não n_times+1
E_full <- matrix(E_debug_at[, 1:n_times], nrow=n_regions)
x_debug_at <- array(rnorm(n_regions * (n_times + 1) * p), dim = c(n_regions, (n_times + 1), p))
# IMPORTANTE: x deve ter dimensão n_times
x_full <- array(x_debug_at[, 1:n_times, ], dim = c(n_regions, n_times, p))
hAI <- matrix(c(1, 0, 0, 0), nrow = n_regions)

# Simular Y_ini (usando a lógica forward do debug_script_v2 para consistência)
epsilon_fixed <- 1 - as.numeric(hAI %*% gamma_ini)
# lambda_at_fwd_true <- matrix(NA, ncol = n_times, nrow = n_regions) # Lambda simulado (não usado diretamente)
y_at <- matrix(NA, nrow = n_regions, ncol = n_times + 1)
y_at[, 1] <- 5 # Valor inicial arbitrário razoável para iniciar a simulação de Y
att_sim <- matrix(NA, nrow = n_regions, ncol = n_times)
btt_sim <- matrix(NA, nrow = n_regions, ncol = n_times)
at_sim  <- matrix(NA, nrow = n_regions, ncol = n_times + 1)
bt_sim  <- matrix(NA, nrow = n_regions, ncol = n_times + 1)
lambda_sim_temp <- matrix(NA, nrow = n_regions, ncol = n_times)

for(i in 1:n_regions) {
  at_sim[i, 1] <- a0
  bt_sim[i, 1] <- b0
  for(t_sim in 1:n_times) {
    # Filtragem
    att_sim[i, t_sim] <- w * at_sim[i, t_sim]
    btt_sim[i, t_sim] <- w * bt_sim[i, t_sim]
    at_sim[i, t_sim + 1] <- att_sim[i, t_sim] + y_at[i, t_sim]
    prod_val_sim <- sum(x_full[i, t_sim, ] * beta_ini) # Usar x_full
    bt_sim[i, t_sim + 1] <- btt_sim[i, t_sim] + E_full[i, t_sim] * epsilon_fixed[i] * exp(prod_val_sim) # Usar E_full
    
    # Amostrar lambda da filtrada para simular Y seguinte
    shape_sim = max(at_sim[i, t_sim + 1], 1e-10)
    rate_sim = max(bt_sim[i, t_sim + 1], 1e-10)
    lambda_sim_temp[i, t_sim] <- rgamma(1, shape = shape_sim, rate = rate_sim)
    
    # Simular Y seguinte
    mu_sim <- lambda_sim_temp[i, t_sim] * exp(prod_val_sim) * epsilon_fixed[i] * E_full[i, t_sim]
    y_at[i, t_sim + 1] <- rpois(1, max(mu_sim, 1e-10)) # Usar Y_{t+1} para próximo passo
  }
}
Y_ini <- round(y_at[, 2:(n_times + 1), drop = FALSE]) # Usar Y maiúsculo para dados, pegar t=1 a n_times

# Definir lista de constantes
constants <- list(
  n_regions = n_regions,
  n_times   = n_times,
  p         = p,
  K         = K,
  a0        = a0,
  b0        = b0,
  w         = w,
  h         = hAI,
  mu_beta   = rep(0, p), # Prior N(0, sd=5) para beta
  a_unif    = a_unif_ini, # Para prior gamma
  b_unif    = b_unif_ini  # Para prior gamma
)

# Definir lista de dados
data <- list(
  Y = Y_ini,
  E = E_full,
  x = x_full
)

# Definir lista de inicializações
lambda_init_guess = matrix(1, nrow = n_regions, ncol = n_times) # Suposição inicial
inits <- list(
  beta = beta_ini,
  gamma = gamma_ini,
  lambda_smooth = lambda_init_guess, # Inicializar o nó suavizado
  lambda_filtrado = lambda_init_guess # Inicializar o nó filtrado (dummy)
  # NIMBLE inicializará 'at', 'bt', 'att', 'btt' deterministicamente se possível,
  # mas pode ser bom fornecer inits se houver problemas.
  # at = matrix(a0, nrow=n_regions, ncol=n_times+1), # Exemplo
  # bt = matrix(b0, nrow=n_regions, ncol=n_times+1)  # Exemplo
)
print("Constantes, dados e inicializações definidos.")

# --- PASSO 4: Criar, configurar e rodar o MCMC ---

# Criar a instância do modelo ANTES de compilar
print("Criando a instância do modelo NIMBLE...")
model <- nimbleModel(code_hibrido,
                     constants = constants,
                     data = data,
                     inits = inits,
                     check = FALSE)
print("Instância do modelo criada.")

# Compilar o modelo primeiro
print("Compilando o modelo...")
# Pode demorar um pouco
Cmodel <- compileNimble(model, resetFunctions = TRUE, showCompilerOutput = FALSE) # showCompilerOutput = TRUE para depurar compilação
print("Modelo compilado.")

# Configurar MCMC
print("Configurando MCMC...")
conf <- configureMCMC(Cmodel, print = FALSE) # Usar Cmodel

# Remover samplers padrão desnecessários ou que serão substituídos
print("Modificando samplers...")
conf$removeSamplers('lambda_smooth')
#conf$removeSamplers('lambda_filtrado')
# Os nós 'at', 'bt', 'att', 'btt' são determinísticos ou intermediários no filtro.
# O NIMBLE pode não atribuir samplers a eles, mas se atribuir, remova-os.
# Verifique com conf$printSamplers() antes de remover se necessário.
# conf$removeSamplers(c('at', 'bt', 'att', 'btt'))


# Adicionar nosso sampler customizado para lambda_smooth
conf$addSampler(target = 'lambda_smooth',
                type = backward_sampler,
                control = list(n_regions = constants$n_regions,
                               n_times = constants$n_times,
                               w = constants$w))
print("Configuração final dos samplers:")
conf$printSamplers(c('beta', 'gamma', 'lambda_smooth'))


# Adicionar monitores
# Adicionar monitores
# conf$monitors <- c('beta', 'gamma', 'lambda_smooth') # Linha original
conf$monitors <- c('beta', 'gamma', 'lambda_smooth',
                   'lambda_filtrado', 'at', 'bt', 'att', 'btt') # Nova linha

# Construir e compilar o MCMC
print("Construindo MCMC...")
Rmcmc <- buildMCMC(conf)
print("Compilando MCMC...")
# Pode demorar um pouco
Cmcmc <- compileNimble(Rmcmc, project = Cmodel, resetFunctions = TRUE, showCompilerOutput = FALSE)
print("MCMC compilado.")

# Rodar o MCMC
print("Iniciando execução do MCMC...")
# Aumentar niter para uma execução real, 1000 é só para teste rápido
n_iter <- 5000
n_burnin <- 1000
samples <- runMCMC(Cmcmc, niter = n_iter, nburnin = n_burnin, nchains = 1, summary = TRUE, samplesAsCodaMCMC = TRUE)
print("Execução do MCMC concluída.")

# --- PASSO 5: Analisar os resultados ---
print("Sumário das amostras:")
print(samples$summary)

# Verificar traceplots (requer pacote coda)
print("Gerando traceplots...")
# plot(samples$samples)

# Exemplo: plotar traceplot para beta[1] e lambda_smooth[1,1]
if (require(coda)) {
  samples_mcmc <- samples$samples
  par(mfrow=c(2,1))
  coda::traceplot(samples_mcmc[, "beta[1]"], main="Traceplot beta[1]")
  coda::traceplot(samples_mcmc[, "lambda_smooth[1, 1]"], main="Traceplot lambda_smooth[1,1]")
  par(mfrow=c(1,1))
} else {
  print("Pacote 'coda' não instalado. Não foi possível gerar traceplots.")
}

print("Script concluído.")