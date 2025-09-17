# Definindo o número de Klusters (K)
K = 4

# Estimando os coeficientes iniciais usando um modelo GLM (Modelo Linear Generalizado)
# O modelo utiliza a variável resposta Y2009_2011 e preditores selecionados da matriz x.
n_regions = dim(E)[1]  # Número de regiões
n_times = dim(E)[2]     # Número de períodos de tempo

# Ajustar o modelo GLM com a função de distribuição Poisson e link log
glm <- glm(Y2009_2011 ~ 1 + x[, 10, 1] + x[, 10, 2] + x[, 10, 3], family = poisson(link = "log"))

# Exibir o resumo do modelo ajustado
summary(glm)

# Extraindo o coeficiente intercepto (lambda0) do modelo
lambda0 = glm$coefficients[1]
print(lambda0)
log(lambda0)
# Extraindo os coeficientes beta iniciais (para as variáveis preditoras)
beta_ini = glm$coefficients[2:4]
print(beta_ini)

# Definindo valores iniciais para os parâmetros a e b
a0 = 1  # Parâmetro a inicial
b0 = 1  # Parâmetro b inicial

# Definindo valores iniciais para os parâmetros gamma
gamma_ini = c(0.05, 0.10, 0.10, 0.15)

# Inicializando os parâmetros a e b para uma distribuição uniforme
a_unif_ini = c(0.00, rep(0, K - 1))  # Parâmetros a uniformes
b_unif_ini = c(0.1, rep(1, K - 1))   # Parâmetros b uniformes

# Estimando epsilon com base em um fator hAI e os parâmetros gamma
epsilon_ini = 1 - hAI %*% gamma_ini
# Expandindo epsilon para uma matriz com 75 regiões e 23 períodos
epsilon_ini = matrix(rep(epsilon_ini, 23), nrow = 75, ncol = 23)

# Estimando lambda inicial utilizando a distribuição gamma
lambda_ini = matrix(rgamma(1725, shape = a0, rate = b0), ncol = 23, nrow = 75)

# Estimando Y com base nos parâmetros lambda, E e epsilon, aplicando a função exponencial nos coeficientes beta
mu_ini = lambda_ini * E * epsilon_ini * exp(beta_ini[1] * x[,, 1] + beta_ini[2] * x[,, 2] + beta_ini[3] * x[,, 3])
y_ini = matrix(rpois(1725,mu_ini),nrow = n_regions,ncol = n_times)
#estimando o theta com os betas iniciais e o x
theta_ini = exp(beta_ini[1] * x[,, 1] + beta_ini[2] * x[,, 2] + beta_ini[3] * x[,, 3])
# Definindo um parâmetro w para ser utilizado em cálculos posteriores

w=0.7
