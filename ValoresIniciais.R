#ajustar o offset globalmente
min_max_scale_matrix <- function(mat, a , b) {
  scaled_mat <- mat  # Cria uma cópia da matriz original
  
  for (j in 1:ncol(mat)) {
    min_val <- min(mat[, j])
    max_val <- max(mat[, j])
    # Aplica Min-Max Scaling para o intervalo [a, b]
    scaled_mat[, j] <- a + ((mat[, j] - min_val) / (max_val - min_val)) * (b - a)
  }
  
  return(scaled_mat)
}

# Exemplo de uso
# Supondo que 'my_matrix' seja a matriz que você deseja escalar
  # Exemplo de matriz
E <- min_max_scale_matrix(E,a=0.5,b=1.5)
##Matriz de lambdas
lambda_at = matrix(, ncol = 23, nrow = 75)
lambda_smooth = matrix(, ncol = 23, nrow = 75)
K = 4 #número de clusters
n_regions = dim(E)[1]  # Número de regiões
n_times = dim(E)[2]     # Número de períodos de tempo
#valores iniciais
y0 = Y1999_2001
glm <- glm(y0~1+x[,1,1]+x[,1,2]+x[,1,3],family = poisson(link = "log"))
summary(glm)
mu_lambda0 = glm$coefficients[1]
var_lambda0 = 2
#calculo a0 e b0 a partir do glm e variancia arbirtraria
b0 = mu_lambda0/var_lambda0
a0= mu_lambda0*(b0^2)
#Definindo valores iniciais de beta a partir do glm
beta_ini = glm$coefficients[2:4]
#número de covariaveis
p = length(beta_ini)
#matriz da média da poisson Y
mu_at <- matrix(nrow = n_regions, ncol = n_times)
mu_smooth <- matrix(nrow = n_regions, ncol = n_times)
y_smooth <- matrix(nrow = n_regions, ncol = n_times)
#matriz de Y
y_at <- matrix(,nrow = n_regions, ncol = n_times+1)
y_at[,1]<-y0
#Matriz de att_ini,btt_ini,at_ini,bt_ini
att_ini <- matrix(,nrow = n_regions, ncol = n_times)
btt_ini <- matrix(,nrow = n_regions, ncol = n_times)
at_ini  <- matrix(,nrow = n_regions, ncol = n_times+1)
bt_ini  <- matrix(,nrow = n_regions, ncol = n_times+1)

#w arbitrário
w = 0.9
# Definindo valores iniciais para os parâmetros gamma
gamma_ini = c(0.05, 0.10, 0.10, 0.15)
#calculo epsilons iniciais
epsilon_ini = 1 - hAI %*% gamma_ini
#CALCULO att_ini,btt_ini,at_ini,bt_ini
for(i in 1:n_regions) {
  at_ini[i,1] <- a0
  bt_ini[i,1] <- b0
  for(t in 2:(n_times+1)) {
    att_ini[i,t-1] <- w * at_ini[i,t-1]
    #cat("Processing region:", i, "at time:", t - 1, "with att:", att_ini[i,t-1], "\n")
    btt_ini[i,t-1] <- w * bt_ini[i,t-1] #/ (epsilon_ini[i] * E[i,t-1] * exp((x[i,t-1,1:p] %*% beta_ini)))
    #cat("Processing region:", i, "at time:", t - 1, "with btt:", btt_ini[i,t-1], "\n")
    at_ini[i,t] <- att_ini[i,t-1] + y_at[i,t-1]
    #cat("Processing region:", i, "at time:", t - 1, "with at:", at_ini[i,t], "\n")
    bt_ini[i,t] <- btt_ini[i,t-1] + (1) * exp((x[i,t-1,1:p] %*% beta_ini)) * epsilon_ini[i] * E[i,t-1]
    #cat("Processing region:", i, "at time:", t - 1, "with bt:", bt_ini[i,t], "\n")
    lambda_at[i,t-1] <- rgamma(1,shape = att_ini[i,t-1],rate = btt_ini[i,t-1])
    cat("Processing region:", i, "at time:", t - 1, "with lambda_at:", lambda_at[i,t-1], "\n")
    mu_at[i,t-1] <- lambda_at[i,t-1] * exp((x[i,t-1,1:p]%*%beta_ini)) * epsilon_ini[i]*E[i,t-1]
    #cat("Processing region:", i, "at time:", t - 1, "with mu:", mu_at[i,t-1], "\n")
    y_at[i,t] <- rpois(1,mu_at[i,t-1])
  }}
lambda_smooth[,n_times] <-lambda_at[,n_times]
y_smooth[,n_times] <- y_at[,n_times]
mu_smooth[,n_times] <- mu_at[,n_times]
for(t in n_times:2){
  lambda_smooth[,t-1] <- rgamma(75,shape = (1-w)*at_ini[,t-1],rate = bt_ini[,t-1])+w*lambda_smooth[,t]
  mu_smooth[,t-1] <- lambda_smooth[,t-1] * exp((x[,t-1,1:p]%*%beta_ini))*epsilon_ini*E[,t-1]
  y_smooth[,t-1] <- rpois(75,mu_smooth[,t-1])
}
y_ini = y_smooth
lambda_ini = lambda_smooth
a_unif_ini = c(0.00, rep(0, K - 1))  # Parâmetros a uniformes
b_unif_ini = c(0.1, rep(0.99, K - 1)) 
