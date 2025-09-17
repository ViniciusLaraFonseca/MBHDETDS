if(!require(dplyr)){ install.packages("dplyr"); require(dplyr)}
if(!require(tidyverse)){ install.packages("tidyverse"); require(tidyverse)}
if(!require(readxl)){ install.packages("readxl"); require(readxl)} 
if(!require(stringr)){ install.packages("stringr"); require(stringr)} 
if(!require(geobr)){ install.packages("geobr"); require(geobr)} 
if(!require(RColorBrewer)){ install.packages("RColorBrewer"); require(RColorBrewer)} 
if(!require(ggplot2)){ install.packages("ggplot2"); require(ggplot2)} 
if(!require(fda)){ install.packages("fda"); require(fda)} 
if(!require(dplyr)){ install.packages("dplyr"); require(dplyr)}   
if(!require(maps)){ install.packages("maps"); require(maps)}   # For adding a map to the plots of Brazil
if(!require(NbClust)){ install.packages("NbClust"); require(NbClust)}   
if(!require(ppclust)){ install.packages("ppclust"); require(ppclust)}   
if(!require(psych)){ install.packages("psych"); require(psych)}   
if(!require(nimble)){ install.packages("nimble"); require(nimble)} # For MCMC computation using NIMBLE.  
if(!require(spdep)){ install.packages("spdep"); require(spdep)}   # For computing the neighbourhood and adjancency objects.
if(!require(coda)){ install.packages("coda"); require(coda)}   
if(!require(viridis)){ install.packages("viridis"); require(viridis)}  # Colour pallette for plots. 
if(!require(ggplot2)){ install.packages("ggplot2"); require(ggplot2)}   
#plotar mapa
#juntar com a microrgião do IA, ordenar e dar seguimento
#DTB = read_xls("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/PesquisaMestrado/RELATORIO_DTB_BRASIL_MUNICIPIO.xls")
#DTB =DTB%>%
#  mutate(COD_Municipio_Reduzido = as.character(COD_Municipio_Reduzido))

#DTB%>%filter(COD_UF =="31")%>%
#  count(COD_Microrregiao_Geografica)

# Ler os dados de códigos de município
cod_micro <- read.table("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/PesquisaMestrado/DATASUS/dados_codMICRO.txt", header = TRUE)

# Renomear a coluna de município e converter para tipo caracter
cod_micro <- cod_micro %>%
  rename(COD_Municipio_Reduzido = COD_MUN) %>%
  mutate(COD_Municipio_Reduzido = as.character(COD_Municipio_Reduzido))

# Importar dados de nascimentos vivos
total = read.table("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/PesquisaMestrado/DATASUS/nascidos_vivos_total.csv", sep = ";", header = TRUE)

# Remover linhas e colunas desnecessárias, renomear e transformar dados
total = total[-5596, -25] %>%
  rename(COD_Municipio_Reduzido = Municipio) %>%
  pivot_longer(cols = 2:24, names_to = "Ano", values_to = "Nascidos_vivos") %>%
  mutate(Ano = str_sub(Ano, 2),  # Remover o primeiro caractere do ano
         Nascidos_vivos = as.integer(Nascidos_vivos), 
         COD_Municipio_Reduzido = as.character(COD_Municipio_Reduzido))

# COVARIÁVEL: PROPORÇÃO DE PRÉ-NATAL INADEQUADO (MENOS DE 5 CONSULTAS)
consultas_ignoradas = read.table("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/PesquisaMestrado/DATASUS/ignorados_consultas.csv", sep = ";", header = TRUE)

# Processar dados de consultas ignoradas
consultas_ignoradas = suppressWarnings({
  consultas_ignoradas[1:5596, -25] %>%
    rename(COD_Municipio_Reduzido = Municipio) %>%
    mutate_if(is.character, as.integer) %>%
    mutate(COD_Municipio_Reduzido = as.character(COD_Municipio_Reduzido)) %>%
    pivot_longer(cols = 2:24, names_to = "Ano", values_to = "Ignorados") %>%
    mutate(Ignorados = ifelse(is.na(Ignorados), 0, Ignorados)) %>%
    filter(COD_Municipio_Reduzido != "NA") %>%
    mutate(Ano = str_sub(Ano, 2)) %>%
    inner_join(total) %>%
    mutate(Nascidos_vivos_corrigido = as.integer(Nascidos_vivos) - Ignorados)
})

# Importar dados de consultas realizadas
consultas = read.table("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/PesquisaMestrado/DATASUS/nascidos_vivos_consultas.csv", sep = ";", header = TRUE)

# Processar dados de consultas
m_consultas = suppressWarnings({
  consultas[, -25] %>%
    rename(COD_Municipio_Reduzido = Municipio) %>%
    mutate_if(is.character, as.integer) %>%
    mutate(COD_Municipio_Reduzido = as.character(COD_Municipio_Reduzido)) %>%
    pivot_longer(cols = 2:24, names_to = "Ano", values_to = "Nascidos_vivos_com_ate_6_consultas") %>%
    mutate(Nascidos_vivos_com_ate_6_consultas = ifelse(is.na(Nascidos_vivos_com_ate_6_consultas), 0, Nascidos_vivos_com_ate_6_consultas)) %>%
    filter(COD_Municipio_Reduzido != "NA") %>%
    mutate(Ano = str_sub(Ano, 2)) %>%
    inner_join(consultas_ignoradas) %>%
    mutate(taxa_consultas = Nascidos_vivos_com_ate_6_consultas / Nascidos_vivos_corrigido) %>%
    inner_join(cod_micro) %>%
    group_by(MICRO_, Ano) %>%
    mutate(Sum_Nascidos_vivos_com_ate_6_consultas = sum(Nascidos_vivos_com_ate_6_consultas),
           SumNascidos_vivos = sum(Nascidos_vivos_corrigido)) %>%
    mutate(taxa_consultas_microrregião = Sum_Nascidos_vivos_com_ate_6_consultas / SumNascidos_vivos) %>%
    arrange(COD_MICRO_VAR) %>%
    ungroup()
})

# COVARIÁVEL: PROPORÇÃO DE MÃES COM ENSINO SUPERIOR
instrucao_ignoradas = read.table("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/PesquisaMestrado/DATASUS/ignorados_instrucao.csv", sep = ";", header = TRUE)

# Processar dados de instrução ignorada
instrucao_ignoradas = suppressWarnings({
  instrucao_ignoradas[, -25] %>%
    rename(COD_Municipio_Reduzido = Municipio) %>%
    mutate_if(is.character, as.integer) %>%
    mutate(COD_Municipio_Reduzido = as.character(COD_Municipio_Reduzido)) %>%
    pivot_longer(cols = 2:24, names_to = "Ano", values_to = "Ignorados") %>%
    mutate(Ignorados = ifelse(is.na(Ignorados), 0, Ignorados)) %>%
    filter(COD_Municipio_Reduzido != "NA") %>%
    mutate(Ano = str_sub(Ano, 2)) %>%
    inner_join(total) %>%
    mutate(Nascidos_vivos_corrigido = as.integer(Nascidos_vivos) - Ignorados)
})

# Importar dados de instrução
instrucao = read.table("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/PesquisaMestrado/DATASUS/nascidos_vivos_instrução2.csv", sep = ";", header = TRUE)

# Processar dados de instrução
m_instrucao = suppressWarnings({
  instrucao[, -25] %>%
    rename(COD_Municipio_Reduzido = Municipio) %>%
    mutate_if(is.character, as.integer) %>%
    mutate(COD_Municipio_Reduzido = as.character(COD_Municipio_Reduzido)) %>%
    pivot_longer(cols = 2:24, names_to = "Ano", values_to = "Nascidos_vivos_com_grau_instrucao_mae") %>%
    mutate(Nascidos_vivos_com_grau_instrucao_mae = ifelse(is.na(Nascidos_vivos_com_grau_instrucao_mae), 0, Nascidos_vivos_com_grau_instrucao_mae)) %>%
    filter(COD_Municipio_Reduzido != "NA") %>%
    mutate(Ano = str_sub(Ano, 2)) %>%
    inner_join(instrucao_ignoradas) %>%
    mutate(taxa_grau0 = Nascidos_vivos_com_grau_instrucao_mae / Nascidos_vivos_corrigido) %>%
    inner_join(cod_micro) %>%
    group_by(MICRO_, Ano) %>%
    mutate(Sum_Nascidos_vivos_com_grau_instrucao_mae = sum(Nascidos_vivos_com_grau_instrucao_mae),
           SumNascidos_vivos = sum(Nascidos_vivos_corrigido)) %>%
    mutate(taxa_instrucao_microrregião = Sum_Nascidos_vivos_com_grau_instrucao_mae / SumNascidos_vivos) %>%
    arrange(COD_MICRO_VAR) %>%
    ungroup()
})

# COVARIÁVEL: PROPORÇÃO DE NASCIDOS VIVOS COM PESO ATÉ 2500g
subnutridos_ignoradas = read.table("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/PesquisaMestrado/DATASUS/ignorados_peso.csv", sep = ";", header = TRUE)

# Processar dados de peso ignorados
subnutridos_ignoradas = suppressWarnings({
  subnutridos_ignoradas[, -25] %>%
    rename(COD_Municipio_Reduzido = Municipio) %>%
    mutate_if(is.character, as.integer) %>%
    mutate(COD_Municipio_Reduzido = as.character(COD_Municipio_Reduzido)) %>%
    pivot_longer(cols = 2:24, names_to = "Ano", values_to = "Ignorados") %>%
    mutate(Ignorados = ifelse(is.na(Ignorados), 0, Ignorados)) %>%
    filter(COD_Municipio_Reduzido != "NA") %>%
    mutate(Ano = str_sub(Ano, 2)) %>%
    inner_join(total) %>%
    mutate(Nascidos_vivos_corrigido = as.integer(Nascidos_vivos) - Ignorados)
})

# Importar dados de nascidos vivos com peso até 2500g
subnutridos = read.table("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/PesquisaMestrado/DATASUS/nascidos_vivos_menos_2500g.csv", sep = ";", header = TRUE)

# Processar dados de subnutridos
m_subnutridos = suppressWarnings({
  subnutridos[, -25] %>%
    rename(COD_Municipio_Reduzido = Municipio) %>%
    mutate_if(is.character, as.integer) %>%
    mutate(COD_Municipio_Reduzido = as.character(COD_Municipio_Reduzido)) %>%
    pivot_longer(cols = 2:24, names_to = "Ano", values_to = "Nascidos_vivos_subnutridos") %>%
    mutate(Nascidos_vivos_subnutridos = ifelse(is.na(Nascidos_vivos_subnutridos), 0, Nascidos_vivos_subnutridos)) %>%
    filter(COD_Municipio_Reduzido != "NA") %>%
    mutate(Ano = str_sub(Ano, 2)) %>%
    inner_join(subnutridos_ignoradas) %>%
    mutate(taxa_subnutrição = Nascidos_vivos_subnutridos / Nascidos_vivos_corrigido) %>%
    inner_join(cod_micro) %>%
    group_by(MICRO_, Ano) %>%
    mutate(Sum_Nascidos_vivos_subnutridos = sum(Nascidos_vivos_subnutridos),
           SumNascidos_vivos = sum(Nascidos_vivos_corrigido)) %>%
    mutate(taxa_subnutrição_microrregiao = Sum_Nascidos_vivos_subnutridos / SumNascidos_vivos) %>%
    arrange(COD_MICRO_VAR) %>%
    ungroup()
})

# Agrupar dados de consultas por microrregião
grouped_consultas_MG = m_consultas %>%
  distinct(MICRO_, Ano, .keep_all = TRUE) %>%
  dplyr::select(MICRO_, Ano, taxa_consultas_microrregião) %>%
  pivot_wider(names_from = Ano, values_from = taxa_consultas_microrregião)

# Agrupar dados de instrução por microrregião
grouped_instrucao_MG = m_instrucao %>%
  distinct(MICRO_, Ano, .keep_all = TRUE) %>%
  dplyr::select(MICRO_, Ano, taxa_instrucao_microrregião) %>%
  pivot_wider(names_from = Ano, values_from = taxa_instrucao_microrregião)

# Agrupar dados de subnutridos por microrregião
grouped_subnutridos_MG = m_subnutridos %>%
  distinct(MICRO_, Ano, .keep_all = TRUE) %>%
  dplyr::select(MICRO_, Ano, taxa_subnutrição_microrregiao) %>%
  pivot_wider(names_from = Ano, values_from = taxa_subnutrição_microrregiao)

# Criar um array para armazenar os dados
dim = c(75, 23, 3)  # Dimensões do array
x = array(data = NA, dim = dim, dimnames = list(grouped_consultas_MG$MICRO_, unique(m_subnutridos$Ano)))

# Preencher o array com dados de consultas, instrução e subnutrição
x[,,1] <- as.matrix(grouped_consultas_MG[, 2:24])  # Pré-natal inadequado
x[,,2] <- as.matrix(grouped_instrucao_MG[, 2:24])  # Baixa escolaridade
x[,,3] <- as.matrix(grouped_subnutridos_MG[, 2:24])  # Peso até suficiente

# Plotar os dados
plot(x[1,,1])  # Gráfico de pré-natal inadequado
plot(x[1,,2])  # Gráfico de baixa escolaridade
plot(x[1,,3])  # Gráfico de peso até suficiente

# Agrupar dados totais de nascimentos
Grouped_total_MG <- total %>%
  inner_join(cod_micro) %>%
  filter(COD_Municipio_Reduzido != "NA") %>%
  arrange(COD_MICRO_VAR) %>%
  group_by(MICRO_, Ano) %>%
  mutate(SumNascidos = sum(Nascidos_vivos)) %>%
  distinct(MICRO_, Ano, .keep_all = TRUE) %>%
  dplyr::select(MICRO_, Ano, SumNascidos) %>%
  pivot_wider(names_from = Ano, values_from = SumNascidos)

# Criar uma matriz a partir dos dados agrupados
E <- as.matrix(Grouped_total_MG[, 2:24])
row.names(E) <- Grouped_total_MG$MICRO_

# Definir variáveis únicas para área e tempo
area = unique(cod_micro$MICRO_)
tempo = unique(total$Ano)

#padronizando as covariáveis
x[,,1] <-scale(x[,,1])
x[,,2] <-scale(x[,,2])
x[,,3] <-scale(x[,,3])
