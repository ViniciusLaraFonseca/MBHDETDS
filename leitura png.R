# Script para organizar e consolidar gráficos de sensibilidade
# Data: 09/12/2025

# Carregar bibliotecas necessárias
library(stringr)
library(fs)

# ==============================================================================
# CONFIGURAÇÕES
# ==============================================================================

# Defina aqui o diretório base onde estão as pastas resultados_sensibilidade_w09 e w07
diretorio_base <- "C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/MBHDETDS"
# Alternativamente, use getwd() se já estiver no diretório correto
# diretorio_base <- getwd()

# Nome da pasta de destino (será criada no diretório base)
pasta_destino <- file.path(diretorio_base, "graficos_sensibilidade_09122025")

# Lista das pastas de origem
pastas_origem <- c(
  file.path(diretorio_base, "resultados_sensibilidade_w09"),
  file.path(diretorio_base, "resultados_sensibilidade_w07")
)

# ==============================================================================
# FUNÇÕES AUXILIARES
# ==============================================================================

# Função para extrair o w da pasta (09 ou 07)
extrair_w_da_pasta <- function(caminho_pasta) {
  # Extrai "w09" ou "w07" do nome da pasta
  if (grepl("w09", basename(caminho_pasta), ignore.case = TRUE)) {
    return("w09")
  } else if (grepl("w07", basename(caminho_pasta), ignore.case = TRUE)) {
    return("w07")
  }
  return(NA)
}

# Função para processar arquivos de uma subpasta
processar_subpasta <- function(caminho_subpasta, w_simulacao) {
  # Extrair nome da subpasta (ex: "C1_T100_w09")
  nome_subpasta <- basename(caminho_subpasta)
  
  # Caminhos para as subpastas internas
  caminho_lambdas <- file.path(caminho_subpasta, "lambdas")
  caminho_traceplots <- file.path(caminho_subpasta, "traceplots")
  
  # Lista para armazenar os arquivos processados
  arquivos_processados <- list()
  
  # ============================================================================
  # PROCESSAR ARQUIVOS EM "lambdas"
  # ============================================================================
  if (dir.exists(caminho_lambdas)) {
    arquivos_lambdas <- list.files(caminho_lambdas, pattern = "\\.png$", full.names = TRUE)
    
    for (arquivo in arquivos_lambdas) {
      nome_arquivo <- basename(arquivo)
      
      # Renomear arquivos de lambdas
      # Ex: painel_C1_T100_w09.png -> painel_w09_C1_T100_w09.png
      novo_nome <- gsub(
        "^([a-zA-Z_]+)_(C\\d+_T\\d+_w\\d{2})\\.png$",
        paste0("\\1_", w_simulacao, "_\\2.png"),
        nome_arquivo
      )
      
      # Verificar se a renomeação foi bem-sucedida
      if (novo_nome != nome_arquivo) {
        arquivos_processados[[length(arquivos_processados) + 1]] <- list(
          origem = arquivo,
          destino_nome = novo_nome,
          tipo = "lambdas"
        )
      } else {
        cat(sprintf("    Aviso: Não foi possível renomear %s\n", nome_arquivo))
      }
    }
  } else {
    cat(sprintf("    Aviso: Pasta 'lambdas' não encontrada em %s\n", nome_subpasta))
  }
  
  # ============================================================================
  # PROCESSAR ARQUIVOS EM "traceplots"
  # ============================================================================
  if (dir.exists(caminho_traceplots)) {
    arquivos_traceplots <- list.files(caminho_traceplots, pattern = "\\.png$", full.names = TRUE)
    
    for (arquivo in arquivos_traceplots) {
      nome_arquivo <- basename(arquivo)
      
      # Renomear arquivos de traceplots
      # Ex: traceplot_gamma_4_.png -> traceplot_gamma_4_w09_C1_T100_w09.png
      # Primeiro, extrair o nome base (sem extensão)
      nome_base <- sub("\\.png$", "", nome_arquivo)
      
      # Remover underscore final se existir
      nome_base <- sub("_$", "", nome_base)
      
      # Construir novo nome
      novo_nome <- paste0(nome_base, "_", w_simulacao, "_", nome_subpasta, ".png")
      
      arquivos_processados[[length(arquivos_processados) + 1]] <- list(
        origem = arquivo,
        destino_nome = novo_nome,
        tipo = "traceplots"
      )
    }
  } else {
    cat(sprintf("    Aviso: Pasta 'traceplots' não encontrada em %s\n", nome_subpasta))
  }
  
  return(arquivos_processados)
}

# ==============================================================================
# FUNÇÃO PRINCIPAL
# ==============================================================================

organizar_graficos_sensibilidade <- function() {
  cat("Iniciando organização dos gráficos de sensibilidade...\n")
  cat("Data: 09/12/2025\n\n")
  
  # Criar pasta de destino se não existir
  if (!dir.exists(pasta_destino)) {
    dir.create(pasta_destino, recursive = TRUE)
    cat(sprintf("Pasta de destino criada: %s\n\n", pasta_destino))
  }
  
  # Contadores para estatísticas
  total_processados <- 0
  total_copiados <- 0
  erros <- 0
  
  # Processar cada pasta de origem
  for (pasta_origem in pastas_origem) {
    if (!dir.exists(pasta_origem)) {
      cat(sprintf("Aviso: Pasta não encontrada: %s\n", pasta_origem))
      next
    }
    
    # Extrair w da simulação (w09 ou w07)
    w_simulacao <- extrair_w_da_pasta(pasta_origem)
    if (is.na(w_simulacao)) {
      cat(sprintf("Aviso: Não foi possível extrair w da pasta: %s\n", pasta_origem))
      next
    }
    
    cat(sprintf("Processando pasta: %s (w = %s)\n", basename(pasta_origem), w_simulacao))
    
    # Listar subpastas (C1_T100_w09, C2_T100_w07, etc.)
    subpastas <- list.dirs(pasta_origem, recursive = FALSE, full.names = TRUE)
    
    for (subpasta in subpastas) {
      nome_subpasta <- basename(subpasta)
      cat(sprintf("  Processando cenário: %s\n", nome_subpasta))
      
      # Processar arquivos na subpasta
      arquivos <- processar_subpasta(subpasta, w_simulacao)
      
      # Copiar arquivos renomeados para a pasta de destino
      for (arquivo_info in arquivos) {
        total_processados <- total_processados + 1
        
        # Verificar se o arquivo de origem existe
        if (!file.exists(arquivo_info$origem)) {
          cat(sprintf("    Erro: Arquivo de origem não encontrado: %s\n", arquivo_info$origem))
          erros <- erros + 1
          next
        }
        
        # Criar caminho de destino
        caminho_destino <- file.path(pasta_destino, arquivo_info$destino_nome)
        
        tryCatch({
          # Copiar arquivo
          file.copy(arquivo_info$origem, caminho_destino, overwrite = FALSE)
          cat(sprintf("    Copiado: %s -> %s\n", basename(arquivo_info$origem), arquivo_info$destino_nome))
          total_copiados <- total_copiados + 1
        }, error = function(e) {
          cat(sprintf("    Erro ao copiar %s: %s\n", 
                      arquivo_info$destino_nome, e$message))
          erros <- erros + 1
        })
      }
    }
    cat("\n")
  }
  
  # ============================================================================
  # RESUMO DA EXECUÇÃO
  # ============================================================================
  linha_separadora <- paste0(rep("=", 60), collapse = "")
  cat(linha_separadora, "\n")
  cat("RESUMO DA ORGANIZAÇÃO\n")
  cat(linha_separadora, "\n")
  cat(sprintf("Pasta de destino: %s\n", pasta_destino))
  cat(sprintf("Arquivos processados: %d\n", total_processados))
  cat(sprintf("Arquivos copiados com sucesso: %d\n", total_copiados))
  cat(sprintf("Erros durante a cópia: %d\n", erros))
  
  # Listar arquivos na pasta de destino
  if (dir.exists(pasta_destino)) {
    arquivos_destino <- list.files(pasta_destino, pattern = "\\.png$", full.names = FALSE)
    cat(sprintf("\nArquivos na pasta de destino (%d arquivos):\n", length(arquivos_destino)))
    
    if (length(arquivos_destino) > 0) {
      for (i in seq_along(arquivos_destino)) {
        cat(sprintf("%3d. %s\n", i, arquivos_destino[i]))
      }
    }
  }
  
  cat("\nProcesso concluído!\n")
}

# ==============================================================================
# EXECUTAR O SCRIPT
# ==============================================================================

# Executar a função principal
organizar_graficos_sensibilidade()

# Função adicional para verificar renomeação
verificar_renomeacao <- function() {
  cat("\n")
  cat(paste0(rep("=", 60), collapse = ""), "\n")
  cat("EXEMPLOS DE RENOMEAÇÃO ESPERADOS:\n")
  cat(paste0(rep("=", 60), collapse = ""), "\n")
  
  cat("Para arquivos em 'lambdas':\n")
  cat("  painel_C1_T100_w09.png -> painel_w09_C1_T100_w09.png\n")
  cat("  painel_mu_samples_C1_T100_w09.png -> painel_mu_samples_w09_C1_T100_w09.png\n")
  cat("  painel_theta_samples_C1_T100_w09.png -> painel_theta_samples_w09_C1_T100_w09.png\n\n")
  
  cat("Para arquivos em 'traceplots':\n")
  cat("  traceplot_gamma_4_.png -> traceplot_gamma_4_w09_C1_T100_w09.png\n")
  cat("  traceplot_beta_1_.png -> traceplot_beta_1_w09_C1_T100_w09.png\n")
  cat("  traceplot_gamma_1_.png -> traceplot_gamma_1_w09_C1_T100_w09.png\n")
}

# Executar verificação
verificar_renomeacao()