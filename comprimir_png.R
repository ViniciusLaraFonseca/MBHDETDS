# ==============================================================================
# SCRIPT SIMPLIFICADO: RECRIAR PASTA E COMBINAR TRACEPLOTS
# ==============================================================================

# Caminhos
base_dir <- "C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/MBHDETDS/graficos_sensibilidade_09122025"
output_dir <- "C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/MBHDETDS/graficos_sensibilidade_comprimidos_09122025"

cat("=== INICIANDO PROCESSAMENTO ===\n")
cat("Pasta original:", base_dir, "\n")
cat("Arquivos na original: 80 PNG\n")
cat("Pasta de saída:", output_dir, "\n\n")

# 1. CRIAR/ESVAZIAR PASTA DE SAÍDA
if (dir.exists(output_dir)) {
  cat("1. Removendo pasta de saída existente...\n")
  unlink(output_dir, recursive = TRUE)
}

cat("2. Criando nova pasta de saída...\n")
dir.create(output_dir, recursive = TRUE)

# 2. COPIAR TODOS OS ARQUIVOS DA ORIGINAL
cat("3. Copiando arquivos da pasta original...\n")
arquivos_originais <- list.files(base_dir, pattern = "\\.png$", full.names = TRUE)

for (i in seq_along(arquivos_originais)) {
  file.copy(arquivos_originais[i], 
            file.path(output_dir, basename(arquivos_originais[i])))
  
  # Mostrar progresso a cada 10 arquivos
  if (i %% 10 == 0 || i == length(arquivos_originais)) {
    cat(sprintf("   [%3d/%3d] arquivos copiados\n", i, length(arquivos_originais)))
  }
}

cat("✓ Cópia concluída: 80 arquivos copiados\n\n")

# 3. VERIFICAR O QUE TEMOS
cat("4. Verificando arquivos copiados...\n")
arquivos_copiados <- list.files(output_dir, pattern = "\\.png$", full.names = FALSE)
cat("   Total na pasta de saída:", length(arquivos_copiados), "\n")

# Contar por tipo
cat("\n   Distribuição inicial:\n")
cat("   - Painéis lambda:", sum(grepl("^painel_w[0-9]{2}_C[1-4]_", arquivos_copiados) & 
                                  !grepl("mu_samples|theta_samples", arquivos_copiados)), "\n")
cat("   - Painéis mu:", sum(grepl("^painel_mu_samples_", arquivos_copiados)), "\n")
cat("   - Painéis theta:", sum(grepl("^painel_theta_samples_", arquivos_copiados)), "\n")
cat("   - Traceplots beta individuais:", sum(grepl("^traceplot_beta_", arquivos_copiados)), "\n")
cat("   - Traceplots gamma individuais:", sum(grepl("^traceplot_gamma_", arquivos_copiados)), "\n")

# 4. INSTALAR PACOTES SE NECESSÁRIO
cat("\n5. Verificando pacotes...\n")
if (!require("png", quietly = TRUE)) {
  install.packages("png")
}
if (!require("grid", quietly = TRUE)) {
  install.packages("grid")
}
if (!require("gridExtra", quietly = TRUE)) {
  install.packages("gridExtra")
}

library(png)
library(grid)
library(gridExtra)

# 5. FUNÇÃO PARA COMBINAR BETAS
cat("\n6. Criando função para combinar betas...\n")

combinar_betas_simples <- function(w_ger, cenario_completo) {
  # Verificar se existem os 3 betas
  betas <- c()
  for (i in 1:3) {
    beta_file <- paste0("traceplot_beta_", i, "_", w_ger, "_", cenario_completo, ".png")
    if (file.exists(file.path(output_dir, beta_file))) {
      betas <- c(betas, beta_file)
    }
  }
  
  if (length(betas) == 3) {
    cat("   Combinando betas para", w_ger, cenario_completo, "...")
    
    tryCatch({
      # Ler as 3 imagens
      img1 <- readPNG(file.path(output_dir, betas[1]))
      img2 <- readPNG(file.path(output_dir, betas[2]))
      img3 <- readPNG(file.path(output_dir, betas[3]))
      
      # Criar grobs
      g1 <- rasterGrob(img1)
      g2 <- rasterGrob(img2)
      g3 <- rasterGrob(img3)
      
      # Combinar horizontalmente
      combinada <- arrangeGrob(g1, g2, g3, nrow = 1, ncol = 3)
      
      # Salvar
      nome_saida <- paste0("beta_combinado_", w_ger, "_", cenario_completo, ".png")
      png(file.path(output_dir, nome_saida), width = 1800, height = 500)
      grid.draw(combinada)
      dev.off()
      
      cat(" OK\n")
      
      # Remover arquivos individuais
      for (beta in betas) {
        file.remove(file.path(output_dir, beta))
      }
      
      return(TRUE)
    }, error = function(e) {
      cat(" ERRO:", e$message, "\n")
      return(FALSE)
    })
  } else {
    cat("   Não encontrou 3 betas para", w_ger, cenario_completo, 
        "(encontrou", length(betas), ")\n")
    return(FALSE)
  }
}

# 6. FUNÇÃO PARA COMBINAR GAMMAS
cat("\n7. Criando função para combinar gammas...\n")

combinar_gammas_simples <- function(w_ger, cenario_completo) {
  # Verificar se existem os 4 gammas
  gammas <- c()
  for (i in 1:4) {
    gamma_file <- paste0("traceplot_gamma_", i, "_", w_ger, "_", cenario_completo, ".png")
    if (file.exists(file.path(output_dir, gamma_file))) {
      gammas <- c(gammas, gamma_file)
    }
  }
  
  if (length(gammas) == 4) {
    cat("   Combinando gammas para", w_ger, cenario_completo, "...")
    
    tryCatch({
      # Ler as 4 imagens
      imagens <- list()
      for (i in 1:4) {
        imagens[[i]] <- rasterGrob(readPNG(file.path(output_dir, gammas[i])))
      }
      
      # Combinar em grid 2x2
      combinada <- arrangeGrob(grobs = imagens, nrow = 2, ncol = 2)
      
      # Salvar
      nome_saida <- paste0("gamma_combinado_", w_ger, "_", cenario_completo, ".png")
      png(file.path(output_dir, nome_saida), width = 1200, height = 800)
      grid.draw(combinada)
      dev.off()
      
      cat(" OK\n")
      
      # Remover arquivos individuais
      for (gamma in gammas) {
        file.remove(file.path(output_dir, gamma))
      }
      
      return(TRUE)
    }, error = function(e) {
      cat(" ERRO:", e$message, "\n")
      return(FALSE)
    })
  } else {
    cat("   Não encontrou 4 gammas para", w_ger, cenario_completo, 
        "(encontrou", length(gammas), ")\n")
    return(FALSE)
  }
}

# 7. LISTA FIXA DOS 8 CENÁRIOS
cat("\n8. Definindo os 8 cenários...\n")

cenarios <- list(
  list(w_ger = "w09", cenario = "C1_T100_w09"),
  list(w_ger = "w09", cenario = "C2_T100_w07"),
  list(w_ger = "w09", cenario = "C3_T23_w09"),
  list(w_ger = "w09", cenario = "C4_T23_w07"),
  list(w_ger = "w07", cenario = "C1_T100_w09"),
  list(w_ger = "w07", cenario = "C2_T100_w07"),
  list(w_ger = "w07", cenario = "C3_T23_w09"),
  list(w_ger = "w07", cenario = "C4_T23_w07")
)

# Mostrar cenários
for (i in seq_along(cenarios)) {
  cat(sprintf("   %d. w_ger=%s, cenario=%s\n", 
              i, cenarios[[i]]$w_ger, cenarios[[i]]$cenario))
}

# 8. PROCESSAR CADA CENÁRIO
cat("\n9. Processando cada cenário...\n")

for (cen in cenarios) {
  w_ger <- cen$w_ger
  cenario_completo <- cen$cenario
  
  cat("\n   --- ", w_ger, cenario_completo, " ---\n")
  
  # Verificar se os arquivos básicos existem
  painel_lambda <- paste0("painel_", w_ger, "_", cenario_completo, ".png")
  painel_mu <- paste0("painel_mu_samples_", w_ger, "_", cenario_completo, ".png")
  painel_theta <- paste0("painel_theta_samples_", w_ger, "_", cenario_completo, ".png")
  
  cat("   Painéis básicos:\n")
  cat("     Lambda:", ifelse(file.exists(file.path(output_dir, painel_lambda)), "✓", "✗"), "\n")
  cat("     Mu:", ifelse(file.exists(file.path(output_dir, painel_mu)), "✓", "✗"), "\n")
  cat("     Theta:", ifelse(file.exists(file.path(output_dir, painel_theta)), "✓", "✗"), "\n")
  
  # Combinar betas
  combinar_betas_simples(w_ger, cenario_completo)
  
  # Combinar gammas
  combinar_gammas_simples(w_ger, cenario_completo)
}

# 9. VERIFICAÇÃO FINAL
cat("\n\n10. VERIFICAÇÃO FINAL:\n")

arquivos_finais <- list.files(output_dir, pattern = "\\.png$")
cat("   Total de arquivos na pasta final:", length(arquivos_finais), "\n\n")

# Contar por tipo
cat("   Distribuição final:\n")
cat("   - Painéis lambda:", sum(grepl("^painel_w[0-9]{2}_C[1-4]_", arquivos_finais) & 
                                  !grepl("mu_samples|theta_samples", arquivos_finais)), "\n")
cat("   - Painéis mu:", sum(grepl("^painel_mu_samples_", arquivos_finais)), "\n")
cat("   - Painéis theta:", sum(grepl("^painel_theta_samples_", arquivos_finais)), "\n")
cat("   - Betas combinados:", sum(grepl("^beta_combinado_", arquivos_finais)), "\n")
cat("   - Gammas combinados:", sum(grepl("^gamma_combinado_", arquivos_finais)), "\n")
cat("   - Traceplots individuais restantes:", sum(grepl("^traceplot_(beta|gamma)_", arquivos_finais)), "\n")

# Verificar cada cenário
cat("\n   Verificação por cenário:\n")
for (cen in cenarios) {
  w_ger <- cen$w_ger
  cenario_completo <- cen$cenario
  
  arquivos_cen <- arquivos_finais[grepl(w_ger, arquivos_finais) & 
                                    grepl(cenario_completo, arquivos_finais)]
  
  n_arquivos <- length(arquivos_cen)
  status <- ifelse(n_arquivos == 5, "✓ COMPLETO", paste("✗ (", n_arquivos, "/5)", sep = ""))
  
  cat(sprintf("   - %s %s: %s\n", w_ger, cenario_completo, status))
  
  if (n_arquivos < 5) {
    cat("     Arquivos encontrados:", paste(arquivos_cen, collapse = ", "), "\n")
  }
}

# 10. CRIAR LISTA DE ARQUIVOS
cat("\n11. Criando lista de arquivos...\n")
lista_arquivos <- data.frame(
  arquivo = arquivos_finais,
  tamanho_kb = round(file.size(file.path(output_dir, arquivos_finais)) / 1024, 1)
)

# Salvar na pasta pai para garantir que conseguimos salvar
caminho_lista <- file.path(dirname(output_dir), "lista_arquivos_finais.csv")
write.csv(lista_arquivos, caminho_lista, row.names = FALSE)
cat("   Lista salva em:", caminho_lista, "\n")

# 11. RESUMO
cat("\n" , rep("=", 60), "\n", sep = "")
cat("RESUMO DO PROCESSAMENTO\n")
cat(rep("=", 60), "\n", sep = "")
cat("Pasta original: ", base_dir, "\n")
cat("   - 80 arquivos PNG\n")
cat("\nPasta de saída: ", output_dir, "\n")
cat("   - ", length(arquivos_finais), " arquivos finais\n", sep = "")
cat("   - Esperado: 40 arquivos (8 cenários × 5 arquivos)\n")

if (length(arquivos_finais) == 40) {
  cat("\n✓ SUCESSO! Todos os arquivos foram processados corretamente.\n")
} else {
  cat("\n⚠ ATENÇÃO: Número de arquivos diferente do esperado.\n")
  cat("   Verifique os erros acima.\n")
}

cat("\nPRÓXIMOS PASSOS:\n")
cat("1. Faça upload da pasta '", basename(output_dir), "' para o Overleaf\n", sep = "")
cat("2. Use o arquivo LaTeX fornecido anteriormente\n")
cat("3. Compile o documento\n")

cat("\nPara verificar manualmente, abra a pasta:\n")
cat(output_dir, "\n")