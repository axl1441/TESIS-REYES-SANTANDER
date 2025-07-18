##=====0. CARGA DE LIBRERÍAS =======
library(DESeq2)
library(AnnotationDbi)
library(ggplot2)
library(org.Hs.eg.db)
library(VennDiagram)
library(ggrepel)
library(reshape2)
library(plotly)
library(RColorBrewer)
library(dplyr)
library(vegan)      
library(data.table)
library(tidyr)
library(AnnotationDbi)
library(msigdbr)
library(DOSE)
library(Biobase)
library(gridExtra)
library(factoextra)
library(EnhancedVolcano)
library(glmnet)
library(caret)
library(e1071)
library(pheatmap)
library(BiocParallel)
library(mixOmics)
library(ggplot2)
library(themis) 
library(genefilter) 
library(ReactomePA) 
library(DO.db)
library(clusterProfiler)
library(enrichplot)
library(MLmetrics)

##===== 1. INTRODUCCIÓN Y PREPROCESAMIENTO DE DATOS INICIALES =======
## Introduciendo los datos
datos <- read.delim("C:/Users/Lenovo/Documents/Accel//Tesis 3/GSE261070_allcounts_PW.txt")
x <- datos[complete.cases(datos),] 
x <- as.data.frame(x)

## Preparando la matriz de conteos
transcriptoma_sin_na <- x[, -1] 
rownames(transcriptoma_sin_na) <- x[, 1] 
colnames(transcriptoma_sin_na) <- colnames(x)[-1] 

## Identificar a qué grupo pertenece cada muestra (fenotipo)
nombres_muestras_conteo <- colnames(transcriptoma_sin_na)
feno <- grepl("IC", nombres_muestras_conteo) 
nombre_feno <- as.factor(ifelse(feno == TRUE, "IC", "PW")) 

## Crear dataframe de metadatos (colData para DESeq2)
midf_sin_na <- data.frame(condicion = nombre_feno)
rownames(midf_sin_na) <- nombres_muestras_conteo

cat("Dimensiones de la matriz de conteos inicial (transcriptoma_sin_na):\n")
print(dim(transcriptoma_sin_na))
cat("Distribución inicial de muestras por condición:\n")
print(table(midf_sin_na$condicion))

##===== 2. ANÁLISIS DE EXPRESIÓN DIFERENCIAL (Con DEseq2) =======

# Crear el objeto DESeqDataSet
transcriptoma_sin_na_integer <- round(as.matrix(transcriptoma_sin_na))
storage.mode(transcriptoma_sin_na_integer) <- "integer"

dds <- DESeqDataSetFromMatrix(countData = transcriptoma_sin_na_integer,
                              colData = midf_sin_na,
                              design = ~ 0+condicion)

# Pre-filtrado para remover genes con pocos conteos
dds <- estimateSizeFactors(dds) 
idx <- rowSums(counts(dds, normalized = TRUE) >= 10) >= 5
dds <- dds[idx,]
cat(paste("Número de genes después del pre-filtrado: ", nrow(dds), "\n"))

# Correr el pipeline de DESeq2
dds <- DESeq(dds) # dds es el "cerebro" del análisis (contiene todo).

# Extraer conteos normalizados para todo el conjunto de datos
contajes_TMM <- counts(dds, normalized = TRUE) # es la versión normalizada para gráficos.

# Obtener los resultados del análisis de expresión diferencial
res <- results(dds, contrast = c("condicion", "PW", "IC"))
res.copy <- res

# Convertir a data.frame y limpiar NAs para EnhancedVolcano
res_df <- as.data.frame(res)
res_df <- na.omit(res_df) 

# Seleccionar los top N genes diferencialmente expresados 
top_genes_padj <- res_df[order(res_df$padj), ]
top_genes_names <- head(rownames(top_genes_padj), 5) 

select_genes <- head(rownames(top_genes_padj), 50)
aux <- reshape2::melt(scale(log2(t(contajes_TMM) + 1)))

ggplot(aux,aes(value))+geom_density()

##===== 3. DETECCIÓN Y REMOCIÓN DE OUTLIERS (Mahalanobis y PCA) =======

# Realizando detección y remoción de outliers, calcular PCA en los datos log-transformados y escalados
pcx <- prcomp(log2(t(contajes_TMM) + 1), scale. = TRUE)
scores_plot <- as.data.frame(pcx$x)

# Crear un DataFrame para visualización de PCA y Mahalanobis
data_of_pca <- data.frame(
  sample_name = rownames(midf_sin_na),
  condicion = midf_sin_na$condicion,
  PC1 = scores_plot$PC1,
  PC2 = scores_plot$PC2
)

# Cálculo y Detección de Outliers por Distancia de Mahalanobis
scores_for_maha <- pcx$x[, 1:2]
df_mahalanobis <- ncol(scores_for_maha)
d2 <- mahalanobis(scores_for_maha, center = colMeans(scores_for_maha), cov = cov(scores_for_maha))
pval <- pchisq(d2, df = df_mahalanobis, lower.tail = FALSE)

# Definir umbral y clasificar outliers
outlier_threshold_pval <- 0.001
is_outlier <- pval < outlier_threshold_pval
data_of_pca$is_outlier <- is_outlier

# Visualización de Outliers en PCA
plot_pca_outliers <- ggplot(data_of_pca, aes(x = PC1, y = PC2, color = condicion, shape = is_outlier)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 4), labels = c("No Outlier", "Outlier")) +
  labs(title = "PCA con Outliers (Mahalanobis)",
       subtitle = paste("Umbral p-valor Mahalanobis:", outlier_threshold_pval)) +
  theme_minimal() +
  geom_text_repel(aes(label = ifelse(is_outlier, sample_name, "")), box.padding = 0.5)

print(plot_pca_outliers)

# Filtrado y Preparación de Datos Limpios sin outliers
muestras_outliers_mahalanobis <- data_of_pca$sample_name[data_of_pca$is_outlier]
if (length(muestras_outliers_mahalanobis) > 0) {
  cat("Muestras identificadas como outliers y removidas:\n")
  print(muestras_outliers_mahalanobis)
} else {
  cat("No se identificaron outliers con el umbral especificado.\n")
}

# Crea un vector con los nombres de las muestras que no son outliers 
muestras_validas <- setdiff(colnames(transcriptoma_sin_na_integer), muestras_outliers_mahalanobis)

#Filtra la matriz de conteos original
transcriptoma_clean_raw <- round(as.matrix(transcriptoma_sin_na_integer[, muestras_validas]))

# Asegura que los conteos sean de tipo entero
storage.mode(transcriptoma_clean_raw) <- "integer"

# Actualiza el data frame de metadatos 
meta_sin_outliers <- midf_sin_na[muestras_validas, , drop = FALSE]

cat("Dimensiones de la matriz de conteos después de remover outliers:\n")
print(dim(transcriptoma_clean_raw))
cat("Distribución de muestras por condición después de remover outliers:\n")
print(table(meta_sin_outliers$condicion))

ggplot(data.frame(Distance = d2, Sample = rownames(scores_for_maha)), 
       aes(x = seq_along(Distance), y = Distance, color = is_outlier)) +
  geom_point() +
  geom_hline(yintercept = qchisq(1 - outlier_threshold_pval, df = df_mahalanobis), 
             linetype = "dashed", color = "red") +
  labs(title = "Distancias de Mahalanobis",
       x = "Índice de muestra", y = "Distancia de Mahalanobis",
       subtitle = paste("Umbral p-valor:", outlier_threshold_pval)) +
  scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
  theme_minimal()

##===== 4. Preparación, Balanceo de Muestras y Actualización de los Metadatos =======

# Balanceo de Muestras y Preparación Final de Datos y transpone la matriz de conteos limpia
csa <- t(transcriptoma_clean_raw)

# Muestra las dimensiones y convierte csa a un data frame y agrega la columna de condición
csa <- as.data.frame(csa)
csa$condicion <- nombre_feno[!is_outlier]

# Separa las muestras en dos data frames diferentes según su condición 
pw <- csa[csa$condicion=="PW",-ncol(csa)]
ic <- csa[csa$condicion=="IC",-ncol(csa)]

# Selecciona aleatoriamente 42 muestras del grupo "PW"
set.seed(123)
idx <- sample(nrow(pw),size = 42,replace = F)
idx

# Crea el data frame final para el grupo "PW" con las muestras submuestreadas.
pw_final <- pw[idx,]

# Combina las muestras submuestreadas de "PW" con todas las muestras de "IC" en un solo data frame final para el análisis.
dats_final <- rbind(pw_final,ic)

# Divide los metadatos limpios según las condiciones.
meta_sin_outliers_pw <- meta_sin_outliers[meta_sin_outliers$condicion=="PW",]
meta_sin_outliers_ic <- meta_sin_outliers[meta_sin_outliers$condicion=="IC",]

# Reconstruye los data frames de metadatos asegurándose de que solo incluyan las muestras que se conservaron en dats_fina
meta_sin_outliers_pw <- data.frame(condicion = meta_sin_outliers_pw[idx],row.names = rownames(pw_final))
meta_sin_outliers_ic <- data.frame(condicion = meta_sin_outliers_ic,row.names = rownames(ic))

# Combina los metadatos actualizados en un único data frame 
meta_final <- rbind(meta_sin_outliers_pw,meta_sin_outliers_ic)

# Crear dds_balanced (objeto DESeq2 balanceado)
dds_balanced <- DESeqDataSetFromMatrix(
  countData = t(dats_final),      
  colData = meta_final,            
  design = ~ condicion             
)

# Pre-filtrado para remover genes con pocos conteos
dds <- estimateSizeFactors(dds_balanced) 
idx <- rowSums(counts(dds, normalized = TRUE) >= 10) >= 5
## FILTRADO POR GENES CON BAJA VARIABLIDAD
## SUMANDO LAS FILAS QUE SON LOS GENES QUE ESTEN EN MAS DE 5 MUESTRAS
dds <- dds[idx,]

# 2. Ejecutar el análisis DEG balanceado
dds <- DESeq(dds)

# 3. Extraer resultados
res_balanced <- results(dds, contrast = c("condicion", "PW", "IC"))

saveRDS(dds_balanced, "dds_balanced.rds")
saveRDS(res_balanced, "res_balanced.rds")

# Considera visualizar la distribución de las distancias
ggplot(data.frame(Distance = d2), aes(x = Distance)) + 
  geom_histogram(bins = 30) +
  geom_vline(xintercept = qchisq(1-outlier_threshold_pval, df = df_mahalanobis), 
             color = "red") +
  labs(title = "Distribución de Distancias de Mahalanobis")

# Extraer conteos normalizados para todo el conjunto de datos
contajes_TMM <- counts(dds, normalized = TRUE)

aux <- reshape2::melt(scale(log2(t(contajes_TMM) + 1)))
ggplot(aux,aes(value))+geom_density()

#============ 5. Visualización PCA de los datos transformados ==================

# Selección de genes significativos
nombres_genes <- rownames(res_balanced)[abs(res_balanced$log2FoldChange) > 2 & res_balanced$padj < 0.01]

# Filtrar genes que existen en dats_final 
genes_validos <- nombres_genes[nombres_genes %in% colnames(dats_final)]
stopifnot(length(genes_validos) > 0)
dats_normalized  <- scale(log2(t(contajes_TMM) + 1))

# PCA (ajustar scale. según necesidad)
pcx_balanced <- prcomp(dats_normalized, scale. = TRUE)

# Scree Plot 
fviz_eig(pcx_balanced, addlabels = TRUE, ylim = c(0, 50)) +
  labs(title = "Scree Plot de PCA Balanceado",
       x = "Componentes Principales",
       y = "% de Varianza Explicada") +
  theme_minimal()

# Gráfico PCA (usar solo componentes con varianza relevante)
percentVar <- round(100 * pcx_balanced$sdev^2 / sum(pcx_balanced$sdev^2), 2)

ggplot(data.frame(PC1 = pcx_balanced$x[,1], 
                  PC2 = pcx_balanced$x[,2],
                  Condicion = meta_final$condicion),
       aes(PC1, PC2, color = Condicion)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +
  labs(x = paste0("PC1 (", percentVar[1], "%)"), 
       y = paste0("PC2 (", percentVar[2], "%)"))

##===== 6. Selección de Genes Basada en PCA Balanceado =======

factoextra::fviz_contrib(pcx_balanced,"var")

# Usar el PCA balanceado (pcx_balanced) y datos balanceados (dats_final_norm)
dd <- facto_summarize(pcx_balanced, element = "var", result = "contrib", axes = 1)
contrib <- dd$contrib
names(contrib) <- rownames(dd)

# Selección por contribución 
qntil_pca <- quantile(dd$contrib, 0.85)
prefinales <- dd$name[dd$contrib > qntil_pca]

# Selección por correlación con PC1 (usando datos balanceados)
cors_mat <- cor(dats_normalized[,prefinales], pcx_balanced$x)  
genes_cor <- rownames(cors_mat[abs(cors_mat[,1]) > 0.5, ])
genes_cor2 <- rownames(cors_mat[abs(cors_mat[,2]) > 0.5, ])

finales_ <- unique(c(genes_cor,genes_cor2))

# Lista única de genes
genes_pca_finales <- unique(c(finales_, as.character(prefinales)))

# Matriz filtrada 
dats_final_norm_f <- dats_normalized[,genes_pca_finales ]
dim(dats_final_norm_f)

##===== 7. Preparación de Datos para Modelado (GLMNET LASSO) =======

# Usar matriz balanceada 
set.seed(123)  # Para reproducibilidad
idx <- sample(1:nrow(dats_final_norm_f), size = nrow(dats_final_norm_f) * 0.8)

# Conjuntos de entrenamiento (80%) y prueba (20%)
X_train <- dats_final_norm_f[idx, ]
X_test <- dats_final_norm_f[-idx, ]

dim(dats_final_norm_f)
# Variables de respuesta (etiquetas balanceadas)
Y_train <- meta_final$condicion[idx]
Y_test <- meta_final$condicion[-idx]

# Dataframes para modelado (opcional, según necesidad)
data_Train <- data.frame(X_train, Y_train = Y_train)
data_test <- data.frame(X_test)

dim(X_train)  
dim(X_test)   
length(Y_train) == nrow(X_train)  
length(Y_test) == nrow(X_test)    
length(genes_pca_finales)  

##===== 8. Ejecución Y Evaluación de GLMNET LASSO =======
cat("\n--- Sección 8: Ejecución y Evaluación de GLMNET LASSO ---\n")

# Inicialización de Variables y Matrices para el Bucle
n_rep <- 500
accuracies <- numeric(n_rep) 
f1_scores <- numeric(n_rep)  # Vector para almacenar F1-scores
coef_freq <- setNames(rep(0, ncol(dats_final_norm_f)), colnames(dats_final_norm_f))
coef_values_sum <- setNames(rep(0, ncol(dats_final_norm_f)), colnames(dats_final_norm_f))

gene_selection_matrix <- matrix(
  0, 
  nrow = ncol(dats_final_norm_f), 
  ncol = n_rep, 
  dimnames = list(colnames(dats_final_norm_f), paste0("Rep_", 1:n_rep)) 
)

# calcular F1-score
calculate_f1 <- function(predicted, actual) {
  conf_matrix <- table(predicted, actual)
  if(nrow(conf_matrix) == 1) return(NA)  
  
  tp <- conf_matrix["PW", "PW"]
  fp <- conf_matrix["PW", "IC"]
  fn <- conf_matrix["IC", "PW"]
  
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  
  if(is.na(precision) || is.na(recall)) return(NA)
  if(precision + recall == 0) return(0)
  
  f1 <- 2 * (precision * recall) / (precision + recall)
  return(f1)
}

# Bucle Principal de GLMNET LASSO con Validación Cruzada Repetida 
cat("\nIniciando el bucle LASSO con", ncol(dats_final_norm_f), "genes y", n_rep, "repeticiones. Esto puede tomar un tiempo...\n")

for(i in 1:n_rep) {
  set.seed(100 + i) 
  
  # Dividir los datos usando los índices de meta_final
  idx <- sample(1:nrow(dats_final_norm_f), nrow(dats_final_norm_f)*0.8)
  
  X_train <- as.matrix(dats_final_norm_f[idx, ])
  X_test  <- as.matrix(dats_final_norm_f[-idx, ])
  
  Y_train <- meta_final$condicion[idx] 
  Y_test  <- meta_final$condicion[-idx] 
  
  y_train_bin <- as.numeric(Y_train == "PW")
  
  # Entrenar el modelo con estandarización interna
  mdl <- cv.glmnet(X_train, y = y_train_bin, family = "binomial", 
                   standardize = TRUE, alpha = 1) # alpha=1 para LASSO
  
  # Predicciones y evaluación
  prob <- predict(mdl, newx = X_test, s = "lambda.1se", type = "response")
  pred <- ifelse(prob > 0.5, "PW", "IC")
  pred <- factor(pred, levels = levels(Y_train))
  
  # Calcular métricas
  accuracies[i] <- mean(pred == Y_test)
  f1_scores[i] <- calculate_f1(pred, Y_test)
  
  # Procesamiento de coeficientes
  coefs <- coef(mdl, s = "lambda.1se")
  coefs_no_intercept <- coefs[-1, 1] 
  nonzero_genes_indices <- which(coefs_no_intercept != 0)
  
  if(length(nonzero_genes_indices) > 0) {
    selected_gene_names <- rownames(coefs)[-1][nonzero_genes_indices]
    coef_values_for_current_iter <- coefs_no_intercept[nonzero_genes_indices] 
    
    coef_freq[selected_gene_names] <- coef_freq[selected_gene_names] + 1
    coef_values_sum[selected_gene_names] <- coef_values_sum[selected_gene_names] + coef_values_for_current_iter
    gene_selection_matrix[selected_gene_names, i] <- 1
  }
  
  if (i %% 10 == 0) cat(paste("  Repetición", i, "de", n_rep, "completada.\n"))
}

cat("Bucle LASSO finalizado.\n")

# Ordenar genes por frecuencia de selección 
selected_sorted_freq <- sort(coef_freq[coef_freq > 50], decreasing = TRUE)
cat("\nTop 15 Variables más frecuentemente seleccionadas por LASSO:\n")
print((selected_sorted_freq)) 
selected_sorted_freq

# Crear dataframe resumen de genes
gene_summary <- data.frame(
  Gene = names(selected_sorted_freq),
  Frequency = selected_sorted_freq,
  Sum_Coef = coef_values_sum[names(selected_sorted_freq)],
  Avg_Coef = coef_values_sum[names(selected_sorted_freq)] / selected_sorted_freq)

genes_finales <- names(selected_sorted_freq)
genes_finales

## Resumen y Análisis de Genes Seleccionados
cat("\n--- Resumen de la Evaluación del Modelo LASSO ---\n")
cat("Accuracy promedio en", n_rep, "repeticiones:", mean(accuracies, na.rm = TRUE), "\n")
cat("Desviación estándar de Accuracy:", sd(accuracies, na.rm = TRUE), "\n")
cat("F1-Score promedio en", n_rep, "repeticiones:", mean(f1_scores, na.rm = TRUE), "\n")
cat("Desviación estándar de F1-Score:", sd(f1_scores, na.rm = TRUE), "\n")

#================ 9. genes importantes =========================================

# Creación del DataFrame final con resultados DEG 
gene_modelo <- sub("\\..*", "", genes_finales)
gene_modelo2 <- clusterProfiler::bitr(gene_modelo, 
                                      fromType = "ENSEMBL", 
                                      toType = "SYMBOL", 
                                      OrgDb = org.Hs.eg.db)
genes <- c(
  "IFITM1", "IFITM2", "IFNA1", "IFNA2", "IFNB1", "IFNG", "IFNAR2", "IFI6", "IFI16", 
  "IFI35", "IFIH1", "IFIT2", "IRF1", "IRF3", "IRF4", "IRF5", "IRF8", "IRGM", 
  "JAK1", "JAK2", "TYK2", "STAT1", "STAT2", "STAT3", "STAT5A", "STAT6", "SOCS1", 
  "MX1", "RSAD2", "OAS3", "OAS", "CLEC1", "HERC5", "IFI44L", "ISG15", "USP18", 
  "CXCL10", "CCL2", "CXCL11", "CX3CL1", "CCL3", "CCL4", "CCL20", "IL1B", "IL2", 
  "IL4", "IL5", "IL6", "IL7", "IL10", "IL12P70", "IL13", "IL17", "IL21", "IL23", 
  "TNF", "GMCSF", "CAV1", "HLAG", "SHB", "TIMP1", "VEGFA", "TCF7", "CD3E", 
  "TRABD2A", "FTH1P20", "MAML1", "LQIRAP1", "CCR7", "FAM102A", "IL7R", "CD6", 
  "PARP9", "CLK1", "DDX58", "TRAF3", "WNT", "CTNNB1", "RUNX3", "ACTB", "B2M", 
  "GAPDH", "HPRT1", "RPLP0", "SCGB3A5", "MUC5AC", "KRT5", "ACTA2"
)

intersect(gene_modelo2$SYMBOL,genes)

res_balanced_ <- res_balanced[which(rownames(res_balanced) %in% genes_finales),]

# Visualización de genes seleccionados
  volcano_data <- data.frame(
  gene = genes_finales,
  logFC = res_balanced_$log2FoldChange,
  pvalue = -log10(res_balanced_$padj))
  
ggplot(volcano_data, aes(x = logFC, y = pvalue)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c("grey", "red")) +
    labs(title = "Genes seleccionados por LASSO y DEG",
         x = "Log2 Fold Change", y = "-Log10 p-value") +
    theme_minimal()

#============================== 10. ORA ===============================================

# Universo
genes_pca_finales  
  
genes_pca_finales_ <- sub("\\..*", "", genes_pca_finales)
entrez_universe <- clusterProfiler::bitr(genes_pca_finales_, 
                                         fromType = "ENSEMBL", 
                                         toType = "ENTREZID", 
                                         OrgDb = org.Hs.eg.db)
    
# cambiar los id de ENSEMBL A ENTREZ
gene_modelo_ensembl_limpio <- sub("\\..*", "", genes_finales)
print(head(gene_modelo_ensembl_limpio))
    
genes_entrezid_df <- bitr(gene_modelo_ensembl_limpio,
                          fromType = "ENSEMBL",   
                          toType = "ENTREZID",   
                          OrgDb = org.Hs.eg.db) 
    
gois <- (genes_entrezid_df$ENTREZID)
  
GO_BP = enrichGO(gois, ont = "BP", universe = entrez_universe$ENTREZID, OrgDb = org.Hs.eg.db,
                pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 2,maxGSSize = 500)
GO_MF = enrichGO(gois, ont = "MF", universe = entrez_universe$ENTREZID, OrgDb = org.Hs.eg.db,
                pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 2,maxGSSize = 500)
Reactome = enrichPathway(gois, universe = entrez_universe$ENTREZID, organism = "human",
                        pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 2,maxGSSize = 500)
DO = enrichDGN(gois, universe = entrez_universe$ENTREZID, 
              pAdjustMethod = "none", pvalueCutoff = 0.1, qvalueCutoff = 0.1,minGSSize = 2,maxGSSize = 500,readable = T)
GO_CC = enrichGO(gois, ont = "CC", universe = entrez_universe$ENTREZID, OrgDb = org.Hs.eg.db,
                pAdjustMethod = "none", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 2,maxGSSize = 500)
KEGG = enrichKEGG(gois, universe = entrez_universe$ENTREZID, organism = 'hsa',
                  pAdjustMethod = "none", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 2,maxGSSize = 500)

    
# Dotplot para GO_BP (Ontología de Procesos Biológicos) ---
if (!is.null(GO_BP) && nrow(GO_BP@result) > 0) {
dotplot(GO_BP, showCategory = 10, title = "GO Biological Process Enrichment")} else {
print("GO_BP no tiene resultados válidos para mostrar un dotplot.")}
    
# Dotplot para GO_MF (Ontología de Funciones Moleculares) ---
if (!is.null(GO_MF) && nrow(GO_MF@result) > 0) {
dotplot(GO_MF, showCategory = 10, title = "GO Molecular Function Enrichment")} else {
print("GO_MF no tiene resultados válidos para mostrar un dotplot.")}
    
# Dotplot para Reactome ---
if (!is.null(Reactome) && nrow(Reactome@result) > 0) {
dotplot(Reactome, showCategory = 10, title = "Reactome Pathway Enrichment")} else {
print("Reactome no tiene resultados válidos para mostrar un dotplot.")}
    
# Dotplot para DO (Disease Ontology) 
if (!is.null(DO) && nrow(DO@result) > 0) {
dotplot(DO, showCategory = 10, title = "Disease Ontology Enrichment")} else {
print("DO no tiene resultados válidos para mostrar un dotplot.")}
  
# Dotplot para GO_CC (Ontología de Componentes Celulares) ---
#if (!is.null(GO_CC) && nrow(GO_CC@result) > 0) {
#  dotplot(GO_CC, showCategory = 5, title = "GO Cellular Component Enrichment")
#} else {
#  print("GO_CC no tiene resultados válidos para mostrar un dotplot.")}
  
# Dotplot para KEGG ---
#if (!is.null(KEGG) && nrow(KEGG@result) > 0) {
#  dotplot(KEGG, showCategory = 10, title = "KEGG Pathway Enrichment",color ="pvalue")
#} else {
#print("KEGG no tiene resultados válidos para mostrar un dotplot.")}
    
#===========================  11. ORA con Genes UP y DOWN ======================================== 
   
# colnames(res_balanced_)
# dim(res_balanced_)
#     
# # Definir los umbrales de significancia y fold change
# padj_threshold <- 0.05    
# logFC_threshold <- 1
#     
# # Identificar genes regulados al alza (UP-regulated)
#   genes_up <- rownames(res_balanced_)[
#   res_balanced_$padj < padj_threshold &
#   res_balanced_$log2FoldChange > logFC_threshold]
#     
# # Identificar genes regulados a la baja (DOWN-regulated)
#   genes_down <- rownames(res_balanced_)[
#   res_balanced_$padj < padj_threshold &
#   res_balanced_$log2FoldChange < -logFC_threshold]
#     
# # Mostrar el número de genes en cada categoría
# cat("Número de genes regulados al alza (UP-regulated):", length(genes_up), "\n")
# cat("Número de genes regulados a la baja (DOWN-regulated):", length(genes_down), "\n")
#     
# # ORA UP
# GO_BP_UP = enrichGO(genes_up, ont = "BP", universe = entrez_universe$ENTREZID, OrgDb = org.Hs.eg.db,
#                    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 2,maxGSSize = 500)
# GO_MF_UP = enrichGO(genes_up, ont = "MF", universe = entrez_universe$ENTREZID, OrgDb = org.Hs.eg.db,
#                    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 2,maxGSSize = 500)
# Reactome_UP = enrichPathway(genes_up, universe = entrez_universe$ENTREZID, organism = "human",
#                            pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 2,maxGSSize = 500)
# DO_UP = enrichDGN(genes_up, universe = entrez_universe$ENTREZID, 
#                  pAdjustMethod = "none", pvalueCutoff = 0.1, qvalueCutoff = 0.1,minGSSize = 2,maxGSSize = 500,readable = T)
# GO_CC_UP = enrichGO(genes_up, ont = "CC", universe = entrez_universe$ENTREZID, OrgDb = org.Hs.eg.db,
#                    pAdjustMethod = "none", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 2,maxGSSize = 500)
# KEGG_UP = enrichKEGG(genes_up, universe = entrez_universe$ENTREZID, organism = 'hsa',
#                     pAdjustMethod = "none", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 2,maxGSSize = 500)
# 
# # ORA DOWN
# GO_BP_down = enrichGO(genes_down, ont = "BP", universe = entrez_universe$ENTREZID, OrgDb = org.Hs.eg.db,
#                      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 2,maxGSSize = 500)
# GO_MF_down = enrichGO(genes_down, ont = "MF", universe = entrez_universe$ENTREZID, OrgDb = org.Hs.eg.db,
#                      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 2,maxGSSize = 500)
# Reactome_down = enrichPathway(genes_down, universe = entrez_universe$ENTREZID, organism = "human",
#                              pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 2,maxGSSize = 500)
# DO_down = enrichDGN(genes_down, universe = entrez_universe$ENTREZID, 
#                    pAdjustMethod = "none", pvalueCutoff = 0.1, qvalueCutoff = 0.1,minGSSize = 2,maxGSSize = 500,readable = T)
# GO_CC_down = enrichGO(genes_down, ont = "CC", universe = entrez_universe$ENTREZID, OrgDb = org.Hs.eg.db,
#                      pAdjustMethod = "none", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 2,maxGSSize = 500)
# KEGG_down = enrichKEGG(genes_down, universe = entrez_universe$ENTREZID, organism = 'hsa',
#                       pAdjustMethod = "none", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 2,maxGSSize = 500)
#     
#===================================== 12. gsea ========================================================================
    
# Obtener estadísticos de prueba ordenados (usando el estadístico Wald)
gene_rank <- res_balanced_$stat
names(gene_rank) <- rownames(res_balanced_)
gene_modelo <- sub("\\..*", "", genes_finales)
    
# Ordenar de mayor a menor (importante para GSEA)
gene_rank <- sort(gene_rank, decreasing = TRUE)
names(gene_rank) <- sub("\\..*", "", names(gene_rank))
    
## Convertir identificadores de genes 
gene_ids <- bitr(names(gene_rank), fromType = "ENSEMBL", 
            toType = "ENTREZID", 
            OrgDb = org.Hs.eg.db)
    
lista_genes <- gene_rank
names(lista_genes) <- gene_ids$ENTREZID
    
lista_genes <- lista_genes[complete.cases(names(lista_genes))]
    
## Realizar el GSEA 
gsea_kegg_result <- gseKEGG(geneList = lista_genes, organism = "hsa", pvalueCutoff = 1, 
                           scoreType = "pos", minGSSize = 2,maxGSSize = 500 )
    
gsea_go_bp_result <- gseGO(geneList = lista_genes,ont = "BP",OrgDb = org.Hs.eg.db, keyType = "ENTREZID", 
                          pvalueCutoff = 1,scoreType = "pos",minGSSize = 2,maxGSSize = 500 )
    
gsea_go_mf_result <- gseGO(geneList = lista_genes, ont = "MF", OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                          pvalueCutoff = 1,scoreType = "pos", minGSSize = 2, maxGSSize = 500)
    
gsea_go_cc_result <- gseGO(geneList = lista_genes,ont = "CC",OrgDb = org.Hs.eg.db,keyType = "ENTREZID",
                          pvalueCutoff = 1,scoreType = "pos",minGSSize = 2,maxGSSize = 500)
    
gsea_reactome_result <- gsePathway(geneList = lista_genes,organism = "human",pvalueCutoff = 1,
                                  scoreType = "pos",minGSSize = 2,maxGSSize = 500 )

# Crear una lista con los resultados de GSEA 
  gsea_results_list <- list(KEGG = gsea_kegg_result,GO_BP = gsea_go_bp_result,GO_MF = gsea_go_mf_result,
  GO_CC = gsea_go_cc_result,Reactome = gsea_reactome_result)
    
# Generar y mostrar los Dot Plots individualmente 
    if (!is.null(gsea_kegg_result) && nrow(gsea_kegg_result@result) > 0) {
      message("Generando Dot Plot para KEGG...")
      p_kegg <- dotplot(gsea_kegg_result,showCategory = 15, orderBy = "X",     
                        font.size = 10,title = "GSEA Dot Plot: KEGG Pathways")
      print(p_kegg)} else {message("No hay resultados de GSEA válidos para KEGG o el objeto es NULL.")}

    if (!is.null(gsea_go_bp_result) && nrow(gsea_go_bp_result@result) > 0) {
      message("Generando Dot Plot para GO Biological Process...")
      p_go_bp <- dotplot(gsea_go_bp_result,showCategory = 15,orderBy = "X",
                         font.size = 10,title = "GSEA Dot Plot: GO Biological Process")
      print(p_go_bp)
    } else {
      message("No hay resultados de GSEA válidos para GO Biological Process o el objeto es NULL.")
    }

    if (!is.null(gsea_go_mf_result) && nrow(gsea_go_mf_result@result) > 0) {
      message("Generando Dot Plot para GO Molecular Function...")
      p_go_mf <- dotplot(gsea_go_mf_result,showCategory = 15,orderBy = "X",
                         font.size = 10,title = "GSEA Dot Plot: GO Molecular Function")
      print(p_go_mf)
    } else {
      message("No hay resultados de GSEA válidos para GO Molecular Function o el objeto es NULL.")
    }
    
    if (!is.null(gsea_go_cc_result) && nrow(gsea_go_cc_result@result) > 0) {
      message("Generando Dot Plot para GO Cellular Component...")
      p_go_cc <- dotplot(gsea_go_cc_result,showCategory = 15,orderBy = "X",
                         font.size = 10,title = "GSEA Dot Plot: GO Cellular Component")
      print(p_go_cc)
    } else {
      message("No hay resultados de GSEA válidos para GO Cellular Component o el objeto es NULL.")
    }

    if (!is.null(gsea_reactome_result) && nrow(gsea_reactome_result@result) > 0) {
      message("Generando Dot Plot para Reactome...")
      p_reactome <- dotplot(gsea_reactome_result,showCategory = 15,orderBy = "X",
                            font.size = 10,title = "GSEA Dot Plot: Reactome Pathways")
      print(p_reactome)
    } else {
      message("No hay resultados de GSEA válidos para Reactome o el objeto es NULL.")
    }
    
    message("Todos los análisis GSEA y dot plots han sido procesados.")