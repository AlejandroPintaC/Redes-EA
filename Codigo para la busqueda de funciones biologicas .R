
 install.packages("gprofiler2")

library(gprofiler2)
 
 #  Analizar la COMUNIDAD 3 (477 genes)
 genes_comunidad <- comm_list[[3]]  # Ajusta el número según la comunidad que quieras
 
 # Verificar los genes
 cat("Analizando", length(genes_comunidad), "genes de la comunidad\n")
 print(head(genes_comunidad, 10))  # Primeros 10 genes
 
 # Correr el análisis completo
 resultados <- gost(
   query = genes_comunidad,     # Tus genes Ensembl
   organism = "hsapiens",       # Organismo (humano)
   sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"),  # Bases de datos
   significant = TRUE,          # Solo resultados significativos
   correction_method = "g_SCS"  # Método de corrección
 )
 
 # Ver resultados en la consola
 print(resultados)
 
 # Ver tabla de resultados
 tabla_resultados <- resultados$result
 head(tabla_resultados[, c("term_name", "p_value", "term_size", "intersection_size")])
 
 # Crear gráfica de los top términos
 gostplot(resultados, capped = TRUE, interactive = TRUE)
 
 # Top 10 términos más significativos
 top_terminos <- head(tabla_resultados[order(tabla_resultados$p_value), ], 10)
 print(top_terminos[, c("term_name", "p_value", "source")])
 
 # Analizar las 3 comunidades principales automáticamente
 for (i in 1:3) {
   cat("\n=== ANALIZANDO COMUNIDAD", i, "===\n")
   
   resultados_com <- gost(
     query = comm_list[[i]],
     organism = "hsapiens",
     sources = c("GO:BP", "KEGG", "REAC"),
     significant = TRUE
   )
   
   # Exportar resultados
   if (!is.null(resultados_com$result)) {
     archivo_nombre <- paste0("enriquecimiento_comunidad_", i, ".csv")
     write.csv(resultados_com$result, archivo_nombre, row.names = FALSE)
     cat("✓", archivo_nombre, "exportado\n")
   }
 }
 
# A partir de aqui, ya terminamos con una comunidad, siguen las otras   
 
# Analizar las comunidades 1, 2 y 3
 comunidades_analizar <- c(1, 2, 3)  # Se ajustan los numeros segun se necesite
 
 for (i in comunidades_analizar) {
   cat("\n Analizando Comunidad", i, "-", length(comm_list[[i]]), "genes\n")
   
   # Ejecutar enriquecimiento funcional
   resultados <- gost(
     query = comm_list[[i]],
     organism = "hsapiens",
     sources = c("GO:BP", "KEGG", "REAC"),
     significant = TRUE
   )
   
   # Guardar resultados en CSV
   if (!is.null(resultados$result)) {
     nombre_archivo <- paste0("enriquecimiento_comunidad_", i, ".csv")
     write.csv(resultados$result, nombre_archivo, row.names = FALSE)
     cat("✓ Resultados guardados en:", nombre_archivo, "\n")
     
     # Mostrar top 3 términos más significativos
     top_terms <- head(resultados$result[order(resultados$result$p_value), ], 3)
     cat("Top 3 términos enriquecidos:\n")
     print(top_terms[, c("term_name", "p_value", "source")])
   }
 }
 
 library(tidyverse)
 
 # Crear dataframe con comunidades a analizar
 comunidades_df <- tibble(
   comunidad_id = c(1, 2, 3),
   genes = comm_list[comunidad_id]
 )
 
 # Analizar y exportar todos
 resultados <- comunidades_df %>%
   mutate(
     analisis = map(genes, ~gost(.x, organism = "hsapiens", 
                                 sources = c("GO:BP", "KEGG", "REAC"))),
     resultados = map(analisis, ~.x$result),
     exportado = map2(resultados, comunidad_id, 
                      ~write.csv(.x, paste0("comunidad_", .y, ".csv"), row.names = FALSE))
   )
 
 # Función para extraer top términos de cada comunidad
 extraer_top_terminos <- function(i) {
   archivo <- paste0("enriquecimiento_comunidad_", i, ".csv")
   if (file.exists(archivo)) {
     datos <- read.csv(archivo)
     datos %>% 
       arrange(p_value) %>%
       head(5) %>%
       select(term_name, p_value, source) %>%
       mutate(comunidad = i)
   }
 }
 
 # Crear comparativa
 comparativa <- map_dfr(c(1, 2, 3), extraer_top_terminos)
 
 # Ver términos más comunes entre comunidades
 comparativa %>%
   group_by(term_name) %>%
   summarise(
     n_comunidades = n_distinct(comunidad),
     min_p_value = min(p_value)
   ) %>%
   arrange(-n_comunidades, min_p_value)
 
 # Análisis completo para todas las comunidades significativas
 analisis_completo <- function(comunidades_ids) {
   # Crear directorio para resultados
   dir.create("resultados_enriquecimiento", showWarnings = FALSE)
   
   # Analizar cada comunidad
   for (i in comunidades_ids) {
     cat("\n", strrep("=", 50))
     cat("\nAnalizando Comunidad", i, "\n")
     cat(strrep("=", 50), "\n")
     
     # Enriquecimiento funcional
     resultados <- gost(comm_list[[i]], 
                        organism = "hsapiens",
                        sources = c("GO:BP", "KEGG", "REAC", "WP"))
     
     if (!is.null(resultados$result)) {
       # Guardar resultados
       archivo <- paste0("resultados_enriquecimiento/comunidad_", i, "_enriquecimiento.csv")
       write.csv(resultados$result, archivo, row.names = FALSE)
       
       # Generar reporte breve
       top_5 <- resultados$result %>%
         arrange(p_value) %>%
         head(5)
       
       cat("Top 5 términos enriquecidos:\n")
       print(top_5[, c("term_name", "p_value")])
       
       # Gráfica de los términos más significativos
       if (nrow(top_5) > 0) {
         p <- gostplot(resultados, capped = FALSE, interactive = FALSE)
         ggsave(paste0("resultados_enriquecimiento/comunidad_", i, "_plot.png"), 
                plot = p, width = 10, height = 6)
       }
     }
   }
 }
 
 # Ejecutar para las 3 comunidades principales
 analisis_completo(c(1, 2, 3))
 
 