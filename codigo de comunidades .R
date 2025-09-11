
library(igraph)

if (exists("g")) {
  cat(" Red cargada exitosamente:\n")
  cat("   Nodos:", vcount(g), "\n")
  cat("   Aristas:", ecount(g), "\n")
  cat("   Atributos de nodos:", vertex_attr_names(g), "\n") 
  cat("   Atributos de aristas:", edge_attr_names(g), "\n")  
}

# Las comunidades ya están calculadas
# Puedes extraerlas directamente:
comunidades <- V(g)$community

# Extraer la membresía de comunidades
membresia <- V(g)$community

# Ver distribución
cat("Distribución de comunidades:\n")
print(table(membresia))

# Número de comunidades únicas
num_comunidades <- length(unique(membresia))
cat("Número de comunidades detectadas:", num_comunidades, "\n")


# Crear comm_list desde los atributos
comm_list <- split(V(g)$name, V(g)$community)

# Ver tamaño de comunidades
cat("Tamaño de las comunidades:\n")
tamanos <- sapply(comm_list, length)
print(sort(tamanos, decreasing = TRUE))

# Top 5 comunidades más grandes
cat("\nTop 5 comunidades más grandes:\n")
top_5 <- head(sort(tamanos, decreasing = TRUE), 5)
print(top_5)


library(igraph)
library(dplyr)

# 1. Verificar la membresía de comunidades
cat("=== ANÁLISIS DE COMUNIDADES EXISTENTES ===\n")
cat("Total de comunidades:", length(unique(V(g)$community)), "\n")
cat("Distribución de comunidades:\n")
print(table(V(g)$community))

# 2. Crear lista de comunidades
comm_list <- split(V(g)$name, V(g)$community)

# 3. Calcular estadísticas por comunidad
comunidad_stats <- data.frame(
  Comunidad = names(comm_list),
  Nodos = sapply(comm_list, length),
  Densidad = sapply(comm_list, function(nodos) {
    subg <- induced_subgraph(g, vids = nodos)
    edge_density(subg)
  })
) %>% arrange(desc(Nodos))

# 4. Mostrar resultados
cat("\n ESTADÍSTICAS POR COMUNIDAD:\n")
print(comunidad_stats)

# 5. Top comunidades
cat("\n TOP 5 COMUNIDADES MÁS GRANDES:\n")
top_5 <- head(comunidad_stats, 5)
print(top_5)

# 6. Exportar resultados
write.csv(comunidad_stats, "estadisticas_comunidades_existentes.csv", row.names = FALSE)

library(igraph)
library(dplyr)


cat("=== DESCRIPTORES DE MÓDULOS CON ENLACES ===\n")

# 1. Calcular estadísticas COMPLETAS por comunidad
comunidad_stats <- data.frame(
  Comunidad = names(comm_list),
  Nodos = sapply(comm_list, length),
  Enlaces = sapply(comm_list, function(nodos) {
    subg <- induced_subgraph(g, vids = nodos)
    ecount(subg)  
  }),
  Densidad = sapply(comm_list, function(nodos) {
    subg <- induced_subgraph(g, vids = nodos)
    edge_density(subg)
  })
) %>% arrange(desc(Nodos))

# 2. Crear tabla 
tabla_tutor <- comunidad_stats %>%
  mutate(
    Modulo_ID = paste0("Módulo_", Comunidad),
    Porcentaje_Nodos = round((Nodos / vcount(g)) * 100, 2)
  ) %>%
  select(Modulo_ID, Nodos, Enlaces, Densidad, Porcentaje_Nodos) %>%
  arrange(desc(Nodos))

# 3. Exportar tabla final
write.csv(tabla_tutor, "descriptores_modulos_COMPLETO.csv", row.names = FALSE)

# 4. Mostrar resumen EJECUTIVO
cat(" ESTADÍSTICAS GENERALES:\n")
cat("--------------------------\n")
cat(" Total de módulos:", nrow(tabla_tutor), "\n")
cat(" Total de nodos en la red:", vcount(g), "\n")
cat(" Total de enlaces en la red:", ecount(g), "\n")
cat(" Modularidad global:", round(modularity(g, V(g)$community), 4), "\n\n")

cat(" TOP 5 MÓDULOS MÁS GRANDES:\n")
cat("-----------------------------\n")
print(head(tabla_tutor, 5))

cat("\n TOP 5 MÓDULOS MÁS CONECTADOS:\n")
cat("---------------------------------\n")
top_conectados <- tabla_tutor %>% arrange(desc(Enlaces)) %>% head(5)
print(top_conectados)

cat("\n TOP 5 MÓDULOS MÁS DENSOS:\n")
cat("-----------------------------\n")
top_densos <- tabla_tutor %>% arrange(desc(Densidad)) %>% head(5)
print(top_densos)
