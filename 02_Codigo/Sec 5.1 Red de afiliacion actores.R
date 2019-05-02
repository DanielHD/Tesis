
library(ggplot2)
library(igraph)
library(qgraph)
source("02_Codigo/multifunction.R")

BDD <- read.csv("01_Datos/02_Datos Coaliciones de eventos/BDD simple LAOMS.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=TRUE)
BDD$C1x1FechaPeriodico <- as.Date(BDD$C1x1FechaPeriodico, "%d/%m/%Y")
BDD$C1x2FechaEvento <- as.Date(BDD$C1x2FechaEvento, "%d/%m/%Y")
actores <- read.csv("01_Datos/02_Datos Coaliciones de eventos/actores.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=TRUE)
actores_eventos <- read.csv("01_Datos/02_Datos Coaliciones de eventos/actores_eventos.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=TRUE)

id_actores <- levels(actores$Nombre.corto)
id_eventos <- levels(BDD$Id)

################################################################################################
### 5.1. Red bipartita: Eventos de coalición y actores colectivos
################################################################################################
dimensiones = c(id_actores, id_eventos)

Adj <- matrix(0L, nrow=length(dimensiones), ncol=length(dimensiones))
rownames(Adj) <- c(dimensiones)
colnames(Adj) <- c(dimensiones)

for(id_ep in id_eventos){
    actores_ep <- as.character(actores_eventos[actores_eventos$Id==id_ep,"Organización"])
    for(actor in actores_ep){ Adj[id_ep,actor] <- Adj[actor,id_ep] <- 1 }
}

g <- graph_from_adjacency_matrix(Adj, mode="undirected")
V(g)$tipo <- c(rep("actor", length(id_actores)), rep("evento", length(id_eventos))) 
V(g)$color <- c(rep("blue", length(id_actores)), rep("red", length(id_eventos))) 
V(g)$type <- c(rep(0, length(id_actores)), rep(1, length(id_eventos))) 
V(g)$grado <- degree(g)
E(g)$color <- "#444444"
V(g)$c_intermediacion <- betweenness(g)
V(g)$c_eigenvector <- eigen_centrality(g)$vector

matriz_distancias <- shortest.paths(g, v=V(g), to=V(g))
diag(matriz_distancias) <- NA
V(g)$c_cercania <- rowSums(1/matriz_distancias, na.rm = TRUE) ### centralidad de cercanía (para grafos desconectados)

e <- get.edgelist(g,names=FALSE)
l <- qgraph.layout.fruchtermanreingold(e, vcount=vcount(g), area=7*(vcount(g)^2), repulse.rad=(vcount(g)^3.3))

plot(g, layout = l, 
     vertex.label = c(V(g)$name[1:length(id_actores)], rep(NA, length(id_eventos))), 
     vertex.label.color = c(rep("#00ccff", length(id_actores)), rep("red", length(id_eventos)) ), 
     vertex.label.cex = 0.2 + V(g)$grado*0.01, 
     vertex.size  = 0.5 + V(g)$grado*0.15, 
     vertex.frame.color = "#888888", edge.color = "#888888", 
     edge.width=0.3)
legend("bottomright", legend=c("Actores", "Eventos"), pch=20, cex=0.8, col=c("blue", "red"), inset = c(0.19,0.12), bty = "n", pt.cex=2, y.intersp=0.4)

###############################################################################################################
#### Medidas de centralidad normalizadas
n1 <- length(id_actores)
n2 <- length(id_eventos)
V(g)$grado_normalizado <- c(V(g)$grado[1:length(id_actores)]/n2, V(g)$grado[length(id_actores)+1:length(id_eventos)]/n1)
V(g)$c_cercania_normalizada <- c(
    V(g)$c_cercania[1:length(id_actores)] / ( n2 + 2*n1 - 2 ), 
    V(g)$c_cercania[length(id_actores)+1:length(id_eventos)] / ( n1 + 2*n2 - 2 ))
V(g)$c_intermediacion_normalizada <- c(
    V(g)$c_intermediacion[1:length(id_actores)] / ( 0.5*n1*(n1-1) + 0.5*(n2-1)*(n2-2) + (n2-1)*(n1-1) ), 
    V(g)$c_intermediacion[length(id_actores)+1:length(id_eventos)] / ( 2*(n2-1)*(n1-1) ))

unidades_sociales <- as_data_frame(g, what="vertices")
unidades_sociales[is.na(unidades_sociales)] <- 0
v_actores <- unidades_sociales[unidades_sociales$tipo=="actor",]
v_eps <- unidades_sociales[unidades_sociales$tipo=="evento",]

###############################################################################################################
####  Descripciones básicas

componentes <- components(g)
componentes$no # 31
componentes$csize
V(g)$componente <- componentes$membership

length(E(g))/(n1*n2) # densidad = 0.007067531
mean_distance(g) # Distancia promedio = 5.8057
length(V(g))*(length(V(g))-1)*(1/sum(1/matriz_distancias, na.rm = TRUE)) # distancia promedio DE ACUERDO A LA ECUACIÓN 1.3 (p.24) = 6.488805

#### Componente principal
g_ppal <- decompose.graph(g)[[2]]
excentricidades <- eccentricity(g_ppal, vids = V(g_ppal))
max(excentricidades) # Diámetro = 16
periferia <- excentricidades[excentricidades==16]  # 21 vértices
min(excentricidades) # Radio = 9
centro <- excentricidades[excentricidades==9] # 21 eventos

mean_distance(g_ppal) # Distancia promedio = 5.808504
hist(excentricidades, main="", ylab="Frecuencia", xlab="Excentrididad", breaks = c(8:16), freq=T)

###############################################################################################################
#### COEFICIENTES DE AGRUPAMIENTO Y CLUSTERING EN LA COMPONENTE PRINCIPAL
Matriz_afiliacion <- Adj[id_actores, id_eventos]
vertices_otras_componentes <- c(names(componentes$membership[componentes$membership!=2]))
Matriz_afiliacion_ppal <- Matriz_afiliacion[!rownames(Matriz_afiliacion) %in% vertices_otras_componentes, !colnames(Matriz_afiliacion) %in% vertices_otras_componentes]

library(tnet)
bip_ppal <- as.tnet(Matriz_afiliacion_ppal, type="binary two-mode tnet")
bip <- as.tnet(Matriz_afiliacion, type="binary two-mode tnet")


# REFORZAMIENTO EMPIRICO Y ALEATORIO
reinforcement_tm(bip_ppal) ## 0.4060121
reinforcement_tm(bip) ## 0.4068178
reforzamiento_aleatorio <- rep(NaN, 10000)
for(i in 1:length(reforzamiento_aleatorio)) {
    red_aleatoria <- rg_tm(ni=432, np=583, ties=1780, weights=1, seed=i)
    reforzamiento_aleatorio[i] <- reinforcement_tm(red_aleatoria)
}
hist(reforzamiento_aleatorio, breaks = "Scott", xlab = "Reforzamiento aleatorio", ylab="Frecuencia", main="")
summary(reforzamiento_aleatorio)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.002216 0.006268 0.007028 0.007059 0.007823 0.011787 


### CLUSTERING EMPIRICO Y ALEATORIO 
clustering_tm(bip, subsample=1) # 0.4024881
clustering_tm(bip_ppal, subsample=1) # 0.4025831

clustering_local_bip <- clustering_local_tm(bip)
rownames(clustering_local_bip) <- rownames(Matriz_afiliacion)
clustering_local_bip[id_eventos,] <- data.frame(node = 999, lc = NaN)
V(g)$clustering_local <- clustering_local_bip[dimensiones,"lc"]
summary(clustering_local_bip[,"lc"])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0000  0.0000  0.2386  0.3085  0.4622  1.0000     838 


clustering_aleatorio <- rep(NaN, 10000)
for(i in 1:length(clustering_aleatorio)) {
    red_aleatoria <- rg_tm(ni=432, np=583, ties=1780, weights=1, seed=i)
    clustering_aleatorio[i] <- clustering_tm(red_aleatoria, subsample=1)
}
hist(clustering_aleatorio, breaks = "Scott", xlab = "Clustering aleatorio", ylab="Frecuencia", main="")
summary(clustering_aleatorio)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.02165 0.02693 0.02832 0.02837 0.02973 0.03695 



##############################################################################################
### GRAFICOS DE DISTRIBUCIONES DE CENTRALIDAD NORMALIZADAS

# Distribuciones de GRADO
p1 <- ggplot(data=NULL) + 
    geom_point( aes(as.double(names(table(v_eps$grado_normalizado))), as.double(table(v_eps$grado_normalizado)/length(v_eps$grado_normalizado)), color="Eventos"), alpha=0.7, size=1.2) + 
    geom_point( aes(as.double(names(table(v_actores$grado_normalizado))), as.double(table(v_actores$grado_normalizado)/length(v_actores$grado_normalizado)), color="Actores"), alpha=0.7, size=1.2) + 
    theme_bw() + xlab("Grado simple normalizado k") + ylab("P(k)") + labs(colour = "Tipo")

plot(g, layout = l,
     vertex.label = c(V(g)$name[1:length(id_actores)], rep(NA, length(id_eventos))), 
     vertex.label.color = c(rep("#00ccff", length(id_actores)), rep("red", length(id_eventos)) ), 
     vertex.label.cex = 0.2 + V(g)$grado_normalizado*5, 
     vertex.size  = 0.5 + (100*V(g)$grado_normalizado), 
     vertex.frame.color = "#888888", edge.color = "#888888", 
     edge.width=0.3)
title("Centralidad de grado normalizada", line = -1.5)
legend("bottomright", legend=c("Actores", "Eventos"), pch=20, cex=0.8, col=c("blue", "red"), inset = c(0.19,0.12), bty = "n", pt.cex=2, y.intersp=0.4)

# Distribuciones de CERCANÍA
p2 <- ggplot(data=NULL) + 
    geom_point( aes(as.double(names(table(v_eps$c_cercania_normalizada))), as.double(table(v_eps$c_cercania_normalizada)/length(v_eps$c_cercania_normalizada)), color="Eventos"), alpha=0.7, size=1.2) + 
    geom_point( aes(as.double(names(table(v_actores$c_cercania_normalizada))), as.double(table(v_actores$c_cercania_normalizada)/length(v_actores$c_cercania_normalizada)), color="Actores"), alpha=0.7, size=1.2) + 
    theme_bw() + xlab("Centralidad de Cercanía normalizada CC") + ylab("P(CC)") + labs(colour = "Tipo")

plot(g, layout = l, 
     vertex.label = c(V(g)$name[1:length(id_actores)], rep(NA, length(id_eventos))), 
     vertex.label.color = c(rep("#00ccff", length(id_actores)), rep("red", length(id_eventos)) ), 
     vertex.label.cex = 0.2 + V(g)$c_cercania_normalizada*2, 
     vertex.size  = 0.5 + (9*V(g)$c_cercania_normalizada)^4, 
     vertex.frame.color = "#888888", edge.color = "#888888", 
     edge.width=0.3)
title("Centralidad de cercanía normalizada", line = -1.5)
legend("bottomright", legend=c("Actores", "Eventos"), pch=20, cex=0.8, col=c("blue", "red"), inset = c(0.19,0.12), bty = "n", pt.cex=2, y.intersp=0.4)

# Distribuciones de INTERMEDIACIÓN
p3 <- ggplot(data=NULL) + 
    geom_point( aes(as.double(names(table(v_eps$c_intermediacion_normalizada))), as.double(table(v_eps$c_intermediacion_normalizada)/length(v_eps$c_intermediacion_normalizada)), color="Eventos"), alpha=0.7, size=1.2) + 
    geom_point( aes(as.double(names(table(v_actores$c_intermediacion_normalizada))), as.double(table(v_actores$c_intermediacion_normalizada)/length(v_actores$c_intermediacion_normalizada)), color="Actores"), alpha=0.7, size=1.2) + 
    theme_bw() + xlab("Centralidad de intermediación normalizada CB") + ylab("P(CB)") + labs(colour = "Tipo")

plot(g, layout = l,
     vertex.label = c(V(g)$name[1:length(id_actores)], rep(NA, length(id_eventos))), 
     vertex.label.color = c(rep("#00ccff", length(id_actores)), rep("red", length(id_eventos)) ), 
     vertex.label.cex = 0.2 + V(g)$c_intermediacion_normalizada*4, 
     vertex.size  = 0.5 + (100*V(g)$c_intermediacion_normalizada), 
     vertex.frame.color = "#888888", edge.color = "#888888", 
     edge.width=0.3)
title("Centralidad de intermediación normalizada", line = -1.5)
legend("bottomright", legend=c("Actores", "Eventos"), pch=20, cex=0.8, col=c("blue", "red"), inset = c(0.19,0.12), bty = "n", pt.cex=2, y.intersp=0.4)

# Distribuciones de VECTOR PROPIO
p4 <- ggplot(data=NULL) + 
    geom_point( aes(as.double(names(table(v_eps$c_eigenvector))), as.double(table(v_eps$c_eigenvector)/length(v_eps$c_eigenvector)), color="Eventos"), alpha=0.7, size=1.2) + 
    geom_point( aes(as.double(names(table(v_actores$c_eigenvector))), as.double(table(v_actores$c_eigenvector)/length(v_actores$c_eigenvector)), color="Actores"), alpha=0.7, size=1.2) + 
    theme_bw() + xlab("Centralidad de vector propio x") + ylab("P(x)") + labs(colour = "Tipo")

plot(g, layout = l, 
     vertex.label = c(V(g)$name[1:length(id_actores)], rep(NA, length(id_eventos))), 
     vertex.label.color = c(rep("#00ccff", length(id_actores)), rep("red", length(id_eventos)) ), 
     vertex.label.cex = 0.2 + V(g)$c_eigenvector, 
     vertex.size  = 0.5 + (60*V(g)$c_eigenvector), 
     vertex.frame.color = "#888888", edge.color = "#888888", 
     edge.width=0.3)
title("Centralidad de vector propio", line = -1.5)
legend("bottomright", legend=c("Actores", "Eventos"), pch=20, cex=0.8, col=c("blue", "red"), inset = c(0.19,0.12), bty = "n", pt.cex=2, y.intersp=0.4)

multiplot(p1, p2, p3, p4, cols=2)


# Distribuciones de CLUSTERING LOCAL
p5 <- ggplot(data=NULL) + 
    geom_point( aes(as.double(names(table(v_actores$clustering_local))), as.double(table(v_actores$clustering_local)/length(v_actores$clustering_local)), color="Actores"), alpha=0.7, size=1.2) + 
    theme_bw() + xlab("Coeficiente de clustering local bipartito Cb") + ylab("P(Cb)") + labs(colour = "Tipo")

V(g)$clustering_local[is.na(V(g)$clustering_local)] <- 0

plot(g, layout = l, 
     vertex.label = c(V(g)$name[1:length(id_actores)], rep(NA, length(id_eventos))), 
     vertex.label.color = c(rep("#00ccff", length(id_actores)), rep("red", length(id_eventos)) ), 
     vertex.label.cex = 0.2 + V(g)$clustering_local*0.7, 
     vertex.size  = 0.5 + (8*V(g)$clustering_local), 
     vertex.frame.color = "#888888", edge.color = "#888888", 
     edge.width=0.3)
legend("bottomright", legend=c("Actores", "Eventos"), pch=20, cex=0.8, col=c("blue", "red"), inset = c(0.19,0.12), bty = "n", pt.cex=2, y.intersp=0.4)

