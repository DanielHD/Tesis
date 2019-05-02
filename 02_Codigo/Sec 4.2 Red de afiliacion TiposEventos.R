
library(ggplot2)
library(igraph)
library(qgraph)

BDD <- read.csv("01_Datos/02_Datos Coaliciones de eventos/BDD simple LAOMS.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=TRUE)
actores <- read.csv("01_Datos/02_Datos Coaliciones de eventos/actores.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=TRUE)
actores_eventos <- read.csv("01_Datos/02_Datos Coaliciones de eventos/actores_eventos.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=TRUE)

tipos_actores <- levels(actores$Tipos)
id_eventos <- levels(BDD$Id)


################################################################################################
### 4.2. Red agregada: Tipos de actor y eventos
################################################################################################

dimensiones = c(tipos_actores, id_eventos)

Adj <- matrix(0L, nrow=length(dimensiones), ncol=length(dimensiones))
rownames(Adj) <- c(dimensiones)
colnames(Adj) <- c(dimensiones)

for(id_ep in id_eventos){
    actores_ep <- as.character(actores_eventos[actores_eventos$Id==id_ep,"Organización"])
    for(actor in actores_ep){
        tipo <- as.character(actores[actores$Nombre.corto==actor,"Tipos"])
        Adj[id_ep,tipo] <- Adj[tipo,id_ep] <- Adj[id_ep,tipo]+1  # proyección ponderada simple
    }
}

g <- graph_from_adjacency_matrix(Adj, mode="undirected", weighted = TRUE)
V(g)$tipo <- c(rep("tipo de actor", length(tipos_actores)), rep("evento", length(id_eventos))) 
V(g)$color <- c(rep("#4da6ff", length(tipos_actores)), rep("red", length(id_eventos))) 
V(g)$shape <- c(rep("circle", length(tipos_actores)), rep("square", length(id_eventos))) 
V(g)$type <- c(rep(0, length(tipos_actores)), rep(1, length(id_eventos))) 

V(g)$c_grado_simple <- apply(Adj, 1, function(x) length(x[x>0]))
V(g)$fortaleza <- apply(Adj, 1, sum)
alfa = 0.5 #### alfa es un parámetro de ajuste a los pesos para medidas de centralidad
V(g)$c_grado_ajustado <- V(g)$fortaleza^(1-alfa)*V(g)$c_grado_simple^(alfa)
V(g)$c_eigen <- eigen_centrality(g)$vector


#######################################################################################################
### Se calculan la centralidad de intermediación y cercanía para grafos ponderados. 
### Como estas medidas son dependientes de las distancias calculadas, se recalculan (también por distancias ponderadas) a partir de la matriz de asyacencia ajustada
matriz_ajustada <- 1/Adj^alfa
matriz_ajustada[is.infinite(matriz_ajustada)] <- 0

g_distancias <- graph_from_adjacency_matrix(matriz_ajustada, weighted=TRUE, mode="undirected")
matriz_distancias <- shortest.paths(g_distancias, v=V(g_distancias), to=V(g_distancias))
diag(matriz_distancias) <- NA
V(g)$c_cercania <- rowSums(1/matriz_distancias, na.rm = TRUE)

V(g)$c_intermediacion <- betweenness(g_distancias)


#######################################################################################################
### Se crea gráfico con layer a partir del algoritmo de Davidson-Harel
set.seed(26)
l <- layout_with_dh(g, coords = NULL, maxiter = 10, fineiter = max(10, log2(vcount(g))), cool.fact = 0.75, weight.node.dist = 1,
                    weight.border = 0, weight.edge.lengths = edge_density(g)/10,
                    weight.edge.crossings = 1 - sqrt(edge_density(g)),
                    weight.node.edge.dist = 0.2 * (1 - edge_density(g)))
## Se fijan manualmente algunas coordenadas para vértices en una componente aislada: 
l[369,] <- c(-135, -75)
l[371,] <- c(-135, -95)
l[49,] <- c(-135, -85)

plot(g, layout = l,
     vertex.size = 1+V(g)$c_grado_ajustado*0.065, 
     vertex.frame.color = "#AAAAAA", 
     edge.width = E(g)$weight/3, edge.color = "#777777",
     vertex.label = c(V(g)$name[1:length(tipos_actores)], rep("",length(id_eventos))),
     vertex.label.cex = c(rep(0.6+V(g)$c_grado_ajustado*0.0014,length(tipos_actores)), rep(0,length(id_eventos))),
     vertex.label.color = "#000066"
)
legend("bottomleft", legend=c("Sectores", "Eventos"), pch=c(16, 15), cex=0.7, col=c("#4da6ff", "red"), inset = c(0.28,0.07), bty = "n", pt.cex=1.2, y.intersp=0.4)



######################################################
### Descripción básica del grafo: 
#############################################
componentes <- components(g)

componentes$csize
# [1] 642   5   6

length(E(g))/(length(tipos_actores)*length(id_eventos)) # densidad = 0.02977211 

componentes$membership[componentes$membership==1] # Componente princial, 642 unidades sociales
componentes$membership[componentes$membership==2] # FE / RP     SPC   C-521   C-523   C-548
componentes$membership[componentes$membership==3] # RL / IN   C-298   C-299   C-300   C-301   C-302

#### Componente principal
g_ppal <- decompose.graph(g)[[1]]
excentricidades <- eccentricity(g_ppal, vids = V(g_ppal))
max(excentricidades) # Diámetro = 10
periferia <- excentricidades[excentricidades==10]  # AV / MIG  PX / RL
min(excentricidades) # Radio = 5
centro <- excentricidades[excentricidades==5] # C-017 C-046 C-061 C-065 C-069 C-076 C-078 C-191 C-193 C-244 C-248 C-250 C-251 C-259 C-275 C-284 C-286 C-288 C-295 C-312 C-320 C-335 C-406 C-442 C-454 C-467 C-483 C-493 C-535

V(g)$c_grado_ajustado[V(g)$name %in% names(periferia)] # Centralidad de grado ajustado en la pariferia= 8, 1
mean_distance(g_ppal) # Distancia promedio = 3.303469
length(V(g))*(length(V(g))-1)*(1/sum(1/matriz_distancias, na.rm = TRUE)) # distancia promedio DE ACUERDO A LA ECUACIÓN 1.3 (p.24) = 2.3698

hist(excentricidades, main="", ylab="Frecuencia", xlab="Excentrididad")



#############################################################################
############## BÚSQUEDA DE COMUNIDADES (DIRTLPAwb+)
#############################################################################
# Se descartan los vértices pertenecientes a las dos componentes aisladas para buscar comunidades solo en la componente principal:
Matriz_afiliacion <- Adj[tipos_actores, id_eventos]
vertices_otras_componentes <- c(names(componentes$membership[componentes$membership==2]), names(componentes$membership[componentes$membership==3]) )
Matriz_afiliacion_ppal <- Matriz_afiliacion[!rownames(Matriz_afiliacion) %in% vertices_otras_componentes, !colnames(Matriz_afiliacion) %in% vertices_otras_componentes]

source("https://raw.githubusercontent.com/sjbeckett/weighted-modularity-LPAwbPLUS/master/code/R/LPA_wb_plus.R")
comunidades_ppal <- DIRT_LPA_wb_plus(Matriz_afiliacion_ppal, 6, 10000) # Es muy lenta la deducción de comunidades!

comunidades_ppal <- list(
    Row_labels = c(1, 8, 5, 5, 3, 6, 1, 3, 8, 8, 6, 7, 1, 1, 1, 7, 6, 1, 6, 8, 5, 8, 8, 6, 6, 2, 7, 5, 5, 6, 6, 3, 1, 6, 7, 8, 6, 8, 8, 1, 5, 1, 8, 3, 8, 3, 7, 8, 5, 6, 8, 7, 5, 5, 5, 4, 5, 5, 5, 4, 7, 5, 5, 4, 3, 8, 7),
    Col_labels = c(5, 6, 6, 8, 6, 6, 6, 6, 6, 1, 5, 6, 6, 5, 5, 5, 5, 6, 5, 5, 5, 8, 8, 5, 3, 4, 4, 5, 1, 1, 1, 1, 1, 4, 4, 4, 4, 5, 6, 6, 5, 4, 5, 5, 7, 5, 7, 7, 7, 4, 4, 4, 4, 4, 4, 4, 4, 7, 7, 2, 5, 7, 7, 7, 5, 7, 2, 5, 5, 5, 4, 7, 5, 5, 6, 7, 8, 5, 4, 2, 3, 8, 6, 5, 1, 3, 5, 4, 5, 8, 6, 5, 1, 6, 7, 2, 6, 8, 2, 7, 8, 3, 2, 2, 4, 2, 2, 4, 2, 5, 2, 4, 2, 4, 4, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 4, 4, 2, 4, 2, 2, 2, 2, 2, 4, 2, 2, 2, 2, 2, 4, 2, 1, 8, 2, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 5, 2, 2, 2, 2, 2, 6, 5, 2, 2, 7, 5, 2, 5, 4, 6, 2, 5, 5, 4, 1, 8, 2, 5, 1, 5, 7, 4, 5, 4, 4, 4, 2, 5, 2, 5, 4, 4, 4, 5, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 6, 1, 5, 5, 2, 1, 1, 1, 1, 1, 4, 2, 4, 6, 1, 5, 1, 1, 4, 3, 8, 6, 6, 5, 1, 1, 7, 8, 7, 3, 1, 1, 7, 4, 8, 5, 1, 3, 7, 5, 7, 7, 2, 5, 1, 1, 8, 5, 7, 7, 7, 5, 7, 5, 7, 7, 7, 5, 2, 7, 7, 7, 5, 8, 5, 7, 5, 1, 1, 1, 7, 7, 7, 7, 7, 8, 7, 7, 5, 1, 5, 1, 1, 7, 5, 5, 7, 7, 8, 1, 7, 8, 5, 5, 4, 4, 7, 1, 1, 5, 5, 2, 7, 4, 4, 2, 5, 5, 5, 2, 3, 5, 5, 5, 4, 1, 7, 2, 2, 2, 2, 2, 6, 1, 8, 8, 4, 5, 8, 8, 1, 3, 6, 2, 2, 4, 5, 5, 8, 1, 2, 5, 7, 6, 5, 5, 2, 4, 5, 5, 4, 4, 4, 4, 5, 4, 4, 3, 4, 3, 5, 8, 5, 2, 5, 5, 2, 4, 6, 5, 5, 8, 5, 1, 5, 7, 7, 7, 7, 5, 8, 8, 7, 2, 2, 2, 1, 2, 5, 5, 1, 1, 1, 3, 8, 2, 2, 2, 5, 8, 7, 2, 5, 5, 5, 2, 6, 2, 1, 2, 2, 4, 1, 6, 1, 5, 4, 2, 1, 2, 1, 2, 4, 7, 1, 1, 1, 1, 5, 1, 1, 3, 2, 2, 2, 2, 2, 2, 5, 2, 4, 1, 7, 2, 4, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2, 4, 4, 2, 2, 2, 2, 2, 3, 2, 4, 4, 5, 2, 2, 5, 2, 5, 2, 2, 4, 2, 2, 2, 2, 2, 2, 2, 8, 4, 2, 1, 4, 2, 4, 4, 2, 2, 7, 2, 2, 1, 4, 2, 2, 4, 2, 2, 2, 2, 5, 2, 7, 2, 2, 4, 4, 2, 2, 2, 2, 2, 2, 5, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 7, 4, 2, 2, 2, 2, 5, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1, 5, 5, 1, 7, 7, 6, 6, 2, 8),
    modularity= c(0.574415)
)

source("https://raw.githubusercontent.com/sjbeckett/weighted-modularity-LPAwbPLUS/master/code/R/MODULARPLOT.R") # Gráfica modular
MODULARPLOT(Matriz_afiliacion_ppal,comunidades_ppal)

source("https://raw.githubusercontent.com/sjbeckett/weighted-modularity-LPAwbPLUS/master/code/R/GetModularInformation.R")
info_comunidades = GetModularInformation(Matriz_afiliacion_ppal,comunidades_ppal)
print(info_comunidades$normalised_modularity)  # 0.7023147
print(info_comunidades$realized_modularity)  # 0.5130533

names(comunidades_ppal$Row_labels) <- rownames(Matriz_afiliacion_ppal)
names(comunidades_ppal$Col_labels) <- colnames(Matriz_afiliacion_ppal)

# Se indican comunidades para los vértices en las componentes aisladas:
comunidad_9 <- rep(0, length(names(componentes$membership[componentes$membership==2])))
names(comunidad_9) <- names(componentes$membership[componentes$membership==2])
comunidad_10 <- rep(0, length(names(componentes$membership[componentes$membership==3])))
names(comunidad_10) <- names(componentes$membership[componentes$membership==3])

comunidades <- c(comunidades_ppal$Row_labels, comunidades_ppal$Col_labels, comunidad_9, comunidad_10)
comunidades <- comunidades[dimensiones]
V(g)$comunidad <- comunidades[dimensiones]


communities <- make_clusters(g, membership = V(g)$comunidad, algorithm = "DIRTLPAwb+", merges = NULL, modularity = F)
sizes(communities) # 0-11;  1-97;  2-165;  3-20;  4-82;  5-122;  6-42;  7-69;  8-45
length(communities) # 8
plot(communities, g, layout = l, 
     vertex.size = 1+V(g)$c_grado_ajustado*0.065, 
     vertex.frame.color = "#AAAAAA", 
     edge.width = E(g)$weight/3, 
     edge.color =c("#555555", "#ff5050")[crossing(communities,g) + 1], 
     vertex.label = c(V(g)$name[1:length(tipos_actores)], rep("",length(id_eventos))),
     vertex.label.cex = c(rep(0.6+V(g)$c_grado_ajustado*0.0014,length(tipos_actores)), rep(0,length(id_eventos))),
     vertex.label.color = "#000066",
     mark.groups = NULL, 
     vertex.color = "white")
legend("bottomleft", legend=c("Sectores", "Eventos"), pch=c(16, 15), cex=0.7, col="gray", inset = c(0.28,0.07), bty = "n", pt.cex=1.2, y.intersp=0.4)



#################################################################################
######### Análisis con tnet ###########################
#################################################################################
library(tnet)

bip_w <- as.tnet(Matriz_afiliacion_ppal, type="weighted two-mode tnet")
bip <- dichotomise_tm(bip_w, GT=0)
bip_w_completa <- as.tnet(Matriz_afiliacion, type="weighted two-mode tnet")

# REFORZAMIENTO EMPIRICO Y ALEATORIO 
reinforcement_tm(bip) # 0.1861475
reforzamiento_aleatorio <- rep(NaN, 10000)
for(i in 1:length(reforzamiento_aleatorio)) {
    red_aleatoria <- rg_tm(ni=70, np=583, ties=1215, weights=1, seed=i)
    reforzamiento_aleatorio[i] <- reinforcement_tm(red_aleatoria)
}
# se puede obtener el promedio o los intervalos de confianza (bajo hipótesis de modelo nulo aleatorio)
hist(reforzamiento_aleatorio, breaks = "Scott", xlab = "reforzamiento aleatorio", ylab="Frecuencia", main="")


### CLUSTERING EMPIRICO
clustering_tm(bip_w, subsample=1)
#     bi        am        gm        ma        mi 
# 0.7682642 0.7763050 0.7742694 0.7825669 0.7688661 
clustering_local_wbip <- clustering_local_tm(bip_w_completa)
clustering_local_wbip[is.na(clustering_local_wbip)] <- 0
V(g)$clustering_local_simple <- 0
V(g)$clustering_local_mediaAritmetica <- 0
V(g)$clustering_local_mediaGeom <- 0
V(g)[V(g)$tipo=="tipo de actor"]$clustering_local_simple <- clustering_local_wbip$lc
V(g)[V(g)$tipo=="tipo de actor"]$clustering_local_mediaAritmetica <- clustering_local_wbip$lc.am
V(g)[V(g)$tipo=="tipo de actor"]$clustering_local_mediaGeom <- clustering_local_wbip$lc.gm

plot(g, layout = l,
     vertex.size = 1+V(g)$clustering_local_mediaGeom*20, 
     vertex.frame.color = "#AAAAAA", 
     edge.width = E(g)$weight/3, edge.color = "#777777",
     vertex.label = c(V(g)$name[1:length(tipos_actores)], rep("",length(id_eventos))),
     vertex.label.cex = c(rep(0.6+V(g)$c_grado_ajustado*0.0014,length(tipos_actores)), rep(0,length(id_eventos))),
     vertex.label.color = "#000066"
)
legend("bottomleft", legend=c("Sectores", "Eventos"), pch=c(16, 15), cex=0.7, col=c("#4da6ff", "red"), inset = c(0.28,0.07), bty = "n", pt.cex=1.2, y.intersp=0.4)

### CLUSTERING EN REDES ALEATORIAS SIMULADAS
clustering_aleatorio <- rep(NaN, 10000)
for(i in 1:length(clustering_aleatorio)){
    red_aleatoria <- rg_tm(ni=70, np=583, ties=1215, weights=1, seed=i)
    clustering_aleatorio[i] <- clustering_tm(red_aleatoria)
}
hist(clustering_aleatorio, breaks = "Scott", xlab = "clustering aleatorio", ylab="Frecuencia", main="")
summary(clustering_aleatorio)
#  Min.   1st Qu. Median   Mean   3rd Qu.   Max. 
# 0.3572  0.3922  0.4001  0.4003  0.4083  0.4470

