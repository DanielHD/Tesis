
library(ggplot2)
library(igraph)
library(qgraph)

BDD <- read.csv("01_Datos/02_Datos Coaliciones de eventos/BDD simple LAOMS.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=TRUE)
BDD$C1x1FechaPeriodico <- as.Date(BDD$C1x1FechaPeriodico, "%d/%m/%Y")
BDD$C1x2FechaEvento <- as.Date(BDD$C1x2FechaEvento, "%d/%m/%Y")
actores <- read.csv("01_Datos/02_Datos Coaliciones de eventos/actores.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=TRUE)
actores_eventos <- read.csv("01_Datos/02_Datos Coaliciones de eventos/actores_eventos.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=TRUE)
demandas_eventos <- read.csv("01_Datos/02_Datos Coaliciones de eventos/demandas_eventos.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=TRUE)


tipos_actores <- levels(actores$Tipos)
id_eventos <- levels(BDD$Id)
catalogo_demandas <- levels(as.factor(demandas_eventos$DemandaFormateada))
clas_demandas <- levels(as.factor(demandas_eventos$Clasificación.en.tesis))


################################################################################################
### 4.3. Agregación simbólica: Tipos de actor y categorías de demandas
################################################################################################
dimensiones = c(tipos_actores, clas_demandas)

Adj <- matrix(0L, nrow=length(dimensiones), ncol=length(dimensiones))
rownames(Adj) <- c(dimensiones)
colnames(Adj) <- c(dimensiones)

for(id_ep in id_eventos){
    actores_ep <- as.character(actores_eventos[actores_eventos$Id==id_ep,"Organización"])
    demandas_ep <- as.character(demandas_eventos[demandas_eventos$Id==id_ep,"Clasificación.en.tesis"])
    for(actor in actores_ep){ 
        tipo <- as.character(actores[actores$Nombre.corto==actor,"Tipos"])
        for(demanda in demandas_ep){
            Adj[tipo,demanda] <- Adj[demanda,tipo] <- Adj[demanda,tipo] + 1
        }
    }
}


g <- graph_from_adjacency_matrix(Adj, mode="undirected", weighted = TRUE)
V(g)$tipo <- c(rep("sector", length(tipos_actores)), rep("clasificación demandas", length(clas_demandas))) 
V(g)$type <- c(rep(0, length(tipos_actores)), rep(1, length(clas_demandas))) 
V(g)$color <- c(rep("#4da6ff", length(tipos_actores)), rep("#ffffcc", length(clas_demandas))) 

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
## GRAFICO

# se fijan las posiciones de forma manual para todos los vétices
l <- matrix(0, ncol=2, nrow=length(dimensiones))
l[,1] <- c(rep(-1, length(tipos_actores)), rep(0, length(clas_demandas)))
l[,2] <- -c(seq(from = -1, to = 1, by=0.0285714286), seq(from = -1.25, to = 1.25, by=0.06097561))

l[c(1:12,14:25,27:50,52:60,63:70),2] <- -seq(from = -1.25, to = 1.25, by=0.03846154)

# Vértices con mayor centralidad de grado:
l[13,] <- c(1, 1) # CA
l[51,] <- c(1, 0.5) # RP
l[26,] <- c(1, 0) # ED / SPB
l[61,] <- c(1, -0.5) # SPB
l[62,] <- c(1, -1) # SPB / ED

plot(g, layout=l, 
     vertex.size=1+V(g)$c_grado_ajustado*0.1, vertex.frame.color="#777777", 
     edge.width=0.05+E(g)$weight/10, edge.color = "#999999", edge.lty = 1,
     vertex.label.cex= 0.4+V(g)$c_grado_ajustado*0.003,
     vertex.label= V(g)$name, rescale=F,
     vertex.label.color= c(rep("#000066", length(tipos_actores)), rep("#660000", length(clas_demandas)))
)
legend("bottomright", legend=c("Sector", "Clasificación de demanda"), pch=16, cex=0.7, col=c("#4da6ff", "#ffffcc"), inset = c(0.05,-0.09), bty = "n", pt.cex=1.2, y.intersp=0.4)

#######################################################################################################
## CARACTERÍSTICAS

length(E(g))/(length(tipos_actores)*length(clas_demandas)) # Densidad = 0.2195122
excentricidades <- eccentricity(g, vids = V(g))
max(excentricidades) # Diámetro = 5
periferia <- excentricidades[excentricidades==5]  # 3 dem y 17 sectores
min(excentricidades) # Radio = 3
centro <- excentricidades[excentricidades==3] # 5 dem y 2 sectores

V(g)$c_grado_ajustado[V(g)$name %in% names(periferia)] # Centralidad de grado ajustado en la pariferia
mean_distance(g) # Distancia promedio = 2.552989
length(V(g))*(length(V(g))-1)*(1/sum(1/matriz_distancias, na.rm = TRUE)) # distancia promedio DE ACUERDO A LA ECUACIÓN 1.3 (p.24) = 0.9034994

hist(excentricidades, main="", ylab="Frecuencia", xlab="Excentrididad")

#############################################################################
############## BÚSQUEDA DE COMUNIDADES (DIRTLPAwb+)
#############################################################################

matriz_afiliacion <- Adj[tipos_actores, clas_demandas]

source("https://raw.githubusercontent.com/sjbeckett/weighted-modularity-LPAwbPLUS/master/code/R/LPA_wb_plus.R")
comunidades <- DIRT_LPA_wb_plus(matriz_afiliacion, 3, 10000) # El proceso puede tardar varias horas!

# comunidades resultantes:
comunidades <- list(
    Row_labels = c(2, 1, 5, 2, 2, 2, 3, 5, 1, 1, 1, 2, 2, 2, 2, 5, 2, 4, 2, 1, 5, 1, 1, 2, 2, 5, 1, 5, 4, 1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 2, 1, 2, 1, 5, 5, 4, 4, 1, 1, 4, 4, 1, 4, 2, 5, 5, 5, 5, 5, 5, 5, 4, 2, 5, 5, 5, 2, 2, 4),
    Col_labels = c(2, 5, 5, 5, 5, 2, 1, 5, 1, 5, 5, 4, 1, 5, 2, 5, 3, 2, 5, 1, 1, 2, 4, 5, 1, 2, 1, 2, 1, 2, 4, 5, 2, 2, 5, 2, 5, 5, 5, 4, 2),
    modularity= c(0.3326181))

source("https://raw.githubusercontent.com/sjbeckett/weighted-modularity-LPAwbPLUS/master/code/R/MODULARPLOT.R")
MODULARPLOT(matriz_afiliacion, comunidades)

source("https://raw.githubusercontent.com/sjbeckett/weighted-modularity-LPAwbPLUS/master/code/R/GetModularInformation.R")
info_comunidades = GetModularInformation(matriz_afiliacion, comunidades)
print(info_comunidades$normalised_modularity)  # 0.508111
print(info_comunidades$realized_modularity)  # 0.3560021

V(g)$comunidad <- c(comunidades$Row_labels, comunidades$Col_labels)
communities <- make_clusters(g, membership = V(g)$comunidad, algorithm = "DIRTLPAwb+", merges = NULL, modularity = F)


sizes(communities)  # 13; 35; 2; 34; 27
length(communities) # 5
plot(communities, g, layout=l, 
     vertex.size=1+V(g)$c_grado_ajustado*0.1, vertex.frame.color="#777777", 
     edge.width=0.05+E(g)$weight/10, edge.lty = 1, edge.color =c("#888888", "#ff5050")[crossing(communities,g) + 1], 
     vertex.label.cex= 0.4+V(g)$c_grado_ajustado*0.003, vertex.label= V(g)$name, rescale=F,  vertex.label.color= "black", 
     mark.groups = NULL)


#################################################################################
######### Análisis con tnet ###########################
#################################################################################
library(tnet)

bip_w <- as.tnet(matriz_afiliacion, type="weighted two-mode tnet")
bip <- dichotomise_tm(bip_w, GT=0)

# REFORZAMIENTO EMPIRICO Y ALEATORIO 
reinforcement_tm(bip) # 0.4876874
reforzamiento_aleatorio <- rep(NaN, 10000)
for(i in 1:length(reforzamiento_aleatorio)) {
    red_aleatoria <- rg_tm(ni=70, np=41, ties=630, weights=1, seed=i)
    reforzamiento_aleatorio[i] <- reinforcement_tm(red_aleatoria)
}
hist(reforzamiento_aleatorio, breaks = "Scott", xlab = "reforzamiento aleatorio", ylab="Frecuencia", main="")
summary(reforzamiento_aleatorio)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2030  0.2154  0.2184  0.2186  0.2216  0.2361


### CLUSTERING EMPIRICO
clustering_tm(bip_w, subsample=1)
# bi        am        gm        ma        mi 
# 0.8827660 0.9009291 0.9072267 0.8987778 0.9024033 

clustering_local_wbip <- clustering_local_tm(bip_w)
rownames(clustering_local_wbip) <- rownames(matriz_afiliacion)
summary(clustering_local_wbip[!is.na(clustering_local_wbip["lc"]),"lc"])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7341  0.8521  0.8814  0.8714  0.8940  0.9826
mean(clustering_local_wbip[!is.na(clustering_local_wbip["lc"]),"lc"])
# [1] 0.8713832
V(g)$clustering_local <- clustering_local_wbip[dimensiones,"lc"]
V(g)$clustering_local_mediaAritmetica <- clustering_local_wbip[dimensiones,"lc.am"]
V(g)$clustering_local_mediaGeometrica <- clustering_local_wbip[dimensiones,"lc.gm"]


### CLUSTERING EN REDES ALEATORIAS SIMULADAS
clustering_aleatorio <- rep(NaN, 1000)
for(i in 1:length(clustering_aleatorio)){
    red_aleatoria <- rg_tm(ni=70, np=41, ties=630, weights=1, seed=i)
    clustering_aleatorio[i] <- clustering_tm(red_aleatoria)
}
hist(clustering_aleatorio, breaks = "Scott", xlab = "clustering aleatorio", ylab="Frecuencia", main="")
summary(clustering_aleatorio)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8370  0.8490  0.8522  0.8525  0.8560  0.8710

