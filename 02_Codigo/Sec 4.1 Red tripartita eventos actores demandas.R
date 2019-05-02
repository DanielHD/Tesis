
library(ggplot2)
library(igraph)
library(qgraph)
source("02_Codigo/multifunction.R")

BDD <- read.csv("01_Datos/02_Datos Coaliciones de eventos/BDD simple LAOMS.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=TRUE)
actores <- read.csv("01_Datos/02_Datos Coaliciones de eventos/actores.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=TRUE)
actores_eventos <- read.csv("01_Datos/02_Datos Coaliciones de eventos/actores_eventos.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=TRUE)
demandas_eventos <- read.csv("01_Datos/02_Datos Coaliciones de eventos/demandas_eventos.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=TRUE)

id_actores <- levels(actores$Nombre.corto)
id_eventos <- levels(BDD$Id)
catalogo_demandas <- levels(as.factor(demandas_eventos$DemandaFormateada))


################################################################################################
### 4.1. Red de eventos, actores y demandas
################################################################################################

dimensiones = c(id_actores, id_eventos, catalogo_demandas)

Adj <- matrix(0L, nrow=length(dimensiones), ncol=length(dimensiones))
rownames(Adj) <- c(dimensiones)
colnames(Adj) <- c(dimensiones)

for(id_ep in id_eventos){
    demandas_ep <- as.character(demandas_eventos[demandas_eventos$Id==id_ep,"DemandaFormateada"])
    for(dem in demandas_ep){ Adj[id_ep,dem] <- Adj[dem,id_ep] <- 1 }
    actores_ep <- as.character(actores_eventos[actores_eventos$Id==id_ep,"Organización"])
    for(actor in actores_ep){ Adj[id_ep,actor] <- Adj[actor,id_ep] <- 1 }
}

g <- graph_from_adjacency_matrix(Adj, mode="undirected")
V(g)$tipo <- c(rep("actor", length(id_actores)), rep("evento", length(id_eventos)), rep("demanda", length(catalogo_demandas))) 
V(g)$color <- c(rep("blue", length(id_actores)), rep("red", length(id_eventos)), rep("yellow", length(catalogo_demandas))) 
V(g)$c_grado <- degree(g)
V(g)$c_eigen <- eigen_centrality(g)$vector
V(g)$c_intermediacion <- betweenness(g)

#######################################################################################################
# Se calcula la centralidad de cercanía, de acuerdo a la ecucación 1.12 (caso para grafo desconcectado)
matriz_distancias <- shortest.paths(g, v=V(g), to=V(g))
diag(matriz_distancias) <- NA
V(g)$c_cercania <- rowSums(1/matriz_distancias, na.rm = TRUE)
#######################################################################################################

## Grafo tripartito simple
e <- get.edgelist(g, names=FALSE)
l <- qgraph.layout.fruchtermanreingold(e, vcount=vcount(g), area=7*(vcount(g)^2), repulse.rad=(vcount(g)^3.1))

plot(g, layout=l, vertex.size=1+V(g)$c_grado*0.06, vertex.frame.color="#888888", edge.width=0.3, 
     vertex.label.cex=0.04+V(g)$c_grado*0.005,  vertex.label.size=0.04, edge.color = "#888888",
     vertex.label.color = c(rep("#00ccff", length(id_actores)), rep("red", length(id_eventos)), rep("#ff00ff", length(catalogo_demandas))) )
legend("bottomright", legend=c("Actores", "Eventos", "Demandas"), pch=20, cex=0.8, col=c("blue", "red", "yellow"), inset = c(0.15,0.05), bty = "n", pt.cex=2, y.intersp=0.4)



###############################################################################################################
#### 4.1.1. Descripciones básicas
########################################

### Componentes
componentes <- components(g)

componentes$no    # número de componentes
componentes$csize # tamaño por componente
V(g)$componente <- componentes$membership

### Descripción de componentes en la red:
componentes$membership[componentes$membership==1] # 4    unidades | 1 EP: C-075 > Queretáro > 2 actores > 1 demanda
componentes$membership[componentes$membership==2] # 1601 unidades | componente principal
componentes$membership[componentes$membership==3] # 12   unidades | 2 EPs: C-245 + C-315 > Edomex & BCS > 8 actores > 2 demandas | *CROC funciona como enlace
componentes$membership[componentes$membership==4] # 4    unidades | 1 EP: C-314 > Hidalgo > 2 actores > 1 demanda
componentes$membership[componentes$membership==5] # 7    unidades | 2 EPs: C-184 + C-041 > SLP > 3 actores > 2 demandas | *#YoSoy132-SLP funciona como enlace
componentes$membership[componentes$membership==6] # 25   unidades | 8 EPs > varios estados > 2 actores | *eventos anuales organizados por el MMM para madres centroamericanas
componentes$membership[componentes$membership==7] # 6    unidades | 1 EP: C-349 > Tabasco > 3 actores > 1 demanda
componentes$membership[componentes$membership==8] # 5    unidades | 1 EP: C-509 > Veracruz > 2 actores > 2 demandas
componentes$membership[componentes$membership==9] # 5    unidades | 1 EP: C-241 > Nuevo León > 3 actores > 1 demanda
componentes$membership[componentes$membership==10] # 12  unidades | 3 EPs: C-521 + C-523 + C-548 > Morelos > 4 actores > 5 demandas | *SUTSPJ y SUTSPLEM funcionan como enlace
componentes$membership[componentes$membership==11] # 5   unidades | 1 EP: C-489 > Guerrero > 2 actores > 2 demandas
componentes$membership[componentes$membership==12] # 5   unidades | 1 EP: C-387 > Sonora > 2 actores > 2 demandas
componentes$membership[componentes$membership==13] # 5   unidades | 1 EP: C-235 > Chiapas > 2 actores > 2 demandas
componentes$membership[componentes$membership==14] # 4   unidades | 1 EP: C-066 > Chiapas > 2 actores > 1 demanda
componentes$membership[componentes$membership==15] # 4   unidades | 1 EP: C-370 > Hidalgo > 2 actores > 1 demanda
componentes$membership[componentes$membership==16] # 7   unidades | 1 EP: C-101 > Chiapas > 2 actores > 4 demandas
componentes$membership[componentes$membership==17] # 4   unidades | 1 EP: C-362 > Morelos > 2 actores > 1 demanda
componentes$membership[componentes$membership==18] # 4   unidades | 1 EP: C-533 > *único Tamaulipas > 2 actores > 1 demanda
componentes$membership[componentes$membership==19] # 6   unidades | 1 EP: C-197 > Sonora > 2 actores > 3 demandas
componentes$membership[componentes$membership==20] # 4   unidades | 1 EP: C-082 > Jalisco > 2 actores > 1 demanda


########################################
### características y propiedades topológicas básicas de la componente principal
g_ppal <- decompose.graph(g)[[2]]

e_ppal <- get.edgelist(g_ppal, names=FALSE)
l_ppal <- qgraph.layout.fruchtermanreingold(e_ppal, vcount=vcount(g_ppal), area=7*(vcount(g_ppal)^2), repulse.rad=(vcount(g_ppal)^3.1))

# gráfico componente principal: 
plot(g_ppal, layout=l_ppal, vertex.label=NA, vertex.size=1+V(g_ppal)$c_grado*0.06, vertex.frame.color="#888888", edge.width=0.3, edge.color = "#888888")
legend("bottomright", legend=c("Actores", "Eventos", "Demandas"), pch=20, cex=0.8, col=c("blue", "red", "yellow"), inset = c(0.15,0.05), bty = "n", pt.cex=2, y.intersp=0.4)

excentricidades <- eccentricity(g_ppal, vids = V(g_ppal))
max(excentricidades) # Diámetro = 12
periferia <- excentricidades[excentricidades==12]  # 47 demandas textuales únicas + 32 actores con poca representación
min(excentricidades) # Radio = 6
centro <- excentricidades[excentricidades==6] # una único centro -> presentación con vida de los 43 de ayotzinapa

V(g)$c_grado[V(g)$name %in% names(periferia)] # la centralidad de grado es igual a uno en casi toda la pariferia
mean_distance(g_ppal) # distancia promedio = 5.732697
length(V(g))*(length(V(g))-1)*(1/sum(1/matriz_distancias, na.rm = TRUE)) # distancia promedio = 5.973377 DE ACUERDO A LA ECUACIÓN 1.3 (p.24), más de una componente

hist(excentricidades, main="", ylab="Frecuencia", xlab="Excentrididad")



###############################################################################################################
#### 4.1.2. Medidas de centralidad
###############################################################################################################
unidades_sociales <- as_data_frame(g, what="vertices")

v_actores <- unidades_sociales[unidades_sociales$tipo=="actor",]
v_demandas <- unidades_sociales[unidades_sociales$tipo=="demanda",]
v_eps <- unidades_sociales[unidades_sociales$tipo=="evento",]

######################################################
######### CENTRALIDAD DE GRADO

### DISTRIBUCIÓN DE GRADO DE ACTORES (escala normal y log-log:
p1 <- ggplot(data=NULL, aes(as.integer(names(table(v_actores$c_grado))), as.double(table(v_actores$c_grado)/length(unidades_sociales$c_grado)))) + 
    geom_point() + theme_bw() + xlab("Grado simple k") + ylab("P(k)")
p2 <- ggplot(data=NULL, aes(as.integer(names(table(v_actores$c_grado))), as.double(table(v_actores$c_grado)/length(unidades_sociales$c_grado)))) + 
    geom_point() + scale_y_log10() + scale_x_log10() + theme_bw() + xlab("log( Grado simple k )") + ylab("log( P(k) )")
multiplot(p1, p2, cols=1)

ggplot(data=NULL) +  # Distribuciones
    geom_point( aes(as.double(names(table(v_eps$c_grado))), as.double(table(v_eps$c_grado)/length(unidades_sociales$c_grado)), color="Eventos"), alpha=0.8, size=1.2) + 
    geom_point( aes(as.double(names(table(v_demandas$c_grado))), as.double(table(v_demandas$c_grado)/length(unidades_sociales$c_grado)), color="Demandas"), alpha=0.8, size=1.2) + 
    geom_point( aes(as.double(names(table(v_actores$c_grado))), as.double(table(v_actores$c_grado)/length(unidades_sociales$c_grado)), color="Actores"), alpha=0.8, size=1.2) + 
    theme_bw() + xlab("Grado simple k") + ylab("P(k)") + labs(colour = "Tipo")

######################################################
######### CENTRALIDAD DE INTERMEDIACIÓN
cor(unidades_sociales$c_grado, unidades_sociales$c_intermediacion, method="pearson")  # 0.7955758

plot(g, layout=l, vertex.size=1+V(g)$c_intermediacion*0.00005, vertex.frame.color="#888888", edge.width=0.3, edge.color = "#888888", 
     vertex.label.cex=0.04+V(g)$c_intermediacion*0.000004,  
     vertex.label.color = c(rep("#00ccff", length(id_actores)), rep("red", length(id_eventos)), rep("#ff00ff", length(catalogo_demandas))) )
legend("bottomright", legend=c("Actores", "Eventos", "Demandas"), pch=20, cex=0.8, col=c("blue", "red", "yellow"), inset = c(0.2,0.06), bty = "n", pt.cex=2, y.intersp=0.4)

ggplot(data=NULL) +  # Distribuciones
    geom_point( aes(as.double(names(table(v_eps$c_intermediacion))), as.double(table(v_eps$c_intermediacion)/length(unidades_sociales$c_intermediacion)), color="Eventos"), alpha=0.9, size=1.3) + 
    geom_point( aes(as.double(names(table(v_demandas$c_intermediacion))), as.double(table(v_demandas$c_intermediacion)/length(unidades_sociales$c_intermediacion)), color="Demandas"), alpha=0.9, size=1.3) + 
    geom_point( aes(as.double(names(table(v_actores$c_intermediacion))), as.double(table(v_actores$c_intermediacion)/length(unidades_sociales$c_intermediacion)), color="Actores"), alpha=0.9, size=1.3) + 
    theme_bw() + xlab("Centralidad de intermediación Cb") + ylab("P(Cb)") + labs(colour = "Tipo")

######################################################
######### CENTRALIDAD DE CERCANÍA
# La correlación con las otras dos centralidades es muy baja dado que el sesgo en las distribuciones no se mantiene
cor(unidades_sociales$c_grado, unidades_sociales$c_cercania, method="pearson") # 0.2929418
cor(unidades_sociales$c_intermediacion, unidades_sociales$c_cercania, method="pearson") # 0.2782752

# se exageran las centralidades, porque el rango de la centralidad en la componente principal es bastante pequeño
plot(g, layout=l, 
     vertex.size= 0.5+(0.0028*V(g)$c_cercania)^8, vertex.frame.color="#888888", edge.width=0.3, edge.color = "#888888",
     vertex.label.cex= (0.002*V(g)$c_cercania)^5, 
     vertex.label.color = c(rep("#00ccff", length(id_actores)), rep("red", length(id_eventos)), rep("#ff00ff", length(catalogo_demandas))) )
legend("bottomright", legend=c("Actores", "Eventos", "Demandas"), pch=20, cex=0.8, col=c("blue", "red", "yellow"), inset = c(0.15,0.05), bty = "n", pt.cex=2, y.intersp=0.4)

ggplot(data=NULL) +  # Distribuciones
    geom_point( aes(as.double(names(table(v_eps$c_cercania))), as.double(table(v_eps$c_cercania)/length(unidades_sociales$c_cercania)), color="Eventos"), alpha=0.8, size=1.3) + 
    geom_point( aes(as.double(names(table(v_demandas$c_cercania))), as.double(table(v_demandas$c_cercania)/length(unidades_sociales$c_cercania)), color="Demandas"), alpha=0.8, size=1.3) + 
    geom_point( aes(as.double(names(table(v_actores$c_cercania))), as.double(table(v_actores$c_cercania)/length(unidades_sociales$c_cercania)), color="Actores"), alpha=0.8, size=1.3) + 
    theme_bw() + xlab("Centralidad de cercanía Cc'") + ylab("P(Cc')") + labs(colour = "Tipo")

######################################################
######### CENTRALIDAD DE VECTOR PROPIO
plot(g, layout=l, vertex.size=0.8+V(g)$c_eigen*15, vertex.frame.color="#888888", edge.width=0.3, edge.color = "#888888", vertex.label=NA,
     vertex.label.color = c(rep("#00ccff", length(id_actores)), rep("red", length(id_eventos)), rep("#ff00ff", length(catalogo_demandas))) )
legend("bottomright", legend=c("Actores", "Eventos", "Demandas"), pch=20, cex=0.8, col=c("blue", "red", "yellow"), inset = c(0.15,0.05), bty = "n", pt.cex=2, y.intersp=0.4)

ggplot(data=NULL) + 
    geom_point( aes(as.double(names(table(v_eps$c_eigen))), as.double(table(v_eps$c_eigen)/length(unidades_sociales$c_eigen)), color="Eventos"), alpha=0.9, size=1.3) + 
    geom_point( aes(as.double(names(table(v_demandas$c_eigen))), as.double(table(v_demandas$c_eigen)/length(unidades_sociales$c_eigen)), color="Demandas"), alpha=0.9, size=1.3) + 
    geom_point( aes(as.double(names(table(v_actores$c_eigen))), as.double(table(v_actores$c_eigen)/length(unidades_sociales$c_eigen)), color="Actores"), alpha=0.9, size=1.3) + 
    theme_bw() + xlab("Centralidad de vector propio x") + ylab("P(x)") + labs(colour = "Tipo")

