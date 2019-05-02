
library(ggplot2)
library(igraph)
library(qgraph)
library(tnet)
source("02_Codigo/multifunction.R")

BDD <- read.csv("01_Datos/02_Datos Coaliciones de eventos/BDD simple LAOMS.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=TRUE)
BDD$C1x1FechaPeriodico <- as.Date(BDD$C1x1FechaPeriodico, "%d/%m/%Y")
BDD$C1x2FechaEvento <- as.Date(BDD$C1x2FechaEvento, "%d/%m/%Y")
actores <- read.csv("01_Datos/02_Datos Coaliciones de eventos/actores.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=TRUE)
actores_eventos <- read.csv("01_Datos/02_Datos Coaliciones de eventos/actores_eventos.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=TRUE)

id_actores <- levels(actores$Nombre.corto)
tipos_actores <- levels(actores$Tipos)
id_eventos <- levels(BDD$Id)


################################################################################################
### 5.2. Red proyectada de coparticipación entre actores
################################################################################################
dimensiones = c(id_actores)

### MATRIZ DE ADYACENCIA - PROYECCIÓN DE NEWMANN (colaborativa)
Adj <- matrix(0L, nrow=length(dimensiones), ncol=length(dimensiones))
rownames(Adj) <- c(dimensiones)
colnames(Adj) <- c(dimensiones)
for(id_ep in id_eventos){
    actores_ep <- as.character(actores_eventos[actores_eventos$Id==id_ep,"Organización"])
    a = 1
    for(actor1 in actores_ep){ 
        if(length(actores_ep) > a){
            a = a+1
            actores_secundarios <- actores_ep[a:length(actores_ep)]
            for(actor2 in actores_secundarios){
                Adj[actor1,actor2] <- Adj[actor2,actor1] <- Adj[actor2,actor1]+( 1/(length(actores_ep)-1) ) 
            }
        }
    }
}

### MATRIZ DE ADYACENCIA - PROYECCIÓN PONDERADA SIMPLE
Adj_ps <- matrix(0L, nrow=length(dimensiones), ncol=length(dimensiones))
rownames(Adj_ps) <- c(dimensiones)
colnames(Adj_ps) <- c(dimensiones)
for(id_ep in id_eventos){
    actores_ep <- as.character(actores_eventos[actores_eventos$Id==id_ep,"Organización"])
    a = 1
    for(actor1 in actores_ep){ 
        if(length(actores_ep) > a){
            a = a+1
            actores_secundarios <- actores_ep[a:length(actores_ep)]
            for(actor2 in actores_secundarios){
                Adj_ps[actor1,actor2] <- Adj_ps[actor2,actor1] <- Adj_ps[actor2,actor1] + 1
            }
        }
    }
}

### MATRIZ DE ADYACENCIA - PROYECCIÓN SIMPLE
Adj_s <- matrix(0L, nrow=length(dimensiones), ncol=length(dimensiones))
rownames(Adj_s) <- c(dimensiones)
colnames(Adj_s) <- c(dimensiones)
for(id_ep in id_eventos){
    actores_ep <- as.character(actores_eventos[actores_eventos$Id==id_ep,"Organización"])
    a = 1
    for(actor1 in actores_ep){ 
        if(length(actores_ep) > a){
            a = a+1
            actores_secundarios <- actores_ep[a:length(actores_ep)]
            for(actor2 in actores_secundarios){
                Adj_s[actor1,actor2] <- Adj_s[actor2,actor1] <- 1
            }
        }
    }
}

g <- graph_from_adjacency_matrix(Adj, mode="undirected", weighted = TRUE)
g_simple <- graph_from_adjacency_matrix(Adj_s, mode="undirected")
V(g_simple)$color <- V(g)$color <- "blue"
V(g_simple)$tipo1 <- V(g)$tipo1 <- as.character(actores[order(actores[,"Nombre.corto"]),]$Tipo1)
V(g_simple)$tipo2 <- V(g)$tipo2 <- as.character(actores[order(actores[,"Nombre.corto"]),]$Tipo2)
V(g_simple)$sector <- V(g)$sector <- as.character(actores[order(actores[,"Nombre.corto"]),]$Tipos)
V(g_simple)$RFOSC <- V(g)$RFOSC <- as.character(actores[order(actores[,"Nombre.corto"]),]$Aparece.en.RFOSC)
V(g_simple)$organizacion <- V(g)$organizacion <- as.character(actores[order(actores[,"Nombre.corto"]),]$Tipo.de.organización)
V(g_simple)$forma_org <- V(g)$forma_org <- as.character(actores[order(actores[,"Nombre.corto"]),]$Forma.organizativa)
V(g_simple)$alcance <- V(g)$alcance <- as.character(actores[order(actores[,"Nombre.corto"]),]$Alcance)
V(g_simple)$alcance_estatal <- V(g)$alcance_estatal <- as.character(actores[order(actores[,"Nombre.corto"]),]$Alcance_Estatal)
V(g_simple)$figura_J <- V(g)$figura_J <- as.character(actores[order(actores[,"Nombre.corto"]),]$Figura.jurídica)
V(g)$num <- c(1:432) #índice simple

V(g)$c_grado_simple <- apply(Adj, 1, function(x) length(x[x>0]))
V(g)$fortaleza <- apply(Adj, 1, sum)
alfa = 0.5 #### alfa es un parámetro de ajuste a los pesos para medidas de centralidad
V(g)$c_grado_ajustado <- (V(g)$c_grado_simple^(1-alfa)) * (V(g)$fortaleza^(alfa))
V(g)$c_eigen <- eigen_centrality(g)$vector


#######################################################################################################
### Se calculan la centralidad de intermediación y cercanía para grafos ponderados. 
### Como estas medidas son dependientes de las distancias calculadas, se recalculan (también por distancias ponderadas) a partir de la matriz de asyacencia ajustada
matriz_ajustada <- 1/(Adj^alfa)
matriz_ajustada[is.infinite(matriz_ajustada)] <- 0

g_distancias <- graph_from_adjacency_matrix(matriz_ajustada, weighted=TRUE, mode="undirected")
matriz_distancias <- shortest.paths(g_distancias, v=V(g_distancias), to=V(g_distancias))
diag(matriz_distancias) <- NA

V(g)$c_cercania <- rowSums(1/matriz_distancias, na.rm = TRUE)
V(g)$c_intermediacion <- betweenness(g_distancias)

#######################################################################################################
### Se calcula la centralidad de Katz para una proyección simple 
V(g)$c_Katz <- 0
1/max(eigen(Adj_s)$values) # inverso del máximo eigenvalor = 0.04305995
alpha <- 0.043
beta <- V(g)$c_Katz+1

sumAnt_Kantz <- 0
repeat{
    V(g)$c_Katz <- ( alpha * Adj_s %*% V(g)$c_Katz) + beta
    sum_Kantz <- sum(V(g)$c_Katz)
    if(abs(sum_Kantz-sumAnt_Kantz)< 10^(-100)){ break }
    else{sumAnt_Kantz <- sum_Kantz}
}

#######################################################################################################
### Se crea gráfico con layer a partir del algoritmo de Davidson-Harel
set.seed(23)
l <- layout_with_dh(g)

# Se ajustan manualmente algunos vértices superpuestos
l[96,] <- c(75, 60.06) # CNTE-40
l[97,] <- c(55, 61.50) # CNTE-40
l[89,] <- c(30.32, 88.89) # CNTE-18
l[81,] <- c(-48.57, -18.23) # CNPA-MLN
l[35,] <- c(-68.97, 25.48) # CAP
l[29,] <- c(-53.19, 0) # Barzón
l[146,] <- c(-51.94, 11.5) # FAC
l[11,] <- c(-45.97, 19.85) # ALCANO
l[166,] <- c(-91.03, -50.99) # FNLS
l[246,] <- c(-85.13, -55.66) # OCEZ-CP
l[351,] <- c(-77.56, -100.92) # SNTMMSSRM-65
l[119,] <- c(84.91, 0.42) # CoordUnivApoyoAyotz
l[339,] <- c(95.16, -5.55) # SNTE-55
l[292,] <- c(99.43, -8.05) # Sind_7_de_Mayo
l[295,] <- c(94.49, -10.17) # Sind_Martires1910
l[280,] <- c(96.59, 8.39) # RNG
l[192,] <- c(11.69, -63.79) # HelechoVerde
l[326,] <- c(92.84, 36.01) # SNTE-31
l[223,] <- c(62.13, -89.81) #MODEVITE
l[139,] <- c(56.93, -92.25) #Diócesis_SnCristóbal
l[108,] <- c(38.41, -33.65) #ColectivoSolecito
l[222,] <- c(0.78, -86) #MOCRI-CNPA-MN
l[203,] <- c(8.27, -89) #MAO-LN
l[38,] <- c(-53.9, 74.37) #CaravanaMadreTierra
l[73,] <- c(-24.16, 68.9) #CmteSolidaridadAyotz
l[372,] <- c(-14.69, 35.29) #STUAdeG
l[56,] <- c(-38.88, -78.84) #CE-UAEH
l[88,] <- c(-40.57, -73.62) #CNTE-15
l[400,] <- c(-76.21, 26.11) #UCEN
l[60,] <- c(-60.50, 14.01) #CIOAC	
l[164,] <- c(-49.67, 19.48) #FIOAC
l[336,] <- c(-14.44, -86.56) #SNTE-5
l[328,] <- c(-17.06, -84.47) #SNTE-35
l[330,] <- c(-15.14, -89.73) #SNTE-38
l[333,] <- c(-19.57, -87.13) #SNTE-44
l[358,] <- c(-11.34, -50.85) #SolidaridadxFam
l[319,] <- c(50.71, 80.61) #SNTE-21
l[98,] <- c(42.35, 63.43) #CNTE-9
l[24,] <- c(18.8, 37.88) #ASPA
l[307,] <- c(-32.59, 38.96) #SME


######################################################################################################
## MEDIDAS DE CENTRALIDAD Y SUS DISTRIBUCIONES
######################################################################################################

# CENTRALIDAD DE GRADO AJUSTADO
p1 <- ggplot(data=NULL) + 
    geom_point( aes(as.double(names(table(V(g)$c_grado_ajustado))), as.double(table(V(g)$c_grado_ajustado)/length(V(g)$c_grado_ajustado))), alpha=0.5, size=1.2) + 
    theme_bw() + xlab(expression("Centralidad de grado ajustado C"["D"]^"wa")) + ylab(expression(paste("P(C"["D"]^"wa",")")))

plot(g, layout = l, vertex.frame.color = "#AAAAAA", edge.color = "#999999", vertex.label.color = "black", vertex.color="#FFD12A", 
     vertex.size = 0.5 + V(g)$c_grado_ajustado*0.35, edge.width = 0.2 + E(g)$weight/3, 
     vertex.label = V(g)$name, vertex.label.cex = 0.2 + V(g)$c_grado_ajustado*0.016
)
title(expression("Centralidad de grado ajustado C"["D"]^"wa"), line = -1)


par(mfrow=c(3,2), oma = c(2, 2, 0, 0), mar = c(1, 1, 0, 0), mgp = c(2, 1, 0), xpd = NA)

# CENTRALIDAD DE GRADO SIMPLE
p2 <- ggplot(data=NULL) + 
    geom_point( aes(as.double(names(table(V(g)$c_grado_simple))), as.double(table(V(g)$c_grado_simple)/length(V(g)$c_grado_simple))), alpha=0.5, size=1.2) + 
    theme_bw() + xlab("Centralidad de grado simple k") + ylab("P(k)") + labs(colour = "Tipo")

plot(g, layout = l, edge.color = "#999999", vertex.frame.color = "#AAAAAA", vertex.label.color = "black", vertex.color="#FFD12A",
     vertex.size = 0.5 + V(g)$c_grado_simple*0.26, edge.width = 0.2 + E(g)$weight/10, 
     vertex.label = V(g)$name, vertex.label.cex = 0.2 + V(g)$c_grado_simple*0.012
)
title(expression("Centralidad de grado simple k"), line = -1)

# FORTALEZA
p3 <- ggplot(data=NULL) + 
    geom_point( aes(as.double(names(table(V(g)$fortaleza))), as.double(table(V(g)$fortaleza)/length(V(g)$fortaleza))), alpha=0.5, size=1.2) + 
    theme_bw() + xlab("Fortaleza s") + ylab("P(s)") + labs(colour = "Tipo")

plot(g, layout = l, edge.color = "#999999", vertex.frame.color = "#AAAAAA", vertex.label.color = "black", vertex.color="#FFD12A",
     vertex.size = 0.5 + V(g)$fortaleza*0.2, edge.width = 0.2 + E(g)$weight/3, 
     vertex.label = V(g)$name, vertex.label.cex = 0.2 + V(g)$fortaleza*0.01
)
title(expression("Fortaleza s"), line = -1)

# CENTRALIDAD DE EIGENVALOR
p4 <- ggplot(data=NULL) + 
    geom_point( aes(as.double(names(table(V(g)$c_eigen))), as.double(table(V(g)$c_eigen)/length(V(g)$c_eigen))), alpha=0.5, size=1.2) + 
    theme_bw() + xlab("Centralidad de vector propio x") + ylab("P(x)") + labs(colour = "Tipo")

plot(g, layout = l, edge.color = "#999999", vertex.frame.color = "#AAAAAA", vertex.label.color = "black", vertex.color="#FFD12A",
     vertex.size = 0.5 + V(g)$c_eigen*40, edge.width = 0.2 + E(g)$weight/3, 
     vertex.label = V(g)$name, vertex.label.cex = 0.2 + V(g)$c_eigen*2.3
)
title(expression("Centralidad de vector propio x"), line = -1)

# CENTRALIDAD DE KATZ
p5 <- ggplot(data=NULL) + 
    geom_point( aes(as.double(names(table(V(g)$c_Katz))), as.double(table(V(g)$c_Katz)/length(V(g)$c_Katz))), alpha=0.5, size=1.2) + 
    theme_bw() + xlab("Centralidad de Katz x'") + ylab("P(x')") + labs(colour = "Tipo")

plot(g, layout = l, edge.color = "#999999", vertex.frame.color = "#AAAAAA", vertex.label.color = "black", vertex.color="#FFD12A",
     vertex.size = 0.5 + V(g)$c_Katz*0.009, edge.width = 0.2 + E(g)$weight/10, 
     vertex.label = V(g)$name, vertex.label.cex = 0.2 + V(g)$c_Katz*0.0006
)
title(expression("Centralidad de Katz x'"), line = -1)

# CENTRALIDAD DE CERCANÍA
p6 <- ggplot(data=NULL) + 
    geom_point( aes(as.double(names(table(V(g)$c_cercania))), as.double(table(V(g)$c_cercania)/length(V(g)$c_cercania))), alpha=0.5, size=1.2) + 
    theme_bw() + xlab(expression("Centralidad de cercanía C'"["C"])) + ylab(expression(paste("P(C'"["C"],")"))) + labs(colour = "Tipo")

plot(g, layout = l, edge.color = "#999999", vertex.frame.color = "#AAAAAA", vertex.label.color = "black", vertex.color="#FFD12A",
     vertex.size = 0.5 + (V(g)$c_cercania*0.02)^3, edge.width = 0.2 + E(g)$weight/3, 
     vertex.label = V(g)$name, vertex.label.cex = 0.2 + (V(g)$c_cercania*0.006)^2
)
title(expression("Centralidad de cercanía C'"["C"]), line = -1)

# CENTRALIDAD DE INTERMEDIACIÓN
p7 <- ggplot(data=NULL) + 
    geom_point( aes(as.double(names(table(V(g)$c_intermediacion))), as.double(table(V(g)$c_intermediacion)/length(V(g)$c_intermediacion))), alpha=0.5, size=1.2) + 
    theme_bw() + xlab(expression("Centralidad de intermediación C"["B"]^"wa")) + ylab(expression(paste("P(C"["B"]^"wa",")"))) + labs(colour = "Tipo")

plot(g, layout = l, edge.color = "#999999", vertex.frame.color = "#AAAAAA", vertex.label.color = "black", vertex.color="#FFD12A",
     vertex.size = 0.5 + V(g)$c_intermediacion*0.002, edge.width = 0.2 + E(g)$weight/3, 
     vertex.label = V(g)$name, vertex.label.cex = 0.2 + V(g)$c_intermediacion*0.0001
)
title(expression("Centralidad de intermediación C"["B"]^"wa"), line = -1)

multiplot(p1, p2, p3, p4, p5, p6, p7, cols=2)
par(mfrow=c(1,1), oma = c(2, 2, 0, 0), mar = c(1, 1, 0, 0), mgp = c(2, 1, 0), xpd = NA)

######################################################################################################
## DESCRIPCIONES BÁSICAS
######################################################################################################
componentes <- components(g)
componentes$no #31
componentes$csize 
#  [1]   2 353   8   4   2   2   2   3   2   3   4   3   2   2   3   4   2   2   2   2   2   2   2   2   2   3   2   4   2   2   2
V(g)$componente <- componentes$membership

# Histograma acotado
hist(E(g)$weight, ylab="Frecuencia", xlab=expression(paste("w"["ij"]," = Nivel de acuerdo e interacción pontencial en coaliciones")), main="", breaks=1000, xlim=c(0,4.6))
length(E(g)[E(g)$weight>4.9]) # Sólo 15 aristas tienen un peso mayor o igual 5

# Histograma completo
hist(E(g)$weight, ylab="Frecuencia", xlab="Contribución marginal en coaliciones", main="", breaks=100)

m <- length(E(g))
n <- length(V(g))
m/(n*(n-1)/2) # DENSIDAD = 0.01721878

#### Componente principal
g_ppal <- decompose.graph(g)[[2]]
excentricidades <- eccentricity(g_ppal, vids = V(g_ppal))
max(excentricidades) # Diámetro = 8
periferia <- excentricidades[excentricidades==8]  # 21 actores (feministas, sindicatos de cultura y otros)
min(excentricidades) # Radio = 5 
centro <- excentricidades[excentricidades==5] # 66 actores (Sindicatos corporativistas, agrupaciones de maestros, campesinos y otros)

mean_distance(g_ppal) # Distancia promedio = 3.344064
hist(excentricidades, main="", ylab="Frecuencia", xlab="Excentricidad")


######################################################################################################
## CLUSTERING
######################################################################################################

g_tnet <- as.tnet(Adj, type="weighted one-mode tnet")

# Medidas de clustering empirico global simple y ponderado
clustering_w(g_tnet, measure="gm")
# am = 0.4036379; gm = 0.3744465; bi = 0.3652564
clustering_aleatorio_gm <- rep(NaN, 1000)
for(i in 1:length(clustering_aleatorio_gm)) {
    red_aleatoria <- rg_reshuffling_tm(bip, option="links", seed=i)
    proyeccion_aleatoria <- projecting_tm(red_aleatoria, method = "Newman")
    clustering_aleatorio_gm[i] <- clustering_w(proyeccion_aleatoria, measure="gm")
}
hist(clustering_aleatorio_gm, breaks = "Scott", xlab = "Clustering aleatorio global ponderado (media geométrica)", ylab="Frecuencia", main="")
summary(clustering_aleatorio_gm)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2214  0.2355  0.2402  0.2401  0.2446  0.2651


# Clustering local ponderado, de acuerdo a Barrat et al. (2004)
V(g)$clustering_local_w <- transitivity(g, type = c("barrat"), isolates = "zero") 
V(g)$clustering_local <- transitivity(g_simple, type = c("local"), isolates = "zero") 
plot(g, layout = l, edge.color = "#999999", vertex.frame.color = "#AAAAAA", vertex.label.color = "black", vertex.color="#FFD12A",
     vertex.size = 0.5 + V(g)$clustering_local_w*6.5, edge.width = 0.2 + E(g)$weight/3, 
     vertex.label = V(g)$name, vertex.label.cex = 0.2 + V(g)$clustering_local_w*0.25
)
hist(V(g)$clustering_local_w, breaks=20, main="", xlab="coeficiente de clustering local ponderado", ylab="Frecuencia")


################################################################################################
###### Correlación entre coeficientes  y centralidad de grado
coeficientes_ponderados <- vector(mode = "numeric", length = 76)
coeficientes_simples <- vector(mode = "numeric", length = 76)
for(i in 1:77){
    coeficientes_ponderados[i] <- mean( V(g)[V(g)$c_grado_simple==i]$clustering_local_w )
    coeficientes_simples[i] <- mean( V(g)[V(g)$c_grado_simple==i]$clustering_local )
}
ggplot(data=NULL) + 
    geom_point( aes(c(1:77), coeficientes_ponderados, color="Cw(k)"), alpha=0.7, size=1.2) + 
    geom_point( aes(c(1:77), coeficientes_simples, color="C(k)"), alpha=0.7, size=1.2) + 
    theme_bw() + xlab("Número de colaboradores (k)") + ylab(expression(paste("Promedio C"["i"], "(k), C"["iw"],"(k)"))) + labs(colour = "Clustering local")


##########################################################################################
##### HOMOFILIA
##########################################################################################

################################################################
# Homofilia por grado

grado_prom <- knn(g, vids = V(g), weights = rep(1, length(E(g)$weight)))
grado_prom_w <- knn(g, vids = V(g), weights = E(g)$weight)

V(g)$grado_prom_vecinos_w <- grado_prom$knn
V(g)$grado_prom_vecinos <- grado_prom_w$knn

# knnk
ggplot(data=NULL) + geom_abline(intercept = 0, slope = 1, linetype=3) + 
    geom_point( aes(c(1:77), grado_prom_w$knnk, color="k^w_nn(k)"), alpha=0.7, size=1.2) + 
    geom_point( aes(c(1:77), grado_prom$knnk, color="k_nn(k)"), alpha=0.7, size=1.2) + 
    theme_bw() + xlab("Número de colaboradores (k)") + ylab(expression(paste("k"["nn"]^"w", "(k), k"["nn"],"(k)"))) + 
    labs(colour = "Conectividad de colaboradores")


# knn
ggplot(data=NULL) + geom_abline(intercept = 0, slope = 1, linetype=3) + 
    geom_point( aes(V(g)$c_grado_simple, V(g)$grado_prom_vecinos, color="k^w_nn,i"), alpha=0.5, size=0.8) +
    geom_point( aes(V(g)$c_grado_simple, V(g)$grado_prom_vecinos_w, color="k_nn,i"), alpha=0.5, size=0.8) + 
    theme_bw() + xlab("Número de colaboradores (k)") + ylab(expression(paste("k"["nn,i"]^"w", ", k"["nn,i"],""))) + 
    labs(colour = "Conectividad de colaboradores") 
cor(V(g)$c_grado_simple, V(g)$grado_prom_vecinos, method = "pearson") # 0.168494
cor(V(g)$c_grado_simple, V(g)$grado_prom_vecinos_w, method = "pearson") # 0.1995411


############################################
# POR ATRIBUTOS CATEGÓRICOS EN RED PROYECTADA SIMPLE
E(g)$h_sector <- "#AAAAAA"
E(g)$h_tipos <- "#AAAAAA"
E(g)$h_organizacion <- "#AAAAAA"
E(g)$h_forma_org <- "#AAAAAA"
E(g)$h_alcance <- "#AAAAAA"
E(g)$h_alcance_estatal <- "#AAAAAA"
E(g)$h_RFOSC <- "#AAAAAA"

aristasXsector <- 0
aristasXtipos <- 0
aristasXorganizacion <- 0
aristasXforma_org <- 0
aristasXalcance <- 0
aristasXalcance_estatal <- 0
aristasXRFOSC <- 0

pesosXsector <- 0
pesosXtipos <- 0
pesosXorganizacion <- 0
pesosXforma_org <- 0
pesosXalcance <- 0
pesosXalcance_estatal <- 0
pesosXRFOSC <- 0
for(i in 1:nrow(Adj_s)){
    for(j in 1:ncol(Adj_s)){
        # sectores identificados
        if(Adj_s[i,j]==1){ 
            if(V(g_simple)$sector[i] == V(g_simple)$sector[j]){ 
                aristasXsector <- aristasXsector + 1 
                pesosXsector <- pesosXsector + Adj[i,j]
                E(g)[inc(V(g)[name==V(g)$name[i]]) & inc(V(g)[name==V(g)$name[j]])]$h_sector <- "#3399ff"
            }
            # organizacion (Discreta | Matriz | Paraguas)
            if(V(g_simple)$organizacion[i] == V(g_simple)$organizacion[j]){ 
                aristasXorganizacion <- aristasXorganizacion + 1 
                pesosXorganizacion <- pesosXorganizacion + Adj[i,j]
                E(g)[inc(V(g)[name==V(g)$name[i]]) & inc(V(g)[name==V(g)$name[j]])]$h_organizacion <- "#3399ff"
            }
            # forma organizativa (Civil | Cultural | Económico-gremial | Indígena | Matriz religiosa | Otra | Urbano-gremial)
            if(V(g_simple)$forma_org[i] == V(g_simple)$forma_org[j]){ 
                aristasXforma_org <- aristasXforma_org + 1 
                pesosXforma_org <- pesosXforma_org + Adj[i,j]
                E(g)[inc(V(g)[name==V(g)$name[i]]) & inc(V(g)[name==V(g)$name[j]])]$h_forma_org <- "#3399ff"
            }
            # alcance (Internacional | Local | Nacional | Regional)
            if(V(g_simple)$alcance[i] == V(g_simple)$alcance[j]){ 
                aristasXalcance <- aristasXalcance + 1 
                pesosXalcance <- pesosXalcance + Adj[i,j]
                E(g)[inc(V(g)[name==V(g)$name[i]]) & inc(V(g)[name==V(g)$name[j]])]$h_alcance <- "#3399ff"
            }
            # alcance estatal (Internacional | Nacional | Regional | ... )
            if(V(g_simple)$alcance_estatal[i] == V(g_simple)$alcance_estatal[j]){ 
                aristasXalcance_estatal <- aristasXalcance_estatal + 1 
                pesosXalcance_estatal <- pesosXalcance_estatal + Adj[i,j]
                E(g)[inc(V(g)[name==V(g)$name[i]]) & inc(V(g)[name==V(g)$name[j]])]$h_alcance_estatal <- "#3399ff"
            }
            # RFSCO (Sí | No)
            if(V(g_simple)$RFOSC[i] == V(g_simple)$RFOSC[j]){ 
                aristasXRFOSC <- aristasXRFOSC + 1 
                pesosXRFOSC <- pesosXRFOSC + Adj[i,j]
                E(g)[inc(V(g)[name==V(g)$name[i]]) & inc(V(g)[name==V(g)$name[j]])]$h_RFOSC <- "#3399ff"
            }
            # coincidencia en al menos un tipo de actor
            if(V(g_simple)$tipo1[i] == V(g_simple)$tipo1[j]){
                aristasXtipos <- aristasXtipos + 1 
                pesosXtipos <- pesosXtipos + Adj[i,j]
                E(g)[inc(V(g)[name==V(g)$name[i]]) & inc(V(g)[name==V(g)$name[j]])]$h_tipos <- "#3399ff"
            } 
            else if( !is.na(V(g_simple)$tipo2[i]) & !is.na(V(g_simple)$tipo2[j]) ){
                if(V(g_simple)$tipo1[i] == V(g_simple)$tipo2[j] | V(g_simple)$tipo2[i] == V(g_simple)$tipo1[j] | V(g_simple)$tipo2[i] == V(g_simple)$tipo2[j] )
                { 
                    aristasXtipos <- aristasXtipos + 1 
                    pesosXtipos <- pesosXtipos + Adj[i,j]
                    E(g)[inc(V(g)[name==V(g)$name[i]]) & inc(V(g)[name==V(g)$name[j]])]$h_tipos <- "#3399ff"
                }
            }
        }
    }
}
homofiliaEmpirica.Sectores <- 1/(2*m)*aristasXsector # 0.2202121
homofiliaEmpirica.Tipos <- 1/(2*m)*aristasXtipos # 0.3643169
homofiliaEmpirica.Organizacion <- 1/(2*m)*aristasXorganizacion # 0.4622583
homofiliaEmpirica.Forma_org <- 1/(2*m)*aristasXforma_org # 0.6556457
homofiliaEmpirica.Alcance <- 1/(2*m)*aristasXalcance # 0.621335
homofiliaEmpirica.Alcance_estatal <- 1/(2*m)*aristasXalcance_estatal # 0.4516532
homofiliaEmpirica.RFOSC <- 1/(2*m)*aristasXRFOSC # 0.834685

w <- sum(Adj)
homofiliaEmpiricaPonderada.Sectores <- 1/(2*w)*pesosXsector # 0.1973679
homofiliaEmpiricaPonderada.Tipos <- 1/(2*w)*pesosXtipos # 0.2503934
homofiliaEmpiricaPonderada.Organizacion <- 1/(2*w)*pesosXorganizacion # 0.2937625
homofiliaEmpiricaPonderada.Forma_org <- 1/(2*w)*pesosXforma_org # 0.3408906
homofiliaEmpiricaPonderada.Alcance <- 1/(2*w)*pesosXalcance # 0.3744652
homofiliaEmpiricaPonderada.Alcance_estatal <- 1/(2*w)*pesosXalcance_estatal # 0.3159425
homofiliaEmpiricaPonderada.RFOSC <- 1/(2*w)*pesosXRFOSC # 0.4328852


# RELACIONES CON HOMOFILIA SECTORIAL
plot(g, layout = l, vertex.frame.color = "#AAAAAA", edge.color = E(g)$h_sector, vertex.label.color = "black", vertex.color="#FFD12A", 
     vertex.size = 0.5 + V(g)$c_grado_ajustado*0.25, edge.width = 0.5 + E(g)$weight/5, 
     vertex.label = V(g)$name, vertex.label.cex = 0.2 + V(g)$c_grado_ajustado*0.010)
# RELACIONES CON HOMOFILIA EN AL MENOS UN TIPO DE ACTOR
plot(g, layout = l, vertex.frame.color = "#AAAAAA", edge.color = E(g)$h_tipos, vertex.label.color = "black", vertex.color="#FFD12A", 
     vertex.size = 0.5 + V(g)$c_grado_ajustado*0.25, edge.width = 0.5 + E(g)$weight/5, 
     vertex.label = V(g)$name, vertex.label.cex = 0.2 + V(g)$c_grado_ajustado*0.010)
title(sub=expression("Relaciones entre actores que comparten (al menos) un mismo tipo"), line = -3)
# RELACIONES CON HOMOFILIA EN ORGANIZACIÓN
plot(g, layout = l, vertex.frame.color = "#AAAAAA", edge.color = E(g)$h_organizacion, vertex.label.color = "black", vertex.color="#FFD12A", 
     vertex.size = 0.5 + V(g)$c_grado_ajustado*0.25, edge.width = 0.5 + E(g)$weight/5, 
     vertex.label = V(g)$name, vertex.label.cex = 0.2 + V(g)$c_grado_ajustado*0.010)
# RELACIONES CON HOMOFILIA EN FORMA ORGANIZATIVA
plot(g, layout = l, vertex.frame.color = "#AAAAAA", edge.color = E(g)$h_forma_org, vertex.label.color = "black", vertex.color="#FFD12A", 
     vertex.size = 0.5 + V(g)$c_grado_ajustado*0.25, edge.width = 0.5 + E(g)$weight/5, 
     vertex.label = V(g)$name, vertex.label.cex = 0.2 + V(g)$c_grado_ajustado*0.010)
# RELACIONES CON HOMOFILIA EN ALCANCE
plot(g, layout = l, vertex.frame.color = "#AAAAAA", edge.color = E(g)$h_alcance, vertex.label.color = "black", vertex.color="#FFD12A", 
     vertex.size = 0.5 + V(g)$c_grado_ajustado*0.25, edge.width = 0.5 + E(g)$weight/5, 
     vertex.label = V(g)$name, vertex.label.cex = 0.2 + V(g)$c_grado_ajustado*0.010)
# RELACIONES CON HOMOFILIA EN ALCANCE ESTATAL
plot(g, layout = l, vertex.frame.color = "#AAAAAA", edge.color = E(g)$h_alcance_estatal, vertex.label.color = "black", vertex.color="#FFD12A", 
     vertex.size = 0.5 + V(g)$c_grado_ajustado*0.25, edge.width = 0.5 + E(g)$weight/5, 
     vertex.label = V(g)$name, vertex.label.cex = 0.2 + V(g)$c_grado_ajustado*0.010)
title(sub=expression("Relaciones entre actores con el mismo alcance estatal"), line = -3)
# RELACIONES CON HOMOFILIA EN RFOSC
plot(g, layout = l, vertex.frame.color = "#AAAAAA", edge.color = E(g)$h_RFOSC, vertex.label.color = "black", vertex.color="#FFD12A", 
     vertex.size = 0.5 + V(g)$c_grado_ajustado*0.35, edge.width = 0.2 + E(g)$weight/3, 
     vertex.label = V(g)$name, vertex.label.cex = 0.2 + V(g)$c_grado_ajustado*0.016)

######################################
# Para evaluar la relevancia de estas medidas de homofilia, se obtienen mil simulaciones de redes aleatorias (por reordenamiento de aristas en red biparita original)
dimensiones_2 = c(id_actores, id_eventos)
Adj_bip <- matrix(0L, nrow=length(dimensiones_2), ncol=length(dimensiones_2))
rownames(Adj_bip) <- colnames(Adj_bip) <- c(dimensiones_2)
for(id_ep in id_eventos){
    actores_ep <- as.character(actores_eventos[actores_eventos$Id==id_ep,"Organización"])
    for(actor in actores_ep){ Adj_bip[id_ep,actor] <- Adj_bip[actor,id_ep] <- 1 }
}
Matriz_afiliacion <- Adj_bip[id_actores, id_eventos]
bip <- as.tnet(Matriz_afiliacion, type="binary two-mode tnet")

homofilia_aleatoria <- list("h_sector" = vector(mode="numeric", length=1000), "h_tipos" = vector(mode="numeric", length=1000),
                            "h_organizacion" = vector(mode="numeric", length=1000), "h_forma_org" = vector(mode="numeric", length=1000),
                            "h_alcance" = vector(mode="numeric", length=1000), "h_alcance_estatal" = vector(mode="numeric", length=1000), 
                            "h_RFOSC" = vector(mode="numeric", length=1000), 
                            
                            "w_RFOSC" = vector(mode="numeric", length=1000), "w_sector" = vector(mode="numeric", length=1000),
                            "w_tipos" = vector(mode="numeric", length=1000), "w_organizacion" = vector(mode="numeric", length=1000),
                            "w_forma_org" = vector(mode="numeric", length=1000), "w_alcance" = vector(mode="numeric", length=1000), 
                            "w_alcance_estatal" = vector(mode="numeric", length=1000), "w_RFOSC" = vector(mode="numeric", length=1000)
                            )
for(i in 1:1000){
    red_aleatoria <- rg_reshuffling_tm(bip, option="links", seed=i)
    proyeccion_aleatoria <- projecting_tm(red_aleatoria, method = "binary")
    proyeccion_aleatoriaPonderada <- projecting_tm(red_aleatoria, method = "Newman")
    aristasXsector <- 0
    aristasXtipos <- 0
    aristasXorganizacion <- 0
    aristasXforma_org <- 0
    aristasXalcance <- 0
    aristasXalcance_estatal <- 0
    aristasXRFOSC <- 0
    
    pesosXsector <- 0
    pesosXtipos <- 0
    pesosXorganizacion <- 0
    pesosXforma_org <- 0
    pesosXalcance <- 0
    pesosXalcance_estatal <- 0
    pesosXRFOSC <- 0
    for(arista in 1:nrow(proyeccion_aleatoria)){
        if(V(g_simple)$sector[ proyeccion_aleatoria[arista,1] ] == V(g_simple)$sector[ proyeccion_aleatoria[arista,2] ]){ aristasXsector <- aristasXsector + 1 }
        if(V(g_simple)$organizacion[ proyeccion_aleatoria[arista,1] ] == V(g_simple)$organizacion[ proyeccion_aleatoria[arista,2] ]){ aristasXorganizacion <- aristasXorganizacion + 1 }
        if(V(g_simple)$forma_org[ proyeccion_aleatoria[arista,1] ] == V(g_simple)$forma_org[ proyeccion_aleatoria[arista,2] ]){ aristasXforma_org <- aristasXforma_org + 1 }
        if(V(g_simple)$alcance[ proyeccion_aleatoria[arista,1] ] == V(g_simple)$alcance[ proyeccion_aleatoria[arista,2] ]){ aristasXalcance <- aristasXalcance + 1 }
        if(V(g_simple)$alcance_estatal[ proyeccion_aleatoria[arista,1] ] == V(g_simple)$alcance_estatal[ proyeccion_aleatoria[arista,2] ]){ aristasXalcance_estatal <- aristasXalcance_estatal + 1 }
        if(V(g_simple)$RFOSC[ proyeccion_aleatoria[arista,1] ] == V(g_simple)$RFOSC[ proyeccion_aleatoria[arista,2] ]){ aristasXRFOSC <- aristasXRFOSC + 1 }
        if(V(g_simple)$tipo1[ proyeccion_aleatoria[arista,1] ] == V(g_simple)$tipo1[ proyeccion_aleatoria[arista,2] ]){ aristasXtipos <- aristasXtipos + 1 }
        else if( !is.na(V(g_simple)$tipo2[ proyeccion_aleatoria[arista,1] ]) & !is.na(V(g_simple)$tipo2[ proyeccion_aleatoria[arista,2] ]) ){
            if(V(g_simple)$tipo1[ proyeccion_aleatoria[arista,1] ] == V(g_simple)$tipo2[ proyeccion_aleatoria[arista,2] ] | 
               V(g_simple)$tipo2[ proyeccion_aleatoria[arista,1] ] == V(g_simple)$tipo1[ proyeccion_aleatoria[arista,2] ] | 
               V(g_simple)$tipo2[ proyeccion_aleatoria[arista,1] ] == V(g_simple)$tipo2[ proyeccion_aleatoria[arista,2] ] )
            { aristasXtipos <- aristasXtipos + 1 }
        }
    }
    for(arista in 1:nrow(proyeccion_aleatoriaPonderada)){
        if(V(g_simple)$sector[ proyeccion_aleatoriaPonderada[arista,1] ] == V(g_simple)$sector[ proyeccion_aleatoriaPonderada[arista,2] ]){ aristasXsector <- aristasXsector + proyeccion_aleatoriaPonderada[arista, "w"] }
        if(V(g_simple)$organizacion[ proyeccion_aleatoriaPonderada[arista,1] ] == V(g_simple)$organizacion[ proyeccion_aleatoriaPonderada[arista,2] ]){ aristasXorganizacion <- aristasXorganizacion + proyeccion_aleatoriaPonderada[arista, "w"] }
        if(V(g_simple)$forma_org[ proyeccion_aleatoriaPonderada[arista,1] ] == V(g_simple)$forma_org[ proyeccion_aleatoriaPonderada[arista,2] ]){ aristasXforma_org <- aristasXforma_org + proyeccion_aleatoriaPonderada[arista, "w"] }
        if(V(g_simple)$alcance[ proyeccion_aleatoriaPonderada[arista,1] ] == V(g_simple)$alcance[ proyeccion_aleatoriaPonderada[arista,2] ]){ aristasXalcance <- aristasXalcance + proyeccion_aleatoriaPonderada[arista, "w"] }
        if(V(g_simple)$alcance_estatal[ proyeccion_aleatoriaPonderada[arista,1] ] == V(g_simple)$alcance_estatal[ proyeccion_aleatoriaPonderada[arista,2] ]){ aristasXalcance_estatal <- aristasXalcance_estatal + proyeccion_aleatoriaPonderada[arista, "w"] }
        if(V(g_simple)$RFOSC[ proyeccion_aleatoriaPonderada[arista,1] ] == V(g_simple)$RFOSC[ proyeccion_aleatoriaPonderada[arista,2] ]){ aristasXRFOSC <- aristasXRFOSC + proyeccion_aleatoriaPonderada[arista, "w"] }
        if(V(g_simple)$tipo1[ proyeccion_aleatoriaPonderada[arista,1] ] == V(g_simple)$tipo1[ proyeccion_aleatoriaPonderada[arista,2] ]){ aristasXtipos <- aristasXtipos + proyeccion_aleatoriaPonderada[arista, "w"] }
        else if( !is.na(V(g_simple)$tipo2[ proyeccion_aleatoriaPonderada[arista,1] ]) & !is.na(V(g_simple)$tipo2[ proyeccion_aleatoriaPonderada[arista,2] ]) ){
            if(V(g_simple)$tipo1[ proyeccion_aleatoriaPonderada[arista,1] ] == V(g_simple)$tipo2[ proyeccion_aleatoriaPonderada[arista,2] ] | 
               V(g_simple)$tipo2[ proyeccion_aleatoriaPonderada[arista,1] ] == V(g_simple)$tipo1[ proyeccion_aleatoriaPonderada[arista,2] ] | 
               V(g_simple)$tipo2[ proyeccion_aleatoriaPonderada[arista,1] ] == V(g_simple)$tipo2[ proyeccion_aleatoriaPonderada[arista,2] ] )
            { aristasXtipos <- aristasXtipos + proyeccion_aleatoriaPonderada[arista, "w"] }
        }
    }
    homofilia_aleatoria[["h_sector"]][i] <- 1/(2*nrow(proyeccion_aleatoria))*aristasXsector
    homofilia_aleatoria[["h_tipos"]][i] <- 1/(2*nrow(proyeccion_aleatoria))*aristasXtipos
    homofilia_aleatoria[["h_organizacion"]][i] <- 1/(2*nrow(proyeccion_aleatoria))*aristasXorganizacion
    homofilia_aleatoria[["h_forma_org"]][i] <- 1/(2*nrow(proyeccion_aleatoria))*aristasXforma_org
    homofilia_aleatoria[["h_alcance"]][i] <- 1/(2*nrow(proyeccion_aleatoria))*aristasXalcance
    homofilia_aleatoria[["h_alcance_estatal"]][i] <- 1/(2*nrow(proyeccion_aleatoria))*aristasXalcance_estatal
    homofilia_aleatoria[["h_RFOSC"]][i] <- 1/(2*nrow(proyeccion_aleatoria))*aristasXRFOSC
    
    homofilia_aleatoria[["w_sector"]][i] <- 1/(2*sum(proyeccion_aleatoriaPonderada$w))*aristasXsector
    homofilia_aleatoria[["w_tipos"]][i] <- 1/(2*sum(proyeccion_aleatoriaPonderada$w))*aristasXtipos
    homofilia_aleatoria[["w_organizacion"]][i] <- 1/(2*sum(proyeccion_aleatoriaPonderada$w))*aristasXorganizacion
    homofilia_aleatoria[["w_forma_org"]][i] <- 1/(2*sum(proyeccion_aleatoriaPonderada$w))*aristasXforma_org
    homofilia_aleatoria[["w_alcance"]][i] <- 1/(2*sum(proyeccion_aleatoriaPonderada$w))*aristasXalcance
    homofilia_aleatoria[["w_alcance_estatal"]][i] <- 1/(2*sum(proyeccion_aleatoriaPonderada$w))*aristasXalcance_estatal
    homofilia_aleatoria[["w_RFOSC"]][i] <- 1/(2*sum(proyeccion_aleatoriaPonderada$w))*aristasXRFOSC
}


homofilia <- data.frame( tipo = c("Sector", "Tipo de actor (1 ó 2)", "Tipo de organización", "Forma organizativa", "Alcance", "Alcance estatal", "Aparición en el RFOSC"), 
                         homofilia_empirica = c(homofiliaEmpirica.Sectores, homofiliaEmpirica.Tipos, homofiliaEmpirica.Organizacion, 
                                                homofiliaEmpirica.Forma_org, homofiliaEmpirica.Alcance, homofiliaEmpirica.Alcance_estatal, 
                                                homofiliaEmpirica.RFOSC), 
                         max_a = c(max(homofilia_aleatoria[["h_sector"]]), max(homofilia_aleatoria[["h_tipos"]]), max(homofilia_aleatoria[["h_organizacion"]]),
                                   max(homofilia_aleatoria[["h_forma_org"]]), max(homofilia_aleatoria[["h_alcance"]]), max(homofilia_aleatoria[["h_alcance_estatal"]]),
                                   max(homofilia_aleatoria[["h_RFOSC"]]) ),
                         min_a = c(min(homofilia_aleatoria[["h_sector"]]), min(homofilia_aleatoria[["h_tipos"]]), min(homofilia_aleatoria[["h_organizacion"]]),
                                   min(homofilia_aleatoria[["h_forma_org"]]), min(homofilia_aleatoria[["h_alcance"]]), min(homofilia_aleatoria[["h_alcance_estatal"]]),
                                   min(homofilia_aleatoria[["h_RFOSC"]]) ), 
                         diferencia_max = c(homofiliaEmpirica.Sectores-max(homofilia_aleatoria[["h_sector"]]),
                                            homofiliaEmpirica.Tipos-max(homofilia_aleatoria[["h_tipos"]]), 
                                            homofiliaEmpirica.Organizacion-max(homofilia_aleatoria[["h_organizacion"]]), 
                                            homofiliaEmpirica.Forma_org-max(homofilia_aleatoria[["h_forma_org"]]), 
                                            homofiliaEmpirica.Alcance-max(homofilia_aleatoria[["h_alcance"]]), 
                                            homofiliaEmpirica.Alcance_estatal-max(homofilia_aleatoria[["h_alcance_estatal"]]), 
                                            homofiliaEmpirica.RFOSC-max(homofilia_aleatoria[["h_RFOSC"]]) ))

homofilia_ponderada <- data.frame( tipo = c("Sector", "Tipo de actor (1 ó 2)", "Tipo de organización", "Forma organizativa", "Alcance", "Alcance estatal", "Aparición en el RFOSC"), 
                         homofilia_empirica = c(homofiliaEmpiricaPonderada.Sectores, homofiliaEmpiricaPonderada.Tipos, homofiliaEmpiricaPonderada.Organizacion, 
                                                homofiliaEmpiricaPonderada.Forma_org, homofiliaEmpiricaPonderada.Alcance, homofiliaEmpiricaPonderada.Alcance_estatal, 
                                                homofiliaEmpiricaPonderada.RFOSC), 
                         max_a = c(max(homofilia_aleatoria[["w_sector"]]), max(homofilia_aleatoria[["w_tipos"]]), max(homofilia_aleatoria[["w_organizacion"]]),
                                   max(homofilia_aleatoria[["w_forma_org"]]), max(homofilia_aleatoria[["w_alcance"]]), max(homofilia_aleatoria[["w_alcance_estatal"]]),
                                   max(homofilia_aleatoria[["w_RFOSC"]]) ),
                         min_a = c(min(homofilia_aleatoria[["w_sector"]]), min(homofilia_aleatoria[["w_tipos"]]), min(homofilia_aleatoria[["w_organizacion"]]),
                                   min(homofilia_aleatoria[["w_forma_org"]]), min(homofilia_aleatoria[["w_alcance"]]), min(homofilia_aleatoria[["w_alcance_estatal"]]),
                                   min(homofilia_aleatoria[["w_RFOSC"]]) ), 
                         diferencia_max = c(homofiliaEmpiricaPonderada.Sectores-max(homofilia_aleatoria[["w_sector"]]),
                                            homofiliaEmpiricaPonderada.Tipos-max(homofilia_aleatoria[["w_tipos"]]), 
                                            homofiliaEmpiricaPonderada.Organizacion-max(homofilia_aleatoria[["w_organizacion"]]), 
                                            homofiliaEmpiricaPonderada.Forma_org-max(homofilia_aleatoria[["w_forma_org"]]), 
                                            homofiliaEmpiricaPonderada.Alcance-max(homofilia_aleatoria[["w_alcance"]]), 
                                            homofiliaEmpiricaPonderada.Alcance_estatal-max(homofilia_aleatoria[["w_alcance_estatal"]]), 
                                            homofiliaEmpiricaPonderada.RFOSC-max(homofilia_aleatoria[["w_RFOSC"]]) )
)

######################################
## Figura 5.12 Homofilia por atributos en proyección simple y ponderada colaborativa.
p1 <- ggplot(homofilia, aes(x=tipo, y=homofilia_empirica) ) + 
    geom_segment(aes(x=tipo, xend=tipo, y = max_a, yend = min_a, color="Homofilia aleatoria"), size=1) + 
    geom_point( aes(tipo, diferencia_max, color="Diferencia"), alpha=0.6) + 
    geom_point( aes(color="Homofilia empirica"), alpha=0.6) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    scale_y_continuous(limits = c(0, 0.85)) + 
    xlab("") + ylab("Homofilia simple") + labs(colour = "Tipo") + 
    scale_color_manual(values=c("black", "red", "#0066ff"))

p2 <- ggplot(homofilia_ponderada, aes(x=tipo, y=homofilia_empirica) ) + 
    geom_segment(aes(x=tipo, xend=tipo, y = max_a, yend = min_a, color="Homofilia aleatoria"), size=1) + 
    geom_point( aes(tipo, diferencia_max, color="Diferencia"), alpha=0.6) + 
    geom_point( aes(color="Homofilia empirica"), alpha=0.6) + 
    scale_y_continuous(limits = c(0, 0.85)) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    xlab("Atributos covariables") + ylab("Homofilia ponderada") + labs(colour = "Tipo") + 
    scale_color_manual(values=c("black", "red", "#0066ff"))

multiplot(p1, p2, cols=1)



####################################
###### Basada en EL EFECTO DEL CLUB RICO. Utiliza "s", que en esta proyección es igual al número de coparticipaciones en coaliciones

# Proyección de Newman (puede tardar varios minutos, se puede probar con menos iteraciones para resultados similares)
club_rico <- weighted_richclub_tm(bip, NR=1000, seed=NULL, projection.method="Newman", nbins=50)

plot(club_rico$x, club_rico$y, xlab = "s = coparticipaciones en coaliciones", ylab=expression(rho^{w}*(s)), pch=20, col="red")
abline(h=1, untf = FALSE, lwd=1, lty=3)
segments(x0 = club_rico[,"x"], y0 = club_rico[,"l95"], x1 = club_rico[,"x"], y1 = club_rico[,"h95"])
par(mfrow=c(1,1))





##########################################################################################
##########################################################################################
##### SIMILARIDADES Y MODELOS DE BLOQUES
##########################################################################################
##########################################################################################

##########################################################################################
#### SIMILARIDADES ESTRUCTURALES (distancia euclidiana)
##########################################################################################

dist_eu <- matrix(0L, nrow=length(dimensiones), ncol=length(dimensiones))
rownames(dist_eu) <- c(dimensiones)
colnames(dist_eu) <- c(dimensiones)

for(i in 1:nrow(Adj)){
    for(j in 1:nrow(Adj)){
        dist_eu[i,j] = sum( (Adj[i,]-Adj[j,])^2 ) / ( sum(Adj[i,])^2 + sum(Adj[j,])^2 )
    }
}

Distancia_euclidiana <- as.dist(dist_eu[ V(g)[V(g)$componente==2]$name,V(g)[V(g)$componente==2]$name ])
hist(Distancia_euclidiana, ylab="Frecuencia", xlab="Distancia euclidiana", breaks=50, main="")
h_clustering <- hclust(Distancia_euclidiana, method="ward.D2")

plot(h_clustering, hang = 0.05,  ann = TRUE, 
     main = "",
     sub = "", 
     xlab = "", ylab = "", 
     cex = 0.19)
rect.hclust(h_clustering, k=51)
abline(h=0.437, untf = FALSE, lwd=1, lty=5, col="red")

plot(h_clustering$height,nrow(dist_eu[ V(g)[V(g)$componente==2]$name,V(g)[V(g)$componente==2]$name ]):2,type="s",xlab="Distancia estructural entre grupos", ylab= "Número de grupos", main="")
abline(v=0.437, untf = FALSE, lwd=1, lty=5, col="red")

# Establece una partición basada en un criterio único y calcula el índice de Dunn
library(clValid)
particion <- cutree(h_clustering, k=51)
clase_equivalencia <- as.factor(particion)

bloques_str <- make_clusters(g_ppal, membership = particion, algorithm = "clustering jerárquico", merges = NULL, modularity = FALSE)
V(g_ppal)$bloque_estr <- clase_equivalencia
V(g_ppal)$color <- as.character(as.factor(rainbow(length(1:51), alpha=0.8))[clase_equivalencia])
V(g_ppal)$color[V(g_ppal)$color=="#FF1E00CC"] <- "#ABCABC"


plot(g_ppal, layout = l[V(g)[V(g)$componente==2]$num,], 
     vertex.size = 0.5 + V(g_ppal)$c_grado_ajustado*0.35, edge.width = 0.2 + E(g_ppal)$weight/3, vertex.color = V(g_ppal)$color,
     vertex.label = V(g_ppal)$name, vertex.label.cex = 0.2 + V(g_ppal)$c_grado_ajustado*0.016, alpha=0.5,  mark.groups=NULL,
     edge.color =c("#444444", "#ff5050")[crossing(bloques_str,g_ppal) + 1], vertex.label.color = "black", vertex.frame.color = "#888888")

dunn(Distancia_euclidiana, particion) # 0.04187764

# Se halla el máximo valor del índice de Dunn
indice_dunn <- 0 # Inf
grupos <- 1 # 352
for(i in 2:431){
    particion <- cutree(h_clustering, k=i)
    clase_equivalencia <- as.factor(particion)
    d <- dunn(Distancia_euclidiana, particion)
    if(d>indice_dunn){
        indice_dunn <- d
        grupos<-i
        print(i)
        print(indice_dunn)
    }
}


###############################################################################################################
###############################################################################################################
###############################################################################################################
#### SIMILARIDADES REGULARES RED PROYECTADA SIMPLE)

1/max(eigen(Adj_s)$values) # inverso del máximo eigenvalor = 0.04305995
D_inv <- diag(length(dimensiones))
diag(D_inv) <- 1/V(g)$c_grado_simple
I <- diag(length(V(g)))
dimensiones = c(id_actores)

alpha <- 0.043
######################################################
######### OPCIÓN 1: SIMILARIDAD DE KATZ
similaridad_KatzSimple <- matrix(0L, nrow=length(dimensiones), ncol=length(dimensiones))
rownames(similaridad_KatzSimple) <- colnames(similaridad_KatzSimple) <- c(dimensiones)

sumAnt_Kantz <- t <- 0
repeat{
    t = t+1 # se requieren 21200 iteraciones para la convergencia de valores... el valor de alpha esta muy cerca del máximo
    similaridad_KatzSimple <- (alpha * (Adj_s %*% similaridad_KatzSimple)) + I
    sum_Kantz <- sum(similaridad_KatzSimple)
    if(abs(sum_Kantz-sumAnt_Kantz)< 10^(-100)){ break }
    else{sumAnt_Kantz <- sum_Kantz}
}
# Se comprueba el sesgo dado por el grado de cada vértice:
hist(similaridad_KatzSimple, ylab="Frecuencia", xlab="Similaridad de Katz", main="", breaks=100) 
max(as.dist(similaridad_KatzSimple)) # 52.09222

###############################################################
######### OPCIÓN 2: SIMILARIDAD DE KATZ NORMALIZADA
similaridad_KatzNorm <- matrix(0L, nrow=length(dimensiones), ncol=length(dimensiones))
rownames(similaridad_KatzNorm) <- c(dimensiones)
colnames(similaridad_KatzNorm) <- c(dimensiones) 

sumAnt_Kantz <- t <- 0
repeat{
    t = t+1 # se requieren 13 iteraciones para la convergencia de valores
    similaridad_KatzNorm <-  ( (alpha * D_inv) %*% (Adj_s %*% similaridad_KatzNorm)) + I
    sum_Kantz <- sum(similaridad_KatzNorm)
    if(abs(sum_Kantz-sumAnt_Kantz)==0){ break }
    else{sumAnt_Kantz <- sum_Kantz}
}
# Dado el pequeño número de iteraciones y valores de alpha todos los valores de similaridad son cercanos a cero
hist(similaridad_KatzNorm, ylab="Frecuencia", xlab="Similaridad de Katz normalizada", main="", breaks=100)
max(as.dist(similaridad_KatzNorm)) # 0.04307965

######################################################
######### OPCIÓN 3: SIMILARIDAD DE KATZ  (alpha= 0.5 Y t=3)
alpha_2 <- 0.5
similaridad_Katz_T3 <- matrix(0L, nrow=length(dimensiones), ncol=length(dimensiones))
rownames(similaridad_Katz_T3) <- colnames(similaridad_Katz_T3) <- c(dimensiones)

for(t in 1:3){
    similaridad_Katz_T3 <- (alpha_2 * (Adj_s %*% similaridad_KatzSimple)) + I
}
hist(as.dist(similaridad_Katz_T3), ylab="Frecuencia", xlab="Similarida de Katz (alpha= 0.5 y t=3)", main="", breaks=100)
max(as.dist(similaridad_Katz_T3)) # 0.5

######################################################
######### OPCIÓN 4: SIMILARIDAD DE KATZ NORMALIZADA (alpha= 0.5 Y t=3)
similaridad_KatzNorm_T3 <- matrix(0L, nrow=length(dimensiones), ncol=length(dimensiones))
rownames(similaridad_KatzNorm_T3) <- colnames(similaridad_KatzNorm_T3) <- c(dimensiones)

for(t in 1:3){
    similaridad_KatzNorm_T3 <- alpha_2 * D_inv %*% Adj_s %*% similaridad_KatzNorm_T3 + I
}
hist(similaridad_KatzNorm_T3, ylab="Frecuencia", xlab="Similarida de Katz normalizada (alpha= 0.5 y t=3)", main="", breaks=100)
max(as.dist(similaridad_KatzNorm_T3)) # 0.5

######################################################
######### OPCIÓN 5: SIMILARIDAD DE KATZ NORMALIZADA, CON NORMALIZACIÓN POR FILAS (descartando la diagonal ppal)

# MAX NORM:
similaridad_KatzNorm_MaxNorm <- similaridad_KatzNorm
diag(similaridad_KatzNorm_MaxNorm) <- 0
similaridad_KatzNorm_MaxNorm <- t(apply(similaridad_KatzNorm_MaxNorm, 1, function(x)(x)/(max(x))))
diag(similaridad_KatzNorm_MaxNorm) <- 1
hist(similaridad_KatzNorm_MaxNorm, ylab="Frecuencia", xlab="Similarida de Katz normalizada (norma max por filas adicinal)", main="", breaks=100)
max(as.dist(similaridad_KatzNorm_MaxNorm)) # 1

# L1 NORM:
similaridad_KatzNorm_L1Norm <- similaridad_KatzNorm
diag(similaridad_KatzNorm_L1Norm) <- 0
similaridad_KatzNorm_L1Norm <- t(apply(similaridad_KatzNorm_L1Norm, 1, function(x)(x)/(sum(x))))
hist(similaridad_KatzNorm_L1Norm, ylab="Frecuencia", xlab="Similarida de Katz normalizada (norma l1 por filas adicinal)", main="", breaks=100)
diag(similaridad_KatzNorm_L1Norm) <- 1
max(as.dist(similaridad_KatzNorm_L1Norm)) # 0.04122723

# L2 NORM:
similaridad_KatzNorm_L2Norm <- similaridad_KatzNorm
diag(similaridad_KatzNorm_L2Norm) <- 0
similaridad_KatzNorm_L2Norm <- t(apply(similaridad_KatzNorm_L2Norm, 1, function(x)(x)/(sqrt(sum(x^2)))))
diag(similaridad_KatzNorm_L2Norm) <- 1
hist(similaridad_KatzNorm_L2Norm, ylab="Frecuencia", xlab="Similarida de Katz normalizada (norma l2 por filas adicinal)", main="", breaks=100)
max(as.dist(similaridad_KatzNorm_L2Norm)) # 0.04299996


######################################################


# Calcula todos los atributos coindicentes entre pares de actores, considerando: 
#  - alcance estatal        (Nacional, regional, [ESTADO])
#  - aparición en el RFOSC  (Sí, No)
#  - forma organizativa     (civil, cultural,... )
#  - tipo de actor (1 ó 2)  (SPB, ED, CA,... )
#  - tipo de organización   (Discreta, madre, paraguas)

similaridad_attr <- matrix(0L, nrow=length(dimensiones), ncol=length(dimensiones))
rownames(similaridad_attr) <- c(dimensiones)
colnames(similaridad_attr) <- c(dimensiones)
for(i in 1:nrow(Adj_s)){
    for(j in 1:ncol(Adj_s)){
        # if(V(g_simple)$sector[i] == V(g_simple)$sector[j]){ similaridad_attr[i,j] <- similaridad_attr[i,j]+1 }
        if(V(g_simple)$organizacion[i] == V(g_simple)$organizacion[j]){ similaridad_attr[i,j] <- similaridad_attr[i,j]+1 }
        if(V(g_simple)$forma_org[i] == V(g_simple)$forma_org[j]){ similaridad_attr[i,j] <- similaridad_attr[i,j]+1 }
        # if(V(g_simple)$alcance[i] == V(g_simple)$alcance[j]){ similaridad_attr[i,j] <- similaridad_attr[i,j]+1 }
        if(V(g_simple)$alcance_estatal[i] == V(g_simple)$alcance_estatal[j]){ similaridad_attr[i,j] <- similaridad_attr[i,j]+1 }
        if(V(g_simple)$RFOSC[i] == V(g_simple)$RFOSC[j]){ similaridad_attr[i,j] <- similaridad_attr[i,j]+1 }
        if(V(g_simple)$tipo1[i] == V(g_simple)$tipo1[j]){ similaridad_attr[i,j] <- similaridad_attr[i,j]+1 }
        else if( !is.na(V(g_simple)$tipo2[i]) & !is.na(V(g_simple)$tipo2[j]) ){
            if(V(g_simple)$tipo1[i] == V(g_simple)$tipo2[j] | V(g_simple)$tipo2[i] == V(g_simple)$tipo1[j] | V(g_simple)$tipo2[i] == V(g_simple)$tipo2[j] )
            { similaridad_attr[i,j] <- similaridad_attr[i,j]+1 }
        }
        similaridad_attr[i,j] <- similaridad_attr[i,j]/5 # Es sólo el promedio de los posibles atriobutos coincidentes
    }
}
hist(as.dist(similaridad_attr), main="", ylab="Frecuencia", xlab="Similaridad por atributos")
I_Sim <- similaridad_attr

###################################################################################################
######### OPCIÓN 7: SIMILARIDAD DE KATZ NORMALIZADA APOYADA EN ATRIBUTOS
similaridad_KatzNorm_Attr <- matrix(0L, nrow=length(dimensiones), ncol=length(dimensiones))
rownames(similaridad_KatzNorm_Attr) <- c(dimensiones)
colnames(similaridad_KatzNorm_Attr) <- c(dimensiones) 

sumAnt_Kantz <- t <- 0
repeat{
    t = t+1 # se requieren 13 iteraciones para la convergencia de valores
    similaridad_KatzNorm_Attr <-  ( (alpha * D_inv) %*% (Adj_s %*% similaridad_KatzNorm_Attr)) + I_Sim
    sum_Kantz <- sum(similaridad_KatzNorm_Attr)
    if(abs(sum_Kantz-sumAnt_Kantz)==0){ break }
    else{sumAnt_Kantz <- sum_Kantz}
}
# Dado el pequeño número de iteraciones y valores de alpha todos los valores de similaridad son cercanos a cero
hist(similaridad_KatzNorm_Attr, ylab="Frecuencia", xlab="Similaridad de Katz normalizada", main="", breaks=200)
max(as.dist(similaridad_KatzNorm_Attr)) # 1.044932


###################################################################################################
######### OPCIÓN 6: SIMILARIDAD DE KATZ NORMALIZADA (alpha= 0.5 Y t=3) APOYADA EN ATRIBUTOS
similaridad_KatzNorm_T3_Attr <- matrix(0L, nrow=length(dimensiones), ncol=length(dimensiones))
rownames(similaridad_KatzNorm_T3_Attr) <- colnames(similaridad_KatzNorm_T3_Attr) <- c(dimensiones)
for(i in 1:9){
    similaridad_KatzNorm_T3_Attr <- 0.5 * D_inv %*% Adj_s %*% similaridad_KatzNorm_T3_Attr + I_Sim
}
# ...se normaliza el rango de valores:
similaridad_KatzNorm_T3_Attr <- t(apply(similaridad_KatzNorm_T3_Attr, 1, function(x)(x)/(max(x))))
# ...y se indica que los elementos en la diagonal principal son iguales a uno (sólo unos cuantos están un poco por debajo):
diag(similaridad_KatzNorm_T3_Attr) <- 1
hist(similaridad_KatzNorm_T3_Attr, ylab="Frecuencia", xlab="Similaridad de Katz", main="", breaks=50)


################################################################################################
# Crea una matriz simétrica, a partir de la anterior, en la que se elige la menor similaridad
Similaridad_Reg <- similaridad_KatzNorm_T3_Attr
for(i in 1:length(Similaridad_Reg[1,])){
    for(j in 1:length(Similaridad_Reg[1,])){
        if(Similaridad_Reg[i,j] > Similaridad_Reg[j,i] ){
            Similaridad_Reg[i,j] <- Similaridad_Reg[j,i]
        }else if(Similaridad_Reg[i,j] < Similaridad_Reg[j,i] ){
            Similaridad_Reg[j,i] <- Similaridad_Reg[i,j]
        }
    }
}
hist(Similaridad_Reg, ylab="Frecuencia", #xlab="Similaridad de Katz normalizada (alpha= 0.5, t=3) y apoyada en atributos", 
     xlab="", main="", breaks=100)

#################################################################################
#################################################################################
#################################################################################
######### PARTICIÓN EN CLASES DE EQUIVALENCIA
h_clustering_reg <- hclust(as.dist(1-Similaridad_Reg), method="ward.D2")

plot(h_clustering_reg, hang = 0.05,  ann = TRUE, 
     main = "",
     sub = "", 
     xlab = "", ylab = "", 
     cex = 0.19)
rect.hclust(h_clustering_reg, k=31)
abline(h=0.853, untf = FALSE, lwd=1, lty=5, col="red")


plot(h_clustering_reg$height,nrow(Adj):2,type="s",xlab="Distancia regular entre grupos", ylab= "Número de grupos", main="")
abline(v=0.853, untf = FALSE, lwd=1, lty=5, col="red")


# Establece una partición basada en un criterio único y calcula el índice de Dunn
particion2 <- cutree(h_clustering_reg, k=31)
clase_equivalencia_reg <- as.factor(particion2)

bloques_reg <- make_clusters(g, membership = particion2, algorithm = "clustering jerárquico", merges = NULL, modularity = FALSE)
V(g)$bloque_reg <- clase_equivalencia_reg
V(g)$color <- as.character(as.factor(rainbow(length(1:max(bloques_reg$membership)), alpha=0.7))[clase_equivalencia_reg])

plot(g, layout = l, 
     vertex.size = 1 + V(g)$c_grado_simple*0.2, edge.width = 0.3, vertex.color = V(g)$color,
     vertex.label = V(g)$name, vertex.label.cex = 0.3 + V(g)$c_grado_simple*0.008, alpha=0.5,  mark.groups=NULL,
     edge.color =c("#444444", "#ff5050")[crossing(bloques_reg,g) + 1], vertex.label.color = "black", vertex.frame.color = "#888888"
)

dunn(1-Similaridad_Reg, particion2) # 0.1108681


# Se halla el máximo valor del índice de Dunn
indice_dunn <- 0 # MAX = Inf
grupos <- 1 # MAX = 343
for(i in 2:431){
    particion_2 <- cutree(h_clustering_reg, k=i)
    clase_equivalencia <- as.factor(particion_2)
    d <- dunn(1-Similaridad_Reg, particion_2)
    if(d>indice_dunn){
        indice_dunn <- d
        grupos <-i
        print(i)
        print(indice_dunn)
    }
}



##########################################################################################
##### INFOMAP (COMUNIDADES DETERMINADAS POR FLUJO)
##########################################################################################
# Los resultados fueron obtenidos por Infomap para 2 niveles en 100 mil intentos, a partir del código descargado en el sitio http://www.mapequation.org/code.html
# Aunque informap sí está implementado en igraph, los resultados varían, además los autores sugieren sólo utilizar el código que ellos proporcionan

V(g_ppal)$comunidades_informap <- c(25, 3, 14, 8, 32, 21, 14, 6, 7, 17, 12, 5, 5, 21, 13, 1, 4, 26, 4, 4, 4, 4, 4, 6, 11, 18, 11, 7, 25, 
                                    18, 1, 13, 26, 11, 1, 5, 8, 8, 8, 1, 4, 11, 29, 1, 4, 3, 14, 1, 11, 4, 17, 7, 8, 4, 8, 9, 4, 12, 5, 8, 
                                    4, 16, 13, 5, 1, 4, 3, 5, 5, 16, 4, 12, 3, 27, 1, 3, 28, 3, 11, 3, 15, 9, 2, 2, 3, 4, 5, 7, 13, 17, 
                                    21, 27, 9, 17, 8, 1, 35, 5, 13, 5, 11, 1, 4, 32, 4, 18, 4, 8, 3, 4, 30, 4, 8, 27, 27, 27, 27, 20, 8, 
                                    5, 14, 8, 7, 3, 8, 17, 4, 18, 3, 4, 6, 11, 8, 34, 25, 1, 21, 10, 5, 7, 8, 13, 3, 1, 13, 28, 5, 18, 15, 
                                    3, 18, 8, 4, 15, 18, 1, 3, 17, 8, 1, 3, 35, 35, 8, 21, 14, 36, 7, 32, 13, 21, 17, 24, 5, 30, 11, 34, 
                                    20, 13, 9, 15, 24, 11, 28, 9, 3, 36, 32, 20, 11, 9, 17, 14, 11, 4, 17, 1, 14, 19, 7, 36, 8, 4, 4, 21, 
                                    4, 6, 8, 14, 5, 36, 17, 29, 29, 3, 11, 21, 28, 18, 4, 6, 17, 8, 14, 21, 5, 14, 14, 29, 10, 19, 4, 11, 
                                    1, 21, 18, 18, 4, 11, 9, 35, 14, 19, 4, 20, 20, 13, 1, 11, 3, 22, 22, 22, 11, 4, 4, 4, 9, 4, 22, 27, 31, 
                                    34, 3, 9, 19, 23, 5, 22, 12, 22, 20, 33, 19, 3, 10, 31, 23, 20, 3, 33, 15, 3, 10, 4, 4, 24, 4, 24, 4, 
                                    24, 15, 4, 6, 14, 15, 3, 4, 5, 3, 1, 5, 4, 1, 4, 4, 15, 4, 4, 3, 1, 1, 1, 1, 3, 1, 4, 4, 4, 15, 4, 29, 
                                    21, 4, 3, 25, 6, 5, 26, 4, 5, 7, 5, 4, 16, 5, 4, 29, 5, 13, 4, 7, 32, 1, 9, 4, 11, 13, 13, 4, 28, 15, 13)
V(g_ppal)$comunidades_flujo <- c(0.000600601, 0.0012012, 0.000600601, 0.0012012, 0.000600601, 0.000600601, 0.0018018, 0.0024024, 0.0012012, 0.000600601, 0.000600601, 
                            0.0012012, 0.0012012, 0.0018018, 0.000600601, 0.0048048, 0.0024024, 0.0018018, 0.0042042, 0.0024024, 0.000600601, 0.0042042, 
                            0.000600601, 0.0222222, 0.0012012, 0.000600601, 0.000600601, 0.0012012, 0.0036036, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 
                            0.000600601, 0.000600601, 0.00600601, 0.000600601, 0.000600601, 0.0012012, 0.0018018, 0.000600601, 0.0018018, 0.000600601, 0.00660661, 
                            0.0012012, 0.0024024, 0.0012012, 0.0558559, 0.000600601, 0.0018018, 0.0018018, 0.0108108, 0.000600601, 0.000600601, 0.000600601, 
                            0.0024024, 0.000600601, 0.0018018, 0.0012012, 0.0024024, 0.000600601, 0.0042042, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 
                            0.000600601, 0.000600601, 0.0312312, 0.003003, 0.000600601, 0.012012, 0.0018018, 0.0012012, 0.0162162, 0.0258258, 0.0018018, 0.03003, 
                            0.003003, 0.00540541, 0.003003, 0.00660661, 0.0726727, 0.0786787, 0.0138138, 0.000600601, 0.0012012, 0.00900901, 0.000600601, 
                            0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.0012012, 0.0012012, 0.000600601, 0.0012012, 0.000600601, 0.000600601, 0.0012012, 
                            0.0012012, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.0018018, 0.0024024, 0.000600601, 0.000600601, 0.000600601, 0.0024024, 
                            0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.0018018, 0.0012012, 
                            0.000600601, 0.0102102, 0.0012012, 0.000600601, 0.000600601, 0.003003, 0.0012012, 0.0048048, 0.000600601, 0.012012, 0.000600601, 
                            0.000600601, 0.000600601, 0.000600601, 0.0138138, 0.000600601, 0.000600601, 0.00900901, 0.0024024, 0.0018018, 0.000600601, 0.000600601, 
                            0.0012012, 0.0042042, 0.000600601, 0.00600601, 0.0018018, 0.0036036, 0.0036036, 0.0012012, 0.000600601, 0.003003, 0.000600601, 
                            0.000600601, 0.000600601, 0.0024024, 0.0012012, 0.000600601, 0.0114114, 0.000600601, 0.000600601, 0.000600601, 0.0036036, 0.000600601, 
                            0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.0012012, 0.0102102, 0.0018018, 0.000600601, 
                            0.0012012, 0.0012012, 0.000600601, 0.000600601, 0.000600601, 0.0018018, 0.0018018, 0.0018018, 0.0012012, 0.0012012, 0.000600601, 
                            0.0012012, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.0012012, 0.000600601, 0.000600601, 0.021021, 0.0012012, 
                            0.000600601, 0.0036036, 0.000600601, 0.003003, 0.0012012, 0.0012012, 0.0012012, 0.000600601, 0.0036036, 0.0012012, 0.0018018, 
                            0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 
                            0.000600601, 0.0018018, 0.0012012, 0.000600601, 0.0012012, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.0018018, 
                            0.000600601, 0.0012012, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.0012012, 0.000600601, 0.000600601, 
                            0.000600601, 0.000600601, 0.0018018, 0.000600601, 0.0012012, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.0018018, 0.0012012, 
                            0.000600601, 0.0018018, 0.000600601, 0.000600601, 0.0012012, 0.0192192, 0.0036036, 0.000600601, 0.000600601, 0.0012012, 0.0024024, 
                            0.0012012, 0.0012012, 0.0018018, 0.003003, 0.0036036, 0.0012012, 0.0012012, 0.0024024, 0.0012012, 0.0024024, 0.0018018, 0.003003, 
                            0.000600601, 0.00840841, 0.0012012, 0.003003, 0.0018018, 0.000600601, 0.0018018, 0.000600601, 0.000600601, 0.00720721, 0.000600601, 
                            0.0012012, 0.000600601, 0.003003, 0.0012012, 0.000600601, 0.0012012, 0.000600601, 0.000600601, 0.0024024, 0.000600601, 0.0012012, 
                            0.000600601, 0.000600601, 0.000600601, 0.0048048, 0.0018018, 0.000600601, 0.0114114, 0.0012012, 0.0018018, 0.000600601, 0.000600601, 
                            0.000600601, 0.0102102, 0.00540541, 0.018018, 0.0018018, 0.000600601, 0.000600601, 0.000600601, 0.0042042, 0.000600601, 0.0018018, 
                            0.0018018, 0.000600601, 0.000600601, 0.0012012, 0.000600601, 0.000600601, 0.0012012, 0.000600601, 0.0132132, 0.000600601, 0.003003, 
                            0.000600601, 0.003003, 0.0012012, 0.0018018, 0.000600601, 0.003003, 0.0036036, 0.000600601, 0.0012012, 0.0018018, 0.000600601, 
                            0.00960961, 0.0114114, 0.000600601, 0.000600601, 0.0024024, 0.0012012, 0.000600601, 0.000600601, 0.0042042, 0.000600601, 0.000600601, 
                            0.000600601, 0.000600601)

hist(V(g_ppal)$comunidades_flujo, breaks=50, xlab="Flujo de actividad en protestas por actor", ylab="Frecuencia", main="")
infomap <- make_clusters(g_ppal, membership = V(g_ppal)$comunidades_informap, algorithm = "infomap", merges = NULL, modularity = FALSE)
V(g_ppal)$color <- as.character(as.factor(rainbow(length(1:max(V(g_ppal)$comunidades_informap)), alpha=0.8))[V(g_ppal)$comunidades_informap])

sizes(infomap) 
length(infomap) # 36
plot(g_ppal, layout = l[V(g)[V(g)$componente==2]$num,], vertex.label.color = "black", vertex.frame.color = "#888888", 
     vertex.size = 1+V(g_ppal)$comunidades_flujo*250, edge.width = 0.2+E(g_ppal)$weight/3, 
     vertex.label = V(g_ppal)$name, vertex.label.cex = 0.2+V(g_ppal)$comunidades_flujo*15,
     edge.color =c("#444444", "#ff5050")[crossing(infomap,g_ppal) + 1],
     vertex.color = V(g_ppal)$color)


########## Red agregada. 
modulos_infomap <- data.frame(name = c('CETEG,...', 'CNTE-7,...', 'CNTE-22,...', 'SME,...', 'CNPA-MLN,...', 'Barzón,...', 'UNTA,...', '¡HastaEncontrarlos!,...',
                             'CNTE-36,...', 'SNTE-42,...', 'CNTE-23,...', 'CNTE,...', '#YoSoy132-CDMX,...', 'AI,...', 'FPLZ,...', 'CNC,...', 'CILAS,...',
                             'CPSCDF,...', 'SNTE-2,...', 'SNTE-31,...', 'ANAD,...', 'SITET,...', 'SNTE-21,...', 'MENA,...', 'CAP,...', 'UCIZONI,...', 
                             'SNTE-10,...', 'MMB,...', 'Tetela,...', 'CUL,...', 'SNTE-13,...', 'MMSP,...', 'SNTE-33,...', 'MBMH,...', 'CONOC,...', 'IPLP,...'), 
                      flujo_dentro =c(0.166366, 0.151351, 0.114114, 0.111712, 0.0864865, 0.0576577, 0.0516517, 0.0234234, 0.0204204, 0.018018, 0.0174174, 
                                      0.0168168, 0.0156156, 0.0126126, 0.012012, 0.0102102, 0.00960961, 0.00960961, 0.00900901, 0.00840841, 0.00840841, 
                                      0.00660661, 0.00660661, 0.00600601, 0.00540541, 0.00540541, 0.00540541, 0.00540541, 0.0048048, 0.0042042, 0.0036036, 
                                      0.0036036, 0.0036036, 0.003003, 0.003003, 0.0024024), 
                      flujo_fuera =c(0.0142643, 0.00892559, 0.0305556, 0.0294504, 0.0285886, 0.0124925, 0.0133734, 0.003003, 0.00643501, 0.0012012, 0.00564136, 
                                     0.00980981, 0.0048048, 0.0024024, 0.00336336, 0.0012012, 0.0032032, 0.003003, 0.000600601, 0.0018018, 0.0018018, 0.000600601, 
                                     0.0012012, 0.0024024, 0.003003, 0.0024024, 0.000600601, 0.0027027, 0.0012012, 0.000600601, 0.0012012, 0.000960961, 0.000600601, 
                                     0.0012012, 0.0012012, 0.000600601), 
                      color = c('#FF000099', '#FF2A0099', '#FF550099', '#FF800099', '#FFAA0099', '#FFD50099', '#FFFF0099', '#D5FF0099', '#AAFF0099', '#80FF0099', 
                                '#55FF0099', '#2BFF0099', '#00FF0099', '#00FF2A99', '#00FF5599', '#00FF8099', '#00FFAA99', '#00FFD499', '#00FFFF99', '#00D4FF99', 
                                '#00AAFF99', '#0080FF99', '#0055FF99', '#002BFF99', '#0000FF99', '#2A00FF99', '#5500FF99', '#8000FF99', '#AA00FF99', '#D400FF99', 
                                '#FF00FF99', '#FF00D599', '#FF00AA99', '#FF008099', '#FF005599', '#FF002B99')
                      )

aristas_infomap <- data.frame(from = c('CETEG,...', 'CETEG,...', 'CETEG,...', 'CETEG,...', 'CETEG,...', 'CETEG,...', 'CETEG,...', 'CETEG,...', 'CETEG,...', 'CETEG,...',
                               'CETEG,...', 'CNTE-7,...', 'CNTE-7,...', 'CNTE-7,...', 'CNTE-7,...', 'CNTE-7,...', 'CNTE-7,...', 'CNTE-7,...', 'CNTE-7,...', 
                               'CNTE-7,...', 'CNTE-22,...', 'CNTE-22,...', 'CNTE-22,...', 'CNTE-22,...', 'CNTE-22,...', 'CNTE-22,...', 'CNTE-22,...', 'CNTE-22,...',
                               'CNTE-22,...', 'CNTE-22,...', 'CNTE-22,...', 'CNTE-22,...', 'CNTE-22,...', 'CNTE-22,...', 'CNTE-22,...', 'CNTE-22,...', 'CNTE-22,...',
                               'CNTE-22,...', 'CNTE-22,...', 'CNTE-22,...', 'SME,...', 'SME,...', 'SME,...', 'SME,...', 'SME,...', 'SME,...', 'SME,...', 'SME,...',
                               'SME,...', 'SME,...', 'SME,...', 'SME,...', 'SME,...', 'SME,...', 'SME,...', 'CNPA-MLN,...', 'CNPA-MLN,...', 'CNPA-MLN,...',
                               'CNPA-MLN,...', 'CNPA-MLN,...', 'CNPA-MLN,...', 'CNPA-MLN,...', 'CNPA-MLN,...', 'CNPA-MLN,...', 'CNPA-MLN,...', 'CNPA-MLN,...',
                               'CNPA-MLN,...', 'CNPA-MLN,...', 'CNPA-MLN,...', 'CNPA-MLN,...', 'CNPA-MLN,...', 'CNPA-MLN,...', 'CNPA-MLN,...', 'Barzón,...',
                               'Barzón,...', 'Barzón,...', 'Barzón,...', 'Barzón,...', 'Barzón,...', 'UNTA,...', 'UNTA,...', 'UNTA,...', 'UNTA,...', 
                               '¡HastaEncontrarlos!,...', '¡HastaEncontrarlos!,...', 'CNTE-36,...', 'CNTE-36,...', 'SNTE-42,...', 'CNTE-23,...', 'CNTE-23,...',
                               'CNTE-23,...', 'CNTE-23,...', 'CNTE-23,...', 'CNTE-23,...', 'CNTE-23,...', 'CNTE-23,...', 'CNTE,...', 'CNTE,...', 'CNTE,...',
                               'CNTE,...', 'CNTE,...', 'CNTE,...', 'CNTE,...', 'CNTE,...', '#YoSoy132-CDMX,...', 'AI,...', 'AI,...', 'FPLZ,...', 'FPLZ,...',
                               'CNC,...', 'CILAS,...', 'CILAS,...', 'SNTE-31,...', 'SNTE-31,...', 'SNTE-31,...', 'ANAD,...', 'SITET,...', 'SNTE-13,...'), 
                      to = c('CNTE-7,...', 'CNTE-22,...', 'SME,...', 'CNPA-MLN,...', 'CNTE-36,...', 'CNTE-23,...', 'AI,...', 'SNTE-31,...', 'MMB,...', 'SNTE-13,...', 
                             'MBMH,...', 'CNTE-22,...', 'SME,...', 'CNTE-36,...', 'SNTE-42,...', 'CNTE-23,...', 'SNTE-31,...', 'MENA,...', 'MMB,...', 'SNTE-13,...', 
                             'SME,...', 'CNPA-MLN,...', '¡HastaEncontrarlos!,...', 'CNTE-36,...', 'SNTE-42,...', 'CNTE-23,...', 'CNTE,...', '#YoSoy132-CDMX,...', 
                             'FPLZ,...', 'SNTE-31,...', 'SITET,...', 'SNTE-21,...', 'MENA,...', 'UCIZONI,...', 'SNTE-10,...', 'MMB,...', 'CUL,...', 'SNTE-13,...', 
                             'SNTE-33,...', 'MBMH,...', 'CNPA-MLN,...', 'Barzón,...', 'UNTA,...', 'CNTE-36,...', 'CNTE-23,...', 'CNTE,...', '#YoSoy132-CDMX,...', 
                             'FPLZ,...', 'CILAS,...', 'CPSCDF,...', 'ANAD,...', 'MENA,...', 'MMB,...', 'CUL,...', 'MMSP,...', 'Barzón,...', 'UNTA,...', 
                             '¡HastaEncontrarlos!,...', 'CNTE-36,...', 'CNTE-23,...', '#YoSoy132-CDMX,...', 'FPLZ,...', 'CNC,...', 'CILAS,...', 'CPSCDF,...', 
                             'ANAD,...', 'SITET,...', 'SNTE-21,...', 'CAP,...', 'Tetela,...', 'MBMH,...', 'CONOC,...', 'IPLP,...', 'UNTA,...', '#YoSoy132-CDMX,...', 
                             'AI,...', 'FPLZ,...', 'CILAS,...', 'CONOC,...', 'CNTE,...', '#YoSoy132-CDMX,...', 'CILAS,...', 'CAP,...', 'CPSCDF,...', 'UCIZONI,...', 
                             '#YoSoy132-CDMX,...', 'Tetela,...', 'CNTE,...', 'FPLZ,...', 'SNTE-31,...', 'SITET,...', 'SNTE-21,...', 'MMB,...', 'Tetela,...', 
                             'SNTE-13,...', 'MBMH,...', '#YoSoy132-CDMX,...', 'AI,...', 'CILAS,...', 'CPSCDF,...', 'SNTE-2,...', 'SNTE-31,...', 'SNTE-21,...', 
                             'MENA,...', 'Tetela,...', 'CILAS,...', 'ANAD,...', 'SITET,...', 'SNTE-21,...', 'CAP,...', 'CPSCDF,...', 'MMSP,...', 'MMB,...', 
                             'SNTE-13,...', 'MBMH,...', 'MENA,...', 'SNTE-21,...', 'MBMH,...'), 
                      flujo = c(0.00129177, 0.00590591, 0.00395812, 0.00015015, 0.000759688, 0.000311025, 0.000600601, 0.000225225, 0.000675676, 0.00023595, 
                                0.00015015, 0.0047047, 0.00113447, 0.000572001, 0.0002002, 0.00023595, 0.00015015, 0.0004004, 0.00015015, 8.58E-05, 0.00610968, 
                                0.0003003, 0.000900901, 0.00296368, 0.0004004, 0.00109395, 0.0003003, 0.000600601, 0.0003003, 0.00045045, 0.0003003, 0.0003003, 
                                0.000600601, 0.0018018, 0.000600601, 0.00112613, 0.0004004, 0.000568426, 0.000600601, 0.000225225, 0.00237237, 0.0012012, 0.00136136, 
                                0.000487988, 0.0018876, 0.0041041, 0.00248248, 0.00048048, 0.00012012, 0.0012012, 0.000533867, 0.000734067, 0.000600601, 0.0002002, 
                                0.00048048, 0.00737738, 0.00802803, 0.00105105, 0.00105105, 0.000975976, 0.000600601, 0.00127628, 0.000600601, 0.0001001, 0.00015015, 
                                0.000600601, 7.51E-05, 7.51E-05, 0.0012012, 0.00035035, 0.000600601, 0.00105105, 0.000600601, 0.00202202, 6.01E-05, 0.0002002, 
                                0.00108108, 0.0004004, 0.00015015, 0.0002002, 6.01E-05, 0.000500501, 0.0012012, 0.00045045, 0.000600601, 0.0004004, 0.0002002, 
                                0.000600601, 7.51E-05, 0.00015015, 7.51E-05, 7.51E-05, 7.51E-05, 0.00045045, 0.000160875, 7.51E-05, 0.0004004, 0.000600601, 
                                0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.000600601, 0.0002002, 0.0004004, 0.000600601, 7.51E-05, 7.51E-05, 
                                0.000600601, 0.000600601, 0.00048048, 7.51E-05, 7.51E-05, 7.51E-05, 6.67E-05, 7.51E-05, 7.51E-05))

g_modulos_infomap <- graph_from_data_frame(aristas_infomap, directed=FALSE, vertices=modulos_infomap)
l_infomap <- as.matrix(data.frame(V1=c(88.350440979, 55, 37.0379295349, -32.59, -48.57, -53.19, -32.6633224487, -92.2706298828, 59.1969223022, 72.7037963867, 
                                       90.9873275757, 29.1660251617, -17.8630561829, 42.4076805115, -47.6372146606, -67.4067687988, 16.9434547424, 25.5657920837, 
                                       44.8087158203, 92.84, 15.9724206924, 16.1179275513, 50.71, 38.6608390808, -68.97, 2.304826498, 71.6457214355, -7.1567931175, 
                                       -30.7502059937, 17.3787746429, 97.0877914429, 6.8339300156, 37.2960739136, -22.276058197, -74.2743301392, -82.7057647705), 
                                  V2=c(13.9084339142, 61.5, 74.5166397095, 38.96, -18.23, 0, 14.2420682907, -26.0223655701, 35.1093215942, 48.4196090698, 
                                       48.6638793945, 27.235490799, 20.9897441864, 0.1934056878, 67.1593322754, 33.283367157, -12.0031938553, -6.89290905, 
                                       -8.9333677292, 36.01, 15.5749177933, 101.9463348389, 80.61, 29.2469329834, 25.48, 68.271697998, -58.0035400391, 78.2148971558, 
                                       -53.3914642334, 82.5672225952, 65.8700485229, -6.0750436783, 77.4263153076, 100.4801330566, -33.0328712463, 66.0332489014)))

plot(g_modulos_infomap, layout = l_infomap, vertex.label.color = "black", vertex.frame.color = "#888888", edge.color = "#0066FF99", 
     vertex.size = 0.3+V(g_modulos_infomap)$flujo_dentro*200, edge.width = 0+E(g_modulos_infomap)$flujo*2000, 
     vertex.label = V(g_modulos_infomap)$name, vertex.label.cex = 0.5+V(g_modulos_infomap)$flujo_dentro*10 )
