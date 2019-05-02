library(ggplot2)
library(dplyr)
library(plyr)
source("02_Codigo/multifunction.R")

BDD <- read.csv("01_Datos/02_Datos Coaliciones de eventos/BDD simple LAOMS.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=FALSE)
BDD$C1x1FechaPeriodico <- as.Date(BDD$C1x1FechaPeriodico, "%d/%m/%Y")
BDD$C1x2FechaEvento <- as.Date(BDD$C1x2FechaEvento, "%d/%m/%Y")
actores <- read.csv("01_Datos/02_Datos Coaliciones de eventos/actores.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=FALSE)

########################################################################################################################
# 3.2.2.Demandas
########################################################################################################################

########################################################################################################################
## Figura 3.17: Eventos de protesta — Tendencias asociadas a demandas.
# 
# Se cuentan EPs únicos -si al menos una de sus demandas pudo ser relacionada a uno de estos temas-
# Como en el caso de las coincidencias por diario, se buscaron coincidencias exactas en la demanda textual.
# Si la demanda textual no indicaba alguno de términos, se añadió la presición de que dicha demanda era en alguno
#  de estos contextos (por ejemplo: una demanda que hacía referencia a algun tema específico descrito en la cronología 
# -Nochixtlán, Paro de labores de la CNTE, Leyes secundarias de la reforma o control del IEEPO-)
########################################################################################################################

ts <- ddply(BDD, .(C1x2FechaEvento), summarize, Ayotzinapa=sum(Ayotzinapa), Reforma.educativa=sum(Reforma.educativa))
ts$Total <- as.vector( table(BDD$C1x2FechaEvento) )
hh <- data.frame(fecha=seq(as.Date("2012-10-16"), as.Date("2016-12-31"), by="days"))
ts <- merge(ts,hh,by.x='C1x2FechaEvento',by.y='fecha',all.x=T,all.y=T)
ts$Total[is.na(ts$Total)] <- 0
ts$Ayotzinapa[is.na(ts$Ayotzinapa)] <- 0
ts$Reforma.educativa[is.na(ts$Reforma.educativa)] <- 0
ts$semana <- format(ts$C1x2FechaEvento, format="%Y-%U")
ts$mes <- format(ts$C1x2FechaEvento, format="%Y-%m")

# Serie de tiempo diaria
p1 <- ggplot(ts, aes(C1x2FechaEvento, Total, group=1, colour=Demandas.presentes.en.EPs )) + 
    xlab("Fecha del evento") + 
    ylab("EPs conjuntos diarios") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_line(aes(y = Total, colour = "Total EPs")) +
    geom_line(aes(y = Ayotzinapa, colour = '"Ayotzinapa"'), alpha=0.8) +
    geom_line(aes(y = Reforma.educativa, colour = '"Reforma educativa"'), alpha=0.8) +
    scale_color_manual(values=c("red", "blue", "#999999")) +
    scale_y_continuous(breaks = seq(0, 13, 2)) +
    scale_x_date(date_breaks = "50 day", date_labels = "%d %b %y", limits = c(as.Date("2012-10-16"), as.Date("2016-12-31")))

# Serie de tiempo semanal
ts_semanal <- ddply(ts, .(semana), summarize, Total=sum(Total), Ayotzinapa=sum(Ayotzinapa), Reforma.educativa=sum(Reforma.educativa))
p2 <- ggplot(ts_semanal, aes(semana, Total, group=1, colour=Demandas.presentes.en.EPs )) + 
    xlab("Año - Semana (evento)") + 
    ylab("EPs conjuntos semanales") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_line(aes(y = Total, colour = "Total EPs"), alpha=0.6) +
    geom_line(aes(y = Ayotzinapa, colour = '"Ayotzinapa"'), alpha=0.8) +
    geom_line(aes(y = Reforma.educativa, colour = '"Reforma educativa"'), alpha=0.8) +
    scale_color_manual(values=c("red", "blue", "#222222")) +
    scale_y_continuous(breaks = seq(0, 20, 2)) +
    scale_x_discrete(breaks= levels(as.factor(ts_semanal$semana))[c(T, rep(F, 6))])

# Serie de tiempo mensual
ts_mensual <- ddply(ts, .(mes), summarize, Total=sum(Total), Ayotzinapa=sum(Ayotzinapa), Reforma.educativa=sum(Reforma.educativa))
p3 <- ggplot(ts_mensual, aes(mes, group=1, colour=Demandas.presentes.en.EPs )) + 
    xlab("Año - Mes (evento)") + 
    ylab("EPs conjuntos mensuales") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_line(aes(y = Ayotzinapa, colour = '"Ayotzinapa"'), alpha=0.8) +
    geom_line(aes(y = Reforma.educativa, colour = '"Reforma educativa"'), alpha=0.8) +
    geom_line(aes(y = Total, colour = "Total EPs"), alpha=0.6) +
    scale_y_continuous(breaks = seq(0, 40, 4)) +
    scale_x_discrete(breaks= levels(as.factor(ts_mensual$mes))[c(T, rep(F, 0))]) + 
    scale_color_manual(values=c("red", "blue", "#222222"))

multiplot(p1, p2, p3,cols=1)




###########################################################################################################################
# Similaridades de Jaccard entre atributos registrados en protestas (clasificación de demandas, repertorios y respuestas
###########################################################################################################################

library(lattice)
library(reshape)
library(scales)
library(corrplot)

### Dada la trivialidad de su cálculo, las matrices de superposición fueron obtenidas por sendas tablas dinámicas en Excel
dem <- read.csv("01_Datos/02_Datos Coaliciones de eventos/superposicion_demandas.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=TRUE)
rep <- read.csv("01_Datos/02_Datos Coaliciones de eventos/superposicion_repertorios.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=TRUE)
res <- read.csv("01_Datos/02_Datos Coaliciones de eventos/superposicion_respuestas.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=TRUE)

####################################################################################
### Figura 3.20: Superposición entre demandas (protestas con más de una categoría).
####################################################################################

dem[is.na(dem)] <- 0
deman <- t(dem[,2:length(colnames(dem))])
sim_demandas <- matrix(0, nrow=nrow(deman), ncol=nrow(deman))
rownames(sim_demandas) <- rownames(deman)
colnames(sim_demandas) <- rownames(deman)

pares <- t(combn(1:nrow(deman), 2))
for (i in 1:nrow(pares)){
    num <- sum(sapply(1:ncol(deman), function(x)(min(deman[pares[i,1],x],deman[pares[i,2],x]))))
    den <- sum(sapply(1:ncol(deman), function(x)(max(deman[pares[i,1],x],deman[pares[i,2],x]))))
    sim_demandas[pares[i,1],pares[i,2]] <- num/den
    sim_demandas[pares[i,2],pares[i,1]] <- num/den  
}
sim_demandas[which(is.na(sim_demandas))] <- 0
diag(sim_demandas) <- 1

tipos <- c("Abusos de policías y militares", "Condiciones de trabajo", "Contra funcionarios electos", "Contratos colectivos", "Corrupción", "Defensa nacional", "Desapariciones forzadas", "Despidos", "Elecciones", "Estabilidad laboral", "Fuentes de empleo", "Homicidios", "Impartición de justicia", "Impuestos", "Indígenas", "Infomación verídica", "Jubilaciones", "Medio ambiente", "Mesa de negociación", "Migrantes", "Mujeres", "Permisos", "Presos políticos", "Prestaciones", "Procuración de justicia", "Recursos (consumidores)", "Recursos (educación)", "Recursos (rurales)", "Recursos (salud)", "Recursos (seguridad)", "Recursos (urbanos)", "Reformas (constitucionales)", "Reformas (no const.)", "Regulación del mercado", "Rendición de cuentas (derechos)", "Rendición de cuentas (informativas)", "Represión de protestas", "Salarios", "Sindicales", "Solidaridad con extranjeros", "Trabajo de autoridades")
rownames(sim_demandas) <- colnames(sim_demandas) <- c(tipos)

sim <- melt(sim_demandas) ## Ejecutar la función al final del código para obtener la imagen de la matriz de Jaccard

####################################################################################
### Figura 3.22: Superposición entre repertorios de protesta.
####################################################################################
rep[is.na(rep)] <- 0
tipos <- colnames(rep)[2:length(colnames(rep))]
totales <- colSums(rep[2:length(colnames(rep))])

sim_repertorios <- matrix(0L, nrow=length(tipos), ncol=length(tipos))
rownames(sim_repertorios) <- colnames(sim_repertorios) <- c(tipos)

for(tipo1 in tipos){
    for(tipo2 in tipos){
        interseccion <- sum(rep[tipo1]*rep[tipo2])
        sim_repertorios[tipo1,tipo2] <- interseccion / (totales[tipo1] + totales[tipo2] - interseccion )
    }
}
tipos <- c("Artístico / simbólico", "Dramático (comedia)", "Dramático (tragedia)", "Bloqueo", "Boicot", "Brigadeo", "Cadena Humana", "Campamento", "Caravana", "Ciberactivismo", "Huelga de hambre", "Huelga", "Marcha", "Marcha conmemorativa", "Mitin", "Paro de labores", "Plantón", "Indeterminado", "Retención pacífica de vehículos", "Rodada", "Toma de casetas", "Toma de instalaciones", "Agresiones contra uniformados", "Ataques a civiles", "Ataques a uniformados", "Destrucción de bienes de los partidos políticos", "Destrucción de bienes privados", "Destrucción de bienes públicos", "Destrucción de infraestructura", "Exhibición de armas blancas", "Exhibición de armas de fuego", "Petardos y molotov contra civiles", "Petardos y molotov contra uniformados", "Proyectiles contra civiles", "Proyectiles contra uniformados", "Retención de civiles", "Retención de funcionarios", "Retención de uniformados", "Retención violenta de vehículos", "Toma u ocupación violenta")
rownames(sim_repertorios) <- colnames(sim_repertorios) <- c(tipos)
sim <- melt(sim_repertorios)
sim$X1 <- factor(as.character(sim$X1), levels=tipos)
sim$X2 <- factor(as.character(sim$X2), levels=tipos)
sim <-sim[with(sim,order(X2, X1)),] ## Ejecutar la función al final del código para obtener la imagen de la matriz de Jaccard


####################################################################################
### Figura 3.24: Superposición entre respuestas a la protesta y uso de violencia.
####################################################################################
res[is.na(res)] <- 0
tipos <- colnames(res)[2:length(colnames(res))]
totales <- colSums(res[2:length(colnames(res))])

sim_respuestas <- matrix(0L, nrow=length(tipos), ncol=length(tipos))
rownames(sim_respuestas) <- colnames(sim_respuestas) <- c(tipos)

for(tipo1 in tipos){
    for(tipo2 in tipos){
        interseccion <- sum(res[tipo1]*res[tipo2])
        sim_respuestas[tipo1,tipo2] <- interseccion / (totales[tipo1] + totales[tipo2] - interseccion )
    }
}
tipos <- c("¿Se abrió interlocución?", "¿Se atendieron las demandas por completo?", "¿Se atendieron las demandas parcialmente?", "¿Intervino la fuerza pública?", "¿Hubo violencia de no uniformados?", "Bloqueos", "Encapsulamientos", "Desalojos", "Agresiones por parte de uniformados", "Gas lacrimógeno y/o pimienta", "Balas de goma", "Disparos", "Toletazos", "Otros", "¿Hubo Manifestantes detenidos?", "¿Hubo Transeúntes detenidos?", "¿Hubo Uniformados retenidos contra su voluntad?", "¿Hubo Transeúntes retenidos contra su voluntad?", "¿Hubo Manifestantes heridos?", "¿Hubo Uniformados heridos?", "¿Hubo Transeúntes heridos?", "¿Hubo Manifestantes muertos?", "¿Hubo destrucción de bienes publicos?", "¿Hubo destrucción de bienes privados?", "¿Repertorios violentos? (dicotomizados)", "¿Repertorios pacíficos? (dicotomizados)")
rownames(sim_respuestas) <- colnames(sim_respuestas) <- c(tipos)
sim <- melt(sim_respuestas)
sim$X1 <- factor(as.character(sim$X1), levels=tipos)
sim$X2 <- factor(as.character(sim$X2), levels=tipos)
sim <-sim[with(sim,order(X2, X1)),] ## Ejecutar la función al final del código para obtener la imagen de la matriz de Jaccard


#######################################################################################################################
##### CONSTRUYE Y EXPORTA LA MATRIZ DE SIMILARIDAD EN UN ARCHIVO PNG
#####  Se incrementa la resolución de las imágenes, ya que por default RStudio no permite tal manipulación. 
#######################################################################################################################
png(filename="MatrizDeSimilaridad.png", width=1600, height=1200, units="px", res=110)

ggplot(sim, aes(X1, X2, fill = value)) + 
    geom_tile(color="#888888") + 
    xlab("") + 
    ylab("") + 
    theme_minimal() +
    scale_fill_gradientn(colours=c("white", "cyan", "#007FFF", "blue","#00007F"),
                         values=rescale(c(0,0.25,0.5,0.75,1)),
                         limits = c(0,1),
                         name = "",
                         breaks = seq(0,1, 0.1),
                         guide = guide_colorbar(nbin = 40,
                                                barwidth = 1,
                                                title.position = "bottom",
                                                title.hjust = 0.5, 
                                                raster = FALSE,
                                                ticks = FALSE)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    geom_text(aes(x=X1,y=X2,label=round(value,2)), size=2.5, color="#ff33cc")

dev.off()
