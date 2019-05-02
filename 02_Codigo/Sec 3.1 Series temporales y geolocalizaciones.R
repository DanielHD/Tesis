library(ggplot2)
library(plyr)
library(dplyr)
source("02_Codigo/multifunction.R")


#########################################################################
### 3.1.1. Cobertura temporal y sobrerrepresentación de eventos
#########################################################################

fuentes <- read.csv("01_Datos/02_Datos Coaliciones de eventos/referencias comunes a eventos.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=FALSE)
fuentes$C1x1FechaPeriodico <- as.Date(fuentes$C1x1FechaPeriodico, "%d/%m/%Y")


##################################################################################################################################################
# Figura 3.1: Comparación entre Referencias periodísticas y Eventos de protesta.

ts_fuentes <-data.frame( fecha = as.Date(rownames( table(fuentes$C1x1FechaPeriodico) ),"%Y-%m-%d" ), protestas = as.vector( table(fuentes$C1x1FechaPeriodico) ))
hh_fuentes <- data.frame(fecha=seq(as.Date("2012-10-16"), as.Date("2016-12-31"), by="days"))
ts_fuentes <- merge(ts_fuentes,hh_fuentes,by.x='fecha',by.y='fecha',all.x=T,all.y=T)
ts_fuentes$protestas[is.na(ts_fuentes$protestas)] <- 0
ts_fuentes$semana <- format(ts_fuentes$fecha, format="%Y-%U")
ts_fuentes$mes <- format(ts_fuentes$fecha, format="%Y-%m")

fuentes_diarias <- ggplot(ts_fuentes, aes(fecha, protestas, group=1)) + 
    xlab("Fecha del periódico") + 
    ylab("Referencias diarias recuperadas") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_line(stat="identity") + 
    scale_y_continuous(breaks = seq(0, 13, 2)) +
    scale_x_date(date_breaks = "50 day", date_labels = "%d %b %y", limits = c(as.Date("2012-10-16"), as.Date("2016-12-31")))

# Serie de tiempo semanal
ts_semanal_fuentes <- ddply(ts_fuentes, .(semana), summarize, protestas=sum(protestas))
fuentes_semanales <- ggplot(ts_semanal_fuentes, aes(semana, protestas, group=1)) + 
    xlab("Año - Semana (periódico)") + 
    ylab("Referencias semanales recuperadas") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_line(stat="identity") + 
    scale_y_continuous(breaks = seq(0, 20, 2)) +
    scale_x_discrete(breaks= levels(as.factor(ts_semanal_fuentes$semana))[c(T, rep(F, 6))]) 

# Serie de tiempo mensual
ts_mensual_fuentes <- ddply(ts_fuentes, .(mes), summarize, protestas=sum(protestas))
fuentes_mensuales <- ggplot(ts_mensual_fuentes, aes(mes, protestas, group=1)) + 
    xlab("Año - Mes (periódico)") + 
    ylab("Referencias mensuales recuperadas") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_point() + 
    geom_line(stat="identity") +
    scale_y_continuous(breaks = seq(0, 40, 5)) +
    scale_x_discrete(breaks= levels(as.factor(ts_mensual_fuentes$mes))[c(T, rep(F, 1))]) 

### SERIE DE TIEMPO DE EPS POR FECHA DE NOTICIAS
BDD <- read.csv("01_Datos/02_Datos Coaliciones de eventos/BDD simple LAOMS.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=FALSE)
BDD$C1x1FechaPeriodico <- as.Date(BDD$C1x1FechaPeriodico, "%d/%m/%Y")

# Serie de tiempo diaria
ts <-data.frame( fecha = as.Date(rownames( table(BDD$C1x1FechaPeriodico) ),"%Y-%m-%d" ), protestas = as.vector( table(BDD$C1x1FechaPeriodico) ))
hh <- data.frame(fecha=seq(as.Date("2012-10-16"), as.Date("2016-12-31"), by="days"))
ts <- merge(ts,hh,by.x='fecha',by.y='fecha',all.x=T,all.y=T)
ts$protestas[is.na(ts$protestas)] <- 0
ts$semana <- format(ts$fecha, format="%Y-%U")
ts$mes <- format(ts$fecha, format="%Y-%m")

p1 <- ggplot(ts, aes(fecha, protestas, group=1)) + 
    xlab("Fecha del periódico") + 
    ylab("EPs conjuntos diarios") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_line(stat="identity") + 
    scale_y_continuous(breaks = seq(0, 13, 2)) +
    scale_x_date(date_breaks = "50 day", date_labels = "%d %b %y", limits = c(as.Date("2012-10-16"), as.Date("2016-12-31")))

# Serie de tiempo semanal
ts_semanal <- ddply(ts, .(semana), summarize, protestas=sum(protestas))
p2 <- ggplot(ts_semanal, aes(semana, protestas, group=1)) + 
    xlab("Año - Semana (periódico)") + 
    ylab("EPs conjuntos semanales") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_line(stat="identity") + 
    scale_y_continuous(breaks = seq(0, 20, 2)) +
    scale_x_discrete(breaks= levels(as.factor(ts_semanal$semana))[c(T, rep(F, 6))]) 

# Serie de tiempo mensual
ts_mensual <- ddply(ts, .(mes), summarize, protestas=sum(protestas))
p3 <- ggplot(ts_mensual, aes(mes, protestas, group=1)) + 
    xlab("Año - Mes (periódico)") + 
    ylab("EPs conjuntos mensuales") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_point() + 
    geom_line(stat="identity") +
    scale_y_continuous(breaks = seq(0, 40, 5)) +
    scale_x_discrete(breaks= levels(as.factor(ts_mensual$mes))[c(T, rep(F, 1))]) 


multiplot(fuentes_diarias, fuentes_semanales, fuentes_mensuales, p1, p2, p3,cols=2)

# Correlaciones
cor(ts_fuentes$protestas, ts$protestas) # 0.7765757
cor(ts_semanal_fuentes$protestas, ts_semanal$protestas) # 0.8198635
cor(ts_mensual_fuentes$protestas, ts_mensual$protestas) # 0.9012282



################################################################
####  3.1.2. Registro de eventos de protesta de larga duración
################################################################

eps_largos <- unique(BDD[!is.na(BDD$id_EPLargo),c("id_EPLargo", "Duración..en.días.", "dias_cobertura")])

##################################################################################################################################################
#  Figura 3.2: Duración aproximada y cobertura diaria de eventos de protesta de larga duración.
p1 <- ggplot(data=eps_largos, aes(x=id_EPLargo, y=Duración..en.días. )) + 
    geom_bar(stat="identity") + 
    scale_y_continuous(breaks = seq(0, 950, 100)) + 
    ylab("Duración en días") + xlab("EPs de larga duración") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
p2 <- ggplot(eps_largos, aes(x=id_EPLargo, y=dias_cobertura) ) + 
    geom_bar(stat="identity") + 
    scale_y_continuous(breaks = seq(0, 70, 10)) + 
    ylab("Cobertura en días") + xlab("EPs de larga duración") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

multiplot(p1, p2, cols=1)




####################################################################################################
####   3.1.3.Cobertura temática
# Se probaron varios términos para otros actores cuyos términos eran relevantes (por ejemplo: 
# Peña nieto, Campesinos, CNTE), pero se eligió sólo "Ayotzinapa" y "Reforma Educativa" al 
# ser términos más contextiales en periodos prolongados.

# Desde un nivel diario y semanal se observa la tendencia, pero se considera un nivel mensual 
# debido a que presenta variaciones más estables y fácilmente identificables
####################################################################################################

noticias <- read.csv("01_Datos/02_Datos Coaliciones de eventos/noticias periodísticas sobre eventos.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=FALSE)
noticias$C1x1FechaPeriodico <- as.Date(noticias$C1x1FechaPeriodico, "%d/%m/%Y")


####################################################################################################
# Figura 3.3: Ocurrencia de términos relevantes en noticias muestreadas.

# Serie de tiempo diaria
ts_noticias <- ddply(noticias, .(C1x1FechaPeriodico), summarize, Ayotzinapa=sum(Referencias.Ayotzinapa), Reforma.educativa=sum(Referencias.Reforma.educativa), 
                     Peña.Nieto=sum(Referencias.Peña.Nieto), Campesinos=sum(Referencias.Campesinos), CNTE=sum(Referencias.CNTE), 
                     Ayotzinapa.rel=sum(Referencias.relevantes.Ayotzinapa), Reforma.educativa.rel=sum(Referencias.relevantes.Reforma.educativa))
ts_noticias$Total <-as.vector( table(noticias$C1x1FechaPeriodico) )

hh_noticias <- data.frame(fecha=seq(as.Date("2012-10-16"), as.Date("2016-12-31"), by="days"))
ts_noticias <- merge(ts_noticias,hh_noticias,by.x='C1x1FechaPeriodico',by.y='fecha',all.x=T,all.y=T)
ts_noticias$Ayotzinapa[is.na(ts_noticias$Ayotzinapa)] <- 0
ts_noticias$Ayotzinapa.rel[is.na(ts_noticias$Ayotzinapa.rel)] <- 0
ts_noticias$Reforma.educativa[is.na(ts_noticias$Reforma.educativa)] <- 0
ts_noticias$Reforma.educativa.rel[is.na(ts_noticias$Reforma.educativa.rel)] <- 0
ts_noticias$Peña.Nieto[is.na(ts_noticias$Peña.Nieto)] <- 0
ts_noticias$Campesinos[is.na(ts_noticias$Campesinos)] <- 0
ts_noticias$CNTE[is.na(ts_noticias$CNTE)] <- 0
ts_noticias$Total[is.na(ts_noticias$Total)] <- 0
ts_noticias$semana <- format(ts_fuentes$fecha, format="%Y-%U")
ts_noticias$mes <- format(ts_fuentes$fecha, format="%Y-%m")

p1 <- ggplot(ts_noticias, aes(C1x1FechaPeriodico, Ayotzinapa, group=1, colour=Tendencia)) + 
    xlab("Fecha del periódico") + 
    ylab("Ocurrencia de términos") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_line(aes(y = Total, colour = "Total")) +
    geom_line(aes(y = Ayotzinapa, colour = '"Ayotzinapa"')) +
    geom_line(aes(y = Reforma.educativa, colour = '"Reforma educativa"')) +
    scale_y_continuous(breaks = seq(0, 13, 2)) +
    scale_x_date(date_breaks = "50 day", date_labels = "%d %b %y", limits = c(as.Date("2012-10-16"), as.Date("2016-12-31"))) + 
    scale_color_manual(values=c("red", "blue", "#999999"))

# Serie de tiempo semanal
ts_noticias_semanal <- ddply(ts_noticias, .(semana), summarize, Ayotzinapa=sum(Ayotzinapa), Reforma.educativa=sum(Reforma.educativa), 
                    Peña.Nieto=sum(Peña.Nieto), Campesinos=sum(Campesinos), CNTE=sum(CNTE), Total=sum(Total), Ayotzinapa.rel=sum(Ayotzinapa.rel), 
                    Reforma.educativa.rel=sum(Reforma.educativa.rel))
p2 <- ggplot(ts_noticias_semanal, aes(semana, group=1, colour=Tendencia )) + 
    xlab("Año - Semana (periódico)") + 
    ylab("Ocurrencia de términos por noticia") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_line(aes(y = Total, colour = "Total")) +
    geom_line(aes(y = Ayotzinapa, colour = '"Ayotzinapa"')) +
    geom_area(aes(y = Ayotzinapa.rel), colour=NA, fill="red", alpha=0.2) +
    geom_line(aes(y = Reforma.educativa, colour = '"Reforma educativa"')) +
    geom_area(aes(y = Reforma.educativa.rel), colour=NA, fill="blue", alpha=0.2) +
    scale_y_continuous(breaks = seq(0, 40, 4)) +
    scale_x_discrete(breaks= levels(as.factor(ts_semanal$semana))[c(T, rep(F, 6))]) + 
    scale_color_manual(values=c("red", "blue", "#999999"))

# Serie de tiempo mensual
ts_noticias_mensual <- ddply(ts_noticias, .(mes), summarize, Ayotzinapa=sum(Ayotzinapa), Reforma.educativa=sum(Reforma.educativa), 
                    Peña.Nieto=sum(Peña.Nieto), Campesinos=sum(Campesinos), CNTE=sum(CNTE), Total=sum(Total), Ayotzinapa.rel=sum(Ayotzinapa.rel), 
                    Reforma.educativa.rel=sum(Reforma.educativa.rel))
p3 <- ggplot(ts_noticias_mensual, aes(mes, group=1, colour=Tendencia )) + 
    xlab("Año - Mes (periódico)") + 
    ylab("Ocurrencia de términos por noticia") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_line(aes(y = Ayotzinapa, colour = '"Ayotzinapa"')) +
    geom_area(aes(y = Ayotzinapa.rel), colour=NA, fill="red", alpha=0.2) +
    geom_line(aes(y = Reforma.educativa, colour = '"Reforma educativa"')) +
    geom_area(aes(y = Reforma.educativa.rel), colour=NA, fill="blue", alpha=0.2) +
    geom_line(aes(y = Total, colour = "Total")) +
    scale_y_continuous(breaks = seq(0, 40, 4)) +
    scale_x_discrete(breaks= levels(as.factor(ts_mensual$mes))[c(T, rep(F, 0))]) + 
    scale_color_manual(values=c("red", "blue", "#999999"))


# comparación de frecuencias diaria, semanal y mensual:
multiplot(p1, p2, p3, cols=1)




####################################################################################################
####   3.1.4.Representación geográfica de eventos
####################################################################################################
library(mxmaps)
library(ggrepel)



estados <- read.csv("01_Datos/02_Datos Coaliciones de eventos/mapas estados.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=FALSE)
df_mxstate$value <- estados$protestas_conjuntas
gg <- mxstate_choropleth(df_mxstate, 
                         num_colors = 1,
                         title = "")

########################################################################################### 
# Figura 3.4: Menciones únicas de estados en noticias consideradas.

df_mxstate$noticias <- estados$noticias
df_mxstate$lat = estados$lat
df_mxstate$lon = estados$lon
df_mxstate$group <- df_mxstate$state_abbr
gg + 
    geom_text_repel(data = df_mxstate, aes(lon, lat, label = noticias), size = 3,
                    box.padding = unit(0.1, 'lines'), force = 0.5)



########################################################################################### 
# Figura 3.5: Eventos de protesta — agregación estatal.

df_mxstate$protestas_conjuntas <- estados$protestas_conjuntas
df_mxstate$lat = estados$lat
df_mxstate$lon = estados$lon
df_mxstate$group <- df_mxstate$state_abbr
gg +
    geom_text_repel(data = df_mxstate, aes(lon, lat, label = protestas_conjuntas), size = 3,
                    box.padding = unit(0.1, 'lines'), force = 0.5)


########################################################################################### 
# Figura 3.6: Eventos de protesta — agregación municipal.

mun <- read.csv("01_Datos/02_Datos Coaliciones de eventos/mapas municipios.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=FALSE)
df_mxmunicipio$value <-  mun$protestas_conjuntas
mxmunicipio_choropleth(df_mxmunicipio, num_colors = 9,
                       zoom = subset(mun, protestas %in% c("S"))$region,
                       title = "") 


########################################################################################### 
# Figura 3.8: Coaliciones de eventos — Geolocalización aproximada.
library(rworldmap)
library(ggmap)
library(ggmapstyles)

geolocalizaciones <- read.csv("01_Datos/02_Datos Coaliciones de eventos/mapas geolocalizaciones.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=FALSE)
geolocalizaciones_lugares <- read.csv("01_Datos/02_Datos Coaliciones de eventos/mapas geolocalizaciones_lugares.csv", header = TRUE, sep = ',', na.strings = "", stringsAsFactors=FALSE)


# Requiere una clave para "Maps Static API", de Google, se puede crear desde: https://console.developers.google.com/apis/api/static_maps_backend?project=_
register_google(key = "...")
mapa <- get_googlemap(center = c(-103.28135795, 23.78539625),
                      zoom = 5,
                      language = "es-ES",
                      sensor = FALSE, 
                      messaging = FALSE,
                      maptype = "roadmap",
                      scale = 2,
                      style="feature:all|element:labels|visibility:off")
plot <- ggmap(mapa) +
    geom_point(aes(x = lon, y = lat, colour=Repertorio_de_referencia_, size=EPs), 
               alpha = 0.5, stroke = 0.5, 
               data = geolocalizaciones) + 
    labs(x = 'Longitud', y = 'Latitud') +
    scale_colour_manual(name="Repertorio de referencia",
                        values = c("Bloqueo"="#ff0000",
                                   "Caravana"="#cc9900", 
                                   "Marcha"="#0000ff",
                                   "Mitin"="#66ffff",
                                   "Otro"="#000000",
                                   "Pacífico indeterminado"="#0066cc", 
                                   "Paro de labores"="#9900cc",
                                   "Plantón"="#ff00ff", 
                                   "Toma de instalaciones"="#33cc33",
                                   "Varios"="#666699"))
print(plot)


########################################################################################### 
# Figura 3.9: Geolocalización de coaliciones de eventos en ciudades con mayor representación mediática.

# (CDMX)
mapa1 <- get_googlemap(center = c(-99.14, 19.38),
                      zoom = 11,
                      language = "es-ES",
                      sensor = FALSE, 
                      messaging = FALSE,
                      maptype = "roadmap",
                      scale = 2, style="feature:all|element:labels|fontSize:6")
local1 <- ggmap(mapa1) +
    geom_point(aes(x = lon, y = lat, size = EPs), 
               data = geolocalizaciones_lugares, colour = "red", alpha= 0.5 ) + 
    labs(x = 'Longitud', y = 'Latitud') +    
    theme(legend.position='none')

# (Tuxtla Gutierrez)
mapa2 <- get_googlemap(center = "Tuxtla Gutiérrez",
                      zoom = 12,
                      language = "es-ES",
                      sensor = FALSE, 
                      messaging = FALSE,
                      maptype = "roadmap",
                      scale = 2, style="feature:all|element:labels|fontSize:6")
local2 <- ggmap(mapa2) +
    geom_point(aes(x = lon, y = lat, size = EPs), 
               data = geolocalizaciones_lugares, colour = "red", alpha= 0.5 ) + 
    labs(x = 'Longitud', y = 'Latitud') +    
    theme(legend.position='none')


# (Chilpancingo de los bravo)
mapa3 <- get_googlemap(center = c(-99.5003046, 17.5495581),
                      zoom = 13,
                      language = "es-ES",
                      sensor = FALSE, 
                      messaging = FALSE,
                      maptype = "roadmap",
                      scale = 2, style="feature:all|element:labels|fontSize:6")
local3 <- ggmap(mapa3) +
    geom_point(aes(x = lon, y = lat, size = EPs), 
               data = geolocalizaciones_lugares, colour = "red", alpha= 0.5 ) + 
    labs(x = 'Longitud', y = 'Latitud') +    
    theme(legend.position='none')

# (Acapulco)
mapa4 <- get_googlemap(center = "Acapulco de Juárez",
                      zoom = 11,
                      language = "es-ES",
                      sensor = FALSE, 
                      messaging = FALSE,
                      maptype = "roadmap",
                      scale = 2)
local4 <- ggmap(mapa4) +
    geom_point(aes(x = lon, y = lat, size = EPs), 
               data = geolocalizaciones_lugares, colour = "red", alpha= 0.5 ) + 
    labs(x = 'Longitud', y = 'Latitud') +    
    theme(legend.position='none')

multiplot(local1, local2, local3, local4, cols=2)

