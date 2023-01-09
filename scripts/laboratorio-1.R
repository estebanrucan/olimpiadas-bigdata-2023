library(tidyverse)
library(pspline)

# Pregunta 1

datos <- read_csv(file.choose())

View(datos)
     # Filas  Columnas
datos[, c("categoria", "grasa_saturada")]

datos %>% # pipe: al resultado de la izquierda, le aplico lo de la derecha 
  group_by(categoria) %>% 
  summarise(promedio = mean(grasa_saturada))

names(datos) # ver nombres de las variables

datos_numericos <- datos[, c("tamano_de_la_porcion", "calorias", "proteina", "sodio", "grasa_saturada", "azucares")]

pairs(datos_numericos)

# Pregunta 2

# install.packages("GGally")
library(GGally)

ggpairs(datos_numericos)
# Calorías tiene la asociación lineal más fuerte respecto a grasa saturada

# Pregunta 3           Y         ~     X
mejor_recta <- lm(grasa_saturada ~ calorias, data = datos_numericos)
mejor_recta

plot(grasa_saturada ~ calorias, data = datos_numericos, col = "red", pch = 17)
abline(mejor_recta)

# Pregunta 4

plot(grasa_saturada ~ sodio, data = datos_numericos, col = "red", pch = 17)

# Debemos encontrar el parámetro suavizamiento más indicado

ancho_de_banda <- seq(50, 200, 10)

120 * 16 # modelos
R <- matrix(nrow = 16, ncol = 120)

# Métrica: Residuo = valor_real - valor_estimado
# Llenaremos R con todos los residuos resultantes

for (i in 1:16) {
  
  for (j in 1:120) {
    
    y_estimado <- ksmooth( # estimador de Nadaraya-Watson
      
      x         = datos_numericos$sodio[-j],
      y         = datos_numericos$grasa_saturada[-j],
      kernel    = "normal",
      bandwidth = ancho_de_banda[i],
      x.points  = datos_numericos$sodio[j] 
      
    )$y
    
    R[i, j] <- datos_numericos$grasa_saturada[j] - y_estimado
    
  }
  
}

CV <- apply(R, 1, crossprod) # 1: Filas, crossprod: eleva al cuadrado y dsp suma.

ancho_de_banda_seleccionado <- ancho_de_banda[which.min(CV)]
ancho_de_banda_seleccionado

ajuste_NW <- ksmooth( # estimador de Nadaraya-Watson
  
  x         = datos_numericos$sodio,
  y         = datos_numericos$grasa_saturada,
  kernel    = "normal",
  bandwidth = ancho_de_banda_seleccionado,
  x.points  = datos_numericos$sodio 
  
)

# Splines cúbicos
library(pspline)

alpha <- 10 ^ seq(-7, 0, len = 100)

spline_cubico <- function(a) {
  
  spline <- sm.spline(
    x      = datos_numericos$sodio,
    y      = datos_numericos$grasa_saturada,
    norder = 2, # Spline cúbico: 2 * n - 1
    spar   = a
  )
  
  cv <- spline$cv
  
  return(cv)
  
}


cv <- unlist(lapply(alpha, spline_cubico))
mejor_alpha <- alpha[which.min(cv)]

ajuste_SC <- sm.spline(
  x      = datos_numericos$sodio,
  y      = datos_numericos$grasa_saturada,
  norder = 2, # Spline cúbico: 2 * n - 1
  spar   = mejor_alpha
)

# Gráfico final
plot(grasa_saturada ~ sodio, data = datos_numericos, pch = 16)
lines(ajuste_NW, col = "blue")
lines(ajuste_SC, col = "red")
abline(lm(grasa_saturada ~ sodio, data = datos_numericos))

# Nos quedamos con estimador NW:
plot(grasa_saturada ~ sodio, data = datos_numericos, pch = 16)
lines(ajuste_NW, col = "blue")

# El nivel de sodio debería estar cerca del valor 3000 para disminuir la grasa_saturada.