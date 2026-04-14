t-SNE para reducción de dimensiones
================
Silvia M. Pierotti
14/04/2026

**Limpieza del set de datos**

``` r
library(readr)
library(patchwork)
#UNificacion de tablas
dat<- read.csv("C:/Users/silvi/Desktop/Master bioinformatica/Febrero/Algoritmos e inteligencia artificial/Actividad 1/data.csv", sep = ",")
niveles<-read.csv("C:/Users/silvi/Desktop/Master bioinformatica/Febrero/Algoritmos e inteligencia artificial/Actividad 1/labels.csv",sep=",")#contienen tipos de cancer (Class)
dat_1<-cbind(niveles,dat)#contiene muestra y tipo de cancer 
datos<-dat_1[,-3]#se elimina la columna repetida para muestra
```

``` r
library(dplyr)
filtrados <- datos %>%
  select(where(~ isTRUE(mean(. > 0) >= 0.5 && (sd(.) / mean(.)) > 0.10))) # se eliminan columnas (genes) con 50% o mas de celdas con ceros y CV menores a 10%. El argumento isTRUE elimina todo lo que es NA para el calculo de CV
datos_limpios<-cbind(datos[,1:2],filtrados)
```

**Metodo t-SNE para la reducción de dimensiones y visualización de
clusters**

``` r
#install.packages("Rtsne")
library(Rtsne)
library(ggplot2)
set.seed(1234)
resultados <- data.frame(Rtsne(scale(sapply(datos_limpios[3:10205], perplexity = 20-30, as.numeric)))$Y)#scale normaliza los datos de expresión génica, perplexity ajusta el  número de vecinos (cada paciente se posiciona considerando ~20–30 pacientes cercanos).

ggplot(resultados, aes(x = X1, y = X2, color = datos_limpios$Class)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "orange", "green", "blue", "purple")) +
  labs(title = "Resultados t-SNE", x = "Dim1", y = "Dim2", color = "Grupo") +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "gray95"),
    plot.title = element_text(size=25,hjust = 0.5),
    legend.title = element_text(size = 20),   
    legend.text = element_text(size = 20)     
  )
```

![](t-SNE_files/figure-gfm/setup_t-SNE-1.png)<!-- --> **Figura 1**:
Resultados del método t-SNE

**Identificacion de los 10 genes con niveles de epresion mas altos en
cada cluster**

``` r
#library(dplyr)
library(tidyr)
#library(ggplot2)

datos <- as.data.frame(scale(sapply(datos_limpios[3:10205], as.numeric)))
datos$Class <- datos_limpios$Class

prom_long <- datos %>%
  pivot_longer(-Class, names_to = "Gene", values_to = "Exp") %>%
  group_by(Class, Gene) %>%
  summarise(media = mean(Exp), .groups = "drop")# promedios de expresion (EXP) de cada gen por class

top_genes <- prom_long %>%
  group_by(Class) %>%
  slice_max(order_by = media, n = 20)# 20 genes mas expresados por class

datos_long <- datos %>%
  pivot_longer(-Class, names_to = "Gene", values_to = "Exp") %>%
  inner_join(top_genes[, c("Class", "Gene")], by = c("Class", "Gene"))# inner_join rescata solo lo 10 genes seleccionados para cada class

ggplot(datos_long, aes(x = Gene, y = Exp, fill = Class)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~Class, scales = "free") +
  coord_flip() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16)
  )
```

![](t-SNE_files/figure-gfm/setup_clusters-1.png)<!-- -->

**Figura 2**: Para cada cluster se identifican los 20 genes con mayor
nivel de expresión.
