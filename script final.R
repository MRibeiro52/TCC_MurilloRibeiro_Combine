rm(list = ls()) # Remover tudo da sessão atual do R
setwd(dirname(rstudioapi::getSourceEditorContext()$path))#Seta o diretório automaticamente na pasta onde o arquivo R está salvo

#Bibliotecas necessárias
library(rio)
library(boot)
library(dplyr)
library(caret)
library(tidyr)
library(purrr)
library(ggplot2)
library(nortest)
library(magrittr)
library(factoextra)
library(data.table)
library(caretEnsemble)

#Funcoes criadas -----------------------------------------------------------------------------------------------------------------------
copy_print <- function(x, rowname = F)
{
  if(class(x)[1] == 'table' & length(dim(x)) == 1){ x = as.data.frame(x); rowname = F}
  else if(class(x)[1] == 'table' & length(dim(x)) == 2){x = as.data.frame.matrix(x); rowname = T}
  
  x <- x %>% mutate_if(is.numeric,~ round(. ,digits = 2)) 
  # %>% 
  #   mutate_if(~is.character(.),~str_replace_all(.,"\\,","."))
  
  DT::datatable(x,
                extensions = 'Buttons', 
                rownames = rowname,
                options = list(dom = 'Bftr',
                               buttons = c('copy'),
                               scrollx = '400px',
                               scrolly = '600px',
                               scrollCollapse = T,
                               pagelength = nrow(x),
                               sep = ";")
  )
} # Funcao para gerar output para copiar direto para o excel de forma pratica
my.summary <- function(x) 
{
  list(Media = mean(x, na.rm = T),
       Mediana = median(x, na.rm = T),
       DesvioPadrao = sd(x, na.rm = T),
       Minimo = min(x, na.rm = T),
       Maximo = max(x, na.rm = T),
       Range = max(x, na.rm = T) - min(x, na.rm = T))
} # Funcao para gerar uma descritiva geral personalizada
teste_modelos <- function(modelos)
{
  lista <- list()
  for (i in 1:length(modelos))
  {
    lista[[i]] = summary(modelos[[i]])
  }
  
  #Desvio escalonado do ajuste para validacao do modelo
  
  desvios = vector()
  for(i in 1:length(modelos))
  {
    desvios[i] = lista[[i]]$deviance/lista[[i]]$dispersion 
  }
  
  #q.quadr com n-p graus de liberdade, onde p eh o numero de parametros do modelo
  
  q.quadr = vector()
  
  for(i in 1:length(modelos))
  {
    q.quadr[i] = qchisq(0.95, lista[[i]]$df.residual)
  }
  
  # Teste global do ajuste do modelo (se Verdadeiro modelo adequado)
  
  v_resp = vector()
  for(i in 1:length(modelos))
  {
    v_resp[i] = q.quadr[i] > desvios[i]
  }
  
  melhor_modelo <- which.min(desvios[which(v_resp)])
  
  lista_final <- list(melhor_modelo = lista[[melhor_modelo]],
                      func_desvio = desvios[melhor_modelo],
                      valor_teste = q.quadr[melhor_modelo])
  lista_final
} # Funcao para testar os modelos através da funcao desvio e de adequacao do ajuste do modelo 
analise_res <- function(modelo, vetor_func_lig)
{
  eta.chapeu = predict(modelo)
  eta.til = vetor_func_lig
  pseudo.R2 = (cor(eta.chapeu,eta.til))^2
  
  #Valor ajustado e desvio residual
  fit = fitted(modelo)
  devres = glm.diag(modelo)$rd
  
  #Grafico de Normalidade
  hist <- hist(devres)
  lt <- lillie.test(devres)
  kst <- ks.test(devres,"pnorm",0,1)
  # qqnorm(devres); qqline(devres, col=2)
  
  #Erros nao-correlacionados
  acf <- acf(devres)
  
  par(mfrow = c(1,2))
  #Histograma do desvio residual
  plot(hist, col = "grey", main = "Histograma do desvio residual")
  #QQPlot
  qqnorm(devres); qqline(devres, col=2)
  
  par(mfrow = c(1,2))
  
  #Verificar a funcao de variancia
  plot(fit,devres,main = "Funcao de variancia")
  #Verificar a funcaoo de Ligacao
  plot(fit,modelo$y,main = "Funcao de ligacao")
  
  par(mfrow = c(1,1))
  #Auto-correlacao
  plot(acf)
  
  printar <- list(pseudoR2 = pseudo.R2,
                  verificao_funcao_ligacao = summary(lm(eta.til~eta.chapeu)),
                  teste_normalidade_lt = lt,
                  teste_normalidade_kst = kst,
                  devres = devres)
} # Funcao que realiza a analise de resíduos, gráficos e testes para verificar normalidade, homocedasticidade, funcao de ligacao e correlacao dos resíduos residuais
analise_diag <- function(x,y,modelo, ref = 3)
{
  #x = desvio residual
  #y = base do modelo
  #Modelo = modelo utilizado
  
  n = nrow(y)
  p = (summary(modelo)$df.null - summary(modelo)$df.residual)
  
  
  plot(x)
  abline(h=-ref, col = 2); abline(h=ref , col=2)
  aberrantes <- which(abs(x) > ref)
  
  #Identificacao de pontos de alavanca
  plot(glm.diag(modelo)$h)
  abline(h=2*(p/n), col = 2)
  alavanca <- which(glm.diag(modelo1)$h > 2*(p/n)) 
  
  #Identificacao de pontos de influencia
  plot(glm.diag(modelo)$cook)
  abline(h=qchisq(0.05,p)/p, col = 2)
  influencia <- which(glm.diag(modelo)$cook > qchisq(0.05,p)/p)
  
  ala_inf <- union(alavanca,influencia) #Uniao dos pontos de alavanca e influencia 
  
  obs_retirar <- aberrantes[which(!(aberrantes %in% ala_inf))] #Pontos aberrantes que nao sao pontos de alavanca e nem pontos de influencia
  
  list(n_aberrantes = length(aberrantes),
       n_alavanca = length(alavanca),
       n_influencia = length(influencia),
       qtd_retirar = length(obs_retirar),
       obs_retirar = unname(obs_retirar))
} # Funcao para realizar a analise de diagnostico para verificar e retirar outliers severos de acordo com um valor de referencia (Padrao = 3)

## Tratamento da base -------------------------------------------------------------------------------------------------------------

base_in <- import("base_consolidada.csv")
base_in %<>% filter(!(Pos %in% c("QB","P","K","FB","LS")))

# Tirar quem nao tem altura ou peso ou participou em 2021
# Converter para altura para cm, peso para kg, e resultados em polegada para cm.
base_tratada = base_in %>% filter(ano != 2021) %>% 
  separate(Ht, c('feet', 'inches'), "-", convert = TRUE) %>% 
  mutate(altura = (12*feet + inches)*2.54,
         peso = Wt*0.453592,
         `Broad Jump` = `Broad Jump`*2.54,
         `Vertical` = `Vertical` * 2.54) %>%
  filter(!is.na(altura)) %>% 
  select(-Wt,-feet,-inches)

## Descritiva ------------------------------------------------------------------

# Geral

aux <- base_tratada %>% select(peso,altura,`40yd`,Bench,`Broad Jump`,`3Cone`)
map_dfr(aux,my.summary) %>% t() %>% data.frame()

# Valores medios e desvio padrão do resultado dos testes fisicos por posicao
  
base_tratada %>% group_by(Pos) %>% 
  summarise(med_40 = mean(`40yd`, na.rm = T),
            dp_40 = sd(`40yd`, na.rm = T),
            med_bj = mean(`Broad Jump`, na.rm = T),
            dp_bj = sd(`Broad Jump`, na.rm = T),
            med_bench = mean(Bench, na.rm = T),
            dp_bench = sd(Bench,na.rm = T),
            med_3cone = mean(`3Cone`, na.rm = T),
            dp_3cone = sd(`3Cone`, na.rm = T)) %>% copy_print()

# Boxplot por posição e teste

ggplot(base_tratada) +
  geom_boxplot(aes(x = Pos, y = `40yd`)) +
  labs(x = "Posições", y = "Tempo (s)", 
       title = "Boxplots por posição do resultado do tiro de 40 jardas no combine no período de 2000 a 2022") +
  theme_minimal()

ggplot(base_tratada) +
  geom_boxplot(aes(x = Pos, y = Bench)) +
  labs(x = "Posições", y = "Repetições", 
       title = "Boxplots por posição do resultado do supino no combine no período de 2000 a 2022") +
  theme_minimal()

ggplot(base_tratada) +
  geom_boxplot(aes(x = Pos, y = `Broad Jump`)) +
  labs(x = "Posições", y = "Distância (cm)", 
       title = "Boxplots por posição do resultado do salto em distância no combine no período de 2000 a 2022") +
  theme_minimal()

ggplot(base_tratada) +
  geom_boxplot(aes(x = Pos, y = `3Cone`)) +
  labs(x = "Posições", y = "Tempo (s)", 
       title = "Boxplots por posição do resultado do teste dos 3 cones no combine no período de 2000 a 2022") +
  theme_minimal()

# Volume testes fisicos

base_tratada %>% filter(!is.na(`40yd`)) %>% nrow() #6.123
base_tratada %>% filter(!is.na(`Broad Jump`)) %>% nrow() #5.025
base_tratada %>% filter(!is.na(`3Cone`)) %>% nrow() #4.039
base_tratada %>% filter(!is.na(Bench)) %>% nrow() #4.612

# Gráficos da media do resultado dos testes fisicos por ano e do peso e altura

basegraf <- base_tratada %>% 
  group_by(ano) %>% 
  summarise(m40 = mean(`40yd`, na.rm = T),
            mb = mean(Bench, na.rm = T),
            mbj = mean(`Broad Jump`, na.rm = T),
            m3c = mean(`3Cone`, na.rm = T),
            alt = mean(altura),
            pes = mean(peso)) 

ggplot(basegraf, aes(x = ano, y = m40)) +
  geom_line(size = 1, color = "darkgreen") +
  labs(x = "Ano", y = "Tempo (s)") +
  theme_minimal()

ggplot(basegraf, aes(x = ano, y = mb)) +
  geom_line(size = 1, color = "darkgreen") +
  labs(x = "Ano", y = "Repetições") +
  theme_minimal()

ggplot(basegraf, aes(x = ano, y = mbj)) +
  geom_line(size = 1, color = "darkgreen") +
  labs(x = "Ano", y = "Distância (cm)") +
  theme_minimal()

ggplot(basegraf, aes(x = ano, y = m3c)) +
  geom_line(size = 1, color = "darkgreen") +
  labs(x = "Ano", y = "Tempo (s)") +
  theme_minimal()

ggplot(basegraf, aes(x = ano, y = alt)) +
  geom_line(size = 1, color = "darkgreen") +
  labs(x = "Ano", y = "Altura (cm)") +
  theme_minimal()

ggplot(basegraf, aes(x = ano, y = pes)) +
  geom_line(size = 1, color = "darkgreen") +
  labs(x = "Ano", y = "Peso (kg)") +
  theme_minimal()

## Agrupamento das posicoes ----------------------------------------------------

# Agrupar posicoes conforme desempenho por teste (Grupos definidos por resultados semelhantes)

# Grupos 40 jardas
g1_40 <- c("OL","C","OG","OT","DT")
g2_40 <- c("ILB","LB","OLB","DL","DE","EDGE","TE")
g3_40 <- c("S","RB","WR","CB","DB")

# Grupos BJ
g1_bj <- c("S","LB","WR","CB","DB")
g2_bj <- c("EDGE","RB","OLB","TE","ILB","DE","DL")
g3_bj <- c("OL","C","OG","OT","DT")

# Grupos Bench
g1_bench <- c("DT","C","OL","OG","DL","OT")
g2_bench <- c("DE","EDGE","OLB","ILB","LB","TE")
g3_bench <- c("RB","S","WR","CB","DB")

# Grupos 3cone
g1_3cone <- c("OG","OT","OL","C","DT","DL")
g2_3cone <- c("DE","TE","ILB","EDGE","OLB","LB","RB")
g3_3cone <- c("S","WR","CB","DB")

base_tratada %<>% mutate(GPos40 = case_when(Pos %in% g1_40 ~ "g1",
                                            Pos %in% g2_40 ~ "g2",
                                            Pos %in% g3_40 ~ "g3"),
                         GPosBench = case_when(Pos %in% g1_bench ~ "g4",
                                               Pos %in% g2_bench ~ "g5",
                                               Pos %in% g3_bench ~ "g6"),
                         GPosBJ = case_when(Pos %in% g1_bj ~ "g7",
                                            Pos %in% g2_bj ~ "g8",
                                            Pos %in% g3_bj ~ "g9"),
                         GPos3CONE = case_when(Pos %in% g1_3cone ~ "g10",
                                               Pos %in% g2_3cone ~ "g11",
                                               Pos %in% g3_3cone ~ "g12"))

rm(g1_40,g2_40,g3_40,
   g1_bj,g2_bj,g3_bj,
   g1_bench,g2_bench,g3_bench,
   g1_3cone,g2_3cone,g3_3cone)

## PCA - Altura e peso -------------------------------------------------------------

cor(base_tratada$peso,base_tratada$altura) #73,75% de correlacao

set.seed(52)

basepca <- base_tratada %>% select(peso,altura)
jogadores <- base_tratada %>% select(Player)

# PCA  padronizando

pca_padronizado <- prcomp(basepca, center = T,scale. = T)
summary(pca_padronizado)
# CP1: 86,88% da variabilidade

summary(pca_padronizado)$rotation
# CP1 peso: -0.7071068, altura: -0.7071068

fviz_eig(pca_padronizado, fill = "lightbue")

fviz_pca_var(pca_padronizado,
             col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800"),
             repel = TRUE,     
             legend.title = "Contribution")

# PCA sem padronizar mas retirando o efeito da media

pca <- prcomp(basepca, center = T)
summary(pca)
# CP1: 95,96% da variabilidade

summary(pca)$rotation
# CP1 peso: -0.9705092, altura: -0.2410642

fviz_eig(pca, fill = "lightbue")

fviz_pca_var(pca,
             col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800"),
             repel = TRUE,     
             legend.title = "Contribution")

#Distribuicao do componente 1 por teste fisico

componente1 <- pca$x[,1]
base_tratada %<>% cbind(.,componente1)
base_tratada %<>% rename(cp1 = componente1)

gr1 <- base_tratada %>% 
  ggplot(., aes(x = cp1,y = `40yd`)) +
  geom_point(colour = "darkolivegreen") +
  labs(x = "Componente 1", y = "Tempo (s)") +
  theme_minimal()

gr2 <- base_tratada %>% 
  ggplot(., aes(x = cp1,y = Bench)) +
  geom_point(colour = "royalblue2") +
  labs(x = "Componente 1", y = "Repetições") +
  theme_minimal()

gr3 <- base_tratada %>% 
  ggplot(., aes(x = cp1,y = `Broad Jump`)) +
  geom_point(colour = "salmon4") +
  labs(x = "Componente 1", y = "Distância (cm)") +
  theme_minimal()

gr4 <- base_tratada %>% 
  ggplot(., aes(x = cp1,y = `3Cone`)) +
  geom_point(colour = "darkorange2") +
  labs(x = "Componente 1", y = "Tempo (s)") +
  theme_minimal()

gridExtra::grid.arrange(gr1,gr2,gr3,gr4, ncol = 2)

### Modelagem ------------------------------------------------------------------

## Modelo 40yds (Tiros de 40 jardas - segundos)---------------------------------
set.seed(52)

# Base 40 jardas
base_40 <- base_tratada %>% select(Player,ano,GPos40,cp1,altura,peso,`40yd`) %>% filter(!is.na(`40yd`))

teste <- sample(1:nrow(base_40), size = round(0.2 * nrow(base_40)), replace = F)
base_teste_40 <- base_40[teste,] # 20% da base para pegar os resultados depois.
base_mod_40 <- base_40[-teste,] # 80% da base para desenvolver o modelo

fwrite(base_teste_40, "base_teste_40.csv")

base_mod_40 %>% 
  ggplot(aes(x = `40yd`)) +
  geom_histogram(fill = 'darkolivegreen4', color = 'darkolivegreen', bins = 20) +
  labs(x = "Tempo (s)", y = "Contagem") +
  theme_minimal()

#Gerando o modelo gamma - funcao log, identidade e inversa
modelo1 <- glm(formula = `40yd` ~ cp1 + GPos40, family = Gamma (link = "identity"), data = base_mod_40)
modelo2 <- glm(formula = `40yd` ~ cp1 + GPos40, family = Gamma (link = "inverse"), data = base_mod_40)
modelo3 <- glm(formula = `40yd` ~ cp1 + GPos40, family = Gamma (link = "log"), data = base_mod_40)

#Testando os modelos
teste_modelos(list(modelo1,modelo2,modelo3)) #Retorna o modelo com menor funcao desvio (Modelo gamma funcao identidade)

# Funcao desvio (4877.53) < valor teste entao modelo e adequado.
# Nao ha necessidade de fazer recategorizacao
# Melhor modelo = Modelo 1, gama com funcao de ligacao identidade

## Analise de residuo

analise_mod_40 <- analise_res(modelo1,vetor_func_lig = base_mod_40$`40yd`) 

#Graficos indicam normalidade, homocedasticidade, um bom ajuste da funcao de ligacao e não correlação dos erros

analise_mod_40$pseudoR2 #Pseudo R2 = 81,63
analise_mod_40$verificao_funcao_ligacao #Coeficiente angular significativo, bom ajuste da funcao de ligacao.
analise_mod_40$teste_normalidade_lt #Segundo teste, desvios seguem uma distribuicao normal (Alpha = 0.01)

devres = analise_mod_40$devres

# Analise de diagnostico
# Pontos Aberrantes, de Alavanca e Influentes

# Pontos aberrantes serao considerados valores onde |DevRes| > 3

analise_diag(devres,base_mod_40,modelo1)
#22 pontos aberrantes, 388 pontos de alavanca, 0 pontos de influencia
#18 observacoes a serem retiradas

obs_retirar <- analise_diag(devres,base_mod_40,modelo1)[[5]]

## Remodelagem excluindo aberrantes nao influentes/alavanca (1)-----------------

base_retirada <- base_mod_40[-obs_retirar,]

modelo_1_1 <- glm(formula = `40yd` ~ cp1 + GPos40, family = Gamma (link = "identity"), data = base_retirada)

#Testando os modelos

teste_modelos(list(modelo_1_1))

# Funcao desvio (4871.17) < valor teste entao modelo e adequado.
# Nao ha necessidade de fazer recategorizacao
# Melhor modelo = Modelo 1, gama com funcao de ligacao identidade

# Analise de residuo

analise_mod_40_1 <- analise_res(modelo_1_1,vetor_func_lig = base_retirada$`40yd`)

#Graficos indicam normalidade, homocedasticidade e um bom ajuste da funcao de ligacao

analise_mod_40_1$pseudoR2 #Pseudo R2 = 82,26
analise_mod_40_1$verificao_funcao_ligacao #Coeficiente angular significativo, bom ajuste da funcao de ligacao.
analise_mod_40_1$teste_normalidade_lt #Segundo teste, desvios seguem uma distribuicao normal (Alpha = 0.01)

devres_1 = analise_mod_40_1$devres

# Analise de diagnostico
# Pontos Aberrantes, de Alavanca e Influentes

# Pontos aberrantes serao considerados valores onde |DevRes| > 3

analise_diag(devres_1,base_retirada,modelo_1_1)
#7 pontos aberrantes, 387 pontos de alavanca, 0 pontos de influencia
#6 observacoes a serem retiradas

obs_retirar_1 <- analise_diag(devres_1,base_retirada,modelo_1_1)[[5]]


## Remodelagem excluindo aberrantes nao influentes/alavanca (2)-----------------

base_retirada_2 <- base_retirada[-obs_retirar_1,]

modelo_1_2 <- glm(formula = `40yd` ~ cp1 + GPos40, family = Gamma (link = "identity"), data = base_retirada_2)

#Testando os modelos

teste_modelos(list(modelo_1_2))

# Funcao desvio (4864.27) < valor teste entao modelo e adequado.
# Nao ha necessidade de fazer recategorizacao
# Melhor modelo = Modelo 1, gama com funcao de ligacao identidade

# Analise de residuo

analise_mod_40_2 <- analise_res(modelo_1_2,vetor_func_lig = base_retirada_2$`40yd`) 
#Graficos a direita.

analise_mod_40_2$pseudoR2 #Pseudo R2 = 82,45
analise_mod_40_2$verificao_funcao_ligacao #Coeficiente angular significativo, bom ajuste da funcao de ligacao.
analise_mod_40_2$teste_normalidade_lt #Segundo teste, desvios seguem uma distribuicao normal (Alpha = 0.01)

devres_2 = analise_mod_40_2$devres

# Analise de diagnostico
# Pontos Aberrantes, de Alavanca e Influentes

# Pontos aberrantes serao considerados valores onde |DevRes| > 3

analise_diag(devres_2,base_retirada_2,modelo_1_2) 
#3 pontos aberrantes, 381 pontos de alavanca, 0 pontos de influencia
#1 observacao a ser retirada (nao vale a pena retirar somente uma observacao)

modelo40 <- modelo_1_2
#Modelo final para tiro de 40, modelo gamma com ligacao identidade

#RMSE do modelo de treino
rmse40_tr <- sqrt(mean((modelo40$y - fitted(modelo40))^2)) #0,1297

saveRDS(modelo40,"modelo_40yds.rds")

# Analise do modelo no conjunto de teste 40------------------------------------------------------------------------------------------------------------------------

base_teste_40 %<>% mutate(predito = predict(modelo40, newdata = base_teste_40)) 

rmse40_tst <- sqrt(mean((base_teste_40$`40yd` - base_teste_40$predito)^2))
#0,1340 segundos

pseudor2_40_tst <- (cor(base_teste_40$`40yd`,base_teste_40$predito))^2
#82,10%

### Modelo broad jump (Salto em distancia - cm) -------------------------------------------------------------------------

set.seed(52)

# Base broad jump
base_bj <- base_tratada %>% 
  select(Player,ano,GPosBJ,cp1,altura,peso,`Broad Jump`) %>% 
  filter(!is.na(`Broad Jump`))

teste <- sample(1:nrow(base_bj), size = round(0.20 * nrow(base_bj)), replace = F)
base_teste_bj <- base_bj[teste,] # 20% da base para pegar os resultados depois.
base_mod_bj <- base_bj[-teste,] # 80% da base para desenvolver o modelo

fwrite(base_teste_bj, "base_teste_bj.csv")

base_mod_bj %>% 
  ggplot(aes(x = `Broad Jump`)) +
  geom_histogram(fill = 'darkolivegreen4', color = 'darkolivegreen', bins = 20) +
  labs(x = "Distancia (cm)", y = "Contagem") +
  theme_minimal()

#Gerando o modelo gamma - Funcao log, identidade e inversa
modelo1 <- glm(formula = `Broad Jump` ~ cp1 + GPosBJ, family = inverse.gaussian (link = "identity"), data = base_mod_bj)
modelo2 <- glm(formula = `Broad Jump` ~ cp1 + GPosBJ, family = inverse.gaussian (link = "inverse"), data = base_mod_bj)
modelo3 <- glm(formula = `Broad Jump` ~ cp1 + GPosBJ, family = inverse.gaussian (link = "log"), data = base_mod_bj)
modelo4 <- glm(formula = `Broad Jump` ~ cp1 + GPosBJ, family = inverse.gaussian (link = "1/mu^2"), data = base_mod_bj)
modelo5 <- glm(formula = `Broad Jump` ~ cp1 + GPosBJ, family = gaussian (link = "identity"), data = base_mod_bj)
modelo6 <- glm(formula = `Broad Jump` ~ cp1 + GPosBJ, family = gaussian (link = "inverse"), data = base_mod_bj)
modelo7 <- glm(formula = `Broad Jump` ~ cp1 + GPosBJ, family = gaussian (link = "log"), data = base_mod_bj)

#Testando os modelos

teste_modelos(list(modelo1,modelo2,modelo3,modelo4,modelo5,modelo6,modelo7))

# Funcao desvio (4016) < valor teste entao modelo e adequado.
# Nao ha necessidade de fazer recategorizacao
# Melhor modelo = Modelo 6, gaussiano com ligacao INVERSA

# Analise de residuo

analise_mod_bj <- analise_res(modelo6,vetor_func_lig = 1/base_mod_bj$`Broad Jump`) 
#graficos a direita.

analise_mod_bj$pseudoR2 #Pseudo R2 = 61,84
analise_mod_bj$verificao_funcao_ligacao #Coeficiente angular significativo, bom ajuste da Funcao de ligacao.
analise_mod_bj$teste_normalidade_lt #Segundo teste, desvios seguem uma distribuicao normal (Alpha = 0.01)

devres = analise_mod_bj$devres

# Analise de diagnostico
# Pontos Aberrantes, de Alavanca e Influentes

# Pontos aberrantes serao considerados valores onde |DevRes| > 3

analise_diag(devres,base_mod_bj,modelo6)
#18 pontos aberrantes, 352 pontos de alavanca, 0 pontos de influencia
#16 observacoes a serem retiradas

obs_retirar <- analise_diag(devres,base_mod_bj,modelo6)[[5]]

# Remodelagem excluindo aberrantes nao influentes (1)-------------------------------------------------------------------------------------------------------

base_retirada <- base_mod_bj[-obs_retirar,]

modelo_1_1 <- glm(formula = `Broad Jump` ~ cp1 + GPosBJ, family = gaussian (link = "inverse"), data = base_retirada)

#Testando os modelos

teste_modelos(list(modelo_1_1))

# Funcao desvio (4000) < valor teste entao modelo e adequado.
# Nao ha necessidade de fazer recategorizacao
# Melhor modelo = Modelo 1, gaussiana com Funcao de ligacao inversa

# Analise de residuo

analise_mod_bj_1 <- analise_res(modelo_1_1,vetor_func_lig = 1/base_retirada$`Broad Jump`) 
#graficos a direita.

analise_mod_bj_1$pseudoR2 #Pseudo R2 = 62,48
analise_mod_bj_1$verificao_funcao_ligacao #Coeficiente angular significativo, bom ajuste da Funcao de ligacao.
analise_mod_bj_1$teste_normalidade_lt #Segundo teste, desvios seguem uma distribuicao normal (Alpha = 0.01)

devres_1 = analise_mod_bj_1$devres

# Analise de diagnostico
# Pontos Aberrantes, de Alavanca e Influentes

# Pontos aberrantes serao considerados valores onde |DevRes| > 3

analise_diag(devres_1,base_retirada,modelo_1_1)
#4 pontos aberrantes, 345 pontos de alavanca, 0 pontos de influencia
#4 observacoes a serem retiradas

obs_retirar_1 <- analise_diag(devres_1,base_retirada,modelo_1_1)[[5]]

# Remodelagem excluindo aberrantes nao influentes (2)-------------------------------------------------------------------------------------------------------

base_retirada_2 <- base_retirada[-obs_retirar_1,]

modelo_1_2 <- glm(formula = `Broad Jump` ~ cp1 + GPosBJ, family = gaussian (link = "inverse"), data = base_retirada_2)

#Testando os modelos

teste_modelos(list(modelo_1_2))

# Funcao desvio (3996) < valor teste entao modelo e adequado.
# Nao ha necessidade de fazer recategorizacao
# Melhor modelo = Modelo 1, gaussiana com Funcao de ligacao inversa

# Analise de residuo

analise_mod_bj_2 <- analise_res(modelo_1_2,vetor_func_lig = 1/base_retirada_2$`Broad Jump`)
#graficos a direita.

analise_mod_bj_2$pseudoR2 #Pseudo R2 = 62,81
analise_mod_bj_2$verificao_funcao_ligacao #Coeficiente angular significativo, bom ajuste da Funcao de ligacao.
analise_mod_bj_2$teste_normalidade_lt #Segundo teste, desvios seguem uma distribuicao normal (Alpha = 0.01)

devres_2 = analise_mod_bj_2$devres

# Analise de diagnostico
# Pontos Aberrantes, de Alavanca e Influentes

# Pontos aberrantes serao considerados valores onde |DevRes| > 3

analise_diag(devres_2,base_retirada_2,modelo_1_2)
#0 pontos aberrantes, 345 pontos de alavanca, 0 pontos de influencia
#0 observacoes a serem retiradas

modeloBJ <- modelo_1_2
#Modelo final para o broad jump, modelo gaussiano com ligacao inversa

#RMSE do modelo de treino
rmseBJ_tr <- sqrt(mean((modeloBJ$y - fitted(modeloBJ))^2)) #14,81

saveRDS(modeloBJ,"modelo_BJ.rds")

# Analise do modelo no conjunto de teste BJ------------------------------------------------------------------------------------------------------------------------

base_teste_bj %<>% mutate(predito = 1/predict(modeloBJ, newdata = base_teste_bj)) 

rmseBJ_tst <- sqrt(mean((base_teste_bj$`Broad Jump` - base_teste_bj$predito)^2))
#16,22 cm

pseudor2_bj_tst <- (cor(base_teste_bj$`Broad Jump`,base_teste_bj$predito))^2
#53,82%

### Modelo Bench Press(Teste forca - Supino QTDE) ----------------------------------------------------------------------------

set.seed(52)

# Base bench
base_bench <- base_tratada %>% 
  select(Player,ano,GPosBench,cp1,altura,peso,Bench) %>% 
  filter(!is.na(Bench))

teste <- sample(1:nrow(base_bench), size = round(0.20 * nrow(base_bench)), replace = F)
base_teste_bench <- base_bench[teste,] # 20% da base para pegar os resultados depois.
base_mod_bench <- base_bench[-teste,] # 80% da base para desenvolver o modelo

fwrite(base_teste_bench, "base_teste_bp.csv")

base_mod_bench$Bench %>% mean() #20,84
base_mod_bench$Bench %>% var() #39,94
#O modelo poisson nao e uma opcao pois a variancia da variavel resposta e o dobro da media.


#Gerando o modelo gamma - Funcao log, identidade e inversa
modelo1 <- glm(formula = Bench ~ cp1 + GPosBench, family = quasipoisson (link = "identity"), data = base_mod_bench)
modelo2 <- glm(formula = Bench ~ cp1 + GPosBench, family = quasipoisson (link = "log"), data = base_mod_bench)
modelo3 <- glm(formula = Bench ~ cp1 + GPosBench, family = quasipoisson (link = "sqrt"), data = base_mod_bench)
modelo4 <- glm(formula = Bench ~ cp1 + GPosBench, family = quasipoisson (link = "inverse"), data = base_mod_bench)

#Testando os modelos

teste_modelos(list(modelo1,modelo2,modelo3,modelo4))

# Funcao desvio (3720.60) < valor teste entao modelo e adequado.
# E necessario reagrupar a variavel de posicao
# Melhor modelo = Modelo 1, quasipoisson com Funcao de ligacao identidade

# Reajustando

base_mod_bench[which(base_mod_bench$GPosBench == "g5"),]$GPosBench = "g4"
modelo1 <- glm(formula = Bench ~ cp1 + GPosBench, family = quasipoisson (link = "identity"), data = base_mod_bench)

teste_modelos(list(modelo1))

# Funcao desvio (3722.02) < valor teste entao modelo e adequado.
# Nao ha necessidade de selecionar variÃ¡veis.
# Melhor modelo = Modelo 1, quasipoisson com Funcao de ligacao identidade

# Analise de residuo

analise_mod_bench <- analise_res(modelo1,vetor_func_lig = base_mod_bench$Bench) 
#graficos a direita.

analise_mod_bench$pseudoR2 #Pseudo R2 = 43,36
analise_mod_bench$verificao_funcao_ligacao #Coeficiente angular significativo, bom ajuste da Funcao de ligacao.
analise_mod_bench$teste_normalidade_lt #Segundo teste, desvios seguem uma distribuicao normal (Alpha = 0.01)

devres = analise_mod_bench$devres

# Analise de diagnostico
# Pontos Aberrantes, de Alavanca e Influentes

# Pontos aberrantes serao considerados valores onde |DevRes| > 3

analise_diag(devres,base_mod_bench,modelo1)
#24 pontos aberrantes, 511 pontos de alavanca, 0 pontos de influencia
#22 observacoes a serem retiradas

obs_retirar <- analise_diag(devres,base_mod_bench,modelo1)[[5]]

# Remodelagem excluindo aberrantes nao influentes (1)-------------------------------------------------------------------------------------------------------

base_retirada <- base_mod_bench[-obs_retirar,]

modelo_1_1 <- glm(formula = Bench ~ cp1 + GPosBench, family = quasipoisson (link = "identity"), data = base_retirada)

#Testando os modelos

teste_modelos(list(modelo_1_1))

# Funcao desvio (3681,79) < valor teste entao modelo e adequado.
# Nao ha necessidade de fazer recategorizacao
# Melhor modelo = Modelo 1, quasipoisson com Funcao de ligacao identidade

# Analise de residuo

analise_mod_bench_1 <- analise_res(modelo_1_1,vetor_func_lig = base_retirada$Bench) 
#graficos a direita.

analise_mod_bench_1$pseudoR2 #Pseudo R2 = 44,38
analise_mod_bench_1$verificao_funcao_ligacao #Coeficiente angular significativo, bom ajuste da Funcao de ligacao.
analise_mod_bench_1$teste_normalidade_lt #Segundo teste, desvios seguem uma distribuicao normal (Alpha = 0.01)

devres_1 = analise_mod_bench_1$devres

# Analise de diagnostico
# Pontos Aberrantes, de Alavanca e Influentes

# Pontos aberrantes serao considerados valores onde |DevRes| > 3

analise_diag(devres_1,base_retirada,modelo_1_1)
#3 pontos aberrantes, 501 pontos de alavanca, 0 pontos de influencia
#2 observacoes a serem retiradas

obs_retirar_1 <- analise_diag(devres_1,base_retirada,modelo_1_1)[[5]]

# Remodelagem excluindo aberrantes nao influentes (2)-------------------------------------------------------------------------------------------------------

base_retirada_2 <- base_retirada[-obs_retirar_1,]

modelo_1_2 <- glm(formula = Bench ~ cp1 + GPosBench, family = quasipoisson (link = "identity"), data = base_retirada_2)

#Testando os modelos

teste_modelos(list(modelo_1_2))

# Funcao desvio (3677,59) < valor teste entao modelo e adequado.
# Nao ha necessidade de fazer recategorizacao
# Melhor modelo = Modelo 1, quasipoisson com Funcao de ligacao identidade

# Analise de residuo

analise_mod_bench_2 <- analise_res(modelo_1_2,vetor_func_lig = base_retirada_2$Bench) 
#graficos a direita.

analise_mod_bench_2$pseudoR2 #Pseudo R2 = 44,38
analise_mod_bench_2$verificao_funcao_ligacao #Coeficiente angular significativo, bom ajuste da Funcao de ligacao.
analise_mod_bench_2$teste_normalidade_lt #Segundo teste, desvios seguem uma distribuicao normal (Alpha = 0.01)

devres_2 = analise_mod_bench_2$devres

# Analise de diagnostico
# Pontos Aberrantes, de Alavanca e Influentes

# Pontos aberrantes serao considerados valores onde |DevRes| > 3

analise_diag(devres_2,base_retirada_2,modelo_1_2)
#1 ponto aberrante, 501 pontos de alavanca, 0 pontos de influencia
#1 observacao a ser retirada (nao vale a pena retirar somente uma observacao)

modeloBENCH <- modelo_1_2
#Modelo final para Bench, modelo quasipoisson com ligacao identidade

#RMSE do modelo de treino
rmseBENCH_tr <- sqrt(mean((modeloBENCH$y - round(fitted(modeloBENCH)))^2)) #4,63 [5] repeticoes

saveRDS(modeloBENCH,"modelo_bench.rds")

# Analise do modelo no conjunto de teste BENCH PRESS------------------------------------------------------------------------------------------------------------------------

base_teste_bench[which(base_teste_bench$GPosBench == "g5"),]$GPosBench = "g4"

base_teste_bench %<>% mutate(predito = round(predict(modeloBENCH, newdata = base_teste_bench))) 

rmseBENCH_tst <- sqrt(mean((base_teste_bench$Bench - base_teste_bench$predito)^2))
#4,83 [5] repeticoes

pseudor2_bj_tst <- (cor(base_teste_bench$Bench,base_teste_bench$predito))^2
#42,27%

### Modelo 3 cone (Teste agilidade - Segundos) ----------------------------------------------------------------------------

set.seed(52)

# Base 3 cone
base_3cone <- base_tratada %>% 
  select(Player,ano,GPos3CONE,cp1,altura,peso,`3Cone`) %>% 
  filter(!is.na(`3Cone`))

teste <- sample(1:nrow(base_3cone), size = round(0.20 * nrow(base_3cone)), replace = F)
base_teste_3cone <- base_3cone[teste,] # 20% da base para pegar os resultados depois.
base_mod_3cone <- base_3cone[-teste,] # 80% da base para desenvolver o modelo

fwrite(base_teste_3cone, "base_teste_3c.csv")

base_mod_3cone %>% 
  ggplot(aes(x = `3Cone`)) +
  geom_histogram(fill = 'darkolivegreen4', color = 'darkolivegreen', bins = 20) +
  labs(x = "Tempo (s)", y = "Contagem") +
  theme_minimal()

#Gerando o modelo gamma - Funcao log, identidade e inversa
modelo1 <- glm(formula = `3Cone` ~ cp1 + GPos3CONE, family = Gamma (link = "identity"), data = base_mod_3cone)
modelo2 <- glm(formula = `3Cone` ~ cp1 + GPos3CONE, family = Gamma (link = "log"),      data = base_mod_3cone)
modelo3 <- glm(formula = `3Cone` ~ cp1 + GPos3CONE, family = Gamma (link = "inverse"),  data = base_mod_3cone)
modelo4 <- glm(formula = `3Cone` ~ cp1 + GPos3CONE, family = inverse.gaussian (link = "identity"),  data = base_mod_3cone)
modelo5 <- glm(formula = `3Cone` ~ cp1 + GPos3CONE, family = inverse.gaussian (link = "inverse"),  data = base_mod_3cone)

#Testando os modelos

teste_modelos(list(modelo1,modelo2,modelo3,modelo4,modelo5))

# Funcao desvio (3198,21) < valor teste entao modelo e adequado.
# nao e necessario retirar/reagrupar nenhuma variavel
# Melhor modelo = Modelo 4, Gaussiana inversa com Funcao de ligacao identidade

# Analise de residuo

analise_mod_3cone <- analise_res(modelo4,vetor_func_lig = base_mod_3cone$`3Cone`) 
#graficos a direita.

analise_mod_3cone$pseudo #Pseudo  = 67,48
analise_mod_3cone$verificao_funcao_ligacao #Coeficiente angular significativo, bom ajuste da Funcao de ligacao.
analise_mod_3cone$teste_normalidade_lt #Segundo teste, desvios nao seguem uma distribuicao normal (Alpha = 0.01)

devres = analise_mod_3cone$devres

# Analise de diagnostico
# Pontos Aberrantes, de Alavanca e Influentes

# Pontos aberrantes serao considerados valores onde |DevRes| > 0.03 (Escala dos residuos diferente)

analise_diag(devres,base_mod_3cone,modelo4,ref = 0.03)
#54 pontos aberrantes, 286 pontos de alavanca, 0 pontos de influencia
#47 observacoes a serem retiradas

obs_retirar <- analise_diag(devres,base_mod_3cone,modelo4,ref = 0.03)[[5]]

# Remodelagem excluindo aberrantes nao influentes (1)-------------------------------------------------------------------------------------------------------

base_retirada <- base_mod_3cone[-obs_retirar,]

modelo_1_1 <- glm(formula = `3Cone` ~ cp1 + GPos3CONE, family = inverse.gaussian (link = "identity"), data = base_retirada)

#Testando os modelos

teste_modelos(list(modelo_1_1))

# Funcao desvio (3170,61) < valor teste entao modelo e adequado.
# Nao ha necessidade de fazer recategorizacao
# Melhor modelo = Modelo 1, gaussiana inversa com Funcao de ligacao identidade

# Analise de residuo

analise_mod_3cone_1 <- analise_res(modelo_1_1,vetor_func_lig = base_retirada$`3Cone`) 
#graficos a direita.

analise_mod_3cone_1$pseudo #Pseudo  = 69,84
analise_mod_3cone_1$verificao_funcao_ligacao #Coeficiente angular significativo, bom ajuste da Funcao de ligacao.
analise_mod_3cone_1$teste_normalidade_lt #Segundo teste, desvios seguem uma distribuicao normal (Alpha = 0.01)

devres_1 = analise_mod_3cone_1$devres

# Analise de diagnostico
# Pontos Aberrantes, de Alavanca e Influentes

# Pontos aberrantes serao considerados valores onde |DevRes| > 3

analise_diag(devres_1,base_retirada,modelo_1_1, ref = 0.03)
#8 pontos aberrantes, 273 pontos de alavanca, 0 pontos de influencia
#7 observacoes a serem retiradas

obs_retirar_1 <- analise_diag(devres_1,base_retirada,modelo_1_1, ref = 0.03)[[5]]

# Remodelagem excluindo aberrantes nao influentes (2)-------------------------------------------------------------------------------------------------------

base_retirada_2 <- base_retirada[-obs_retirar_1,]

modelo_1_2 <- glm(formula = `3Cone` ~ cp1 + GPos3CONE, family = inverse.gaussian (link = "identity"), data = base_retirada_2)

#Testando os modelos

teste_modelos(list(modelo_1_2))

# Funcao desvio (3163,1) < valor teste entao modelo e adequado.
# Nao ha necessidade de fazer recategorizacao
# Melhor modelo = Modelo 1, gaussiana inversa com Funcao de ligacao identidade

# Analise de residuo

analise_mod_3cone_2 <- analise_res(modelo_1_2,vetor_func_lig = base_retirada_2$`3Cone`) 
#graficos a direita.

analise_mod_3cone_2$pseudo #Pseudo  = 70,31
analise_mod_3cone_2$verificao_funcao_ligacao #Coeficiente angular significativo, bom ajuste da Funcao de ligacao.
analise_mod_3cone_2$teste_normalidade_lt #Segundo teste, desvios seguem uma distribuicao normal (Alpha = 0.01)

devres_2 = analise_mod_3cone_2$devres

# Analise de diagnostico
# Pontos Aberrantes, de Alavanca e Influentes

# Pontos aberrantes serao considerados valores onde |DevRes| > 0.03

analise_diag(devres_2,base_retirada_2,modelo_1_2, ref = 0.03)
#1 pontos aberrantes, 273 pontos de alavanca, 0 pontos de influencia
#1 observacao a ser retirada, serÃ¡ mantida.

modelo3Cone <- modelo_1_2
#Modelo final para Bench, modelo quasipoisson com ligacao identidade


#RMSE do modelo de treino
rmse3CONE_tr <- sqrt(mean((modelo3Cone$y - round(fitted(modelo3Cone)))^2)) #0,288 segundos

saveRDS(modelo3Cone,"modelo_3Cone.rds")

# Analise do modelo no conjunto de teste 3 Cones------------------------------------------------------------------------------------------------------------------------

base_teste_3cone %<>% mutate(predito = predict(modelo3Cone, newdata = base_teste_3cone)) 

rmse3CONE_tst <- sqrt(mean((base_teste_3cone$`3Cone` - base_teste_3cone$predito)^2))
#0,244 segundos

pseudor2_3CONE_tst <- (cor(base_teste_3cone$`3Cone`,base_teste_3cone$predito))^2
#67,46%

# Exportar base tratada

fwrite(base_tratada, "base_tratada.csv")

## Analises finais -------------------------------------------------------------

rm(list = ls())

#Bases -------------------------------------------------------------------------

base_tratada <- data.table::fread("base_tratada.csv")

modelo40 <- readRDS("modelo_40yds.rds") #Gamma ligacao identidade, n = 4.874
modelobj <- readRDS("modelo_BJ.rds")  #Normal ligacao inversa, n = 4.000
modelobp <- readRDS("modelo_bench.rds") #Quasipoisson ligacao identidade, n = 3.666
modelo3c <- readRDS("modelo_3Cone.rds") #Guassiana inversa funcao identidade, n = 3.177

# Graficos media estimado x real

#40
base_teste_40 <- fread("base_teste_40.csv") 
base_teste_40 %<>% mutate(predito_40 = predict(modelo40, newdata = base_teste_40))

graf_40_realxpred <- base_teste_40 %>% 
  select(ano,`Valor real` = `40yd`,`Valor predito` = predito_40) %>% 
  gather(key = "Origem", value = "Resultado", -ano) %>% 
  group_by(ano,Origem) %>% 
  summarise(media = mean(Resultado)) %>% 
  ggplot(aes(x = ano, y = media, color = Origem)) +
  geom_line(size = 1.2) +
  labs(x = "Ano", y = "Tempo (s)", color = "") +
  theme_minimal()

#Bench
base_teste_bp <- data.table::fread("base_teste_bp.csv") 
base_teste_bp[which(base_teste_bp$GPosBench == "g5"),]$GPosBench = "g4"

base_teste_bp %<>% mutate(predito_bp = round(predict(modelobp, newdata = base_teste_bp)))

graf_bp_realxpred <- base_teste_bp %>% 
  select(ano,`Valor real` = Bench,`Valor predito` = predito_bp) %>% 
  gather(key = "Origem", value = "Resultado", -ano) %>% 
  group_by(ano,Origem) %>% 
  summarise(media = mean(Resultado)) %>% 
  ggplot(aes(x = ano, y = media, color = Origem)) +
  geom_line(size = 1.2) +
  labs(x = "Ano", y = "Repetições", color = "") +
  theme_minimal()

#Broad jump
base_teste_bj <- data.table::fread("base_teste_bj.csv") 

base_teste_bj %<>% mutate(predito_bj = 1/predict(modelobj, newdata = base_teste_bj))

graf_bj_realxpred <- base_teste_bj %>% 
  select(ano,`Valor real` = `Broad Jump`,`Valor predito` = predito_bj) %>% 
  gather(key = "Origem", value = "Resultado", -ano) %>% 
  group_by(ano,Origem) %>% 
  summarise(media = mean(Resultado)) %>% 
  ggplot(aes(x = ano, y = media, color = Origem)) +
  geom_line(size = 1.2) +
  labs(x = "Ano", y = "Distância (cm)", color = "") +
  theme_minimal()

#3 Cone
base_teste_3c <- data.table::fread("base_teste_3c.csv") 

base_teste_3c %<>% mutate(predito_3c = predict(modelo3c, newdata = base_teste_3c))

graf_3c_realxpred <- base_teste_3c %>% 
  select(ano,`Valor real` = `3Cone`,`Valor predito` = predito_3c) %>% 
  gather(key = "Origem", value = "Resultado", -ano) %>% 
  group_by(ano,Origem) %>% 
  summarise(media = mean(Resultado)) %>% 
  ggplot(aes(x = ano, y = media, color = Origem)) +
  geom_line(size = 1.2) +
  labs(x = "Ano", y = "Tempo (s)", color = "") +
  theme_minimal()


# Base de teste global

set.seed(52)

base_filtrada <- base_tratada %>% filter(!is.na(`40yd`), !is.na(`3Cone`),!is.na(`Broad Jump`),!is.na(Bench))

# Gerando resultados a partir do modelo

base_filtrada[which(base_filtrada$GPosBench == "g5"),]$GPosBench = "g4"

base_filtrada %<>% mutate(predito_40 = predict(modelo40, newdata = base_filtrada),
                     indice_40 = predito_40 - `40yd`,
                     ind_pad_40 = ((indice_40 - min(indice_40))/(max(indice_40) - min(indice_40)))*100,
                     
                     predito_bj = 1/predict(modelobj, newdata = base_filtrada),
                     indice_bj = `Broad Jump` - predito_bj,
                     ind_pad_bj = ((indice_bj - min(indice_bj))/(max(indice_bj) - min(indice_bj)))*100,
                     
                     predito_bp = round(predict(modelobp, newdata = base_filtrada)),
                     indice_bp = Bench - predito_bp,
                     ind_pad_bp = ((indice_bp - min(indice_bp))/(max(indice_bp) - min(indice_bp)))*100,
                     
                     predito_3c = predict(modelo3c, newdata = base_filtrada),
                     indice_3c = predito_3c - `3Cone`,
                     ind_pad_3c = ((indice_3c - min(indice_3c))/(max(indice_3c) - min(indice_3c)))*100)


base_filtrada %>% 
  group_by(ano) %>% 
  summarise(media = mean(ind_pad_40)) %>% 
  ggplot(aes(x = ano, y = media)) +
  geom_line(size = 1, color = "darkolivegreen") +
  labs(x = "Ano", y = "Indíce padronizado", color = "") +
  theme_minimal()

base_filtrada %>% 
  group_by(ano) %>% 
  summarise(media = mean(ind_pad_bp)) %>% 
  ggplot(aes(x = ano, y = media)) +
  geom_line(size = 1, color = "darkolivegreen") +
  labs(x = "Ano", y = "Indíce padronizado", color = "") +
  theme_minimal()

base_filtrada %>% 
  group_by(ano) %>% 
  summarise(media = mean(ind_pad_bj)) %>% 
  ggplot(aes(x = ano, y = media)) +
  geom_line(size = 1, color = "darkolivegreen") +
  labs(x = "Ano", y = "Indíce padronizado", color = "") +
  theme_minimal()

base_filtrada %>% 
  group_by(ano) %>% 
  summarise(media = mean(ind_pad_3c)) %>% 
  ggplot(aes(x = ano, y = media)) +
  geom_line(size = 1, color = "darkolivegreen") +
  labs(x = "Ano", y = "Indíce padronizado", color = "") +
  theme_minimal()


## Pca dos indices ----------------------------------------------------------------------------------------------------------

basepca <- base_filtrada %>% select(ind_pad_40,ind_pad_bj,ind_pad_bp,ind_pad_3c)

basepca %>% summarise(m_ind_40 = mean(ind_pad_40),
                      d_ind_40 = sd(ind_pad_40),
                      m_ind_bp = mean(ind_pad_bp),
                      d_ind_bp = sd(ind_pad_bp),
                      m_ind_bj = mean(ind_pad_bj),
                      d_ind_bj = sd(ind_pad_bj),
                      m_ind_3c = mean(ind_pad_3c),
                      d_ind_3c = sd(ind_pad_3c))

#Media dos indices variando entre 41,82 ate 62,75. Desvio padroes entre 10,50 e 13,52

pca <- prcomp(basepca, center = F)
summary(pca)
# CP1 97% da variabilidade

summary(pca)$rotation
# CP1 
#40 = -0.5859289
#bj = -0.4508889
#bp = -0.3916531
#3c = -0.5477175

fviz_eig(pca, fill = "lightbue")

fviz_pca_var(pca,
             col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800"),
             repel = TRUE,     
             legend.title = "Contribution")

#Distribuicao do componente 1 

componente1 <- pca$x[,1]
base_indice <- base_filtrada %>% cbind(.,componente1) %>% rename(indice_global = componente1) %>% 
  mutate(indice_global = -1*indice_global,
         indice_global =((indice_global - min(indice_global))/(max(indice_global) - min(indice_global)))*100)

base_indice %>% 
  group_by(ano) %>% 
  summarise(indice_global_md = mean(indice_global)) %>% 
  ggplot(., aes(x = ano,y = indice_global_md)) +
  geom_line(size = 1, color = "darkolivegreen") +
  labs(x = "Ano", y = "Indíce padronizado") +
  theme_minimal()


# medida APC

pw <- prais::prais_winsten(indice_global ~ ano, data = base_indice, index = "ano")
summary(pw)

ep = 0.03396
b= pw$coeff[2]
gl =pw$df.res
t = qt(0.97598,gl)
minb = b -t*ep
maxb = b + t*ep
APC = -1+10^b
minAPC = -1+10^minb
maxAPC = -1+10^maxb
ICAPC95 = c(minAPC, maxAPC)
ICAPC95
APCperc = 100*APC
ICAPC95perc = 100*c(minAPC, maxAPC)
ICAPC95perc

