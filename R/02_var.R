

## =========================================================
## 02_var.R
## Monografia: dinamica inflacionaria brasileira
##
## Estima o VAR e produz os resultados do Capitulo 4.
##
## ENTRADA : painel_var.csv   (gerado por 01_dados.R)
##
## MAPA DOS RESULTADOS
##   1. Selecao de defasagens ......... Secao 3.3.1
##   2. Diagnosticos do VAR ........... Secao 3.3.1 / Apendice A
##   3. Causalidade de Granger ........ Secao 3.4.3
##   4. Funcoes de resposta a impulso .. Secao 3.4.1 / Tabela 4.1
##   5. Papel do cambio ............... Secao 3.4.2 / Tabela 4.2
##   6. Decomposicao da variancia ..... Secao 3.4.3 / Tabela 4.3
##   7. Subamostras ................... Secao 3.5.1 / Tabela 4.4
##   8. Assimetria .................... Secao 3.5.2 / Tabela 4.5
##   9. Robustez da ordenacao ......... Secao 4.6 / Tabela 4.7
## =========================================================

## ---------------------------------------------------------
## Preparacao do ambiente
## ---------------------------------------------------------
## Instala os pacotes que faltarem, para o script rodar em qualquer
## computador sem preparacao previa.

pacotes <- c("vars", "dynlm", "lmtest")
faltando <- pacotes[!pacotes %in% rownames(installed.packages())]
if (length(faltando) > 0) install.packages(faltando, repos = "https://cloud.r-project.org")

## Coloca a pasta de trabalho na pasta deste arquivo. Sem isso, o R
## procura os dados em qualquer pasta que estiver aberta, e nao acha.
definir_pasta <- function() {
  # 1) rodando pelo RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    caminho <- try(rstudioapi::getSourceEditorContext()$path, silent = TRUE)
    if (!inherits(caminho, "try-error") && nzchar(caminho)) {
      setwd(dirname(caminho)); return(invisible(TRUE))
    }
  }
  # 2) rodando pelo terminal com Rscript
  args <- commandArgs(trailingOnly = FALSE)
  arq  <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(arq) > 0) { setwd(dirname(normalizePath(arq))); return(invisible(TRUE)) }
  invisible(FALSE)
}
try(definir_pasta(), silent = TRUE)
cat("Pasta de trabalho:", getwd(), "\n")

library(vars)
library(dynlm)
library(lmtest)

set.seed(20260819)     # fixa o bootstrap das bandas de confianca

## ---------------------------------------------------------
## Onde estao os arquivos
## ---------------------------------------------------------
## Funciona tanto com as pastas do projeto quanto com tudo junto
## numa pasta so, como no repositorio do GitHub.

achar <- function(nomes, pastas = c(".", "Paineis", "../Paineis", "..")) {
  for (n in nomes) for (p in pastas) {
    caminho <- file.path(p, n)
    if (file.exists(caminho)) return(caminho)
  }
  stop("Nao encontrei: ", paste(nomes, collapse = " nem "),
       "\nProcurei em: ", paste(pastas, collapse = ", "),
       "\nPasta de trabalho atual: ", getwd(),
       "\nRode o 01_dados.R antes deste script.")
}

pasta_saida <- function(nome) {
  for (p in c(nome, file.path("..", nome))) if (dir.exists(p)) return(p)
  dir.create(nome, showWarnings = FALSE)
  nome
}

PASTA_TABS <- pasta_saida("Tabs")
PASTA_GRAF <- pasta_saida("graficos")

## ---------------------------------------------------------
## Graficos
## ---------------------------------------------------------
## Cada grafico ocupa uma pagina inteira do PDF, em tamanho grande.
## O padrao do pacote vars empilha as sete respostas numa figura so,
## o que fica ilegivel; por isso os graficos principais sao desenhados
## aqui, um de cada vez.

# Grafico de uma funcao de resposta a impulso (FRI) acumulada,
# com banda de confianca.
plot_resposta <- function(objeto, impulso, resposta, escala, titulo, ylab) {
  m  <- objeto$irf[[impulso]][,   resposta] / escala
  lo <- objeto$Lower[[impulso]][, resposta] / escala
  hi <- objeto$Upper[[impulso]][, resposta] / escala
  h  <- 0:(length(m) - 1)
  plot(h, m, type = "n", ylim = range(c(lo, hi)),
       main = titulo, xlab = "meses apos o choque", ylab = ylab,
       cex.main = 1.4, cex.axis = 1.2, cex.lab = 1.2)
  polygon(c(h, rev(h)), c(lo, rev(hi)), col = "grey88", border = NA)
  grid(col = "grey92")
  abline(h = 0, col = "red", lty = 2)
  lines(h, m, lwd = 2.5, col = "grey10")
  legend("topleft", c("resposta estimada", "intervalo de 90%"),
         col = c("grey10", "grey88"), lwd = c(2.5, 8), bty = "n", cex = 1.1)
}

## ---------------------------------------------------------
## Montando a base
## ---------------------------------------------------------

dados <- read.csv(achar("painel_var.csv"))

## Confere se as colunas esperadas vieram todas, para o erro aparecer
## aqui, com a causa explicita, e nao dentro da funcao seguinte.
esperadas <- c("dlicbr", "dlpimp", "lvix", "hiato", "infl", "dselic", "dlcambio")
faltando <- setdiff(esperadas, names(dados))
if (length(faltando) > 0)
  stop("Faltam colunas em painel_var.csv: ", paste(faltando, collapse = ", "),
       "\nForam lidas: ", paste(names(dados), collapse = ", "),
       "\nRode o 01_dados.R novamente para regerar o arquivo.")

y <- ts(as.matrix(dados[, esperadas]), start = c(2003, 2), frequency = 12)

# A ordem das colunas E a ordenacao de Cholesky da Secao 3.3.2:
# bloco externo primeiro, cambio por ultimo.
colnames(y) <- c("comm", "pimp", "vix", "hiato", "ipca", "selic", "cambio")

cat("Painel:", nrow(y), "observacoes,", ncol(y), "variaveis\n")
head(y)
tail(y)

## Uma pagina por variavel, em graficos/02_series.pdf
rotulos <- c("Commodities (variacao % mensal)",
             "Precos de importacao (variacao % mensal)",
             "VIX (logaritmo)", "Hiato do produto (%)",
             "IPCA dessazonalizado (% a.m.)",
             "Selic (variacao mensal, p.p.)",
             "Cambio (variacao % mensal)")

pdf(file.path(PASTA_GRAF, "02_series.pdf"), width = 11, height = 6.5)
par(mar = c(4, 5, 4, 2))
for (j in 1:ncol(y)) {
  plot(y[, j], type = "l", lwd = 1.6, col = "grey20",
       main = rotulos[j], xlab = "", ylab = "",
       cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2)
  abline(h = 0, col = "red", lty = 2)
  grid(col = "grey88")
}
invisible(dev.off())

## ---------------------------------------------------------
## 1. Selecao de defasagens (Secao 3.3.1)
## ---------------------------------------------------------

VARselect(y, lag.max = 12, type = "const")

# Os criterios apontam ordens curtas, mas elas deixam autocorrelacao
# nos residuos. Testa-se de p = 1 a p = 8 pelo teste de Breusch-Godfrey
# e adota-se o menor p sem autocorrelacao a 5%.
#
# ATENCAO ao numero de defasagens do teste. O BG multivariado tem
# K^2 x lags.bg graus de liberdade. Com K = 7 variaveis e 280 observacoes,
# lags.bg = 1 usa 49 gl e lags.bg = 2 usa 98 gl, que cabem na amostra.
# Valores maiores nao cabem: lags.bg = 6 ja pediria 294 gl.

for (p in 1:8) {
  v  <- VAR(y, p = p, type = "const")
  b1 <- serial.test(v, lags.bg = 1, type = "BG")$serial$p.value
  b2 <- serial.test(v, lags.bg = 2, type = "BG")$serial$p.value
  cat(sprintf("p = %d  |  BG(1): p = %.4f   BG(2): p = %.4f   %s\n",
              p, b1, b2,
              ifelse(b1 >= 0.05 & b2 >= 0.05, "<-- sem autocorrelacao", "")))
}

# Escolha a menor ordem marcada acima e use-a na linha P <- ... abaixo.

## ---------------------------------------------------------
## 2. Estimacao e diagnosticos (Secao 3.3.1)
## ---------------------------------------------------------

P <- 5     # ordem escolhida na saida do bloco anterior
var1 <- VAR(y, p = P, type = "const")

summary(var1)

# Estabilidade: todas as raizes devem estar dentro do circulo unitario
roots(var1)
cat("\nMaior raiz:", round(max(roots(var1)), 4), "\n")

# Estabilidade estrutural dos parametros, uma equacao por pagina
# O plot do CUSUM monta sozinho um painel por equacao, entao a pagina
# precisa ser maior e as margens ficam por conta dele.
pdf(file.path(PASTA_GRAF, "03_cusum.pdf"), width = 11, height = 9)
ok <- try(plot(stability(var1, type = "OLS-CUSUM")), silent = TRUE)
invisible(dev.off())
if (inherits(ok, "try-error"))
  cat("Aviso: o grafico do CUSUM falhou. O teste em si nao e afetado.\n")

# Autocorrelacao dos residuos do modelo escolhido.
# lags.bg pequeno, pelo mesmo motivo explicado no bloco de selecao acima.
serial.test(var1, lags.bg = 1, type = "BG")
serial.test(var1, lags.bg = 2, type = "BG")

## ---------------------------------------------------------
## Reamostragem para os intervalos de confianca
## ---------------------------------------------------------
## Gera uma serie artificial a partir do VAR estimado, sorteando os
## residuos com reposicao. Cada replicacao do bootstrap reestima o
## modelo nessa serie e refaz a conta de interesse do zero.

simular <- function(v) {
  K  <- v$K
  p  <- v$p
  A  <- Bcoef(v)
  U  <- resid(v)
  Y  <- as.matrix(v$y)
  Tn <- nrow(Y)
  Ys <- matrix(0, Tn, K, dimnames = list(NULL, colnames(Y)))
  Ys[1:p, ] <- Y[1:p, ]
  sorteio <- sample(nrow(U), Tn - p, replace = TRUE)
  for (t in (p + 1):Tn) {
    x <- c(as.vector(t(Ys[(t - 1):(t - p), , drop = FALSE])), 1)
    Ys[t, ] <- A %*% x + U[sorteio[t - p], ]
  }
  ts(Ys, start = start(v$y), frequency = 12)
}

## ---------------------------------------------------------
## 3. Causalidade de Granger (Secao 3.4.3)
## ---------------------------------------------------------
## H0: a variavel indicada nao causa, no sentido de Granger,
## as demais variaveis do sistema.

## A funcao devolve uma lista; a estatistica e o valor-p do teste de
## Granger ficam em $Granger. A tabela abaixo reune os tres testes.

linhas_granger <- list()
for (v_ in c("comm", "pimp", "cambio")) {
  g <- causality(var1, cause = v_)$Granger
  linhas_granger[[v_]] <- data.frame(
    variavel  = v_,
    F         = round(as.numeric(g$statistic), 3),
    gl1       = g$parameter[1],
    gl2       = g$parameter[2],
    p_valor   = round(g$p.value, 4))
}
tab_granger <- do.call(rbind, linhas_granger)

cat("\n=== CAUSALIDADE DE GRANGER ===\n")
cat("H0: a variavel indicada nao causa, no sentido de Granger, as demais.\n")
print(tab_granger, row.names = FALSE)
write.csv(tab_granger, file.path(PASTA_TABS, "tab4_6_granger.csv"),
          row.names = FALSE)

## ---------------------------------------------------------
## 4. Funcoes de resposta a impulso, FRI (Secao 3.4.1, Tabela 4.1)
## ---------------------------------------------------------
## irf() com ortho = TRUE usa a decomposicao de Cholesky, e portanto
## depende da ordenacao das variaveis definida acima.
##
## O irf() devolve a resposta a um choque de um desvio padrao. Os objetos
## criados aqui servem aos graficos; o coeficiente de repasse da Tabela
## 4.1 e calculado logo adiante, dividindo a resposta acumulada do IPCA
## pela variacao acumulada da propria variavel de origem.

H <- 24

irf_comm <- irf(var1, impulse = "comm", n.ahead = H, ortho = TRUE,
                cumulative = TRUE, boot = TRUE, runs = 1000, ci = 0.90)
irf_pimp <- irf(var1, impulse = "pimp", n.ahead = H, ortho = TRUE,
                cumulative = TRUE, boot = TRUE, runs = 1000, ci = 0.90)
irf_camb <- irf(var1, impulse = "cambio", n.ahead = H, ortho = TRUE,
                cumulative = TRUE, boot = TRUE, runs = 1000, ci = 0.90)

# impacto do choque na propria variavel, usado so nos graficos
esc_comm <- irf_comm$irf$comm[1, "comm"]
esc_pimp <- irf_pimp$irf$pimp[1, "pimp"]
esc_camb <- irf_camb$irf$cambio[1, "cambio"]

## Os quatro graficos que entram no Capitulo 4, um por pagina.
pdf(file.path(PASTA_GRAF, "04_respostas.pdf"), width = 10, height = 6)
par(mar = c(5, 5, 4, 2))

plot_resposta(irf_comm, "comm", "ipca", esc_comm,
              "Resposta do IPCA a um choque de 1% nas commodities",
              "resposta acumulada (%)")
plot_resposta(irf_pimp, "pimp", "ipca", esc_pimp,
              "Resposta do IPCA a um choque de 1% nos precos de importacao",
              "resposta acumulada (%)")
plot_resposta(irf_camb, "cambio", "ipca", esc_camb,
              "Resposta do IPCA a uma depreciacao cambial de 1%",
              "resposta acumulada (%)")
plot_resposta(irf_comm, "comm", "cambio", esc_comm,
              "Resposta do cambio a um choque de 1% nas commodities",
              "resposta acumulada (%)")

invisible(dev.off())

## Figuras do Capitulo 4: as tres respostas do IPCA numa pagina so,
## primeiro pontuais e depois acumuladas. A versao pontual mostra em que
## mes o efeito aparece; a acumulada e a que se le como repasse.

irf_pt_comm <- irf(var1, impulse = "comm",   n.ahead = H, ortho = TRUE,
                   cumulative = FALSE, boot = TRUE, runs = 1000, ci = 0.90)
irf_pt_pimp <- irf(var1, impulse = "pimp",   n.ahead = H, ortho = TRUE,
                   cumulative = FALSE, boot = TRUE, runs = 1000, ci = 0.90)
irf_pt_camb <- irf(var1, impulse = "cambio", n.ahead = H, ortho = TRUE,
                   cumulative = FALSE, boot = TRUE, runs = 1000, ci = 0.90)

painel_tres <- function(arquivo, o1, o2, o3, sufixo) {
  pdf(file.path(PASTA_GRAF, arquivo), width = 15, height = 5)
  par(mfrow = c(1, 3), mar = c(5, 5, 4, 1))
  plot_resposta(o1, "comm",   "ipca", esc_comm,
                paste("Choque de 1% nas commodities", sufixo), "% do IPCA")
  plot_resposta(o2, "pimp",   "ipca", esc_pimp,
                paste("Choque de 1% nos precos de importacao", sufixo), "% do IPCA")
  plot_resposta(o3, "cambio", "ipca", esc_camb,
                paste("Depreciacao cambial de 1%", sufixo), "% do IPCA")
  invisible(dev.off())
}

painel_tres("08_respostas_pontuais.pdf",
            irf_pt_comm, irf_pt_pimp, irf_pt_camb, "(resposta mensal)")
painel_tres("09_respostas_acumuladas.pdf",
            irf_comm, irf_pimp, irf_camb, "(resposta acumulada)")

## Todas as respostas do sistema, uma por pagina, para o Apendice.
## Uma pagina por par choque-resposta, com o mesmo desenho dos graficos
## principais. Nao se usa o plot() do pacote vars aqui porque ele monta
## as sete respostas numa figura so e falha em paginas deste tamanho.

nomes_var <- c(comm = "commodities", pimp = "precos de importacao",
               vix = "VIX", hiato = "hiato do produto", ipca = "IPCA",
               selic = "Selic", cambio = "cambio")

choques <- list(list(obj = irf_comm, imp = "comm",   esc = esc_comm),
                list(obj = irf_pimp, imp = "pimp",   esc = esc_pimp),
                list(obj = irf_camb, imp = "cambio", esc = esc_camb))

pdf(file.path(PASTA_GRAF, "05_respostas_completas.pdf"), width = 10, height = 6)
par(mar = c(5, 5, 4, 2))
for (ch in choques) {
  for (resp in colnames(y)) {
    plot_resposta(ch$obj, ch$imp, resp, ch$esc,
                  paste0("Resposta de ", nomes_var[[resp]],
                         " a um choque de 1% em ", nomes_var[[ch$imp]]),
                  "resposta acumulada")
  }
}
invisible(dev.off())

## ---------------------------------------------------------
## Coeficiente de repasse (Tabela 4.1)
## ---------------------------------------------------------
## O repasse NAO e a resposta do IPCA dividida pelo impacto inicial do
## choque. A variavel de origem continua se movendo depois do primeiro
## mes, e o que interessa e quanto o nivel de precos subiu por cada 1%
## de variacao ACUMULADA dessa variavel. Divide-se, portanto, a resposta
## acumulada do IPCA pela resposta acumulada da propria variavel de
## origem, no mesmo horizonte.

repasse_pt <- function(v, H = 24, hh = c(6, 12, 24)) {
  sapply(c("comm", "pimp", "cambio"), function(sh) {
    ir <- irf(v, impulse = sh, response = c("ipca", sh), n.ahead = H,
              ortho = TRUE, cumulative = TRUE, boot = FALSE)
    (ir$irf[[sh]][, "ipca"] / ir$irf[[sh]][, sh])[hh + 1]
  })
}

pt_irf <- repasse_pt(var1, H)

## Como o repasse e uma razao entre duas respostas estimadas, o intervalo
## de confianca precisa vir de um bootstrap da propria razao. Usar as
## bandas que o irf() devolve para o numerador subestimaria a incerteza,
## por ignorar a do denominador.

REPS_IRF <- 1000
boot_irf <- array(NA_real_, c(REPS_IRF, 3, 3),
                  dimnames = list(NULL, c("h6", "h12", "h24"),
                                  c("comm", "pimp", "cambio")))
for (r in 1:REPS_IRF) {
  x <- try(repasse_pt(VAR(simular(var1), p = P, type = "const"), H),
           silent = TRUE)
  if (!inherits(x, "try-error")) boot_irf[r, , ] <- x
  if (r %% 100 == 0) cat("bootstrap Tabela 4.1:", r, "de", REPS_IRF, "\n")
}

nomes_choque <- c(comm = "Commodities (1%)", pimp = "P. importacao (1%)",
                  cambio = "Cambio (1%)")

tab_irf <- do.call(rbind, lapply(c("comm", "pimp", "cambio"), function(sh)
  data.frame(choque    = nomes_choque[[sh]],
             horizonte = c(6, 12, 24),
             repasse   = round(pt_irf[, sh], 3),
             ic_inf    = round(apply(boot_irf[, , sh], 2, quantile, 0.05, na.rm = TRUE), 3),
             ic_sup    = round(apply(boot_irf[, , sh], 2, quantile, 0.95, na.rm = TRUE), 3))))

cat("\n=== TABELA 4.1: REPASSE ACUMULADO AO IPCA (%) ===\n")
cat("resposta acumulada do IPCA / variacao acumulada da variavel de origem\n")
print(tab_irf, row.names = FALSE)
write.csv(tab_irf, file.path(PASTA_TABS, "tab4_1_irf.csv"), row.names = FALSE)

## ---------------------------------------------------------
## 5. O papel do cambio: contrafactual Sims-Zha (Secao 3.4.2, Tabela 4.2)
## ---------------------------------------------------------
## Decompoe a resposta da inflacao ao choque de commodities em
##
##     resposta total = resposta sem o canal do cambio + canal do cambio
##
## O contrafactual mantem o cambio parado na trajetoria de base. Para
## isso, a cada mes do horizonte injeta-se um choque cambial do tamanho
## exato do movimento que o cambio faria, com sinal contrario. Como a
## dinamica do sistema volta a mover o cambio no mes seguinte, a
## compensacao precisa ser refeita mes a mes.
##
## O procedimento e o de Sims e Zha (1995), tal como aplicado por
## Bernanke, Gertler e Watson (1997, p. 106).
##
## A diferenca entre as duas trajetorias recolhe TODOS os caminhos que
## passam pelo cambio, inclusive as realimentacoes (commodities -> cambio
## -> Selic -> inflacao, por exemplo), e nao apenas o caminho direto.

## Funcao de decomposicao, escrita na forma companheira do VAR.
## shock       : variavel que recebe o choque (commodities)
## mediator    : variavel a ser mantida parada (cambio)
## target      : variavel cuja resposta se quer decompor (IPCA)
## compensator : choque usado para segurar o mediador (o proprio cambio)
## shut_impact : TRUE  neutraliza o cambio desde o mes 0
##               FALSE neutraliza apenas a partir do mes 1

bgw_channel <- function(A_list, B, shock, mediator, target,
                        H = 24, compensator = mediator,
                        shut_impact = FALSE) {
  K <- nrow(B); p <- length(A_list); Kp <- K * p
  
  F <- matrix(0, Kp, Kp)
  F[1:K, ] <- do.call(cbind, A_list)
  if (p > 1)
    F[(K + 1):Kp, 1:(Kp - K)] <- diag(Kp - K)
  G <- rbind(B, matrix(0, Kp - K, K))
  J <- cbind(diag(K), matrix(0, K, Kp - K))
  
  q <- as.numeric(J[mediator, ])
  d <- as.numeric(G[, compensator])
  kappa <- sum(q * d)
  if (abs(kappa) < 1e-10)
    stop("O choque compensatorio nao move o mediador no impacto.")
  
  Pm <- diag(Kp) - outer(d, q) / kappa
  x_base <- as.numeric(G[, shock])
  x_cf <- if (shut_impact) Pm %*% x_base else x_base
  
  base <- cf <- matrix(NA_real_, H + 1, K)
  base[1, ] <- J %*% x_base
  cf[1, ]   <- J %*% x_cf
  
  for (h in 1:H) {
    x_base    <- F %*% x_base
    x_cf_pred <- F %*% x_cf
    nu_h      <- -sum(q * x_cf_pred) / kappa
    x_cf      <- x_cf_pred + d * nu_h
    base[h + 1, ] <- J %*% x_base
    cf[h + 1, ]   <- J %*% x_cf
  }
  
  data.frame(h = 0:H,
             total             = base[, target],
             sem_canal         = cf[, target],
             canal             = base[, target] - cf[, target],
             resposta_choque   = base[, shock],
             resposta_mediador = base[, mediator],
             mediador_cf       = cf[, mediator])
}

## ---------------------------------------------------------
## Aplicacao ao VAR estimado
## ---------------------------------------------------------

i_comm   <- which(colnames(y) == "comm")
i_cambio <- which(colnames(y) == "cambio")
i_ipca   <- which(colnames(y) == "ipca")

## A matriz de impacto B e o fator de Cholesky da matriz de covariancia
## dos residuos, na mesma ordenacao da Secao 3.3.2. E calculada aqui do
## mesmo jeito que o pacote vars faz por dentro do irf(), com correcao de
## graus de liberdade, para que o contrafactual e as FRIs da Tabela 4.1
## fiquem exatamente na mesma escala.

matriz_B <- function(v) {
  U     <- resid(v)
  npar  <- ncol(v$datamat) - v$K      # regressores por equacao: K*p + constante
  Sigma <- crossprod(U) / (v$obs - npar)
  B     <- t(chol(Sigma))
  dimnames(B) <- list(colnames(U), colnames(U))
  B
}

## Conferencia: a matriz B calculada aqui tem que dar o mesmo impacto
## que o irf() do pacote, senao as duas tabelas ficariam em escalas
## diferentes.
cat("\nConferencia da escala do choque de commodities\n")
cat("  pela matriz B :", round(matriz_B(var1)[i_comm, i_comm], 6), "\n")
cat("  pelo irf()    :", round(esc_comm, 6), "\n")

## Rotina completa: roda a decomposicao, acumula as respostas e
## normaliza para um choque de 1% nas commodities.
##
## As respostas saem mes a mes. Como o IPCA entra no VAR em variacao
## mensal, a soma das respostas ate h e o efeito acumulado sobre o
## NIVEL de precos, que e o que se quer reportar.

decompor <- function(v, H = 24, shut_impact = TRUE) {
  B  <- matriz_B(v)
  d  <- bgw_channel(Acoef(v), B, i_comm, i_cambio, i_ipca,
                    H = H, compensator = i_cambio, shut_impact = shut_impact)
  # Mesma normalizacao da Tabela 4.1: divide-se pela variacao ACUMULADA
  # das commodities no horizonte, e nao pelo impacto do primeiro mes.
  esc <- cumsum(d$resposta_choque)
  data.frame(h         = d$h,
             total     = cumsum(d$total)     / esc,
             sem_canal = cumsum(d$sem_canal) / esc,
             canal     = cumsum(d$canal)     / esc,
             cambio    = cumsum(d$resposta_mediador) / esc,
             cambio_cf = cumsum(d$mediador_cf)       / esc,
             erro      = max(abs(d$mediador_cf[if (shut_impact) TRUE else -1])))
}

cf <- decompor(var1, H, shut_impact = TRUE)

cat("\nConferencia: maior desvio do cambio no contrafactual =",
    format(cf$erro[1], digits = 3), "(deve ser praticamente zero)\n")

hh <- c(6, 12, 24)
tab_cambio <- data.frame(
  horizonte     = hh,
  total         = round(cf$total[hh + 1],     3),
  sem_canal     = round(cf$sem_canal[hh + 1], 3),
  canal_cambio  = round(cf$canal[hh + 1],     3))

cat("\n=== TABELA 4.2: DECOMPOSICAO DA RESPOSTA DO IPCA ===\n")
print(tab_cambio, row.names = FALSE)
cat("total       : resposta acumulada do IPCA ao choque de 1% nas commodities\n")
cat("sem_canal   : a mesma resposta com o cambio mantido parado\n")
cat("canal_cambio: diferenca entre as duas, negativa se o cambio amortece\n")

## ---------------------------------------------------------
## Robustez da convencao de impacto (nota metodologica, secao 3.3)
## ---------------------------------------------------------
## Na ordenacao de Cholesky adotada o cambio vem depois do IPCA, entao
## um choque cambial nao afeta a inflacao no mes do impacto. Por isso as
## duas convencoes devem dar resultados proximos. A comparacao abaixo
## verifica isso.

cf_h1 <- decompor(var1, H, shut_impact = FALSE)
cat("\nRobustez, cambio neutralizado so a partir do mes 1:\n")
print(data.frame(horizonte = hh,
                 canal_desde_mes0 = round(cf$canal[hh + 1],    3),
                 canal_desde_mes1 = round(cf_h1$canal[hh + 1], 3)),
      row.names = FALSE)

## ---------------------------------------------------------
## Grafico das duas trajetorias
## ---------------------------------------------------------

pdf(file.path(PASTA_GRAF, "07_contrafactual.pdf"), width = 10, height = 6)
par(mar = c(5, 5, 4, 2))
# folga no topo para a legenda nao cobrir a curva contrafactual
lim <- range(c(cf$total, cf$sem_canal, cf$canal))
lim[2] <- lim[2] + 0.45 * diff(lim)

plot(cf$h, cf$sem_canal, type = "l", lwd = 2.5, lty = 2, col = "grey35",
     ylim = lim,
     main = "Resposta do IPCA a um choque de 1% nas commodities",
     xlab = "meses apos o choque", ylab = "resposta acumulada (%)",
     cex.main = 1.4, cex.axis = 1.2, cex.lab = 1.2)
grid(col = "grey92")
abline(h = 0, col = "red", lty = 3)
lines(cf$h, cf$total, lwd = 2.5, col = "grey10")
lines(cf$h, cf$canal, lwd = 2,   col = "grey60")
legend("topleft",
       c("total (cambio responde)", "sem o canal do cambio",
         "canal do cambio (diferenca)"),
       col = c("grey10", "grey35", "grey60"), lwd = c(2.5, 2.5, 2),
       lty = c(1, 2, 1), bty = "n", cex = 1.1)
invisible(dev.off())

## ---------------------------------------------------------
## Intervalo de confianca por bootstrap
## ---------------------------------------------------------
## A nota metodologica e explicita: a cada replicacao o VAR precisa ser
## reestimado e reidentificado, e o contrafactual refeito. Nao se pode
## manter A e B fixos. Demora alguns minutos.

REPS <- 500
boot_total <- boot_sem <- boot_canal <- matrix(NA_real_, REPS, length(hh))

for (r in 1:REPS) {
  tentativa <- try({
    vb <- VAR(simular(var1), p = P, type = "const")
    decompor(vb, H, shut_impact = TRUE)
  }, silent = TRUE)
  if (!inherits(tentativa, "try-error")) {
    boot_total[r, ] <- tentativa$total[hh + 1]
    boot_sem[r, ]   <- tentativa$sem_canal[hh + 1]
    boot_canal[r, ] <- tentativa$canal[hh + 1]
  }
  if (r %% 50 == 0) cat("bootstrap:", r, "de", REPS, "\n")
}

banda <- function(m) round(t(apply(m, 2, quantile, c(0.05, 0.95), na.rm = TRUE)), 3)

tab_cambio$total_inf    <- banda(boot_total)[, 1]
tab_cambio$total_sup    <- banda(boot_total)[, 2]
tab_cambio$sem_inf      <- banda(boot_sem)[, 1]
tab_cambio$sem_sup      <- banda(boot_sem)[, 2]
tab_cambio$canal_inf    <- banda(boot_canal)[, 1]
tab_cambio$canal_sup    <- banda(boot_canal)[, 2]

cat("\n=== TABELA 4.2 com intervalos de 90% ===\n")
print(tab_cambio, row.names = FALSE)
write.csv(tab_cambio, file.path(PASTA_TABS, "tab4_2_cambio.csv"),
          row.names = FALSE)

## ---------------------------------------------------------
## 6. Decomposicao da variancia (Secao 3.4.3, Tabela 4.3)
## ---------------------------------------------------------

fevd1 <- fevd(var1, n.ahead = H)

## Grafico proprio: participacao de cada choque na variancia do IPCA.
pdf(file.path(PASTA_GRAF, "06_decomposicao_variancia.pdf"), width = 10, height = 7)
par(mar = c(7.5, 5, 4, 2), xpd = TRUE)
cores <- grey.colors(ncol(y), start = 0.15, end = 0.9)
barplot(t(100 * fevd1$ipca[c(3, 6, 12, 24), ]),
        names.arg = paste(c(3, 6, 12, 24), "meses"),
        col = cores, border = "white", ylim = c(0, 100),
        ylab = "% da variancia do IPCA",
        main = "Decomposicao da variancia do erro de previsao do IPCA",
        cex.main = 1.4, cex.axis = 1.2, cex.lab = 1.2, cex.names = 1.2)
legend("bottom", inset = c(0, -0.26), legend = colnames(y),
       fill = cores, border = "white", bty = "n", cex = 1.15, ncol = 4)
invisible(dev.off())

tab_fevd <- round(100 * fevd1$ipca[c(3, 12, 24), ], 1)
rownames(tab_fevd) <- paste0("h=", c(3, 12, 24))

cat("\n=== TABELA 4.3: DECOMPOSICAO DA VARIANCIA DO IPCA (%) ===\n")
print(tab_fevd)
write.csv(tab_fevd, file.path(PASTA_TABS, "tab4_3_fevd.csv"))

## ---------------------------------------------------------
## 7. Subamostras (Secao 3.5.1, Tabela 4.4)
## ---------------------------------------------------------
## O repasse cambial em 12 meses e comparado em duas janelas.

## Mesma definicao de repasse da Tabela 4.1: resposta acumulada do IPCA
## dividida pela variacao acumulada do cambio em doze meses.

repasse12 <- function(dados_ts) {
  v  <- VAR(dados_ts, p = P, type = "const")
  ir <- irf(v, impulse = "cambio", response = c("ipca", "cambio"),
            n.ahead = 12, ortho = TRUE, cumulative = TRUE, boot = FALSE)
  ir$irf$cambio[13, "ipca"] / ir$irf$cambio[13, "cambio"]
}

repasse12_ic <- function(dados_ts, reps = 500) {
  pt <- repasse12(dados_ts)
  v  <- VAR(dados_ts, p = P, type = "const")
  bs <- rep(NA_real_, reps)
  for (r in 1:reps) {
    x <- try(repasse12(simular(v)), silent = TRUE)
    if (!inherits(x, "try-error")) bs[r] <- x
  }
  c(resposta = pt,
    inf = as.numeric(quantile(bs, 0.05, na.rm = TRUE)),
    sup = as.numeric(quantile(bs, 0.95, na.rm = TRUE)))
}

j1 <- window(y, end = c(2013, 12))
j2 <- window(y, start = c(2014, 1))

tab_sub <- rbind(
  data.frame(janela = "2003-2013", n = nrow(j1), t(round(repasse12_ic(j1), 3))),
  data.frame(janela = "2014-2026", n = nrow(j2), t(round(repasse12_ic(j2), 3))))

cat("\n=== TABELA 4.4: REPASSE CAMBIAL POR SUBAMOSTRA (12 meses) ===\n")
print(tab_sub, row.names = FALSE)
cat("Comparar os intervalos nao e um teste formal de igualdade:\n")
cat("intervalos que se sobrepoem apenas nao permitem afirmar a diferenca.\n")
write.csv(tab_sub, file.path(PASTA_TABS, "tab4_4_subamostras.csv"), row.names = FALSE)

## ---------------------------------------------------------
## 8. Assimetria: depreciacoes x apreciacoes (Secao 3.5.2)
## ---------------------------------------------------------
## Separa a variacao cambial em parcela positiva e negativa e testa
## se a soma dos coeficientes e igual nos dois casos.

infl   <- y[, "ipca"]
cambio <- y[, "cambio"]
dep <- pmax(cambio, 0)        # depreciacoes
apr <- pmin(cambio, 0)        # apreciacoes
dep <- ts(dep, start = start(y), frequency = 12)
apr <- ts(apr, start = start(y), frequency = 12)

eq_livre <- dynlm(infl ~ L(dep, 0:6) + L(apr, 0:6) + L(infl, 1:2))
summary(eq_livre)

# H0 de simetria: soma dos coeficientes das depreciacoes igual a
# soma dos coeficientes das apreciacoes.
b <- coef(eq_livre)
soma_dep <- sum(b[2:8])
soma_apr <- sum(b[9:15])

R <- rep(0, length(b)); R[2:8] <- 1; R[9:15] <- -1
V <- vcov(eq_livre)
F_est <- (soma_dep - soma_apr)^2 / as.numeric(t(R) %*% V %*% R)
p_val <- 1 - pf(F_est, 1, df.residual(eq_livre))

tab_assim <- data.frame(
  soma_depreciacoes = round(soma_dep, 3),
  soma_apreciacoes  = round(soma_apr, 3),
  diferenca = round(soma_dep - soma_apr, 3),
  F = round(F_est, 2),
  p_valor = round(p_val, 3),
  rejeita_simetria = ifelse(p_val < 0.05, "sim", "nao"))

cat("\n=== TABELA 4.5: ASSIMETRIA CAMBIAL ===\n")
print(tab_assim, row.names = FALSE)
write.csv(tab_assim, file.path(PASTA_TABS, "tab4_5_assimetria.csv"), row.names = FALSE)

## ---------------------------------------------------------
## 9. Robustez da ordenacao de Cholesky (Secao 4.6)
## ---------------------------------------------------------
## A ordenacao adotada tem justificativa economica, mas e sempre em
## alguma medida arbitraria. O exercicio abaixo reestima o modelo
## mudando a posicao das variaveis e recalcula o repasse de 12 meses.
## Se os numeros pouco se alteram, a leitura nao depende da ordem.

ordenacoes <- list(
  "baseline"                          = c("comm","pimp","vix","hiato","ipca","selic","cambio"),
  "cambio antes do bloco domestico"   = c("comm","pimp","vix","cambio","hiato","ipca","selic"),
  "commodities por ultimo no externo" = c("pimp","vix","comm","hiato","ipca","selic","cambio"),
  "bloco domestico invertido"         = c("comm","pimp","vix","selic","ipca","hiato","cambio"))

linhas_ord <- list()
for (nome in names(ordenacoes)) {
  o  <- ordenacoes[[nome]]
  yo <- y[, o]
  vo <- VAR(yo, p = P, type = "const")
  r12 <- repasse_pt(vo, H)[2, ]          # linha 2 = horizonte de 12 meses
  # canal do cambio na mesma ordenacao
  io <- match(c("comm","cambio","ipca"), o)
  dd <- bgw_channel(Acoef(vo), matriz_B(vo), io[1], io[2], io[3],
                    H = H, compensator = io[2], shut_impact = TRUE)
  esc <- cumsum(dd$resposta_choque)
  canal12 <- (cumsum(dd$canal) / esc)[13]
  linhas_ord[[nome]] <- data.frame(
    ordenacao = nome,
    comm      = round(r12[["comm"]],   3),
    pimp      = round(r12[["pimp"]],   3),
    cambio    = round(r12[["cambio"]], 3),
    canal     = round(canal12,         3))
}
tab_ordem <- do.call(rbind, linhas_ord)

cat("\n=== TABELA 4.6: REPASSE EM 12 MESES SOB ORDENACOES ALTERNATIVAS ===\n")
print(tab_ordem, row.names = FALSE)
write.csv(tab_ordem, file.path(PASTA_TABS, "tab4_7_ordenacao.csv"),
          row.names = FALSE)

cat("\nTabelas gravadas em", PASTA_TABS, "\n")
cat("Graficos gravados em", PASTA_GRAF, "\n")
cat("  02_series.pdf                series do VAR, uma por pagina\n")
cat("  03_cusum.pdf                 teste CUSUM de estabilidade\n")
cat("  04_respostas.pdf             as quatro respostas do Capitulo 4\n")
cat("  05_respostas_completas.pdf   todas as respostas, para o Apendice\n")
cat("  06_decomposicao_variancia.pdf\n")
cat("  07_contrafactual.pdf          decomposicao Sims-Zha do canal cambial\n")
cat("  08_respostas_pontuais.pdf     respostas mensais do IPCA aos tres choques\n")
cat("  09_respostas_acumuladas.pdf   as mesmas respostas, acumuladas\n")
