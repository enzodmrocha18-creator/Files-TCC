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
## 3. Causalidade de Granger (Secao 3.4.3)
## ---------------------------------------------------------
## H0: a variavel indicada nao causa, no sentido de Granger,
## as demais variaveis do sistema.

causality(var1, cause = "comm")
causality(var1, cause = "pimp")
causality(var1, cause = "cambio")

## ---------------------------------------------------------
## 4. Funcoes de resposta a impulso, FRI (Secao 3.4.1, Tabela 4.1)
## ---------------------------------------------------------
## irf() com ortho = TRUE usa a decomposicao de Cholesky, e portanto
## depende da ordenacao das variaveis definida acima.
##
## O irf() devolve a resposta a um choque de um desvio padrao.
## Para ler o resultado como coeficiente de repasse, divide-se pela
## resposta da propria variavel de origem no momento do choque, o que
## equivale a normalizar o choque para 1%.

H <- 24

irf_comm <- irf(var1, impulse = "comm", n.ahead = H, ortho = TRUE,
                cumulative = TRUE, boot = TRUE, runs = 1000, ci = 0.90)
irf_pimp <- irf(var1, impulse = "pimp", n.ahead = H, ortho = TRUE,
                cumulative = TRUE, boot = TRUE, runs = 1000, ci = 0.90)
irf_camb <- irf(var1, impulse = "cambio", n.ahead = H, ortho = TRUE,
                cumulative = TRUE, boot = TRUE, runs = 1000, ci = 0.90)

# fatores de normalizacao: impacto do choque na propria variavel
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

## Confere a estrutura devolvida pelo irf() antes de montar a tabela.
cat("Estrutura da resposta a impulso:",
    paste(dim(irf_comm$irf$comm), collapse = " linhas x "), "colunas\n")
cat("Colunas:", paste(colnames(irf_comm$irf$comm), collapse = ", "), "\n")
cat("Bandas de confianca disponiveis:",
    !is.null(irf_comm$Lower), "\n\n")

linha <- function(objeto, impulso, escala, nome) {
  h <- c(7, 13, 25)                      # horizontes 6, 12 e 24 meses
  m  <- objeto$irf[[impulso]][h, "ipca"] / escala
  # Se por algum motivo as bandas nao vierem, a tabela sai so com o
  # ponto estimado, em vez de o script parar.
  lo <- if (!is.null(objeto$Lower)) objeto$Lower[[impulso]][h, "ipca"] / escala else NA_real_
  hi <- if (!is.null(objeto$Upper)) objeto$Upper[[impulso]][h, "ipca"] / escala else NA_real_
  data.frame(choque = nome, horizonte = c(6, 12, 24),
             resposta = round(m, 3),
             ic_inf   = round(lo, 3),
             ic_sup   = round(hi, 3))
}

tab_irf <- try(rbind(
  linha(irf_comm, "comm",   esc_comm, "Commodities (1%)"),
  linha(irf_pimp, "pimp",   esc_pimp, "P. importacao (1%)"),
  linha(irf_camb, "cambio", esc_camb, "Cambio (1%)")), silent = TRUE)

if (inherits(tab_irf, "try-error"))
  stop("Falhou ao montar a Tabela 4.1. Motivo:\n", tab_irf,
       "\nConfira a saida de estrutura impressa logo acima.")

cat("\n=== TABELA 4.1: RESPOSTA ACUMULADA DO IPCA (%) ===\n")
print(tab_irf, row.names = FALSE)
write.csv(tab_irf, file.path(PASTA_TABS, "tab4_1_irf.csv"), row.names = FALSE)

## ---------------------------------------------------------
## 5. O papel do cambio (Secao 3.4.2, Tabela 4.2)
## ---------------------------------------------------------
## A_n = - (resposta do cambio ao choque de commodities)
##         x (resposta do IPCA ao choque cambial)
##
## Como a alta das commodities tende a apreciar o real, a primeira
## parcela e negativa e o sinal invertido torna A_n uma medida
## positiva de amortecimento.

h <- c(7, 13, 25)
theta <- irf_comm$irf$comm[h, "cambio"] / esc_comm     # resposta do cambio
beta  <- irf_camb$irf$cambio[h, "ipca"] / esc_camb     # repasse cambial
An    <- -theta * beta

tab_cambio <- data.frame(
  horizonte = c(6, 12, 24),
  resposta_cambio = round(theta, 3),
  repasse_cambial = round(beta, 3),
  amortecimento = round(An, 3),
  resposta_ipca_observada = round(irf_comm$irf$comm[h, "ipca"] / esc_comm, 3))

cat("\n=== TABELA 4.2: O PAPEL DO CAMBIO ===\n")
print(tab_cambio, row.names = FALSE)
cat("resposta_cambio negativa indica apreciacao do real.\n")
write.csv(tab_cambio, file.path(PASTA_TABS, "tab4_2_cambio.csv"), row.names = FALSE)

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

repasse12 <- function(dados_ts) {
  v <- VAR(dados_ts, p = P, type = "const")
  # Uma unica chamada, pedindo as duas respostas de que precisamos.
  # A primeira linha da resposta do cambio ao proprio choque e o impacto,
  # usado para normalizar o choque para 1%.
  ir <- irf(v, impulse = "cambio", response = c("ipca", "cambio"),
            n.ahead = 12, ortho = TRUE, cumulative = TRUE,
            boot = TRUE, runs = 500, ci = 0.90)
  esc <- ir$irf$cambio[1, "cambio"]
  c(resposta = ir$irf$cambio[13, "ipca"] / esc,
    inf = ir$Lower$cambio[13, "ipca"] / esc,
    sup = ir$Upper$cambio[13, "ipca"] / esc)
}

j1 <- window(y, end = c(2013, 12))
j2 <- window(y, start = c(2014, 1))

tab_sub <- rbind(
  data.frame(janela = "2003-2013", n = nrow(j1), t(round(repasse12(j1), 3))),
  data.frame(janela = "2014-2026", n = nrow(j2), t(round(repasse12(j2), 3))))

cat("\n=== TABELA 4.4: REPASSE CAMBIAL POR SUBAMOSTRA (12 meses) ===\n")
print(tab_sub, row.names = FALSE)
cat("So ha evidencia de mudanca se os intervalos nao se sobrepuserem.\n")
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

cat("\nTabelas gravadas em", PASTA_TABS, "\n")
cat("Graficos gravados em", PASTA_GRAF, "\n")
cat("  02_series.pdf                series do VAR, uma por pagina\n")
cat("  03_cusum.pdf                 teste CUSUM de estabilidade\n")
cat("  04_respostas.pdf             as quatro respostas do Capitulo 4\n")
cat("  05_respostas_completas.pdf   todas as respostas, para o Apendice\n")
cat("  06_decomposicao_variancia.pdf\n")
