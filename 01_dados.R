## =========================================================
## 01_dados.R
## Monografia: dinamica inflacionaria brasileira
##
## Monta o painel descrito na Secao 3.1 e roda os testes de
## raiz unitaria da Secao 3.2 (Tabela 3.1).
##
## ENTRADA : painel_completo.csv   (mesma pasta deste script)
## SAIDAS  : painel_var.csv        painel pronto para o 02_var.R
##           tab_adf.csv           Tabela 3.1
## =========================================================

library(seasonal)
library(mFilter)
library(urca)

## ---------------------------------------------------------
## Graficos
## ---------------------------------------------------------
## Todos os graficos sao gravados em graficos/01_series.pdf, um por
## pagina, em tamanho grande. Isso evita o problema de varias series
## espremidas num painel unico, ilegiveis no RStudio.

dir.create("graficos", showWarnings = FALSE)

serie_plot <- function(x, titulo, unidade = "", zero = FALSE) {
  plot(x, type = "l", lwd = 1.6, col = "grey20",
       main = titulo, xlab = "", ylab = unidade,
       cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2)
  if (zero) abline(h = 0, col = "red", lty = 2)
  grid(col = "grey85")
}

## ---------------------------------------------------------
## 1. Leitura e construcao das series como objetos ts
## ---------------------------------------------------------

## O arquivo pode estar em dois formatos. O original usa virgula como
## separador e ponto decimal. Se ele for aberto e salvo pelo Excel em
## portugues, passa a usar ponto e virgula como separador e virgula
## decimal. A funcao abaixo reconhece os dois casos.

ler_painel <- function(arquivo) {
  primeira <- readLines(arquivo, n = 1)
  if (grepl(";", primeira)) read.csv2(arquivo) else read.csv(arquivo)
}

dados <- ler_painel("painel_completo.csv")

## Confere se as colunas necessarias vieram todas. Se o arquivo tiver sido
## salvo em outro formato, o erro aparece aqui, com a causa explicita.
necessarias <- c("ipca", "cambio", "selic", "ibcbr", "icbr_usd", "pimp", "vix")
faltando <- setdiff(necessarias, names(dados))
if (length(faltando) > 0) {
  stop("Nao encontrei as colunas: ", paste(faltando, collapse = ", "),
       "\nForam lidas ", ncol(dados), " colunas. ",
       "Verifique se o painel_completo.csv nao foi salvo pelo Excel ",
       "em formato diferente do original.")
}

cat("Painel lido:", nrow(dados), "meses,", ncol(dados), "colunas\n")

# a amostra comeca em janeiro de 2003 e as linhas ja estao em ordem
inicio <- c(2003, 1)

icbr   <- ts(dados$icbr_usd, start = inicio, frequency = 12)   # IC-Br em US$
pimp   <- ts(dados$pimp,     start = inicio, frequency = 12)   # precos de importacao
cambio <- ts(dados$cambio,   start = inicio, frequency = 12)   # R$/US$
ipca   <- ts(dados$ipca,     start = inicio, frequency = 12)   # IPCA, % a.m.
selic  <- ts(dados$selic,    start = inicio, frequency = 12)   # Selic, % a.a.
ibcbr  <- ts(dados$ibcbr,    start = inicio, frequency = 12)   # IBC-Br dessaz.
vix    <- ts(dados$vix,      start = inicio, frequency = 12)   # VIX

pdf("graficos/01_series.pdf", width = 11, height = 6.5)
par(mar = c(4, 5, 4, 2))

serie_plot(icbr,   "Indice de Commodities Brasil (IC-Br), em US$", "indice")
serie_plot(pimp,   "Indice de precos de importacao, em US$", "indice")
serie_plot(cambio, "Taxa de cambio nominal R$/US$", "R$/US$")
serie_plot(ipca,   "IPCA, variacao mensal", "% a.m.", zero = TRUE)
serie_plot(selic,  "Taxa Selic acumulada no mes", "% a.a.")
serie_plot(ibcbr,  "IBC-Br dessazonalizado", "indice")
serie_plot(vix,    "Indice VIX", "pontos")

## ---------------------------------------------------------
## 2. Ajustamento sazonal do IPCA (Secao 3.1.3)
## ---------------------------------------------------------
## O IBC-Br do SGS 24364 ja e divulgado dessazonalizado.

infl <- final(seas(ipca))

# comparacao entre a serie observada e a dessazonalizada
plot(ipca, type = "l", lwd = 1.4, col = "grey65",
     main = "IPCA: observado e dessazonalizado", xlab = "", ylab = "% a.m.",
     cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2)
lines(infl, lwd = 2, col = "black")
abline(h = 0, col = "red", lty = 2)
grid(col = "grey85")
legend("topleft", c("observado", "dessazonalizado"),
       col = c("grey65", "black"), lwd = c(1.4, 2), bty = "n", cex = 1.2)

## ---------------------------------------------------------
## 3. Hiato do produto pelo filtro HP (Secao 3.1.2)
## ---------------------------------------------------------

hp    <- hpfilter(log(ibcbr), freq = 14400, type = "lambda")
hiato <- 100 * hp$cycle

serie_plot(hiato, "Hiato do produto (filtro HP)", "% do produto potencial",
           zero = TRUE)

## ---------------------------------------------------------
## 4. Transformacoes (Secao 3.1.3)
## ---------------------------------------------------------

dlicbr   <- 100 * diff(log(icbr))     # variacao % mensal
dlpimp   <- 100 * diff(log(pimp))
dlcambio <- 100 * diff(log(cambio))
dselic   <- diff(selic)               # Selic em primeira diferenca
lvix     <- log(vix)                  # VIX em logaritmo

## ---------------------------------------------------------
## 5. Testes de raiz unitaria (Secao 3.2, Tabela 3.1)
## ---------------------------------------------------------
## Series em nivel: constante e tendencia.
## Series transformadas: apenas constante.
## Em ambos os casos, ate 12 defasagens, escolhidas por Schwarz.

adf <- function(x, tipo) {
  teste <- ur.df(na.omit(x), type = tipo, lags = 12, selectlags = "BIC")
  crit  <- if (tipo == "trend") teste@cval[1, "5pct"] else teste@cval[1, "5pct"]
  data.frame(estatistica = round(as.numeric(teste@teststat[1]), 2),
             critico_5pct = crit,
             conclusao = ifelse(as.numeric(teste@teststat[1]) < crit, "I(0)", "I(1)"))
}

tab_adf <- rbind(
  cbind(serie = "ln IC-Br",            adf(log(icbr),   "trend")),
  cbind(serie = "ln P. importacao",    adf(log(pimp),   "trend")),
  cbind(serie = "ln Cambio",           adf(log(cambio), "trend")),
  # A Selic e testada apenas com constante: uma taxa de juros de politica
  # nao tem tendencia deterministica.
  cbind(serie = "Selic",               adf(selic,       "drift")),
  cbind(serie = "D.ln IC-Br",          adf(dlicbr,      "drift")),
  cbind(serie = "D.ln P. importacao",  adf(dlpimp,      "drift")),
  cbind(serie = "D.ln Cambio",         adf(dlcambio,    "drift")),
  cbind(serie = "D.Selic",             adf(dselic,      "drift")),
  cbind(serie = "ln VIX",              adf(lvix,        "drift")),
  cbind(serie = "Hiato",               adf(hiato,       "drift")),
  cbind(serie = "IPCA dessazonalizado",adf(infl,        "drift")))

cat("\n=== TABELA 3.1: TESTES DE RAIZ UNITARIA (ADF) ===\n")
print(tab_adf, row.names = FALSE)
write.csv(tab_adf, "tab_adf.csv", row.names = FALSE)

## ---------------------------------------------------------
## 6. Painel final
## ---------------------------------------------------------
## As series em primeira diferenca perdem a primeira observacao,
## entao o painel comeca em fevereiro de 2003.

painel <- na.omit(cbind(dlicbr, dlpimp, lvix, hiato, infl, dselic, dlcambio))
colnames(painel) <- c("dlicbr", "dlpimp", "lvix", "hiato", "infl", "dselic", "dlcambio")

cat("\nPainel:", start(painel)[1], "/", start(painel)[2],
    "a", end(painel)[1], "/", end(painel)[2],
    "|", nrow(painel), "observacoes\n")

# uma pagina por variavel do VAR, ja transformada
rotulos <- c("Commodities (variacao % mensal)",
             "Precos de importacao (variacao % mensal)",
             "VIX (logaritmo)",
             "Hiato do produto (%)",
             "IPCA dessazonalizado (% a.m.)",
             "Selic (variacao mensal, p.p.)",
             "Cambio (variacao % mensal)")
for (j in 1:ncol(painel))
  serie_plot(painel[, j], rotulos[j], "", zero = TRUE)

dev.off()
cat("\nGraficos gravados em graficos/01_series.pdf\n")

write.csv(data.frame(data = time(painel), painel), "painel_var.csv", row.names = FALSE)
cat("\nArquivos gravados: painel_var.csv e tab_adf.csv\n")