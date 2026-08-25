

pacotes <- c("seasonal", "mFilter", "urca")
faltando <- pacotes[!pacotes %in% rownames(installed.packages())]
if (length(faltando) > 0) install.packages(faltando, repos = "https://cloud.r-project.org")

definir_pasta <- function() {
  
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    caminho <- try(rstudioapi::getSourceEditorContext()$path, silent = TRUE)
    if (!inherits(caminho, "try-error") && nzchar(caminho)) {
      setwd(dirname(caminho)); return(invisible(TRUE))
    }
  }
 
  args <- commandArgs(trailingOnly = FALSE)
  arq  <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(arq) > 0) { setwd(dirname(normalizePath(arq))); return(invisible(TRUE)) }
  invisible(FALSE)
}
try(definir_pasta(), silent = TRUE)
cat("Pasta de trabalho:", getwd(), "\n")

library(seasonal)
library(mFilter)
library(urca)


achar <- function(nomes, pastas = c(".", "Paineis", "../Paineis", "..")) {
  for (n in nomes) for (p in pastas) {
    caminho <- file.path(p, n)
    if (file.exists(caminho)) return(caminho)
  }
  stop("Nao encontrei: ", paste(nomes, collapse = " nem "),
       "\nProcurei em: ", paste(pastas, collapse = ", "),
       "\nPasta de trabalho atual: ", getwd())
}

pasta_saida <- function(nome) {
  for (p in c(nome, file.path("..", nome))) if (dir.exists(p)) return(p)
  dir.create(nome, showWarnings = FALSE)
  nome
}

PASTA_PAINEIS <- pasta_saida("Paineis")
PASTA_TABS    <- pasta_saida("Tabs")
PASTA_GRAF    <- pasta_saida("graficos")



serie_plot <- function(x, titulo, unidade = "", zero = FALSE) {
  plot(x, type = "l", lwd = 1.6, col = "grey20",
       main = titulo, xlab = "", ylab = unidade,
       cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2)
  if (zero) abline(h = 0, col = "red", lty = 2)
  grid(col = "grey85")
}



ler_painel <- function(arquivo) {
  primeira <- readLines(arquivo, n = 1)
  if (grepl(";", primeira)) read.csv2(arquivo) else read.csv(arquivo)
}

dados <- ler_painel(achar(c("painel_completo.csv", "painel_completo_1.csv")))


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

pdf(file.path(PASTA_GRAF, "01_series.pdf"), width = 11, height = 6.5)
par(mar = c(4, 5, 4, 2))

serie_plot(icbr,   "Indice de Commodities Brasil (IC-Br), em US$", "indice")
serie_plot(pimp,   "Indice de precos de importacao, em US$", "indice")
serie_plot(cambio, "Taxa de cambio nominal R$/US$", "R$/US$")
serie_plot(ipca,   "IPCA, variacao mensal", "% a.m.", zero = TRUE)
serie_plot(selic,  "Taxa Selic acumulada no mes", "% a.a.")
serie_plot(ibcbr,  "IBC-Br dessazonalizado", "indice")
serie_plot(vix,    "Indice VIX", "pontos")



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



hp    <- hpfilter(log(ibcbr), freq = 14400, type = "lambda")
hiato <- 100 * hp$cycle

serie_plot(hiato, "Hiato do produto (filtro HP)", "% do produto potencial",
           zero = TRUE)


dlicbr   <- 100 * diff(log(icbr))     # variacao % mensal
dlpimp   <- 100 * diff(log(pimp))
dlcambio <- 100 * diff(log(cambio))
dselic   <- diff(selic)               # Selic em primeira diferenca
lvix     <- log(vix)                  # VIX em logaritmo



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
write.csv(tab_adf, file.path(PASTA_TABS, "tab_adf.csv"), row.names = FALSE)



painel <- na.omit(cbind(dlicbr, dlpimp, lvix, hiato, infl, dselic, dlcambio))
colnames(painel) <- c("dlicbr", "dlpimp", "lvix", "hiato", "infl", "dselic", "dlcambio")

cat("\nPainel:", start(painel)[1], "/", start(painel)[2],
    "a", end(painel)[1], "/", end(painel)[2],
    "|", nrow(painel), "observacoes\n")


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
cat("\nGraficos gravados em", file.path(PASTA_GRAF, "01_series.pdf"), "\n")

write.csv(data.frame(data = time(painel), painel),
          file.path(PASTA_PAINEIS, "painel_var.csv"), row.names = FALSE)
cat("\nGravados:", file.path(PASTA_PAINEIS, "painel_var.csv"),
    "e", file.path(PASTA_TABS, "tab_adf.csv"), "\n")
