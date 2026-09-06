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

## Conferencia da escala: o valor abaixo tem que ser igual ao esc_comm
## usado na Tabela 4.1.
cat("\nEscala do choque de commodities\n")
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
  esc <- B[i_comm, i_comm]          # impacto do choque nas proprias commodities
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
plot(cf$h, cf$sem_canal, type = "l", lwd = 2.5, lty = 2, col = "grey35",
     ylim = range(c(cf$total, cf$sem_canal, cf$canal)),
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
