# Files-TCC
R scripts worth running
Dinâmica inflacionária brasileira: choques de commodities,
 inflação importada e o papel do câmbio

Monografia de final de curso, PUC-Rio, 2026.
Autor: Enzo De Martino da Rocha
Orientador: Marco Cavalcanti

 Arquivos
- painel_completo.csv : séries mensais, 2003:01 a 2026:05
- 01_dados.R : monta o painel e roda os testes de raiz unitária
- 02_var.R : estima o VAR e produz os resultados

Rode 01_dados.R antes de 02_var.R.
Pacotes: seasonal, mFilter, urca, vars, dynlm, lmtest.
