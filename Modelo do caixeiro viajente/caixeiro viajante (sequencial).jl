# Título do modelo: problema do caixeiro viajante (modelo sequencial).

# Autor da modelagem: Carlos Werner Ribeiro de Castro.

# Descrição do modelo
#
# Problemas de caixeiro-viajante (CV) envolvem um conjunto de cidades, em que o caixeiro sai de
# uma cidade base ou um depósito, visita todas as cidades ou um subconjunto delas somente uma
# vez, e retorna à cidade base de modo a otimizar um ou mais objetivos.
#
# Iniciamos pelo problema clássico e mais importante, objeto de grande número de trabalhos na
# literatura, que trata da minimização da distância da rota percorrida. Considere um grafo não
# orientado G = (N, E), em que o conjunto N consiste de n cidades e E representa o conjunto de
# arestas entre cidades. Suponha que G é um grafo completo, isto é, para qualquer par de cidades
# i, j ∈ N , i ≠ j , existe uma aresta (i,j). A distância entre as cidades i e j é c(i,j),
# e quando c(i,j) = c(j,i), o problema é dito simétrico; caso contrário, é chamado de assimétrico.
# Um caixeiro deve visitar n cidades, passando por cada cidade somente uma vez, e retornar à
# cidade de partida. Este percurso é denominado ciclo hamiltoniano do grafo G, e o problema
# consiste em determinar o ciclo hamiltoniano, ou rota, de distância mínima. Este é um dos
# problemas combinatórios mais conhecidos e pesquisados devido à sua aplicação em diversas áreas,
# tais como manufatura de circuitos, programação da produção, telecomunicações.

# Índices
# | n = número total de cidades para visitar
# | i = 1..n
# | j = 1..n

# Parâmetros
# | c(i,j) = custos totais de transporte da cidade i para cidade j

# Variáveis
# | x(i,j) = 1 se a rota passa da cidade i para cidade j e 0 caso contrário
# | u(i) = variável auxiliar

# Modelo matemático
#
# Função objetivo
# |          n   n
# | Min z =  ∑   ∑   cij * xij      minimizar os custos totais de transporte ou distância
# |         i=1 j=1
#
# Restrições
# |  n
# |  ∑ xij = 1      ∀ j ∈ N         garante que o vendedor só chega a cada cidade j uma única vez
# |i = 1,
# |i ≠ j
# |
# |  n
# |  ∑ xij = 1      ∀ i ∈ N         garante que cada vendedor só sai de cada cidade i uma única vez
# |j = 1,
# |i ≠ j
# |
# | ui - uj + n * xij ⪯ n - 1       restrição para bloqueio de sub-rotas.
# |i > 1,
# |j > 1,
# |i ≠ j
# |
# | u(i) são variáveis auxiliares
# |
# | x(i,j) ∈ {0,1}   (i = 1....n, j = 1...n, i ≠ j)
# |
# | u(i) >= 0 (i = 1....n)

using JuMP, Cbc

# importando o arquivo com o módulo para cálculo da matriz de distâncias, rota e gráfico
include("src.jl")

# pcv = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_THREADS" => 4, "CPX_PARAM_TILIM" => 800))
pcv = Model(optimizer_with_attributes(Cbc.Optimizer, "threads" => 4, "seconds" => 800))

# índices
# # número de cidades/clientes
n = 24

# parâmetros
# gerando coordenadas das cidades
α = rand(-40.0:40.0, n) # coordenadas das cidades em x

β = rand(-40.0:40.0, n) # coordenadas das cidades em y

# α = [10.0 , 20.0 , 50.0 , 70.0 , 90.0]
# #
# β = [30.0 , 50.0 , 90.0 , 30.0 ,50.0]
# #
# n = length(α)

# calculando a matriz de distâncias
c = distances(α,β)

# variáveis
@variable(pcv, x[i = 1:n, j = 1:n], Bin)

@variable(pcv, u[i = 1:n] >= 0)

# função objetivo
@objective(pcv, Min, sum(c[i,j] * x[i,j] for i in 1:n, j in 1:n))

# restrições
@constraint(pcv, [j in 1:n], sum(x[i,j] for i in 1:n if i != j) == 1)

@constraint(pcv, [i in 1:n], sum(x[i,j] for j in 1:n if i != j) == 1)

@constraint(pcv, u[1] == 1)

@constraint(pcv, [i = 1:n, j = 2:n; i != j], u[j] >= u[i] + 1 - n * (1-x[i,j]))

# resolvendo o modelo
optimize!(pcv)

# validando o status da solução e gerando o relatório de saída
outputmodel(pcv,α,β)
