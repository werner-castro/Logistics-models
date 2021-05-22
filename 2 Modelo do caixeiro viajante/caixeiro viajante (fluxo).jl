# Título do modelo: problema do caixeiro viajante (modelo de fluxo).

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
# |
# |
# | x(i,j) ∈ {0,1}   (i = 1....n, j = 1...n, i ≠ j)
# |

using JuMP, Cbc

include("src.jl")

# pcv = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_THREADS" => 4, "CPX_PARAM_TILIM" => 800))
pcv = Model(optimizer_with_attributes(Cbc.Optimizer, "threads" => 5, "seconds" => 800))

# índices
# número de cidades/clientes
n = 24

# # # parâmetros
# # gerando coordenadas das cidades
α = rand(1.0:45.0,n) # coordenadas das cidades em x

β = rand(1.0:45.0,n) # coordenadas das cidades em y

# α = [10.0 , 3.0 , 55.0 , 73.0 , 20.0]
# #
# β = [28.0 , 40.0 , 16.0 , 18.0 ,33.0]
# #
# n = length(α)

# calculando a matriz de distâncias
c = distances(α,β)

# variáveis
@variable(pcv, x[i = 1:n, j = 1:n], Bin)

@variable(pcv, y[i = 1:n, j = 1:n] >= 0)

# função objetivo
@objective(pcv, Min, sum(c[i,j] * x[i,j] for i in 1:n, j in 1:n))

# restrições
@constraint(pcv, [j in 1:n], sum(x[i,j] for i in 1:n if i != j) == 1)

@constraint(pcv, [i in 1:n], sum(x[i,j] for j in 1:n if i != j) == 1)

@constraint(pcv, [i = 1:n, j = 1:n; i != j], y[i,j] <= (n - 1) * x[i,j])

@constraint(pcv, sum(y[1,j] for j = 2:n) == n - 1)

@constraint(pcv, [j = 2:n], sum(y[i,j] for i = 1:n if i ≠ j) - sum(y[j,i] for i = 1:n if j ≠ i) == 1)

# resolvendo o modelo
optimize!(pcv)

# validando o status da solução e gerando o relatório de saída
outputmodel(pcv, α, β)
