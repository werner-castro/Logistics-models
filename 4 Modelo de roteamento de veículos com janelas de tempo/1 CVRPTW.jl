# problema de roteamento de veículos com janelas de tempo, restrição de capacidade e frota homogênea (CVRPTW)

# formulaçao matemática proposta por: Larsen (1999) e Ombuki et al. (2006)

# Instância retirada do livro: Arenales, do mesmo problema do reteamento de veículos.

# implementação da modelagem: Carlos Werner Ribeiro de Castro

# Parâmetros do modelo:

# n = 1:nc, número de clientes á serem visitados

# N = [0, 1:n, n + 1], em que 0 e n+1 representam o depósito

# E = {(i, j) : i, j ∈ N, i ≠ j, i ≠ n + 1, j ≠ 0}

# nv = número total de veículos

# K = [1:nv], conjunto de veículos utilizados

# Cap = capacidade máxima dos veículos

# D(i) = demanda do cliente i

# C(i,j) = custo de transporte de i para j

# t(i,j) = tempo de viagem de i para j

# e(i) = instante inicial da janela de tempo do cliente i

# l(i) = instante final da janela de tempo do cliente i

# s(i) = tempo de serviço do cliente i

# M = número grande.

# Variaveis do modelo:

# x(i,j,k) = 1 se o veículo k vai passa pelo arco (i,j) ou pelos cliente (i,j)

# b(i,k) = instante de inicio do atendimento do cliente i pelo veículo k

# Modelo matemático:

#            N   N   N
# F.o = Min: ∑   ∑   ∑ c(i,j) * x(i,j,k)
#           i=0 j=0 k=1

# S.t:

#  K   N
#  ∑   ∑ x(i,j,k) = 1                       ∀ i = 1:n
# k=1 j=1

#  N       N
#  ∑ d(i)  ∑ x(i,j,k)   <= Cap              ∀ k = 1:nv
# i=1     j=1


#  N
#  ∑ x(0,j,k) = 1                           ∀ k = 1:nv
# j=0


#  N              N
#  ∑ x(i,h,k)  -  ∑ x(h,j,k) = 0            ∀ k = 1:nv, h = 1:n
# i=0            j=0


#  N
#  ∑ x(i,n+1,k) = 1                         ∀ k = 1:nv
# i=0


# b(i,k) + s(i) + t(i,j) - M * (1 - x(i,j,k)) <= b(j,k)       ∀ i,j = 1:n, k = 1:nv


# e(i) <= b(i,k) <= l(i)                    ∀ i = 1:n, k = 1:nv


# x(i,j,k) = {0,1}                          ∀ i,j = 1:n, k = 1:nv


# b(i,k)   >= 0                             ∀ i = 1:n, k = 1:nv

using JuMP, Cbc, OffsetArrays, Plots, PlotThemes

theme(:juno)

include("src.jl")

cvrpwt = Model(optimizer_with_attributes(Cbc.Optimizer, "threads" => 4, "seconds" => 3600.0))

set_silent(cvrpwt) # omitindo a saída do solver

# Parâmetros:

# número de clientes
n = 9

# coordenadas dos clientes
vx = [50.0, 16.0, 23.0, 40.0, 9.0, 97.0, 78.0, 20.0, 71.0, 64.0, 50.0]
vy = [50.0, 32.0, 1.0, 65.0, 77.0, 71.0, 24.0, 26.0, 98.0, 55.0, 50.0]

# conjunto de clientes com o depósito
N = [0:n...]

# número total de veículos
carros = 3

# conjunto de veículos utilizados
K = [1:carros...]

# capacidade dos veículos
Cap = fill(60.0, length(K))

# demanda dos clientes
D = [11, 35, 2, 9, 3, 18, 8, 10, 11]

# custo de transporte do cliente i para o cliente j
custo = OffsetArray(distances(vx,vy), 0:n+1, 0:n+1)

# tempo de viagem de i para j
t = custo

# instante inicial
itime = OffsetVector([0.0, 45.0, 11.0, 25.0, 20.0, 15.0, 50.0, 10.0, 40.0, 10.0, 0.0], 0:n+1)

# instante final
ftime = OffsetVector([0.0, 70.0, 145.0, 40.0, 100.0, 80.0, 190.0, 110.0, 190.0, 45.0, 400.0], 0:n+1)

# número de clientes
C = [1:n...]

# tempo de serviço no cliente
s = fill(0.0,n)

# número suficientemente grande
M = 500

# vetor de domínio dos índices
E = [(i,j) for i in N, j in N if j > 0 && i ≠ j]

# variáveis
@variable(cvrpwt, x[i in N, j in N, k in K], Bin)

@variable(cvrpwt, b[i in N, k in K] ≥ 0)

# função objetivo
@objective(cvrpwt, Min, sum(x[i,j,k] * custo[i,j] for (i,j) in E, k in K))

# restrições
@constraint(cvrpwt, [i in C], sum(x[i,j,k] for j in N, k in K if i ≠ j) == 1)

@constraint(cvrpwt, [k in K], sum(D[i] * sum(x[i,j,k] for j in N) for i in C) ≤ Cap[k])

@constraint(cvrpwt, [h in C, k in K], sum(x[i,h,k] for i in N) - sum(x[h,j,k] for j in N) == 0)

@constraint(cvrpwt, [k in K], sum(x[0,j,k] for j in N) == 1)

@constraint(cvrpwt, [k in K], sum(x[i,0,k] for i in N) == 1)

@constraint(cvrpwt, [(i,j) in E, k in K], b[i,k] + t[i,j] ≤ b[j,k] + (1 - x[i,j,k]) * M)

@constraint(cvrpwt, [i in N, k in K; i > 0], itime[i] ≤ b[i,k] ≤ ftime[i])

# otimizando
optimize!(cvrpwt)

# gerando o relatório de saída
report(cvrpwt, value.(x), value.(b), itime, ftime, vx, vy, custo)
