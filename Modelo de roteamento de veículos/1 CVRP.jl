
# Roteamento de veículos com retrição de capacidade (CVRP)

# Modelagem: Carlos Werner Ribeiro de Castro

# Parâmetros:
#
# 1      = depot
# m      = number of vehicles
# n      = number of clients
# N      = set of clients, with N = {1...n}
# V      = set of vertices {1,N+2}
# c(i,j) = cost of travel over arc (i,j) ∈ V
# Q      = capacity of vehicle
# D(j)   = demand of client j

# variáveis
# x(i,j) = 1 se o i vai para j e 0 caso contrário
# u(i) = variável auxiliar

# Modelo matemático
#
# | Min z =   ∑ c(i,j) * x(i,j)
# |         (i,j) ∈ V
# |
# | sujeito a:
# |
# |  ∑  x(i,j) = 1                         i ∈ V
# |j ∈ V
# |
# |  ∑  x(i,j) = 1                         j ∈ V
# |i ∈ V
# |
# | u(i) - u(j) + c(i,j) < C - D(j)        i,j ∈ V, i ≠ j, such that  D(i) + D(j) <= C
# |
# | D(i) <= u(i) <= Q                      i ∈ V

using JuMP, Cbc, Plots, DelimitedFiles, DataStructures, PlotThemes

theme(:juno)

include("src.jl")

# criando o modelo
cvrp = Model(optimizer_with_attributes(Cbc.Optimizer, "threads" => 5, "seconds" => 1200))
# cvrp = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_THREADS" => 5, "CPX_PARAM_TILIM" => 600))

# omitindo a saída do solver
# set_silent(cvrp)

# n = 8                                          # número de clientes
# m = 4                                          # número máximo de rotas/carros
# D = rand(10:20,1,n+1)                          # demanda dos clientes
# D[1] = 0                                       # demanda fictícia do depósito
# Q = 10 * m                                     # capacidade total
# V = [1:n+1...]                                 # vetor dos clientes com o depósito

D = [0.0, 11.0, 35.0, 2.0, 9.0, 3.0, 18.0, 8.0, 10.0, 11.0]
V = [1:10...]
Q = 60
m = 3
n = 9

# coordenadas do depósito
cdx = 50.0
cdy = 50.0
# cdx = rand(-100.0:100.0)
# cdy = rand(-100.0:100.0)

# coordenadas dos clientes
# cx = rand(-100.0:100.0,n,1)
# cy = rand(-100.0:100.0,n,1)

# coordenadas dos clientes
cx = [16.0, 23.0, 40.0, 9.0, 97.0, 78.0, 20.0, 71.0, 64.0]
cy = [32.0, 1.0, 65.0, 77.0, 71.0, 24.0, 26.0, 98.0, 55.0]

deposito = [cdx cdy]
clientes = [cx cy]

# adicionando todas as coordenadas
coord = [deposito; clientes; deposito]
coordx = coord[:,1]
coordy = coord[:,2]

c = distances(coordx, coordy)

# variáveis
@variable(cvrp, x[i in V, j in V], Bin)

@variable(cvrp, u[i in V] ≥ 0, Int)

# função objetivo
@objective(cvrp, Min, sum(x[i,j] * c[i,j] for i in V, j in V))

# restrições
@constraint(cvrp, [j in V; j > 1], sum(x[i,j] for i ∈ V if i ≠ j) == 1)

@constraint(cvrp, [i in V; i > 1], sum(x[i,j] for j ∈ V if i ≠ j) == 1)

@constraint(cvrp, [i in V; i > 1], u[i] ≤ Q)

@constraint(cvrp, [i in V; i > 1], u[i] ≥ D[i])

@constraint(cvrp, [i in V, j in V; i > 1 && i ≠ j && D[i] + D[j] ≤ Q], u[i] - u[j] + Q * x[i,j] ≤ Q - D[j])

@constraint(cvrp, [j in V], sum(x[1,j]) ≤ m)

@constraint(cvrp, [i in V], sum(x[i,1]) ≤ m)

# resolvendo o modelo
optimize!(cvrp)

# validando o status da solução e gerando o relatório de saída
if termination_status(cvrp) == MOI.OPTIMAL
    println("Solução ótima encontrada !")

    obj = objective_value(cvrp)

    rota = tour(value.(x))

    # relatorio de saída
    relatorio(cvrp, obj, rota, m, n)

    # gráfico das rotas
    pltour(coordx, coordy, rota, animated = 1)
elseif termination_status(cvrp) == MOI.TIME_LIMIT && has_values(cvrp)
    println("Solução sub-ótima encontrada !")

    obj = objective_value(cvrp)

    rota = tour(value.(x))

    # relatorio de saída
    relatorio(cvrp, obj, rota, m, n)

    # gráfico das rotas
    pltour(coordx, coordy, rota, animated = 1)
else
    println("Solução não encontrada !")
end
