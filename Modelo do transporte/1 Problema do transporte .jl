# problema do transporte

# Índices
# i = 1...m = número de origens
# j = 1...n = número de destinos/clientes

# Parâmetros
# s(i)      = disponibilidade/capacidade da origem i
# d(j)      = demanda do destino j
# c(i,j)    = custo do transporte de i para j

# Variáveis
# x(i,j)    = quantidade transportada de i para j

# Modelo matemático

# função objetivo
# |
# |         m    n
# | Min z = ∑    ∑ c(i,j) * x(i,j)
# |        i=1  j=1
# |
# restrições
# |
# |  n
# |  ∑  x(i,j) <= oferta(i)             ∀ i = 1..m
# | j=1
# |
# |  m
# |  ∑  x(i,j) = d(j)              ∀ j = 1..n
# | i=1
# |
# |
# | x(i,j) >= 0                    ∀ i = 1..m, j = 1...n

using JuMP, Cbc

prt = Model(Cbc.Optimizer)

# matriz de custos de i para j
c = [
    6  5 8
    13 12 1
     7 9 5
    10 6 4
]

# vetor de demanda
demanda = [8 32 15]

# vetor de oferta
oferta = [10 20 12 13]

# índices
m,n = size(c)

# variáveis
@variable(prt, x[i = 1:m, j = 1:n] >= 0)

# função objetivo
@objective(prt, Min, sum(c[i,j] * x[i,j] for i = 1:m, j = 1:n))

# restrições
@constraint(prt, [i = 1:m], sum(x[i,j] for j = 1:n) <= oferta[i])

@constraint(prt, [j = 1:n], sum(x[i,j] for i = 1:m) == demanda[j])

# resolvendo o modelo
optimize!(prt)

# validando o status da solução e gerando um relatório de saída
if termination_status(prt) == MOI.OPTIMAL
    clearconsole()
    println("Solução ótima encontrada ! \n")
    println("Custos totais com transporte foi = $(objective_value(prt)) reais \n")
    for i = 1:m, j = 1:n
        if value.(x[i,j]) > 0
            println("Foi transportado da origem $(i) a quantidade de $(value.(x[i,j])) unidades para o destino $(j)")
        end
    end
else
    println("Solução não encontrada !")
end
