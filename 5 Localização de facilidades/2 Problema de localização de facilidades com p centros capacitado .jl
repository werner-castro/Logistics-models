# Problema de localização de facilidades com P-centros capacitado
# Modelagem em julia: Carlos Werner Ribeiro de Castro
# Modelo retirado do livro; Pesquisa Operacional para os cursos de engenharia, Arenales 2007, pág = 201

# Descrição do problema e modelo;
# Este problema envolve a localização de p facilidades e a designação de clientes a facilidades, de modo
# a minimizar a distância máxima de clientes a facilidades. Este problema admite variações do mo-
# delo básico. O problema de p-centros-nós restringe os nós de facilidades aos nós de clientes, en-
# quanto o problema de p-centros-absolutos permite que os nós de facilidades estejam em qualquer
# lugar dos arcos que ligam nós de clientes.
# Para formular este problema, considere as variáveis do problema de p-medianas e a seguinte
# variável adicional:
#
# Indíces:
# I = número de candidatos a facilidades
# J = numero de clientes
#
# Parâmetros:
# p = número de facilidades que poderão ser abertas
# i = 1........I número de facilidades
# j = 1........J número de clientes
# c(i,j) = custo do cliente i ser alocado na facilidade j
# d(i,j) = distância entre o clente j e a facilidade i
# q(j)   = demanda do cliente j
# Q(i)   = capacidade da facilidade i
# Variáveis:
# x(i,j) = se a facilidade i recebe o cliente j
# y(i) = 1 se a facilidade i é aberta e 0 caso contrário
# r = distância máxima de cliente quando designado a uma facilidade
#
# Modelo matemático:
#
# Função objetivo:
#
# Min z = r
#
# Sujeito a:
#  I
#  ∑  x(i,j) = 1                          ∀ j = 1.....J. Cada cliente j só vai para uma facilidade
# i=1
#
# x(i,j) <= y(i)                          ∀ i = 1.....I, j = 1.....J. Os clientes são irão para as facilidades abertas
#
#  I
#  ∑  y(i) = p                            ∀ i = 1.....I. O número de facilidades tem que ser menor ou igual a quantidade permitida
# i=1
#
# |     I
# |r >= ∑ d(i,j) * x(i,j)                 ∀ j = 1......J. Minimiza a distância máxima do cliente a facilidade
# |    i=1
#
#  J
#  ∑ q(j) * x(i,j) <= Q(i) * y(i)         ∀ i = 1.....I
# j=1
#
# x(i,j), y(i) ∈ {0,1}                    ∀ i = 1.....I, j = 1.....J

using JuMP, Cbc, Plots, PlotThemes; theme(:juno)

# construindo o modelo
plf = Model(with_optimizer(Cbc.Optimizer, threads = 5, seconds = 3600.0))

# Índices;
# número de facilidades
I = 5

# número de clientes
J = 20

# número de facilidades que podem ser abertas
p = I

# distância do cliente j para a facilidade i
αc = rand(-10:100,J,1) # coletando as coordenadas dos clientes e das facilidades
βc = rand(-10:100,J,1) # coletando as coordenadas dos clientes e das facilidades

αf = rand(-100:100,I,1) # coletando as coordenadas das facilidades
βf = rand(-100:100,I,1) # coletando as coordenadas das facilidades

# calculando a distância heuclidiana
d = fill(0,I,J)
for i in 1:I
   for j in 1:J
      d[i,j] = sqrt((αc[j] - αf[i])^2) + sqrt((βc[j] - βf[i])^2)
   end
end

# custo fixo por unidade métrica percorrida do cliente j ser alocado na facilidade i
c = 5.67 .* d

# demanda do cliente j
q = rand(1:100,J,1)

# capacidade da facilidade i
Q = rand(100:1000,I,1)

# construindo as variáveis
@variable(plf, x[i = 1:I, j = 1:J], Bin)

@variable(plf, γ[i = 1:I], Bin)

@variable(plf, r >= 0)

# função objetivo
@objective(plf, Min, r)

# restrições
@constraint(plf, [j in 1:J], sum(x[i,j] for i in 1:I) == 1)

@constraint(plf, [i in 1:I, j in 1:J], x[i,j] <= γ[i])

@constraint(plf, sum(γ[i] for i in 1:I) == p)

@constraint(plf, [j in 1:J], r >= sum(d[i,j] * x[i,j] for i in 1:I))

@constraint(plf, [i in 1:I], sum(q[j] * x[i,j] for j in 1:J) <= Q[i] * γ[i])

# resolvendo o modelo
optimize!(plf)

# validando o status da solução encontrada
if termination_status(plf) == MOI.OPTIMAL
   obj = round.(Int, objective_value(plf))
   b = fill(0,J,2)
   for i in 1:I, j in 1:J
      if value.(x[i,j]) > 0.99
         b[j,:] = [i,j]
      end
   end

   # gerando o relatório de saída
   print("\nSolução ótima encontrada !\n\n")
   println("Custos totais de alocação = $(objective_value(plf)) u.m")
   println(" ")
   f = [i for i in 1:I if sum(value.(x[i,1:end])) == 0]
   if f == []
      println("Todas as facilidades foram abertas")
   else
      println("Facilidades que não foram abertas: $(f) \n")
   end
   for i in 1:I
      if sum(value.(x[i,1:end])) > 0
         print("A facilidade $(i) recebeu os clientes:")
         for j in 1:J
            if value.(x[i,j]) > 0
               print(" $([j])")
            end
         end
      else
         continue
      end
      println(" ")
   end

   # gerando o gráfico de com as alocações
   ac = Dict{Any,Any}() # criando um dicionário para ordenação das coordenadas
   for j in 1:J
      ac[j] = (αc[j], βc[j])
   end
   af = Dict{Any,Any}() # criando um dicionário para ordenação das coordenadas das facilidades
   for i in 1:I
      af[i] = (αf[i],βf[i])
   end
   ac = sort(ac) # ordenando as coordenadas
   b = sortslices(b, dims = 1)
   V = collect(1:J)
   # plotando os clientes
   Plots.scatter(αc,βc,label = "Clientes", color = "black")
   # annotate!(αc[1:end],βc[1:end],V[1:end], font(7, :white))  # numeração dos clientes

   # plotando todas as facilidades
   Plots.scatter!(αf, βf, legend = false,
      marker=(:square, 6),
      color = "black"
   )
   Plots.scatter!(αf[b[1:end,1]], βf[b[1:end,1]],
      marker=(:square,6),
      color = b[1:end,1]
   )
   annotate!(αf[1:end], βf[1:end],1:I, font(7, :white))  # numeração das facilidades
   function gloc()
      plt = []
      for j in 1:J
         plt = Plots.scatter!([af[b[j,1]], ac[b[j,2]]], color = b[j,1])
         plt = plot!([af[b[j,1]], ac[b[j,2]]],
            title = "Localização de facilidades com p-centros capacitado",
            color = b[j,1]
         )
         display(plt)
      end
      # return plt
   end
   # plotando o gráfico de localizações
   gloc()

   # plotando as facilidades abertas
   Plots.scatter!(αf[b[1:end,1]], βf[b[1:end,1]],
      marker=(:square,6),
      color = b[1:end,1]
   )
else
   println(" ")
   println("       Modelo executado !       \n")
   println("==================================")
   println("| Solução ótima não encontrada ! |")
   println("==================================")
   println(" ")
end
