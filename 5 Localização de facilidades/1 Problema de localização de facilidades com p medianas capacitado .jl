# Problema de localização de facilidades com P-medianas capacitado
# Modelagem: Carlos Werner Ribeiro de Castro
# Modelo retirado do livro; Pesquisa Operacional para os cursos de engenharia, Arenales 2007, pág = 200

# Descrição do problema e modelo;
# A localização de facilidades é um aspecto crítico do planejamento estratégico de empresas pri-
# vadas e públicas. Exemplos típicos no setor público envolvem decisões de localização de centros
# de saúde, escolas e estações de bombeiros, enquanto no setor privado tem-se a localização de fá-
# bricas, armazéns e centros de distribuição. Em diversas situações, tais como em sistemas de dis-
# tribuição, as decisões da localização de facilidades e de designação de clientes a facilidades
# são feitas simultaneamente. A seguir, apresentam-se modelos matemáticos de alguns problemas im-
# portantes de localização. Para tal, considere os seguintes parâmetros;
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
# q(j)   = demanda do cliente j
# Q(i)   = capacidade da facilidade i
#
# Variáveis:
# x(i,j) = se a facilidade i rebece o cliente j
# y(i) = se a facilidade i é aberta
#
# Modelo matemático:
#
# Função objetivo:
#
#  I   J
#  ∑   ∑  c(i,j) * x(i,j)
# i=1 j=1
#
# Sujeito a:
#
#  I
#  ∑  x(i,j) = 1                         ∀ j = 1.....J. Cada cliente j só vai para uma facilidade i
# i=1
#
# x(i,j) <= y(i)                         ∀ i = 1.....I, j = 1.....J. Os clientes são irão para as facilidades abertas
#
#  I
#  ∑  y(i) = p                           ∀ i = 1.....I. O número de facilidades tem que ser igual a quantidade permitida
# i=1
#
# J
# ∑ q(j) * x(i,j) <= Q(i) * y(i)         ∀ i = 1.....I.
# j=1
#
# x(i,j), y(i) ∈ {0,1}      ∀ i = 1.....I, j = 1.....J.


using JuMP, Cbc, Plots, PlotThemes; theme(:juno)

# construindo o modelo
plf = Model(with_optimizer(Cbc.Optimizer, threads = 5, seconds = 3600))

# Índices;

# número de facilidades
I = 10

# número de clientes
J = 40

# número de facilidades que podem ser abertas
p = 4

αc = rand(-100:100,J,1) # coletando as coordenadas dos clientes
βc = rand(-100:100,J,1) # coletando as coordenadas dos clientes
αf = rand(-90:90,I,1)   # coletando as coordenadas das facilidades
βf = rand(-90:90,I,1)   # coletando as coordenadas das facilidades

# custo do cliente i ser alocado na facilidade j
c = fill(0,I,J)
for i in 1:I
   for j in 1:J
      c[i,j] = sqrt((αc[j] - αf[i])^2) + sqrt((βc[j] - βf[i])^2)
   end
end

# demanda do cliente j
q = rand(1:100,J,1)

# capacidade da facilidade i
Q = rand(100:1000,I,1)

# construindo as variáveis
@variable(plf, x[i = 1:I, j = 1:J], Bin)

@variable(plf, y[i = 1:I], Bin)

# função objetivo
@objective(plf, Min, sum(c[i,j] * x[i,j] for i in 1:I, j in 1:J))

# restrições
@constraint(plf,[j in 1:J], sum(x[i,j] for i in 1:I) == 1)

@constraint(plf,[i in 1:I, j in 1:J], x[i,j] <= y[i])

@constraint(plf, sum(y[i] for i in 1:I) == p)

@constraint(plf, [i in 1:I], sum(q[j] * x[i,j] for j in 1:J) <= Q[i] * y[i])

# resolvendo o modelo
optimize!(plf)

# validando o status da solução encontrada
if termination_status(plf) == MOI.OPTIMAL
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
   println("Facilidades que não foram abertas: $(f) \n")
   for i in 1:I
      if sum(value.(x[i,1:end])) > 0.9
         print("A facilidade $(i) recebeu os clientes: ")
         for j in 1:J
            if value.(x[i,j]) > 0.9
               print(" $([j])")
            end
         end
      else
         continue
      end
      println(" ")
   end

   # gerando o gráfico de com as alocações
   ac = Dict{Any,Any}()   # criando um dicionário para ordenação das coordenadas
   for j in 1:J
      ac[j] = (αc[j], βc[j])
   end
   af = Dict{Any,Any}()  # criando um dicionário para ordenação das coordenadas das facilidades
   for i in 1:I
      af[i] = (αf[i], βf[i])
   end
   ac = sort(ac) # ordenando as coordenadas
   b = sortslices(b, dims = 1)
   V = collect(1:J)

   Plots.scatter(αc, βc, color = "black", legend = false, size = (800, 600))
   annotate!(αc[1:end] .+ 2, βc[1:end] .+ 1, V[1:end], font(7, :white)) # numeração dos clientes

   # plotando todas as facilidades
   Plots.scatter!(
      αf, βf,
      marker = (:square, 6),
      color = "black"
   )

   # identificando das facilidades
   annotate!(αf[1:end], βf[1:end], 1:I, font(7, :white))

   function gloc()
      plt = []
      for j in 1:J
         plt = Plots.scatter!([af[b[j,1]], ac[b[j,2]]], color = b[j,1])
         plt = plot!([af[b[j,1]], ac[b[j,2]]],
         title = "Localização de facilidades capacitado com p-medianas ",
         xlabel = "Modelo com $(I) candidatos, $(J) clientes e $(p) facilidades demandadas",
         color = b[j,1])
         # display(plt)
      end
      return plt
   end

   Plots.scatter!(
      αf[b[1:end,1]], βf[b[1:end,1]],
      legend = false,
      marker=(:square,6),
      color = b[1:end,1]
   )

   # plotando o gráfico de localizações
   gloc()
   Plots.scatter!(
      αf[b[1:end,1]], βf[b[1:end,1]],
      legend = false,
      marker=(:square,6),
      color = b[1:end,1]
   )
else
   println("Solução ótima não encontrada !")
end
