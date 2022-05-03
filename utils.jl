
using JuMP, Cbc, Ipopt, LightGraphs, Combinatorics, LinearAlgebra

#====================================================================================================================#
# ESTRUTURAS DE DADOS                                                                                                # 
#====================================================================================================================#

abstract type Data end        # tipo abstrato 

# estrutura para entrada de dados
struct ModelData <: Data
    n::Int64                  # número de clientes 
    coord_x::Vector{Float64}  # latitude dos clientes
    coord_y::Vector{Float64}  # longitude dos clientes
    demand::Vector{Float64}   # demanda dos clientes
    capacity::Int64           # capacidade dos veículos
    cars::Int64               # número de carros
    distance::Matrix{Float64} # matriz de distâncias
    obj::Int64                # valor da solução objetivo pelo método exato
end

# estrutura para saída dos dados
struct ReportData <: Data
    custo::Float64               # custo da roteirização
    nome::String                 # nome da heurística utilizada
    rotas::Vector{Vector{Int64}} # rotas geradas pela heurística
    imp::Bool                    # melhoria nas rotas geradas
end

struct SolverParameters <: Data
    gap::Float64
    timelimit::Float64
    threads::Int64
end

# método que recebe um diretório e retorna os dados para o problema do cvrp
function readfile(path::String, index_file::Int)
    file = readdir(path)[index_file]
    instance = CSV.File(path * file) |> DataFrame
    coord_init = 7                                    # número de linhas do cabeçalho
    n = (parse(Int64, instance[3, 2]))
    coord_end = coord_init + (n-1)
    dem_init  = coord_init + n+1
    dem_end   = dem_init + n-1
    cx   = [parse(Int64, split.(instance[i, 1], " ")[3]) for i = coord_init:coord_end] 
    cy   = [parse(Int64, split.(instance[i, 1], " ")[4]) for i = coord_init:coord_end] 
    dem  = [parse(Int64, split.(instance[i, 1], " ")[2]) for i = dem_init:dem_end] 
    cap  = parse(Int, split(instance[5,2], " ")[2])
    car  = parse(Int, split(instance[1,3], ",")[1])
    obj  = parse(Int, split(instance[1,4], ")")[1])
    return n, cx, cy, dem, cap, car, obj
end

# método para geração aleatória de instancias
# t = é o número de vezes que a maior demanda gerada é multiplicada (capacidade de cada veículo)
function createfile(num_clients::Int, t::Int)
    cx  = rand(-100:100, num_clients)
    cy  = rand(-100:100, num_clients)
    dem = rand(1:10, num_clients)
    dem[1] = 0 # demanda do depósito
    cap = maximum(dem) * t
    car = round(Int, sum(dem) / num_clients)
    obj = 99999
    return num_clients, cx, cy, dem, cap, car, obj
end

#====================================================================================================================#
# MÉTODOS AUXILIARES                                                                                                 # 
#====================================================================================================================#

# método que remove rotas que possuem somente o depósito
function remove_dep_rote(rotas::Vector{Vector{Int64}})
    for i = 1:length(rotas)
        if length(rotas[i]) > 1 && rotas[i][1] ≠ 1
            rotas[i] = rotas[i]
        else
            continue
        end
    end 
    return rotas
end 

# método para cálculo matriz de distâncias pelo método euclidiano
function distances(coord_x::Vector{Float64}, coord_y::Vector{Float64})
    n = length(coord_x)
    d = [sqrt((coord_x[j] - coord_x[i])^2 + (coord_y[j] - coord_y[i])^2) for i = 1:n, j = 1:n]
    return d
end

# esse método recebe a matriz de variáveis binárias do modelo tsp e retorna todas as rotas geradas para o problema de roteamento de veículos
function routesmodel(data::Data, x::Matrix{Float64})
    inicios = [j for j = 2:data.n if x[1,j] > 0.99]
    k = inicios[1]
    rota = [1,k]
    rotas = []
    while inicios != []
        for i = 1:data.n
            j = argmax(x[k,:])
            if j != 1          # se o veículo não voltar para o deposíto
                push!(rota, j) # cliente entra na rota
                k = j
            else
                push!(rota, j)
                inicios = deleteat!(inicios, 1)
                push!(rotas, rota)
                if inicios ≠ []
                    k = inicios[1]
                    rota = [1,k]
                else
                    break
                end
            end
        end
    end
    return rotas
end

# método de clusterização pelo problema de alocação de tarefas
function alocation(data::ModelData, opt::SolverParameters)
    
    rotas::Vector{Vector{Int64}} = []
    
    clt = Model(
        optimizer_with_attributes(
            Cbc.Optimizer,
            "seconds"      => opt.timelimit, 
            "allowableGap" => opt.gap, 
            "Threads"      => opt.threads
        )
    )

    set_silent(clt)

    @variable(clt, x[i = 1:data.cars, j = 1:data.n], Bin)

    @objective(clt, Min, sum(data.distance[i,j] * x[i,j] for i = 1:data.cars, j = 1:data.n))

    @constraint(clt, [j = 1:data.n], sum(x[i,j] for i = 1:data.cars) == 1)
    @constraint(clt, [i = 1:data.cars], sum(data.demand[j] * x[i,j] for j = 1:data.n) ≤ data.capacity)
 
    optimize!(clt)

    if termination_status(clt) == MOI.OPTIMAL
    println("Modelo de alocação: solução ótima encontrada !")
        for i = 1:data.cars
            rota::Vector{Int64} = []
            for j = 1:data.n 
                if value.(x[i,j]) > 0.99
                    push!(rota, j)
                end
            end
            push!(rotas, rota)
        end
        return rotas
    else
        println("Modelo de alocação: sem solução encontrada !")
        return 0
    end
end

# função que gera a rota do tsp e o valor da função função objetivo
function route_generate(model::Model, x::AbstractArray, n::Int)
    rota::Vector{Int64} = []
    k = 1
    for i = 1:n
        j = argmax(value.(x[k,:]))
        push!(rota, j)
        k = j
    end
    return rota
end

# método que recebe a matriz de distâncias do tsp e a variável binária x e retorna a rota gerada
function lazy_constraint(model::Model, x::AbstractArray)

    n = size(x, 1)

    done = false
    
    while !done

        optimize!(model) # otimizando o modelo
        
        vx = value.(x)
        
        g = SimpleDiGraph(n)

        for i = 1:n, j = 1:n
            if vx[i,j] > 0.99
                add_edge!(g,i,j)
            end
        end
        c = connected_components(g)
        if length(c) == 1
            done = true
        else
            for s ∈ c
                len_s = length(s)
                @constraint(model, sum(x[i,j] for i ∈ s, j ∈ s) ≤ len_s - 1)
            end
        end
    end

    # gerando a rota
    rota = route_generate(model, x, n)
    return rota
end

# modelo do caixeiro viajante com a técnica de geração linhas (lazy constraints) para o problema de roteirização
function tsp_lazyconstraint(data::Data, rota::Vector{Int64})
    
    # criando o modelo
    model = Model(Cbc.Optimizer)

    set_silent(model)

    # adicionando o depósito a rota
    if 1 ∉ rota
        rota = [1;rota]
    end
    
    # número de clientes
    n = length(rota)

    coordx = data.coord_x[rota]
    coordy = data.coord_y[rota]

    # matriz de custos dos clientes selecionados
    c = distances(coordx, coordy)

    # variáveis: x(i,j) = 1 se o carro vai do cliente i para o cliente j
    @variable(model, x[i = 1:n, j = 1:n], Bin)

    # função objetivo: minimizar o custo total de transporte
    @objective(model, Min, sum(x .* c))

    # restrições
    @constraint(model, [j = 1:n], sum(x[i,j] for i = 1:n if i ≠ j) == 1)

    @constraint(model, [i = 1:n], sum(x[i,j] for j = 1:n if i ≠ j) == 1)

    @constraint(model, [i = 1:n], x[i,i] == 0)

    # adicionando as restrições violadas e reexecutando a otimização
    r = lazy_constraint(model, x)
    if r[end] ≠ 1
        r = [r;1]
    end 
    if r[1] ≠ 1
        r = [1;r]
    end 
    return rota[r]
end

# esse método gera o modelo TSP-MTZ para o problema do caixeiro viajante
function tsp_mtz(data::ModelData, rota::Vector{Int64}, opt::SolverParameters, ativo::Bool)

    # construção do modelo
    tsp = Model(Cbc.Optimizer)

    set_silent(tsp)

    # adicionando o depósito na rota
    if 1 ∉ rota
        rota = [1;rota]
    end
    
    # omitindo a execução do modelo
    set_silent(tsp)

    # setando as dimensões
    n = length(rota)

    cx = data.coord_x[rota]
    cy = data.coord_y[rota]
    
    dist = distances(cx, cy)

    # variáveis
    @variable(tsp, x[i = 1:n, j = 1:n], Bin)
    @variable(tsp, u[i = 1:n] >= 0)

    # função objetivo
    @objective(tsp, Min, sum(dist[i,j] * x[i,j] for i in 1:n, j in 1:n))

    # construção das restrições
    @constraint(tsp, [j = 1:n], sum(x[i,j] for i = 1:n if i != j) == 1)
    @constraint(tsp, [i = 1:n], sum(x[i,j] for j = 1:n if i != j) == 1)
    @constraint(tsp, u[1] == 1)
    @constraint(tsp, [i = 1:n, j = 2:n; i != j], u[j] >= u[i] + 1 - n * (1-x[i,j]))
    
    # resolvendo o modelo
    optimize!(tsp)

    # validando o status da solução
    if termination_status(tsp) == MOI.OPTIMAL
        x = value.(x)
        r = route_generate(tsp, x, n)   
        if r[end] ≠ 1
            r = [r;1]
        end 
        if r[1] ≠ 1
            r = [1;r]
        end         
        custo = objective_value(tsp)
        return rota[r]
    elseif termination_status(tsp) == MOI.TIME_LIMIT && has_values(tsp)
        x = value.(x)
        r = route_generate(tsp, x, n)
        if r[end] ≠ 1
            r = [r;1]
        end
        if r[1] ≠ 1
            r = [1; r]
        end    
        custo = objective_value(tsp)
        return rota[r]
    else
        println("Solução não encontrada !")
    end
end

# método para geração da matriz de rotas para a heurística master_problem
function generate_routes_matrix(data::Data, resultado::ReportData)
    a = zeros(Int, data.n) # matriz de rotas
    t = zeros(Int, data.n) # vetor para armazenamento de rotas
    for i = 1:length(resultado.rotas)
        for j = 2:length(resultado.rotas[i])  # removendo o depósito da rota
            if i == 1
                a[resultado.rotas[i][j],1] = 1
            else
                t[resultado.rotas[i][j],1] = 1
            end
        end
        if i > 1
            a = [a t]
        end
        t .= 0
    end
    return a
end

# método para calcular os custos de cada rota gerada
function routes_cost(data::ModelData, resultado::ReportData)
    # calculando os custos das r rotas
    custo::Vector{Float64} = []
    for i = 1:length(resultado.rotas)
        push!(custo, totalcost(resultado.rotas[i], data.distance))
    end
    return custo
end

# método com o problema master para a heurística geração de colunas
function master_problem(data::Data, a::Matrix{Int64}, resultado::ReportData)
    
    # Indices:
    # N = conjunto de clientes / linhas
    # Ω = conjunto de rotas    / colunas
    
    # Parâmetros:
    # a(i,r) = matriz de alocação igual a 1 se o cliente i está na rota r, 0 caso contrário
    # c(r)   = custo da rota r (calculado inicialmente pelo método construtivo selecionado)
    
    # variáveis:
    # x(r)   = 1 se a rota r é selecionada, 0 caso contrário

    # Modelo:

    # Min  ∑ c(r) * x(r)
    #    r ∈ Ω 

    # s.t:

    # 
    #   ∑ a(i,r) * x(r) >= 1                       ∀ i ∈ N  
    # r ∈ Ω

    # x(r) >= 0                                    ∀ r ∈ Ω

    # acessando o número de linhas e colunas da matriz
    N,Ω = size(a)

    custo = routes_cost(data, resultado) 

    lsp = Model(Cbc.Optimizer)

    set_silent(lsp)

    @variable(lsp, x[r = 1:Ω], Bin)

    @objective(lsp, Min, sum(custo * x))

    @constraint(lsp, alocation[i = 1:N], sum(a[i,r] * x[r] for r = 1:Ω) == 1)

    optimize!(lsp)

    u = dual.(alocation) # acessando os valores duais da restrção de alocação

    return a, u, custo
end

# método com o sub-problema para a heurística geração de colunas
function pricing_problem(data::Data, a::Matrix{Int64}, custo::Vector{Float64}, u::Vector{Float64}, resultado::ReportData)

end

# método para cálculo da melhor inserção do cliente em uma rota
function best_insection(data::Data, rota, v::Int64)
    c = copy(data.distance)
    min_arg = data.n
    min_val = c[rota[end], v] + c[v, rota[1]] - c[rota[end], rota[1]]
    for i = 2:length(rota)
        d = c[rota[i - 1], v] + c[v, rota[i]] - c[rota[i-1], rota[i]]
        if d < min_val
            min_val = d
            min_arg = i
        end
    end
    min_arg
    min_val
    return min_arg, min_val
end

# método que calcula os savings da heurística de clark and wright
function saving(data::Data)
    s = zeros(data.n, data.n)
    c = data.distance
    for i = 1:data.n
        for j = 2:i
            s[i, j] = c[i, 1] + c[1, j] - c[i,j]
            s[j, i] = c[j, 1] + c[1, i] - c[j,i]
        end
    end 
    return c, s
end

# esse método aloca o cliente de acordo com a posição na rota respeitando a capacidade disponível
function firstfit(data::Data, rota::Vector{Int64})
    rotas = []
    r = [1]
    carga::Float64 = 0.0
    for i = 1:length(rota) 
        if carga + data.demand[rota[i]] <= data.capacity
            push!(r, rota[i])
            carga += data.demand[rota[i]]
        else
            push!(r,1)
            push!(rotas, r)
            r = [1, rota[i]]
            carga = data.demand[rota[i]]
        end
    end
    push!(r,1)
    push!(rotas,r)
    return rotas
end

# esse método destroi e reconstroi a rota de forma aleatória
function routedestroy(data::Data)
    rota = shuffle([2:data.n...])
    rotas = firstfit(data, rota)
    custo = totalcost(rotas, data.distance)
    nome = " "
    imp = false
    return ReportData(custo, nome, rotas, imp)
end

# método que calcula a melhoria na troca de pontos no algorítmo 2-opt
function prv(distmat::AbstractMatrix{T} , path::AbstractVector{S}, revL::Int, revH::Int) where {T <: Real, S <: Integer}
    
    cost_delta = zero(eltype(distmat))
    
    # if there's an initial unreversed section
    if revL > 1
        # new step onto the reversed section
        cost_delta += distmat[path[revL - 1], path[revH]]
        # no longer pay the cost of old step onto the reversed section
        cost_delta -= distmat[path[revL - 1], path[revL]]
    end
    
    # The actual cost of walking along the reversed section doesn't change
    
    # because the distance matrix is symmetric.
    # if there's an unreversed section after the reversed bit
    if revH < length(path)
        # new step out of the reversed section
        cost_delta += distmat[path[revL], path[revH + 1]]
        # no longer pay the old cost of stepping out of the reversed section
        cost_delta -= distmat[path[revH], path[revH + 1]]
    end
    return cost_delta
end

# método para a conversão de um vetor de vetores (rotas) em uma matriz quadrada binária para utilização no modelo como MIP-START
function convertmatrix(data::Data, resultado::ReportData)
    vx = zeros(eltype(resultado.rotas[1][1]), data.n, data.n)
    for i = 1:length(resultado.rotas), j = 1:length(resultado.rotas[i])-1
        vx[resultado.rotas[i][j], resultado.rotas[i][j+1]] = 1
    end
    return vx
end

# método que recebe como dados de entrada, a variável de decisão x do modelo e o resultado da heurística
function warmstart(data::Data, x,  resultado::ReportData)
    vx = convertmatrix(data, resultado)
    for i = 1:data.n
        println(vx[i, :])
    end
    # setando o valor da variável binária x com o valor da matriz advinda da heurística
    for i = data.n, j = 1:data.n
        JuMP.set_start_value(x[i, j], vx[i,j])
    end
end

# método de melhoria 2-opt
function tsp_two_opt(data::ModelData, rota::AbstractVector{S}) where {S <: Integer}
    
    # size checks
    n = length(rota)
      
    # how much must each swap improve the cost ?
    thresh = 0.02

    # main loop
    # check every possible switch until no 2-swaps reduce objective
    # if the path passed in is a loop (first/last nodes are the same)
    # then we must keep these the endpoints of the path the same
    # ie just keep it a loop, and therefore it doesn't matter which node is at the end
    # if the path is not a cycle, we should respect the endpoints
    
    switchLow = 2
    switchHigh = n - 1
    need_to_loop = true # always look for swaps at least once
    while need_to_loop
        need_to_loop = false
        # we can't change the first
        for i in switchLow:(switchHigh-1)
            for j in switchHigh:-1:(i+1)
                cost_change = prv(data.distance, rota, i, j)
                if cost_change + thresh <= 0
                    need_to_loop = true
                    reverse!(rota, i, j)
                end
            end
        end
    end
    return rota
end

#====================================================================================================================#
# MÉTODOS PARA RESULTADOS                                                                                            # 
#====================================================================================================================#

# método para custo total das rotas
function totalcost(rotas, c) 
    r = length(rotas)
    custo::Float64 = 0.0
    for i = 1:r
        rota = rotas[i]
        if length(rota) == 1 && rota[1] != 1
            custo += c[1, rota[1]]
        else
            for j = 1:length(rota)-1
                custo += c[rota[j], rota[j+1]]
            end
        end
        custo += c[rota[end], rota[1]] 
    end
    return custo
end

# exibição no terminal do resultado da roteirização
function output(resultado::ReportData)
    println(" ")
    println("Algorítmo   : ", resultado.nome)
    println("Custo total : ", resultado.custo)
    println(" ")
    k = 1
    for i = 1:length(resultado.rotas)
        rota = resultado.rotas[i] .- 1
        if rota[end] != 0
            push!(rota, 0)
        end
        if length(rota) >= 3 && sum(rota) != 0
            println("Rota $k: ", rota)
            k += 1
        end
    end
end

# método para geração de gráfico das rotas
function graph(data::Data, resultado::ReportData, g::Bool)

    if g == true
        # plotando o depósito e clientes
        scatter([data.coord_x[1]], [data.coord_y[1]], shape =:rect, markersize = 10, color = :white, legend = false, size=(900, 750))
        scatter!(data.coord_x[2:end], data.coord_y[2:end], legend = false)
        annotate!(data.coord_x, data.coord_y .+ 3, [1:length(data.coord_y)...] .- 1, font(8, :white))
        
        title!("CVRP - " * resultado.nome)
        
        # calculando o quão próximo da solução ótima, está a solução da heurística
        porcent = (1 - (data.obj / resultado.custo)) * 100
        
        if data.obj != 99999
            xlabel!("Custo obtido: $(round(Int, resultado.custo)) | Custo ótimo: $(data.obj) | $(round(Int, porcent))% do ótimo")
        else
            xlabel!("Custo obtido: $(round(Int, resultado.custo))")
        end
        
        # plotando as rotas
        r = length(resultado.rotas)
        plt = []
        for i = 1:r
            rota = resultado.rotas[i]
            if length(resultado.rotas[i]) > 1
                if resultado.rotas[i][end] != 1
                    rota = push!(resultado.rotas[i], 1)
                end
            end
            for j = 1:length(rota)-1
                plt = scatter!(data.coord_x[rota], data.coord_y[rota], color = i)
                plt = plot!( ([data.coord_x[rota[j]], data.coord_x[rota[j+1]]], [data.coord_y[rota[j]], data.coord_y[rota[j+1]]] ), arrow=true,  color = i)
            end
        end
        display(plt)
    end
end
