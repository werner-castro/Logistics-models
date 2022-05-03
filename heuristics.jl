
export ModelData, readfile, distances, totalcost, two_opt, output, graph, model, clark_wright, next_fit, best_fit, two_opt_star, two_opt_star_best_improvement, cluster_router, rms, column_generation

module Heuristic

    include("utils.jl")

    using Plots, CSV, JuMP, Cbc, DataFrames, Distributed, Random, DelimitedFiles, PlotThemes; theme(:juno)
    
    using Base.Threads: @threads, @spawn

    #====================================================================================================================#
    # MODELO MATEMÁTICO                                                                                                  # 
    #====================================================================================================================#
    function model(data::Data, g::Bool, resultado::Union{Nothing, ReportData}, opt::SolverParameters)

        # modelo baseado em fluxo
        cvrp = Model(
            optimizer_with_attributes(
                Cbc.Optimizer, 
                "Threads" => opt.threads
            )
        )

        # setando o tempo do solver
        set_time_limit_sec(cvrp, opt.timelimit)

        # setando o valor de GAP
        set_optimizer_attribute(cvrp, "allowableGap", opt.gap)

        # variáveis
        @variable(cvrp, x[i = 1:data.n, j = 1:data.n], Bin)                                         # se o carro passa do cliente i para o cliente j
        @variable(cvrp, f[i = 1:data.n, j = 1:data.n] >= 0)                                         # quantidade transportada de i para j
            
        # aplicando a técnica de warmstart
        if typeof(resultado) == ReportData
            warmstart(data, x, resultado)
        end

        # função objetivo
        # @objective(cvrp, Min, sum(data.distance[i,j] * x[i,j] for i = 1:data.n, j = 1:data.n))    # minimizar a distância total percorrida
        @objective(cvrp, Min, sum(data.distance[i...] * x[i] for i in eachindex(x)))

        # restrições
        @constraint(cvrp, [i = 2:data.n], sum(x[i,j] for j = 1:data.n if i != j) == 1)
        @constraint(cvrp, [j = 2:data.n], sum(x[i,j] for i = 1:data.n if j != i) == 1)
        @constraint(cvrp, sum(x[1,j] for j = 2:data.n) == sum(x[j,1] for j = 2:data.n))
        @constraint(cvrp, [j = 2:data.n], sum(f[j,i] for i = 1:data.n) - sum(f[i,j] for i = 1:data.n) == data.demand[j])
        @constraint(cvrp, [i = 1:data.n, j = 1:data.n; i != j], f[i,j] <= data.capacity * x[i,j])
        @constraint(cvrp, [i = 1:data.n], x[i,i] == 0)
        
        # resolvendo o modelo
        optimize!(cvrp)

        # validando o status da solução e gerando relarório de saída
        if termination_status(cvrp) == MOI.OPTIMAL
            x = value.(x)
            rotas = routesmodel(data, x)
            custo = objective_value(cvrp)
            nome  = "CVRP | Status: Optimal [GAP = ($(relative_gap(cvrp)))%]"
            imp   = true
            resultado = ReportData(custo, nome, rotas, imp)
            output(resultado)
            graph(data, resultado, g)
            return resultado
        elseif termination_status(cvrp) == MOI.TIME_LIMIT && has_values(cvrp)
            x = value.(x)
            rotas = routesmodel(data, x)
            custo = objective_value(cvrp)
            nome  = "CVRP | Status: Suboptimal [GAP = ($(relative_gap(cvrp)))%]"
            imp   = true
            resultado = ReportData(custo, nome, rotas, imp)
            output(resultado)
            graph(data, resultado, g)
            return resultado
        else
            println("Solução não encontrada !")
        end
    end

    #====================================================================================================================#
    # HEURÍSTICAS MIP                                                                                                    # 
    #====================================================================================================================#

    # esse método resolve a restrição de capacidade pelo modelo de alocação de tarefas e roteiriza pelo modelo do caixeiro viajante
    # se o parametro 'resultado' for do tipo 'Reportdata' o modelo recebe o resultado da heurística e tenta a melhoria pelo tsp, senão, gera a solução inicial pelo modelo de alocação de tarefas
    function cluster_router(data::Data, parallel::Bool, resultado::Union{ReportData, Bool}, opt::SolverParameters, g::Bool)
        
        println(" ")
        println(" ")
                                    
        if parallel == true
            if resultado == false                          
                rotas = alocation(data, opt)               # fase 1 : alocação de clientes pelo modelo de alocação
                n = length(rotas)
            else
                rotas = resultado.rotas                    # fase 1 : recebendo a solução inicial de uma heurística para tentar gerar uma melhoria
                n = length(rotas)
            end
            println("Iniciando a roteirização:")
            # println("Modelo rodando com: ",Threads.nthreads()," threads")
            routes = rotas
            @spawn begin                                   # executando o loop em multiplas threads
                for i = 1:n                                # fase 2 : roteirização de clientes em pararelo
                    println("Rota [$i] executada na thread [$(Threads.threadid())]")
                    rotas[i] = tsp_lazyconstraint(data, routes[i])
                end
            end
        else
            if resultado == false                          
                rotas = alocation(data, opt)                # fase 1 : alocação de clientes pelo modelo de alocação
                n = length(rotas)
            else
                rotas = resultado.rotas                     # fase 1 : recebendo a solução inicial de uma heurística para tentar gerar uma melhoria
                n = length(rotas)
            end
            for i = 1:n 
                rota = rotas[i]                             # fase 2 : roteirização de clientes em série
                # rotas[i] = tsp_mtz(data, rota, opt, true)
                rotas[i] = tsp_lazyconstraint(data, rota)
            end
        end
        nome = "cluster and router"
        imp = true
        custo = totalcost(rotas, data.distance)
        resultado = ReportData(custo, nome, rotas, imp)
        output(resultado)
        graph(data, resultado, g)
        return resultado
    end
    
    # esse método aplica a heurística de geração de colunas
    function column_generation(data::Data, resultado::Union{ReportData,Bool}, g::Bool)

        # Algorítmo: column generation

        # passo 1: gere um conjunto de rotas factíveis (pode utilizar um método contrutivo / heurística)
        # passo 2: resolva o problema master
        # passo 3: recupere o valor dual (u i valores) no problema master.
        # passo 4: resolva o sub-problema tentando encontrar custo negativo reduzido, 
        #   se a rota adicional tem o custo negativo, pare, a solução ótima foi encontrada, 
        #   senão adicione a nova rota no conjunto de rotas factíveis e volte ao passo 2.

        # gerando as rotas factíveis
        a = generate_routes_matrix(data, resultado)

        negative_cost = false

        while negative_cost == false
            a, u, custo = master_problem(data, a, resultado)
            y, resultado, corrent_cost = pricing_problem(data, a, custo, u, resultado)
            if corrent_cost == false
                a = [a y]  # adicionando uma nova rota na matriz de rotas 'a'
                negative_cost = corrent_cost
            else
                a = [a y]
                negative_cost = corrent_cost
                break
            end
        end
        return a
    end

    #====================================================================================================================#
    # HEURÍSTICAS CONSTRUTIVAS                                                                                           # 
    #====================================================================================================================#
    
    # HEURÍSTICA DE CLARK AND WRIGHT
    function clark_wright(data::Data, g::Bool)

        load = [data.demand[i] for i = 1:data.n]
        rotas = []
        
        # criando as rotas triviais: indo do depósito até um cliente
        for i = 2:data.n
            push!(rotas, [1, i])
        end

        # calculando os savings
        c, s = saving(data)
    
        # concatenando as rotas
        while true
            argmax = nothing
            maxval = 0
            for (k, rk) in enumerate(rotas)
                for (l, rl) in enumerate(rotas)
                    if k != l && maxval < s[rk[end], rl[2]] && load[k] + load[l] <= data.capacity
                        argmax = [k, l]
                        maxval = s[rk[end], rl[2]]
                    end
                end
            end
            if typeof(argmax) != Nothing
                k, l = argmax
                rotas[k] = [rotas[k]; rotas[l][2:end]]
                load[k] += load[l]
                deleteat!(rotas, l)
                deleteat!(load, l)
            else
                break
            end
        end
        imp = true
        nome = "Heuristic: clark and wright"
        custo = totalcost(rotas, c)
        resultado = ReportData(custo, nome, rotas, imp)
        output(resultado)
        graph(data, resultado, g)
        return resultado
    end

    # HEURÍSTICA TSP FIT (heurística de roteirização e alocação)
    # roteiriza como um problema do caixeiro viajante e depois aloca respeitando capaciade dos veículos
    function tspfit(data::Data, g::Bool)
        c = copy(data.distance)
        cliente = 1
        j = 0
        rota = [cliente]
        # roteirização pelo vizinho mais próximo
        for i = 1:data.n
            c[cliente,cliente] = 1000.0
            j = argmin(c[cliente,:])
            c[cliente, :] .= 1000.0
            c[:, cliente] .= 1000.0
            cliente = j
            push!(rota, cliente)
        end
        rotas = firstfit(data, rota)
        custo = totalcost(rotas, data.distance)
        nome = "Heuristic: tsp fit"
        imp = true
        resultado = ReportData(custo, nome, rotas, imp)      
        output(resultado)
        graph(data, resultado, g)
        return resultado
    end

    # HEURÍSTICA ANGULAR FIT
    function angularfit(data::Data, g::Bool)

        # quando x < 0, calcula o arco tangente (y / x) + 180
        # quando y < 0, calcula o arco tangente (y / x) + 360

        angulos = []

        # trazendo as coordenadas para o centro do plano cartesiano
        coord_x = data.coord_x .- data.coord_x[1]
        coord_y = data.coord_y .- data.coord_y[1]

        # calculando os angulos entre os pontos
        for i = 1:data.n
            if coord_x[i] > 0.0 && coord_y[i] > 0.0
                push!(angulos, atand(coord_y[i] / coord_x[i]))
            elseif coord_x[i] > 0.0 && coord_y[i] < 0.0
                push!(angulos, atand(coord_y[i] / coord_x[i]) + 360)
            else
                push!(angulos, atand(coord_y[i] / coord_x[i]) + 180)
            end
        end

        df = [[1:data.n...] angulos]

        # ordenando a matriz com base na coluna dos angulos
        rota = sortslices(df, dims=1, by= x-> (x[2]), rev=false)[:,1]
        
        rotas = firstfit(data, rota)
        custo = totalcost(rotas, data.distance)
        nome = "Heuristic: angular fit"
        imp = true
        resultado = ReportData(custo, nome, rotas, imp)
        output(resultado)
        graph(data, resultado, g)
        return resultado
    end

    #====================================================================================================================#
    # HEURÍSTICAS DE MELHORIA                                                                                            # 
    #====================================================================================================================#
    
    # MELHORIA INTRA ROTA (TROCA CLIENTES NA MESMA ROTA)
    #====================================================================================================================#

    # 2-OPT
    # esse método faz o 2-opt para todas as rotas
    function two_opt(data::Data, resultado::ReportData, g::Bool, parallel::Bool, report::Bool)
        initial_cost = resultado.custo
        if parallel == true
            rotas = resultado.rotas
            @threads for k = 1:length(resultado.rotas)
                resultado.rotas[k] = tsp_two_opt(data, rotas[k])
            end
        else
            for k = 1:length(resultado.rotas)
                resultado.rotas[k] = tsp_two_opt(data, resultado.rotas[k])
            end
        end
        nome = "Heuristic: 2-opt"
        custo = totalcost(resultado.rotas, data.distance)
        if custo ≤ initial_cost
            imp = true
        else
            imp = false
        end
        resultado = ReportData(custo, nome, resultado.rotas, imp)
        if report == true
            output(resultado)
        end
        graph(data, resultado, g)
        return resultado
    end
    
    # MELHORIA INTER ROTAS (TROCA DE CLIENTES EM ROTAS DIFERENTES)
    #====================================================================================================================#

    # 2-OPT STAR FIRST IMPROVEMENT
    # esse método executa o 2-opt em rotas diferentes, pegando o primeiro resultado encontrado
    function two_opt_star(data::Data, resultado::ReportData, g::Bool, report::Bool)
        initial_cost = resultado.custo
        cost = 0.0
        chg = false
        imp = true
        c = copy(data.distance)
        d = copy(data.demand)
        q = copy(data.capacity)
        while imp
            imp = false
            for a = 1:length(resultado.rotas)
                ra = resultado.rotas[a]
                if length(ra) < 3
                    continue
                end
                for i = 1:length(ra)-1
                    vi = ra[i]
                    if (i+1) % length(ra) == 0
                        vni = 1
                    else
                        vni = (i+1) % length(ra)
                    end
                    for b = 1:length(a)
                        rb = resultado.rotas[b]
                        if length(rb) < 3
                            continue
                        end
                        for j = 1:b
                            vj = rb[j]
                            if (j + 1) % length(rb) == 0
                                rvj = 1
                            else
                                rvj = (j + 1) % length(rb)
                            end
                            vnj = rb[rvj]
                            delta = c[vj, vni] + c[vi, vnj] - c[vi, vni] - c[vj, vnj]
                            if delta < -1e-3
                                if sum(d[ra[1:i + 1]]) + sum(d[rb[j + 1:end]]) <= q && sum(d[rb[1:j + 1]]) + sum(d[ra[i + 1:end]]) <= q
                                    na = [ra[1:i + 1]; rb[j + 1:end]]
                                    nb = [rb[1:j + 1]; ra[i + 1:end]]
                                    empty!(ra)
                                    ra = [ra; na]
                                    empty!(rb)
                                    rb = [rb; nb]
                                    chg = imp = true
                                    cost += delta
                                    resultado.rotas[a] = ra
                                    break
                                end
                            end
                            delta = c[vnj, vni] + c[vi, vj] - c[vi, vni] - c[vj, vnj]
                            if delta < -1e-3
                                if sum(d[ra[1:i + 1]]) + sum(d[rb[1:j + 1]]) <= q && sum(d[rb[j + 1:end]]) + sum(d[ra[i + 1:end]]) <= q
                                    na = [ra[1:i + 1]; rb[j:-1:1]]
                                    nb = [1; ra[i:end]; rb[1:j + 1]]
                                    empty!(ra)
                                    ra = [ra; na]
                                    empty!(rb)
                                    rb = [rb; nb]
                                    chg = imp = true
                                    cost += delta
                                    resultado.rotas[a] = ra
                                    break
                                end
                            end
                        end
                    end
                    if imp == false
                        break
                    end
                end
                if imp == false
                    break
                end
            end
            if imp == false
                break
            end
        end
        rotas = remove_dep_rote(resultado.rotas)
        custo = totalcost(rotas, data.distance)
        if custo ≤ initial_cost
            imp = true
        else
            imp = false
        end
        nome = "Heuristic: 2-opt star"
        resultado = ReportData(custo, nome, rotas, imp)
        if report == true
            output(resultado)
        end
        graph(data, resultado, g)
        return resultado
    end

    # 2-OPT STAR BEST IMPROVEMENT
    # esse método executa o 2-opt em rotas diferentes, pegando o melhor resultado obtido
    function two_opt_star_best_improvement(data::Data, resultado::ReportData, g::Bool, report::Bool)
        initial_cost = resultado.custo
        cost = 0.0
        chg = false
        imp = true
        c = copy(data.distance)
        d = copy(data.demand)
        q = copy(data.capacity)
        load = [sum(data.demand[r]) for r in resultado.rotas]
        while imp
            imp = false
            min_delta = 0
            arg = nothing
            for a = 1:length(resultado.rotas)
                ra = resultado.rotas[a]
                for i = 1:length(ra)-1
                    vi = ra[i]
                    rvi = (i + 1) % length(ra)
                    if rvi == 0
                        rvi = 1
                    end
                    vip = ra[rvi]
                    for b = 1:a
                        rb = resultado.rotas[b]
                        for j = 1:length(rb)-1
                            vj = rb[j]
                            rvj = (j + 1) % length(rb)
                            if rvj == 0
                                rvj = 1
                            end
                            vjp = rb[rvj]
                            delta = c[vj, vip] + c[vi, vjp] - c[vi, vip] - c[vj, vjp]
                            if delta < -1e-3
                                if sum(d[ra[1:i + 1]]) + sum(d[rb[j + 1:end]]) <= q && sum(d[rb[1:j + 1]]) + sum(d[ra[i + 1:end]]) <= q
                                    if delta < min_delta
                                        min_delta = delta
                                        arg = a, ra, b, rb, i, j
                                    end
                                end
                            end
                        end
                    end
                end
            end
            if typeof(arg) != Nothing
                a, ra, b, rb, i, j = arg
                na = [ra[1:i + 1]; rb[j + 1:end]]
                nb = [rb[1:j + 1]; ra[i + 1:end]]
                empty!(ra)
                ra = [ra; na]
                empty!(rb)
                rb = [rb; nb]
                load[a] = sum(d[ra])
                load[b] = sum(d[rb])
                chg = imp = true
                cost += min_delta
                resultado.rotas[a] = ra
                if imp == false
                    break
                end
            end             
        end
        resultado.rotas = remove_dep_rote(resultado.rotas)
        nome = "Heuristic: 2-opt star best improvement"
        custo = totalcost(resultado.rotas, data.distance)
        if custo ≤ initial_cost
            imp = true
        else
            imp = false
        end
        resultado = ReportData(custo, nome, resultado.rotas, imp)
        if report == true
            output(resultado)
        end
        graph(data, resultado, g)
        return resultado
    end

    # REPLACE
    # esse método tira o cliente de uma rota e o adiciona em outra rota caso haja melhora da solução e sua alocação não exceda a capacidade
    function replace(data::Data, resultado::ReportData, g::Bool, report::Bool)
        initial_cost = resultado.custo
        cost = 0.0
        chg = false
        imp = true
        c = copy(data.distance)
        d = copy(data.demand)
        q = copy(data.capacity)
        load = [sum(data.demand[r]) for r in resultado.rotas]
        while imp
            imp = false
            for (a, ra) in enumerate(resultado.rotas)
                for (i, vi) in enumerate(ra)
                    if i == 1
                        continue
                    end
                    vra = (i + 1) % length(ra)
                    if vra == 0
                        vra = 1
                    end
                    rem_cost = c[ra[i - 1], ra[vra]] - c[ra[i - 1], ra[i]] - c[ra[i], ra[vra]]
                    if rem_cost > -1e-3
                        continue
                    end
                    min_val = Inf
                    min_arg = nothing
                    for (b, rb) in enumerate(resultado.rotas)
                        if load[b] + d[vi] <= q && a != b
                            insert_pos, add_cost = best_insection(data, rb, vi)
                            if add_cost < min_val && add_cost + rem_cost < -1e-3
                                min_val = add_cost
                                min_arg = b, insert_pos
                                if min_val < 1e-3
                                    break
                                end
                            end
                        end
                    end
                    if typeof(min_arg) ≠ Nothing && min_val + rem_cost < -1e-3 
                        deleteat!(ra, i)
                        load[a] -= d[vi]
                        insert!(resultado.rotas[min_arg[1]], min_arg[2], vi) # insert!(collection, index, item)
                        load[min_arg[1]] += d[vi]
                        chg = imp = true
                        cost += min_val + rem_cost
                        break
                    end 
                end
            end
        end
        rotas = remove_dep_rote(resultado.rotas)
        nome  = "Heuristic: replace"
        custo = totalcost(rotas, data.distance)
        if custo ≤ initial_cost
            imp = true
        else
            imp = false
        end
        resultado = ReportData(custo, nome, rotas, imp)
        if report == true
            output(resultado)
        end
        graph(data, resultado, g)
        return resultado
    end

    # SWAP
    # 
    function swap(data::Data, resultado::ReportData, g::Bool, report::Bool)
        initial_cost = resultado.custo
        cost = 0.0
        q = copy(data.capacity)
        c = copy(data.distance)
        d = copy(data.demand)
        imp = true
        chg = false
        load = [sum(data.demand[r]) for r in resultado.rotas]
        while imp
            imp = false
            for a = 1:length(resultado.rotas)
                ra = resultado.rotas[a]
                for i = 2:length(resultado.rotas[a])
                    vi = ra[i]
                    for b = 1:a
                        rb = resultado.rotas[b]
                        for j = 2:length(rb)
                            vj = rb[j]
                            if load[a] + d[vj] - d[vi] <= q && load[b] + d[vi] - d[vj] <= q
                                if (i + 1) % length(ra) == 0
                                    rva = 1
                                else
                                    rva = (i + 1) % length(ra)
                                end
                                if (j + 1) % length(rb) == 0
                                    rvb = 1
                                else
                                    rvb = (j + 1) % length(rb)
                                end
                                delta = c[ra[i - 1], vj] + c[vj, ra[rva]] - c[ra[i - 1], vi] - c[vi, ra[rva]] + c[rb[j - 1], vi] + c[vi, rb[rvb]] - c[rb[j - 1], vj] - c[vj, rb[rvb]]
                                if delta < -1e-3
                                    ra[i] = vj
                                    rb[j] = vi
                                    load[a] += d[vj] - d[vi]
                                    load[b] += d[vi] - d[vj]
                                    chg = imp = true
                                    vi, vj = vj, vi
                                    cost += delta
                                    resultado.rotas[a] = ra

                                    # adpatação para o tabu search
                                    # 
                                    #
                                    #
                                    ##############################  
                                else
                                    imp = false
                                end
                            end
                        end
                    end
                end
            end
        end
        rotas = remove_dep_rote(resultado.rotas)
        nome  = "Heuristic: swap"
        custo = totalcost(rotas, data.distance)
        if custo ≤ initial_cost
            imp = true
        else
            imp = false
        end
        resultado = ReportData(custo, nome, rotas, imp)
        if report == true
            output(resultado)
        end
        graph(data, resultado, g)
        return resultado
    end

    #====================================================================================================================#
    # METAHEURÍSTICAS                                                                                                    #
    #====================================================================================================================#

    # VARIABLE NEIGHBORHOOD SEARCH
    # esse método executa uma série de buscas locais e selecionando a melhor visinhança
    function vnd(data::Data, resultado::ReportData, g::Bool)
        iter = 10
        imp = true
        count = 0
        while count <= iter
            imp = false
            if imp != true
                count += 1
                resultado = replace(data, resultado, false, false)
                imp = resultado.imp
            end
            if imp != true
                count += 1
                resultado = swap(data, resultado, false, false)
                imp = resultado.imp
            end
            if imp != true
                count += 1
                resultado = two_opt_star(data, resultado, false, false)
                imp = resultado.imp
            end
            if imp != true
                count += 1
                resultado = two_opt(data, resultado, false, false, false)
                imp = resultado.imp
            end
            if imp != true
                if count == iter
                    nome = "Metaheuristic: VND"
                    custo = totalcost(resultado.rotas, data.distance)
                    resultado = ReportData(custo, nome, resultado.rotas, imp) 
                    output(resultado)
                    graph(data, resultado, g)
                    return resultado
                    break
                end
            end
        end
        nome = "Metaheuristic: VND"
        custo = totalcost(resultado.rotas, data.distance)
        resultado = ReportData(custo, nome, resultado.rotas, imp) 
        output(resultado)
        graph(data, resultado, g)
        return resultado
    end

    # RANDOM MULTI START
    # esse método executa
    function rms(data::Data, parallel::Bool, iter::Int64, g::Bool)
        bestcost = Inf
        bestsolution = nothing

        if parallel == true
            @spawn for i = 1:iter
                println("iteração $i excutada na thread $(Threads.threadid())")
                initial_sol = routedestroy(data)
                current_sol = vnd(data, initial_sol, false)
                if bestcost > current_sol.custo
                    bestcost = current_sol.custo
                    bestsolution = current_sol.rotas
                end
            end
        else
            @inbounds for i = 1:iter
                println("iteração $i excutada na thread $(Threads.threadid())")
                initial_sol = routedestroy(data)
                current_sol = vnd(data, initial_sol, false)
                if bestcost > current_sol.custo
                    bestcost = current_sol.custo
                    bestsolution = current_sol.rotas
                end
            end
        end
        nome = "Metaheuristic: RMS"
        imp = true
        custo = bestcost
        resultado = ReportData(custo, nome, bestsolution, imp)
        output(resultado)
        graph(data, resultado, g)
        return resultado
    end

    # ITERATED LOCAL SEARCH
    #
    function ils(data::Data, resultado::ReportData)
        
    end

    # GREEDY RANDOMIZED ADAPTATIVE SEARCH PROCEDURE
    #
    function grasp(data::Data, resultado::ReportData)
        
    end

    # TABU-SEARCH
    #
    function tabu_search(data::Data, resultado::ReportData)
        
    end

    # SCATTER-SEARCH
    #
    function scatter_search(data::Data, resultado::ReportData)
        
    end

    # SIMULATED ANNEALING
    #
    function simulated_annealing(data::Data, resultado::ReportData)
        
    end

    # ANT-COLONY
    # 
    function ant_colony(data::Data, resultado::ReportData)

    end
end # end module
