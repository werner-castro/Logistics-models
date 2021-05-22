
# matriz de distâncias
function distances(x::Vector{Float64}, y::Vector{Float64})
    ζ = [sqrt( (x[j] - x[i])^2 + (y[j] - y[i])^2 ) for i = 1:length(x), j = 1:length(x)]
    return ζ
end

# coletando e ordenando os pontos para rota
function tour(mt)
    inicios = [j for j = 1:n+1 if value.(x[1,j]) > 0.99]
    mt = value.(x).data
    k = inicios[1]
    rota = [1,k]
    while inicios != []
        for i = 1:size(mt,1)
            j = argmax(mt[k,:])
            if j != 1 # se o veículo não voltar para o deposíto
                push!(rota,j) # cliente entra na rota
                k = j
            else
                push!(rota,j)
                inicios = deleteat!(inicios,1)
                if inicios != []
                    k = inicios[1]
                    push!(rota,k)
                else
                    break
                end
            end
        end
    end
    return rota
end

# Gráfico com o Plots
function pltour(α::Vector{Float64}, β::Vector{Float64}, rota::Vector{Int64}; animated::Int64)

    # gerando as cores em função da rota
    i = 1
    cor::Vector{Int64} = []
    for j = 1:length(rota)
        if rota[j] != 1
            push!(cor, i)
        else
            i += 3
            push!(cor, i)
        end
    end

    plt = []

    # criando os nomes das cidades
    cidades = "C " .* string.([1:length(α)-2]...)

    # plotando os clientes
    plt = Plots.scatter(α, β, color = "white")

    # plotando a identificação dos clientes
    plt = Plots.annotate!(α[2:end-1], β[2:end-1].+ 4, cidades, font(7, color = "white"))

    if animated != 1

        plt = Plots.plot!(
            arrow = :arrow,
            α[rota],
            β[rota],
            color = cor,
            opacity = .7,
            title = "Roteamento de veículos capacitado - CVRP",
            legend = false
        )

        # colorindo os clientes de acordo com a sua respectiva rota
        plt = scatter!(α[rota], β[rota], color = cor)

        # plotando o ponto de partida (depósito)
        plt = Plots.scatter!(
            [α[rota[1]]],
            [β[rota[1]]],
            legend = false,
            marker=(:square, 6),
            color = "white"
        )
        return plt
    else
        # plotando a rota com animação
        route = [rota[1:end-1] rota[2:end]]
        city = SortedDict{Any,Any}()
        for i in 1:length(α)
            city[i] = (α[i], β[i])
        end
        for i = 1:size(route, 1)
            plt = Plots.plot!(
                arrow = :arrow,
                [city[route[i,1]], city[route[i,2]]],
                color = cor[i],
                opacity =.7,
                size = (700, 550),
                title = "Roteamento de veículos capacitado - CVRP",
                legend = false
            )

            # colorindo das cidades
            plt = Plots.scatter!([city[route[i,1]], city[route[i,2]]], color = cor[i])

            # plotando o ponto de partida (depósito)
            plt = Plots.scatter!(
                [α[rota[1]]], [β[rota[1]]],
                legend = false,
                marker=(:square, 6),
                color = "white"
            )
            display(plt)
        end
    end
end

function relatorio(model::Model, obj::Float64, rota::Vector{Int64}, m::Int64, n::Int64)

    # clearconsole() # essa função ainda não funciona no vscode. 

    rotas = split(rota[2:end-1], 1)

    # relátório de saída
    println(" ")
    println("| Modelo: Roteamento de veículos capacitado - (CVRP)")
    println("| ")
    println("| Número de veículos: $m")
    println("| ")
    println("| Número de rotas: $(length(rotas))")
    println("| ")
    println("| Número de clientes: $n")
    println("| ")
    println("| Custo total da rota = $(obj)")
    println(" ")
    for i = 1:length(rotas)
        rotas[i] = [0; rotas[i] .- 1; 0]
        println("| Rota $i = $(rotas[i])")
    end
end
