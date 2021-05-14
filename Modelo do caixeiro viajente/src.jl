
# Arquivo com recursos para o problema do caixeiro viajante com:

# 1 função para calculo da matriz de distâncias
# 2 função para a geração da rota
# 3 funçoes para geração dos gráficos com pacotes de plotagem distintos

using LinearAlgebra, Plots, PlotlyJS, DataStructures, PlotThemes; theme(:juno)

# função usando array comprehension para o cálculo da matriz de distâncias
function distances(α::Vector{Float64}, β::Vector{Float64})
    n = length(α)
    ζ = [sqrt((α[j] - α[i])^2 + (β[j] - β[i])^2) for i = 1:n, j = 1:n]
    return ζ
end

# gerando a rota
function tour(mt::Array{Float64,2})
    rota = []
    k = 1
    for i = 1:size(mt,1)
        j = argmax(mt[k,:])
        push!(rota,j)
        k = j
    end
    rota = [rota[1:end]; rota[1]]
    return rota
end

# gráfico da rota com plotlyJS e mapbox
function plt_map(α::AbstractArray, β::AbstractArray, rota::AbstractArray)

    # criando os nomes das cidades
    city = "C " .* string.([1:length(α)]...)

    mapbox_access_token = "pk.eyJ1Ijoid2VybmVyMTk4NCIsImEiOiJjazNqOTQxcHYwZ2tmM29ybXRidnZiM2d6In0.pyqB8vW8ljbi_WhBXOz1mg"

    # dados das coordenadas e nomes do locais
    data = scattermapbox(
        lat = α[rota],
        lon = β[rota],
        text = city,
        mode = "markers+lines",
        marker_size = 8,
    )

    layout = Layout(
        autosize = true,
        title = "Gráfico da rota",
        mapbox = attr(
            accesstoken = mapbox_access_token,
            bearing = 0,
            center_lat = α[rota[1]],
            center_lon = β[rota[1]],
            pitch = 0,
            zoom = 4,
            style = "dark"
        )
    )

    # plotando o gráfico
    return PlotlyJS.plot(data, layout)
end

# gráfico somente com plotlyJS
function plt_pltly(α::AbstractArray, β::AbstractArray, rota::AbstractArray)

    # criando os nomes das cidades
    city = "C " .* string.([1:length(α)]...)

    plt1 = PlotlyJS.scatter(
        x = α[rota],      # rodando a coordenada na posição da rota gerada
        y = β[rota],      # rodando a coordenada na posição da rota gerada
        mode = "markers", # definindo a localização das cidades
        name = "Cidade",
        text = city       # nomeando as localidades
    )

    plt2 = PlotlyJS.scatter(
        x = α[rota],    # rodando a coordenada na posição da rota gerada
        y = β[rota],    # rodando a coordenada na posição da rota gerada
        mode = "lines", # definindo a rota
        name = "Rota",
        text = city
    )

    layout = PlotlyJS.Layout( # definindo o layout do gráfico
        title = "Modelo: Caixeiro viajante com $n cidades",
        xaxis=attr(title="Latitude", showgrid=true, zeroline=false),
        yaxis=attr(title="Longitude", showgrid=true, zeroline=false)
    )

    data = [plt1, plt2] # agrupando os gráficos

    return PlotlyJS.plot(data, layout) # plotando
end

# Gráfico com o Plots
function plt_plots(α::Array{Float64,1}, β::Array{Float64,1}, rota::Array{Any,1})
    ENV["GKSwstype"] = "100"
    plt = []
    n = length(α)
    # criando os nomes das cidades
    city = "C " .* string.([1:length(α)]...)

    plt = Plots.scatter(α, β, size = (800, 500))                      # plotando os clientes
    plt = Plots.annotate!(α, β.+ 1.5, city, font(6, color = "white")) # plotando a identiicação dos clientes

    # plotando o ponto de partida
    plt = Plots.scatter!(
        arrow = :arrow,
        [α[rota[1]]], [β[rota[1]]],
        legend = false,
        marker=(:square, 6),
        color = "green"
    )

    # plotando a rota
    plt = Plots.plot!(
        α[rota], β[rota],
        # xlabel = "Distância total percorrida = $(obj) km",
        color  = "white",
        opacity =.7,
        title = "Caixeiro viajante - modelo matemático com $(n) cidades",
        legend = false
    )
    return plt
end

function outputmodel(modelo::Model, α::AbstractArray, β::AbstractArray)

    # validando o status da solução e gerando o relatório de saída
    if termination_status(modelo) == MOI.OPTIMAL

        # limpando o terminal após a otimização
        #clearconsole()

        # pegando o resultado da função objetivo
        obj = round.(Int, objective_value(modelo))

        # gerando a rota
        rota = tour(value.(x))

        # exibindo resultado e rota e a distância total percorrida
        println(" ")
        println("======= Solução ótima encontrada ========")
        println(" ")
        println("Distância total percorrida = $(obj) km ")
        println(" ")
        println("Rota gerada = ", [i for i in rota])
        println(" ")

        # gráfico da rota
        plt_plots(α,β,rota)
    else
        println("Solução ótima não encontrada !")
    end
end
