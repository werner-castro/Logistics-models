
# Arquivo com recursos para o problema do caixeiro viajante com:

# 1 função para calculo da matriz de distâncias
# 2 função para a geração da rota

using Plots, DataStructures, PlotThemes; theme(:juno)

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
