
# matriz de distâncias
function distances(x::Vector{Float64}, y::Vector{Float64})
    m = length(y)
    ζ = [sqrt((x[j] - x[i])^2 + (y[j] - y[i])^2) for i = 1:m, j = 1:m]
    return ζ
end

# gerando a rota
function tour(mt::AbstractArray)
    k = 1
    rota = [k]
    for i = 1:size(mt,1)
        j = argmax(mt[k,:])
        if j != 1         # se o veículo não voltar para o deposíto
            push!(rota,j) # cliente entra na rota
            k = j         # o carro vai para o próximo cliente k
        else
            push!(rota,j) # cliente entra na rota e a rota é encerrada
            break
        end
    end
    return rota
end

# função de saída de dados e plotagem da rota
function report(
    modelo::Model,
    varx::AbstractArray,
    varb::AbstractArray,
    inicio::OffsetArray{Float64,1,Array{Float64,1}},
    fim::OffsetArray{Float64,1,Array{Float64,1}},
    vx::AbstractArray,
    vy::AbstractArray,
    cost::OffsetArray{Float64,2,Array{Float64,2}})

    if termination_status(modelo) == MOI.OPTIMAL
        K = [1:size(varx)[3]...]
        obj = 0.0
        plt = []
        plt = Plots.scatter(vx, vy, c =:white, legend = false, size = (900, 600))
        plt = Plots.scatter!([vx[1]], [vy[1]], c =:white, marker =:rect, markersize = 5)
        println(" ")
        println("| Modelo: Roteamento de veículos capacitado com janelas de tempo - (CVRPTW)")
        println("|  ")
        println("| Número de veículos: $(size(varx)[3])")
        println("| ")
        println("| Número de rotas: $(size(varx)[3])")
        println("| ")
        println("| Número de clientes: $(size(varx)[2] - 1)")
        println(" ")
        for k in K
            if sum(varx[:,:,k]) > 1
                rota = tour(varx[:,:,k].data)
                I = rota[1:end-1] .- 1
                J = rota[2:end] .- 1
                for i in I, j in J
                    obj += varx[i,j,k] * cost[i,j]
                end
                print("| Rota $k = ")
                print("$([i-1 for i in rota])")
                println(" ")
                plt = Plots.plot!(
                    vx[rota],
                    vy[rota],
                    c = k,
                    w = 1, # largura da linha
                    arrow = :arrow
                )
                tw = varb[:,k].data # janela de tempo do carro k
                rota = rota[2:end]
                plt = Plots.annotate!(vx .- 4, vy .+ 4, inicio[0:end-1], font(6, color=:white))
                plt = Plots.annotate!(vx .+ 4, vy .+ 4, fim[0:end-1], font(6, color=:white))
                plt = Plots.annotate!(vx[rota] .+ 0, vy[rota] .+ 4, round.(Int, tw[rota]), font(6, color=:steelblue))
                plt = Plots.scatter!(
                    title = "Modelo: CVRPTW",
                    vx[rota],
                    vy[rota],
                    markersize = 9,
                    msc = k, # cor da borda
                    msw = 3, # espessura da borda
                    c = 3    # cor do ponto / marcador
                )
            else
                continue
            end
        end
        println(" ")
        println("| Custo total da rota: $(obj)")
        plt = annotate!(vx[2:end-1], vy[2:end-1], C, font(6, color=:white))
        plt = Plots.scatter!(
            [vx[1]],
            [vy[1]],
            c =:white,
            marker =:rect,
            markersize = 8
        )
        return plt
    else
        println("Sem solução encontrada !")
    end
end
