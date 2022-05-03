
# meta pacote para o problema de roteamento de veículos capacitado com heurísticas e modelos
module Cvrp
    include("heuristics.jl")
    include("utils.jl")
    using .Heuristic
end # end module