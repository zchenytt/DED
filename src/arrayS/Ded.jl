module Ded
    export Settings, Build, Train
    include("Settings.jl")
    include("LoadGen.jl")
    include("WindGen.jl")
    include("WindCur.jl")
    include("DR.jl")
    include("Generators.jl")
    include("Reserve1.jl")
    include("Reserve2.jl")
    include("Balance.jl")
    include("Costs.jl")
    include("Build.jl")
    include("Train.jl")
end
