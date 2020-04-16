import Pkg
Pkg.activate("/home/arnau/code/Qaze")
using Qaze

#wind = initialize(string(@__DIR__,"/config_test.toml"))
#include("radiations_tests.jl")
#include("integration_tests.jl")
include("tau_tests.jl")
