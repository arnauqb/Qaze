using Qaze

wind = initialize(string(@__DIR__,"/config_test.toml"))
include("radiations_tests.jl")
include("integration_tests.jl")