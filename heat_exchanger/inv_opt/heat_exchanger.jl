using JuMP, Gurobi, Ipopt, LinearAlgebra, Random, PyPlot
include("fop.jl")
include("main_algorithm.jl")
include("iop_blocks.jl")
include("vars_initialization.jl")
rng = MersenneTwister(12345)

I = 10
u = rand(rng, 0:200, I)/10 .+ 573
x = [forward_problem(u[i]) for i in 1:I]
p = init_constraint_params()
Q, λ = init_objective_params()
s = [ones(Float64, 2, I), ones(Float64, 5, I), ones(Float64, I)]
vars = IOPVars(x, λ, s, Q, p)

penalty_coeffs = deepcopy(s)
penalty_coeffs = 1.0.*penalty_coeffs

params = AlgParams(1, 1e6, penalty_coeffs, 10, 5)
solve_bcd(vars, params)