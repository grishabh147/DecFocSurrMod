using JuMP, Gurobi, Ipopt, LinearAlgebra, Random, PyPlot
import("fop.jl")

I = 100

T5 = rand(rng, 10:40, I)/10 .+ 586
x = [forward_problem(T5[i]) for i in 1:10]


