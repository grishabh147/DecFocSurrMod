mutable struct IOPVars
    primal::AbstractArray
    dual::AbstractArray
    slacks::AbstractArray
    obj_params::AbstractArray
    constraint_params::AbstractArray
end

mutable struct AlgParams
    η::Float64
    γ::Float64
    penalty_params::AbstractArray
    max_outer_iter::Int64
    max_inner_iter::Int64
end

function solve_bcd()

vars = IOPVars(x, λ, s, Q, p)
alg_params = AlgParams(10, 1e-6, 1, 100, 10)

    for iter_count in 1:params.max_outer_iter
        for _ in 1:params.max_inner_iter
            subproblem_1!()
            subproblem_2!()
            subproblem_3!()
        end
    end
end

