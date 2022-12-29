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

function solve_bcd(vars::IOPVars, params::AlgParams)

    # update variables
    init_constraint_params!()
    init_objective_params!()

    for iter_count in 1:params.max_outer_iter
        for _ in 1:params.max_inner_iter
            subproblem_1!()
            subproblem_2!()
            subproblem_3!()
        end
    end
end
