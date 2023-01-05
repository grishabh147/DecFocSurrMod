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
    for iter_count in 1:params.max_outer_iter
        println(iter_count)
        for _ in 1:params.max_inner_iter
            subproblem_1!(vars, params)
            subproblem_2!(vars, params)
            subproblem_3!(vars, params)
            println(norm(x - vars.primal))
            println(norm(vars.slacks, 1))
            println(vars.obj_params)
            println(vars.constraint_params)
        end
        params.penalty_params .+= params.η .* vars.slacks
    end
end

