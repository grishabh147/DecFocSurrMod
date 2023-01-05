function subproblem_1!(vars::IOPVars, params::AlgParams)
    T3 = 388
    NoS = length(u)
    IOP = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
    @JuMP.variables(IOP, begin
        hat_Qc[1:NoS]
        hat_FH[1:NoS]
        a[1:NoS]
        b[1:NoS]
        slack1[1:2, 1:NoS] >= 0
        slack2[1:5, 1:NoS] >= 0
        slack3[1:NoS] <= 0
    end)
    @constraint(IOP, [s = 1:NoS], 2*vars.obj_params[1]*hat_Qc[s]+vars.obj_params[2]-
    0.5*vars.dual[s][1]-(a[s]-1)*vars.dual[s][2]+vars.dual[s][3]+vars.dual[s][4]-vars.dual[s][5] <= slack1[1, s])
    @constraint(IOP, [s = 1:NoS], -(2*vars.obj_params[1]*hat_Qc[s]+vars.obj_params[2]-
    0.5*vars.dual[s][1]-(a[s]-1)*vars.dual[s][2]+vars.dual[s][3]+vars.dual[s][4]-vars.dual[s][5]) <= slack1[1, s])

    @constraint(IOP, [s = 1:NoS], 2*vars.obj_params[3]*hat_FH[s]+vars.obj_params[4]-
    (u[s]-T3-170+b[s])*vars.dual[s][2]-(u[s]-393)*vars.dual[s][3]-
    (u[s]-313)*vars.dual[s][4]+(u[s]-323)*vars.dual[s][5] <= slack1[2, s])
    @constraint(IOP, [s = 1:NoS], -(2*vars.obj_params[3]*hat_FH[s]+vars.obj_params[4]-
    (u[s]-T3-170+b[s])*vars.dual[s][2]-(u[s]-393)*vars.dual[s][3]-
    (u[s]-313)*vars.dual[s][4]+(u[s]-323)*vars.dual[s][5]) <= slack1[2, s])

    @constraint(IOP, [s = 1:NoS], hat_Qc[s]/2 + 553 - T3 >= 0)
    @constraint(IOP, [s = 1:NoS], -10 - hat_Qc[s] + (u[s]- T3 - 170)*hat_FH[s] + 
    a[s]*hat_Qc[s] + b[s]*hat_FH[s] >= slack3[s])
    @constraint(IOP, [s = 1:NoS], 2*T3 - 786 - hat_Qc[s] + (u[s] - 393)*hat_FH[s] >= 0)
    @constraint(IOP, [s = 1:NoS], 2*T3 - 1026 - hat_Qc[s] + (u[s] - 313)*hat_FH[s] >= 0)
    @constraint(IOP, [s = 1:NoS], 2*T3 - 1026 - hat_Qc[s] + (u[s] - 323)*hat_FH[s] <= 0)

    @constraint(IOP, [s = 1:NoS], a[s] == sum(vars.constraint_params[1, j]*u[s]^(j-1) for j = 1:4))
    @constraint(IOP, [s = 1:NoS], b[s] == sum(vars.constraint_params[2, j]*u[s]^(j-1) for j = 1:4))

    @constraint(IOP, [s = 1:NoS], (hat_Qc[s]/2 + 553 - T3)*vars.dual[s][1] <= slack2[1, s])
    @constraint(IOP, [s = 1:NoS], -((hat_Qc[s]/2 + 553 - T3)*vars.dual[s][1]) <= slack2[1, s])

    @constraint(IOP, [s = 1:NoS], (-vars.dual[s][2]*10 - vars.dual[s][2]*hat_Qc[s] + 
    (u[s]- T3 - 170)*vars.dual[s][2]*hat_FH[s] + a[s]*vars.dual[s][2]*hat_Qc[s] + b[s]*vars.dual[s][2]*hat_FH[s]) <= slack2[2, s])
    @constraint(IOP, [s = 1:NoS], -(-vars.dual[s][2]*10 - vars.dual[s][2]*hat_Qc[s] + 
    (u[s]- T3 - 170)*vars.dual[s][2]*hat_FH[s] + a[s]*vars.dual[s][2]*hat_Qc[s] + b[s]*vars.dual[s][2]*hat_FH[s]) <= slack2[2, s])

    @constraint(IOP, [s = 1:NoS], (2*T3 - 786 - hat_Qc[s] + (u[s] - 393)*hat_FH[s])*vars.dual[s][3] <= slack2[3, s])
    @constraint(IOP, [s = 1:NoS], -((2*T3 - 786 - hat_Qc[s] + (u[s] - 393)*hat_FH[s])*vars.dual[s][3]) <= slack2[3, s])

    @constraint(IOP, [s = 1:NoS], ((2*T3 - 1026 - hat_Qc[s] + (u[s] - 313)*hat_FH[s])*vars.dual[s][4]) <= slack2[4, s])
    @constraint(IOP, [s = 1:NoS], -((2*T3 - 1026 - hat_Qc[s] + (u[s] - 313)*hat_FH[s])*vars.dual[s][4]) <= slack2[4, s])

    @constraint(IOP, [s = 1:NoS], ((2*T3 - 1026 - hat_Qc[s] + (u[s] - 323)*hat_FH[s])*vars.dual[s][5]) <= slack2[5, s])
    @constraint(IOP, [s = 1:NoS], -((2*T3 - 1026 - hat_Qc[s] + (u[s] - 323)*hat_FH[s])*vars.dual[s][5]) <= slack2[5, s])

    @objective(IOP, Min, sum((hat_Qc[s] - x[s][1])^2 for s=1:NoS) + sum((hat_FH[s] - x[s][2])^2 for s=1:NoS)
    + (sum(params.penalty_params[1].*slack1) + sum(params.penalty_params[2].*slack2) - 
    sum(params.penalty_params[3].*slack3)))
    optimize!(IOP)

    for s in 1:NoS
        vars.primal[s][1] = value(hat_Qc[s])
        vars.primal[s][2] = value(hat_FH[s])
    end

    vars.slacks = [value.(slack1), value.(slack2), value.(slack3)]
    println(params.penalty_params[1])
    println(vars.slacks)
end

function subproblem_2!(vars::IOPVars, params::AlgParams)
    T3 = 388
    NoS = length(u)
    IOP = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
    @JuMP.variables(IOP, begin
        λ[1:NoS, 1:5] >= 0
        Q[i = 1:4]
        a[1:NoS]
        b[1:NoS]
        slack1[1:2, 1:NoS] >= 0
        slack2[1:5, 1:NoS] >= 0
    end)

    @constraint(IOP, [s = 1:NoS], 2*Q[1]*vars.primal[s][1]+Q[2]-0.5*λ[s, 1]-
    (a[s]-1)*λ[s, 2]+λ[s, 3]+λ[s, 4]-λ[s, 5] <= slack1[1, s])
    @constraint(IOP, [s = 1:NoS], -(2*Q[1]*vars.primal[s][1]+Q[2]-0.5*λ[s, 1]-
    (a[s]-1)*λ[s, 2]+λ[s, 3]+λ[s, 4]-λ[s, 5]) <= slack1[1, s])

    @constraint(IOP, [s = 1:NoS], 2*Q[3]*vars.primal[s][2]+Q[4]-(u[s]-T3-170+
    b[s])*λ[s, 2]-(u[s]-393)*λ[s, 3]-(u[s]-313)*λ[s, 4]+(u[s]-323)*λ[s, 5] <= slack1[2, s])
    @constraint(IOP, [s = 1:NoS], -(2*Q[3]*vars.primal[s][2]+Q[4]-(u[s]-T3-170+
    b[s])*λ[s, 2]-(u[s]-393)*λ[s, 3]-(u[s]-313)*λ[s, 4]+(u[s]-323)*λ[s, 5]) <= slack1[2, s])

    @constraint(IOP, [s = 1:NoS], a[s] == sum(vars.constraint_params[1, j]*u[s]^(j-1) for j = 1:4))
    @constraint(IOP, [s = 1:NoS], b[s] == sum(vars.constraint_params[2, j]*u[s]^(j-1) for j = 1:4))

    @constraint(IOP, [s = 1:NoS], (vars.primal[s][1]/2 + 553 - T3)*λ[s, 1] <= slack2[1, s])
    @constraint(IOP, [s = 1:NoS], -((vars.primal[s][1]/2 + 553 - T3)*λ[s, 1]) <= slack2[1, s])

    @constraint(IOP, [s = 1:NoS], (-λ[s, 2]*10 - λ[s, 2]*vars.primal[s][1] + 
    (u[s]- T3 - 170)*λ[s, 2]*vars.primal[s][2] + a[s]*λ[s, 2]*vars.primal[s][1] + b[s]*λ[s, 2]*vars.primal[s][2]) <= slack2[2, s])
    @constraint(IOP, [s = 1:NoS], -(-λ[s, 2]*10 - λ[s, 2]*vars.primal[s][1] + 
    (u[s]- T3 - 170)*λ[s, 2]*vars.primal[s][2] + a[s]*λ[s, 2]*vars.primal[s][1] + b[s]*λ[s, 2]*vars.primal[s][2]) <= slack2[2, s])

    @constraint(IOP, [s = 1:NoS], (2*T3 - 786 - vars.primal[s][1] + 
    (u[s] - 393)*vars.primal[s][2])*λ[s, 3] <= slack2[3, s])
    @constraint(IOP, [s = 1:NoS], -((2*T3 - 786 - vars.primal[s][1] + 
    (u[s] - 393)*vars.primal[s][2])*λ[s, 3]) <= slack2[3, s])

    @constraint(IOP, [s = 1:NoS], ((2*T3 - 1026 - vars.primal[s][1] + 
    (u[s] - 313)*vars.primal[s][2])*λ[s, 4]) <= slack2[4, s])
    @constraint(IOP, [s = 1:NoS], -((2*T3 - 1026 - vars.primal[s][1] + 
    (u[s] - 313)*vars.primal[s][2])*λ[s, 4]) <= slack2[4, s])

    @constraint(IOP, [s = 1:NoS], ((2*T3 - 1026 - vars.primal[s][1] + 
    (u[s] - 323)*vars.primal[s][2])*λ[s, 5]) <= slack2[5, s])
    @constraint(IOP, [s = 1:NoS], -((2*T3 - 1026 - vars.primal[s][1] + 
    (u[s] - 323)*vars.primal[s][2])*λ[s, 5]) <= slack2[5, s])
    @constraint(IOP, Q[3] == 1)

    @objective(IOP, Min, sum(params.penalty_params[1].*slack1) + sum(params.penalty_params[2].*slack2) + 
    (1/(2*params.γ)) * (dot(Q - vars.obj_params, Q - vars.obj_params) + 
    sum((λ[s, i] - vars.dual[s][i])^2 for s in 1:NoS, i in 1:5)))
    optimize!(IOP)
    
    vars.obj_params = value.(Q)
    for s in 1:NoS
        for i in 1:5
            vars.dual[s][i] = value(λ[s, i])
        end
    end
end

function subproblem_3!(vars::IOPVars, params::AlgParams)
    T3 = 388
    NoS = length(u)
    IOP = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
    @JuMP.variables(IOP, begin
        a[1:NoS]
        b[1:NoS]
        slack1[1:2, 1:NoS] >= 0
        slack2[1:5, 1:NoS] >= 0
        slack3[1:NoS] <= 0
        p[1:2, 1:4]
        t[1:2, 1:4] >= 0
    end)
    @constraint(IOP, [s = 1:NoS], 2*vars.obj_params[1]*vars.primal[s][1]+vars.obj_params[2]-
    0.5*vars.dual[s][1]-(a[s]-1)*vars.dual[s][2]+vars.dual[s][3]+vars.dual[s][4]-vars.dual[s][5] <= slack1[1, s])
    @constraint(IOP, [s = 1:NoS], -(2*vars.obj_params[1]*vars.primal[s][1]+vars.obj_params[2]-
    0.5*vars.dual[s][1]-(a[s]-1)*vars.dual[s][2]+vars.dual[s][3]+vars.dual[s][4]-vars.dual[s][5]) <= slack1[1, s])

    @constraint(IOP, [s = 1:NoS], 2*vars.obj_params[3]*vars.primal[s][2]+vars.obj_params[4]-
    (u[s]-T3-170+b[s])*vars.dual[s][2]-(u[s]-393)*vars.dual[s][3]-(u[s]-313)*vars.dual[s][4]+
    (u[s]-323)*vars.dual[s][5] <= slack1[2, s])
    @constraint(IOP, [s = 1:NoS], -(2*vars.obj_params[3]*vars.primal[s][2]+vars.obj_params[4]-
    (u[s]-T3-170+b[s])*vars.dual[s][2]-(u[s]-393)*vars.dual[s][3]-(u[s]-313)*vars.dual[s][4]+
    (u[s]-323)*vars.dual[s][5]) <= slack1[2, s])

    @constraint(IOP, [s = 1:NoS], -10 - vars.primal[s][1] + (u[s]- T3 - 170)*vars.primal[s][2] + 
    a[s]*vars.primal[s][1] + b[s]*vars.primal[s][2] >= slack3[s])

    @constraint(IOP, [s = 1:NoS], a[s] == sum(p[1, j]*u[s]^(j-1) for j = 1:4))
    @constraint(IOP, [s = 1:NoS], b[s] == sum(p[2, j]*u[s]^(j-1) for j = 1:4))

    @constraint(IOP, [s = 1:NoS], (-vars.dual[s][2]*10 - vars.dual[s][2]*vars.primal[s][1] + 
    (u[s]- T3 - 170)*vars.dual[s][2]*vars.primal[s][2] + a[s]*vars.dual[s][2]*vars.primal[s][1] + b[s]*vars.dual[s][2]*vars.primal[s][2]) <= slack2[2, s])
    @constraint(IOP, [s = 1:NoS], -(-vars.dual[s][2]*10 - vars.dual[s][2]*vars.primal[s][1] + 
    (u[s]- T3 - 170)*vars.dual[s][2]*vars.primal[s][2] + a[s]*vars.dual[s][2]*vars.primal[s][1] + b[s]*vars.dual[s][2]*vars.primal[s][2]) <= slack2[2, s])
    @constraint(IOP, [i = 1:2, j = 1:4], p[i, j] <= t[i, j])
    @constraint(IOP, [i = 1:2, j = 1:4], -p[i, j] <= t[i, j])

    @objective(IOP, Min, (sum(params.penalty_params[1].*slack1) + sum(params.penalty_params[2].*slack2) - 
    sum(params.penalty_params[3].*slack3)) + 
    1e-3*sum(t) + 
    (1/(2*params.γ))*sum((p - vars.constraint_params).*(p - vars.constraint_params)))
    optimize!(IOP)

    vars.constraint_params = value.(p)
end