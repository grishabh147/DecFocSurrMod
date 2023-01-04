function subproblem_1()
    T3 = 388
    NoS = length(u)
    IOP = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
    @JuMP.variables(IOP, begin
        hat_Qc[i = 1:NoS]
        hat_FH[i = 1:NoS]
        a[1:NoS]
        b[1:NoS]
        slack1[1:2, 1:NoS] >= 0
        slack2[1:5, 1:NoS] >= 0
        slack3[1:NoS] <= 0
        w1[1:NoS]
        w2[1:NoS]
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

    @constraint(IOP, [s = 1:NoS], w1[s] == a[s]*vars.dual[s][2])
    @constraint(IOP, [s = 1:NoS], w2[s] == b[s]*vars.dual[s][2])

    @constraint(IOP, [s = 1:NoS], (hat_Qc[s]/2 + 553 - T3)*vars.dual[s][1] <= slack2[1, s])
    @constraint(IOP, [s = 1:NoS], -((hat_Qc[s]/2 + 553 - T3)*vars.dual[s][1]) <= slack2[1, s])

    @constraint(IOP, [s = 1:NoS], (-vars.dual[s][2]*10 - vars.dual[s][2]*hat_Qc[s] + 
    (u[s]- T3 - 170)*vars.dual[s][2]*hat_FH[s] + w1[s]*hat_Qc[s] + w2[s]*hat_FH[s]) <= slack2[2, s])
    @constraint(IOP, [s = 1:NoS], -(-vars.dual[s][2]*10 - vars.dual[s][2]*hat_Qc[s] + 
    (u[s]- T3 - 170)*vars.dual[s][2]*hat_FH[s] + w1[s]*hat_Qc[s] + w2[s]*hat_FH[s]) <= slack2[2, s])

    @constraint(IOP, [s = 1:NoS], (2*T3 - 786 - hat_Qc[s] + (u[s] - 393)*hat_FH[s])*vars.dual[s][3] <= slack2[3, s])
    @constraint(IOP, [s = 1:NoS], -((2*T3 - 786 - hat_Qc[s] + (u[s] - 393)*hat_FH[s])*vars.dual[s][3]) <= slack2[3, s])

    @constraint(IOP, [s = 1:NoS], ((2*T3 - 1026 - hat_Qc[s] + (u[s] - 313)*hat_FH[s])*vars.dual[s][4]) <= slack2[4, s])
    @constraint(IOP, [s = 1:NoS], -((2*T3 - 1026 - hat_Qc[s] + (u[s] - 313)*hat_FH[s])*vars.dual[s][4]) <= slack2[4, s])

    @constraint(IOP, [s = 1:NoS], ((2*T3 - 1026 - hat_Qc[s] + (u[s] - 323)*hat_FH[s])*vars.dual[s][5]) <= slack2[5, s])
    @constraint(IOP, [s = 1:NoS], -((2*T3 - 1026 - hat_Qc[s] + (u[s] - 323)*hat_FH[s])*vars.dual[s][5]) <= slack2[5, s])

    @objective(IOP, Min, 1000*sum((hat_Qc[s] - x[s][1])^2 for s=1:NoS) + sum((hat_FH[s] - x[s][2])^2 for s=1:NoS)
    + (sum(slack1) + sum(slack2) - sum(slack3) - sum(dual4)) + 1e-5*sum(t))
    optimize!(IOP)
    return value.(Q), value.(p), value.(dual1), value.(dual2), value.(dual3)
end

function subproblem_2()
    T3 = 388
    NoS = length(u)
    IOP = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
    @JuMP.variables(IOP, begin
        hat_Qc[i = 1:NoS]
        hat_FH[i = 1:NoS]
        a[1:NoS]
        b[1:NoS]
        slack1[1:2, 1:NoS] >= 0
        slack2[1:5, 1:NoS] >= 0
        slack3[1:NoS] <= 0
        w1[1:NoS]
        w2[1:NoS]
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

    @constraint(IOP, [s = 1:NoS], w1[s] == a[s]*vars.dual[s][2])
    @constraint(IOP, [s = 1:NoS], w2[s] == b[s]*vars.dual[s][2])

    @constraint(IOP, [s = 1:NoS], (hat_Qc[s]/2 + 553 - T3)*vars.dual[s][1] <= slack2[1, s])
    @constraint(IOP, [s = 1:NoS], -((hat_Qc[s]/2 + 553 - T3)*vars.dual[s][1]) <= slack2[1, s])

    @constraint(IOP, [s = 1:NoS], (-vars.dual[s][2]*10 - vars.dual[s][2]*hat_Qc[s] + 
    (u[s]- T3 - 170)*vars.dual[s][2]*hat_FH[s] + w1[s]*hat_Qc[s] + w2[s]*hat_FH[s]) <= slack2[2, s])
    @constraint(IOP, [s = 1:NoS], -(-vars.dual[s][2]*10 - vars.dual[s][2]*hat_Qc[s] + 
    (u[s]- T3 - 170)*vars.dual[s][2]*hat_FH[s] + w1[s]*hat_Qc[s] + w2[s]*hat_FH[s]) <= slack2[2, s])

    @constraint(IOP, [s = 1:NoS], (2*T3 - 786 - hat_Qc[s] + (u[s] - 393)*hat_FH[s])*vars.dual[s][3] <= slack2[3, s])
    @constraint(IOP, [s = 1:NoS], -((2*T3 - 786 - hat_Qc[s] + (u[s] - 393)*hat_FH[s])*vars.dual[s][3]) <= slack2[3, s])

    @constraint(IOP, [s = 1:NoS], ((2*T3 - 1026 - hat_Qc[s] + (u[s] - 313)*hat_FH[s])*vars.dual[s][4]) <= slack2[4, s])
    @constraint(IOP, [s = 1:NoS], -((2*T3 - 1026 - hat_Qc[s] + (u[s] - 313)*hat_FH[s])*vars.dual[s][4]) <= slack2[4, s])

    @constraint(IOP, [s = 1:NoS], ((2*T3 - 1026 - hat_Qc[s] + (u[s] - 323)*hat_FH[s])*vars.dual[s][5]) <= slack2[5, s])
    @constraint(IOP, [s = 1:NoS], -((2*T3 - 1026 - hat_Qc[s] + (u[s] - 323)*hat_FH[s])*vars.dual[s][5]) <= slack2[5, s])

    @objective(IOP, Min, 1000*sum((hat_Qc[s] - x[s][1])^2 for s=1:NoS) + sum((hat_FH[s] - x[s][2])^2 for s=1:NoS)
    + (sum(slack1) + sum(slack2) - sum(slack3) - sum(dual4)) + 1e-5*sum(t))
    optimize!(IOP)
end

function subproblem_3()
end