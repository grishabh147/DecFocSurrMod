function init_objective_params()
    # This would also initialize the dual vars
    T3 = 387
    NoS = length(u)
    IOP = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
    # IOP = Model(optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => 2, "TimeLimit" => 200, "OutputFlag" => 1))
    @JuMP.variables(IOP, begin
        λ[1:NoS, 1:5] >= 0
        0 <= t[i = 1:2, j = 1:4] <= 10000
        Q[i = 1:4]
        dual1[1:2, 1:NoS] >= 0
        dual2[1:5, 1:NoS] >= 0
        dual3[1:NoS] <= 0
        dual4[1:NoS] <= 0
        w1[1:NoS]
        w2[1:NoS]
    end)
    # hat_Qc = Qc
    # hat_FH = FH
    @constraint(IOP, [s = 1:NoS], 2*Q[1]*hat_Qc[s]+Q[2]-0.5*λ[s, 1]-
    (a[s]-1)*λ[s, 2]+λ[s, 3]+λ[s, 4]-λ[s, 5] <= dual1[1, s])
    @constraint(IOP, [s = 1:NoS], -(2*Q[1]*hat_Qc[s]+Q[2]-0.5*λ[s, 1]-
    (a[s]-1)*λ[s, 2]+λ[s, 3]+λ[s, 4]-λ[s, 5]) <= dual1[1, s])

    @constraint(IOP, [s = 1:NoS], 2*Q[3]*hat_FH[s]+Q[4]-(u[s]-T3-170+
    b[s])*λ[s, 2]-(u[s]-393)*λ[s, 3]-(u[s]-313)*λ[s, 4]+(u[s]-323)*λ[s, 5] 
    <= dual1[2, s])
    @constraint(IOP, [s = 1:NoS], -(2*Q[3]*hat_FH[s]+Q[4]-(u[s]-T3-170+
    b[s])*λ[s, 2]-(u[s]-393)*λ[s, 3]-(u[s]-313)*λ[s, 4]+(u[s]-323)*λ[s, 5]) 
    <= dual1[2, s])

    @constraint(IOP, [s = 1:NoS], w1[s] == a[s]*λ[s, 2])
    @constraint(IOP, [s = 1:NoS], w2[s] == b[s]*λ[s, 2])

    @constraint(IOP, [s = 1:NoS], (hat_Qc[s]/2 + 553 - T3)*λ[s, 1] 
    <= dual2[1, s])
    @constraint(IOP, [s = 1:NoS], -((hat_Qc[s]/2 + 553 - T3)*λ[s, 1]) 
    <= dual2[1, s])

    @constraint(IOP, [s = 1:NoS], (-λ[s, 2]*10 - λ[s, 2]*hat_Qc[s] + 
    (u[s]- T3 - 170)*λ[s, 2]*hat_FH[s] + w1[s]*hat_Qc[s] + w2[s]*hat_FH[s]) 
    <= dual2[2, s])
    @constraint(IOP, [s = 1:NoS], -(-λ[s, 2]*10 - λ[s, 2]*hat_Qc[s] + 
    (u[s]- T3 - 170)*λ[s, 2]*hat_FH[s] + w1[s]*hat_Qc[s] + w2[s]*hat_FH[s]) 
    <= dual2[2, s])

    @constraint(IOP, [s = 1:NoS], (2*T3 - 786 - hat_Qc[s] + 
    (u[s] - 393)*hat_FH[s])*λ[s, 3] <= dual2[3, s])
    @constraint(IOP, [s = 1:NoS], -((2*T3 - 786 - hat_Qc[s] + 
    (u[s] - 393)*hat_FH[s])*λ[s, 3]) <= dual2[3, s])

    @constraint(IOP, [s = 1:NoS], ((2*T3 - 1026 - hat_Qc[s] + 
    (u[s] - 313)*hat_FH[s])*λ[s, 4]) <= dual2[4, s])
    @constraint(IOP, [s = 1:NoS], -((2*T3 - 1026 - hat_Qc[s] + 
    (u[s] - 313)*hat_FH[s])*λ[s, 4]) <= dual2[4, s])

    @constraint(IOP, [s = 1:NoS], ((2*T3 - 1026 - hat_Qc[s] + 
    (u[s] - 323)*hat_FH[s])*λ[s, 5]) <= dual2[5, s])
    @constraint(IOP, [s = 1:NoS], -((2*T3 - 1026 - hat_Qc[s] + 
    (u[s] - 323)*hat_FH[s])*λ[s, 5]) <= dual2[5, s])
    @constraint(IOP, Q[3] == 1)

    @objective(IOP, Min, sum(dual1) + sum(dual2) - sum(dual3) - sum(dual4))
    optimize!(IOP)
end

function init_constraint_params()
    T3 = 388
    Init = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
    @JuMP.variables(Init, begin
    p[1:2, 1:4]
    a[1:NoS]
    b[1:NoS]
    slack[1:NoS] <= 0
    end)

    @constraint(Init, [s = 1:NoS], -10 - hat_Qc[s] + (u[s]- T3 - 170)*hat_FH[s] + a[s]*hat_Qc[s] + b[s]*hat_FH[s] >= slack[s])
    @constraint(Init, [s = 1:NoS], a[s] == sum(p[1, j]*u[s]^(j-1) for j = 1:4))
    @constraint(Init, [s = 1:NoS], b[s] == sum(p[2, j]*u[s]^(j-1) for j = 1:4))

    @objective(Init, Max, sum(slack))
    optimize!(Init)
    return value.(p)
end