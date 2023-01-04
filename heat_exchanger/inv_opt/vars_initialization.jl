function init_objective_params()
    # This would also initialize the dual vars
    T3 = 388
    NoS = length(u)
    IOP = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
    # IOP = Model(optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => 2, "TimeLimit" => 200, "OutputFlag" => 1))
    @JuMP.variables(IOP, begin
        λ[1:NoS, 1:5] >= 0
        Q[i = 1:4]
        slack1[1:2, 1:NoS] >= 0
        slack2[1:5, 1:NoS] >= 0
        w1[1:NoS]
        w2[1:NoS]
    end)
    @constraint(IOP, [s = 1:NoS], 2*Q[1]*x[s][1]+Q[2]-0.5*λ[s, 1]-
    (a[s]-1)*λ[s, 2]+λ[s, 3]+λ[s, 4]-λ[s, 5] <= slack1[1, s])
    @constraint(IOP, [s = 1:NoS], -(2*Q[1]*x[s][1]+Q[2]-0.5*λ[s, 1]-
    (a[s]-1)*λ[s, 2]+λ[s, 3]+λ[s, 4]-λ[s, 5]) <= slack1[1, s])

    @constraint(IOP, [s = 1:NoS], 2*Q[3]*x[s][2]+Q[4]-(u[s]-T3-170+
    b[s])*λ[s, 2]-(u[s]-393)*λ[s, 3]-(u[s]-313)*λ[s, 4]+(u[s]-323)*λ[s, 5] 
    <= slack1[2, s])
    @constraint(IOP, [s = 1:NoS], -(2*Q[3]*x[s][2]+Q[4]-(u[s]-T3-170+
    b[s])*λ[s, 2]-(u[s]-393)*λ[s, 3]-(u[s]-313)*λ[s, 4]+(u[s]-323)*λ[s, 5]) 
    <= slack1[2, s])
    @constraint(IOP, [s = 1:NoS], a[s] == sum(p[1, j]*u[s]^(j-1) for j = 1:4))
    @constraint(IOP, [s = 1:NoS], b[s] == sum(p[2, j]*u[s]^(j-1) for j = 1:4))

    @constraint(IOP, [s = 1:NoS], w1[s] == a[s]*λ[s, 2])
    @constraint(IOP, [s = 1:NoS], w2[s] == b[s]*λ[s, 2])

    @constraint(IOP, [s = 1:NoS], (x[s][1]/2 + 553 - T3)*λ[s, 1] 
    <= slack2[1, s])
    @constraint(IOP, [s = 1:NoS], -((x[s][1]/2 + 553 - T3)*λ[s, 1]) 
    <= slack2[1, s])

    @constraint(IOP, [s = 1:NoS], (-λ[s, 2]*10 - λ[s, 2]*x[s][1] + 
    (u[s]- T3 - 170)*λ[s, 2]*x[s][2] + w1[s]*x[s][1] + w2[s]*x[s][2]) 
    <= slack2[2, s])
    @constraint(IOP, [s = 1:NoS], -(-λ[s, 2]*10 - λ[s, 2]*x[s][1] + 
    (u[s]- T3 - 170)*λ[s, 2]*x[s][2] + w1[s]*x[s][1] + w2[s]*x[s][2]) 
    <= slack2[2, s])

    @constraint(IOP, [s = 1:NoS], (2*T3 - 786 - x[s][1] + 
    (u[s] - 393)*x[s][2])*λ[s, 3] <= slack2[3, s])
    @constraint(IOP, [s = 1:NoS], -((2*T3 - 786 - x[s][1] + 
    (u[s] - 393)*x[s][2])*λ[s, 3]) <= slack2[3, s])

    @constraint(IOP, [s = 1:NoS], ((2*T3 - 1026 - x[s][1] + 
    (u[s] - 313)*x[s][2])*λ[s, 4]) <= slack2[4, s])
    @constraint(IOP, [s = 1:NoS], -((2*T3 - 1026 - x[s][1] + 
    (u[s] - 313)*x[s][2])*λ[s, 4]) <= slack2[4, s])

    @constraint(IOP, [s = 1:NoS], ((2*T3 - 1026 - x[s][1] + 
    (u[s] - 323)*x[s][2])*λ[s, 5]) <= slack2[5, s])
    @constraint(IOP, [s = 1:NoS], -((2*T3 - 1026 - x[s][1] + 
    (u[s] - 323)*x[s][2])*λ[s, 5]) <= slack2[5, s])
    @constraint(IOP, Q[3] == 1)

    @objective(IOP, Min, sum(slack1) + sum(slack2))
    optimize!(IOP)
    return value.(Q), value.(λ)
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

    @constraint(Init, [s = 1:NoS], -10 - x[s][1] + (u[s]- T3 - 170)*x[s][2] + a[s]*x[s][1] + b[s]*x[s][2] >= slack[s])
    @constraint(Init, [s = 1:NoS], a[s] == sum(p[1, j]*u[s]^(j-1) for j = 1:4))
    @constraint(Init, [s = 1:NoS], b[s] == sum(p[2, j]*u[s]^(j-1) for j = 1:4))

    @objective(Init, Max, sum(slack))
    optimize!(Init)
    return value.(p)
end