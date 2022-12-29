using JuMP, Ipopt, LinearAlgebra, Random, BARON, Gurobi, PyPlot

function ForwardProblem(T3, T5)
    # FP = Model(optimizer_with_attributes(BARON.Optimizer, "PrLevel" => 0, 
    # "CplexLibName" => "C:\\Program Files\\IBM\\ILOG\\CPLEX_Studio128\\cplex\\bin\\x64_win64\\cplex1280.dll"))
    FP = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0, "NonConvex" => 2))
    @JuMP.variables(FP, begin
    Qc
    FH
    end)

    @constraint(FP, Qc/2 + 553 - T3 >= 0)
    @constraint(FP, -10 - Qc + (T5- T3 - 170 + 0.5*Qc)*FH >= 0)
    @constraint(FP, 2*T3 - 786 - Qc + (T5 - 393)*FH >= 0)
    @constraint(FP, 2*T3 - 1026 - Qc + (T5 - 313)*FH >= 0)
    @constraint(FP, 2*T3 - 1026 - Qc + (T5 - 323)*FH <= 0)
    @objective(FP, Min, 1e-2*Qc + 4*(FH - 1.7)^2)
    optimize!(FP)
    return value.(Qc), value.(FH)
end

function BoundingProblem(T3, T5, guess)
    # FP = Model(optimizer_with_attributes(BARON.Optimizer, "PrLevel" => 0))
    FP = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
    @JuMP.variables(FP, begin
    Qc >= 0, (start = guess)
    FH >= 0
    end)

    @constraint(FP, Qc/2 + 553 - T3 >= 0)
    @NLconstraint(FP, -10 - Qc + (T5- T3 - 170 + 0.5*Qc)*FH >= 0)
    @constraint(FP, 2*T3 - 786 - Qc + (T5 - 393)*FH >= 0)
    @constraint(FP, 2*T3 - 1026 - Qc + (T5 - 313)*FH >= 0)
    @constraint(FP, 2*T3 - 1026 - Qc + (T5 - 323)*FH <= 0)
    @objective(FP, Min, Qc)
    optimize!(FP)
    lb = objective_value(FP)
    @objective(FP, Max, Qc)
    optimize!(FP)
    ub = objective_value(FP)
    return lb, ub
end

function SurrogateBoundingProblem(Q, p, T5)
    T3 = 388
    SFP = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
    a = sum(p[1, j]*(T5)^(j-1) for j = 1:4)
    b = sum(p[2, j]*(T5)^(j-1) for j = 1:4)

    @JuMP.variables(SFP, begin
    Qc >= 0
    FH >= 0
    end)

    @constraint(SFP, Qc/2 + 553 - T3 >= 0)
    @constraint(SFP, -10 - Qc + (T5- T3 - 170)*FH + a*Qc + b*FH >= 0)
    @constraint(SFP, 2*T3 - 786 - Qc + (T5 - 393)*FH >= 0)
    @constraint(SFP, 2*T3 - 1026 - Qc + (T5 - 313)*FH >= 0)
    @constraint(SFP, 2*T3 - 1026 - Qc + (T5 - 323)*FH <= 0)
    @objective(SFP, Min, Qc)
    optimize!(SFP)
    lb = objective_value(SFP)
    @objective(SFP, Max, Qc)
    optimize!(SFP)
    ub = objective_value(SFP)
    return lb, ub
end

function SurrogateFP(Q, p, T5)
    SFP = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
    T3 = 388
    a = sum(p[1, j]*(T5)^(j-1) for j = 1:4)
    b = sum(p[2, j]*(T5)^(j-1) for j = 1:4)

    @JuMP.variables(SFP, begin
    Qc >= 0
    FH >= 0
    end)

    @constraint(SFP, Qc/2 + 553 - T3 >= 0)
    @constraint(SFP, -10 - Qc + (T5- T3 - 170)*FH + a*Qc + b*FH >= 0)
    @constraint(SFP, 2*T3 - 786 - Qc + (T5 - 393)*FH >= 0)
    @constraint(SFP, 2*T3 - 1026 - Qc + (T5 - 313)*FH >= 0)
    @constraint(SFP, 2*T3 - 1026 - Qc + (T5 - 323)*FH <= 0)
    @objective(SFP, Min, Q[1]*Qc^2 + Q[2]*Qc + Q[3]*FH^2 + Q[4]*FH)
    optimize!(SFP)
    return value.(Qc), value.(FH)

end

function InitializeP(u, hat_Qc, hat_FH)
    T3 = 388
    NoS = length(u)
    Init = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
    @JuMP.variables(Init, begin
    p[1:2, 1:4]
    a[1:NoS]
    b[1:NoS]
    dual3[1:NoS] <= 0
    end)

    @constraint(Init, [s = 1:NoS], -10 - hat_Qc[s] + (u[s]- T3 - 170)*hat_FH[s] + a[s]*hat_Qc[s] + b[s]*hat_FH[s] >= dual3[s])
    @constraint(Init, [s = 1:NoS], a[s] == sum(p[1, j]*u[s]^(j-1) for j = 1:4))
    @constraint(Init, [s = 1:NoS], b[s] == sum(p[2, j]*u[s]^(j-1) for j = 1:4))

    @objective(Init, Max, sum(dual3))
    optimize!(Init)
    return value.(p)
end

function InitializeQ(T5, hat_Qc, hat_FH, p)
    T3 = 388
    NoS = length(u)
    Init = Model(optimizer_with_attributes(Ipopt.Optimizer))
    @JuMP.variables(Init, begin
    p[1:2, 1:4]
    a[1:NoS]
    b[1:NoS]
    dual3[1:NoS] <= 0
    end)

end

function SeparationProblem(Q, p)
    T3 = 388
    SP = Model(optimizer_with_attributes(BARON.Optimizer, "PrLevel" => 1, 
    "CplexLibName" => "C:\\Program Files\\IBM\\ILOG\\CPLEX_Studio128\\cplex\\bin\\x64_win64\\cplex1280.dll", "MaxTime" => 200))
    @JuMP.variables(SP, begin
        λ[1, 1:5] >= 0
        700>= hat_Qc[1] >=0
        6 >= hat_FH[1] >= 0
        a[1]
        b[1]
        593>=u[1]>=586, (start = 592)
        -1000 <= delta <= 0
    end)

    @NLconstraint(SP, [s = 1], 2*Q[1]*hat_Qc[s]+Q[2]-0.5*λ[s, 1]-(a[s]-1)*λ[s, 2]+λ[s, 3]+λ[s, 4]-λ[s, 5] == 0)
    # @NLconstraint(IOP, [s = 1:NoS], -(2*Q[1]*hat_Qc[s]+Q[2]-0.5*λ[s, 1]-(a[s]-1)*λ[s, 2]+λ[s, 3]+λ[s, 4]-λ[s, 5]) <= dual1[1, s])

    @NLconstraint(SP, [s = 1], 2*Q[3]*hat_FH[s]+Q[4]-(u[s]-T3-170+b[s])*λ[s, 2]-(u[s]-393)*λ[s, 3]-(u[s]-313)*λ[s, 4]+(u[s]-323)*λ[s, 5] == 0)
    # @NLconstraint(IOP, [s = 1:NoS], -(2*Q[3]*hat_FH[s]+Q[4]-(u[s]-T3-170+b[s])*λ[s, 2]-(u[s]-393)*λ[s, 3]-(u[s]-313)*λ[s, 4]+(u[s]-323)*λ[s, 5]) <= dual1[2, s])

    @constraint(SP, [s = 1], hat_Qc[s]/2 + 553 - T3 >= 0)
    @NLconstraint(SP, [s = 1], -10 - hat_Qc[s] + (u[s]- T3 - 170)*hat_FH[s] + a[s]*hat_Qc[s] + b[s]*hat_FH[s] >= 0)
    @constraint(SP, [s = 1], 2*T3 - 786 - hat_Qc[s] + (u[s] - 393)*hat_FH[s] >= 0)
    @constraint(SP, [s = 1], 2*T3 - 1026 - hat_Qc[s] + (u[s] - 313)*hat_FH[s] >= 0)
    @constraint(SP, [s = 1], 2*T3 - 1026 - hat_Qc[s] + (u[s] - 323)*hat_FH[s] <= 0)

    @NLconstraint(SP, [s = 1], (hat_Qc[s]/2 + 553 - T3)*λ[s, 1] == 0)
    # @NLconstraint(IOP, [s = 1], -((hat_Qc[s]/2 + 553 - T3)*λ[s, 1]) <= dual2[1, s])

    @NLconstraint(SP, [s = 1], (λ[s, 2]*(-10 - hat_Qc[s] + (u[s]- T3 - 170)*hat_FH[s] + a[s]*hat_Qc[s] + b[s]*hat_FH[s])) == 0)
    # @NLconstraint(IOP, [s = 1:NoS], -(λ[s, 2]*(-10 - hat_Qc[s] + (u[s]- T3 - 170)*hat_FH[s] + a[s]*hat_Qc[s] + b[s]*hat_FH[s])) <= dual2[2, s])

    @NLconstraint(SP, [s = 1], (2*T3 - 786 - hat_Qc[s] + (u[s] - 393)*hat_FH[s])*λ[s, 3] == 0)
    # @NLconstraint(IOP, [s = 1:NoS], -((2*T3 - 786 - hat_Qc[s] + (u[s] - 393)*hat_FH[s])*λ[s, 3]) <= dual2[3, s])

    @NLconstraint(SP, [s = 1], ((2*T3 - 1026 - hat_Qc[s] + (u[s] - 313)*hat_FH[s])*λ[s, 4]) == 0)
    # @NLconstraint(IOP, [s = 1:NoS], -((2*T3 - 1026 - hat_Qc[s] + (u[s] - 313)*hat_FH[s])*λ[s, 4]) <= dual2[4, s])

    @NLconstraint(SP, [s = 1], ((2*T3 - 1026 - hat_Qc[s] + (u[s] - 323)*hat_FH[s])*λ[s, 5]) == 0)
    # @NLconstraint(IOP, [s = 1:NoS], -((2*T3 - 1026 - hat_Qc[s] + (u[s] - 323)*hat_FH[s])*λ[s, 5]) <= dual2[5, s])

    @NLconstraint(SP, [s = 1], a[s] == sum(p[1, j]*u[s]^(j-1) for j = 1:4))
    @NLconstraint(SP, [s = 1], b[s] == sum(p[2, j]*u[s]^(j-1) for j = 1:4))

    # @NLconstraint(SP,  >= delta)

    @objective(SP, Min, -10 - hat_Qc[1] + (u[1] - T3 - 170 + 0.5*hat_Qc[1])*hat_FH[1])
    optimize!(SP)
    return value.(u[1]), objective_value(SP)

end

function InverseProblem(u, Qc, FH, init_p, count, init_Q)
    T3 = 388
    NoS = length(u)
    IOP = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
    # IOP = Model(optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => 2, "TimeLimit" => 200, "OutputFlag" => 1))
    if count <= 10
        @JuMP.variables(IOP, begin
            λ[1:NoS, 1:5] >= 0
            hat_Qc[i = 1:NoS], (start = Qc[i])
            hat_FH[i = 1:NoS], (start = FH[i])
            a[1:NoS]
            b[1:NoS]
            p[i = 1:2, j = 1:4]
            0 <= t[i = 1:2, j = 1:4] <= 10000
            Q[i = 1:4]
            dual1[1:2, 1:NoS] >= 0
            dual2[1:5, 1:NoS] >= 0
            dual3[1:NoS] <= 0
            dual4[1:NoS] <= 0
            w1[1:NoS]
            w2[1:NoS]
        end)
    else
        @JuMP.variables(IOP, begin
        λ[1:NoS, 1:5] >= 0
        hat_Qc[i = 1:NoS], (start = Qc[i])
        hat_FH[i = 1:NoS], (start = FH[i])
        a[1:NoS]
        b[1:NoS]
        p[i = 1:2, j = 1:4] , (start = init_p[i, j])
        0 <= t[i = 1:2, j = 1:4] <= 10000
        Q[i = 1:4], (start = init_Q[i])
        dual1[1:2, 1:NoS] >= 0
        dual2[1:5, 1:NoS] >= 0
        dual3[1:NoS] <= 0
        dual4[1:NoS] <= 0
        w1[1:NoS]
        w2[1:NoS]
        end)
    end
    # hat_Qc = Qc
    # hat_FH = FH
    @constraint(IOP, [s = 1:NoS], 2*Q[1]*hat_Qc[s]+Q[2]-0.5*λ[s, 1]-(a[s]-1)*λ[s, 2]+λ[s, 3]+λ[s, 4]-λ[s, 5] <= dual1[1, s])
    @constraint(IOP, [s = 1:NoS], -(2*Q[1]*hat_Qc[s]+Q[2]-0.5*λ[s, 1]-(a[s]-1)*λ[s, 2]+λ[s, 3]+λ[s, 4]-λ[s, 5]) <= dual1[1, s])

    @constraint(IOP, [s = 1:NoS], 2*Q[3]*hat_FH[s]+Q[4]-(u[s]-T3-170+b[s])*λ[s, 2]-(u[s]-393)*λ[s, 3]-(u[s]-313)*λ[s, 4]+(u[s]-323)*λ[s, 5] <= dual1[2, s])
    @constraint(IOP, [s = 1:NoS], -(2*Q[3]*hat_FH[s]+Q[4]-(u[s]-T3-170+b[s])*λ[s, 2]-(u[s]-393)*λ[s, 3]-(u[s]-313)*λ[s, 4]+(u[s]-323)*λ[s, 5]) <= dual1[2, s])

    @constraint(IOP, [s = 1:NoS], hat_Qc[s]/2 + 553 - T3 >= 0)
    @constraint(IOP, [s = 1:NoS], -10 - hat_Qc[s] + (u[s]- T3 - 170)*hat_FH[s] + a[s]*hat_Qc[s] + b[s]*hat_FH[s] >= dual3[s])
    @constraint(IOP, [s = 1:NoS], 2*T3 - 786 - hat_Qc[s] + (u[s] - 393)*hat_FH[s] >= 0)
    @constraint(IOP, [s = 1:NoS], 2*T3 - 1026 - hat_Qc[s] + (u[s] - 313)*hat_FH[s] >= 0)
    @constraint(IOP, [s = 1:NoS], 2*T3 - 1026 - hat_Qc[s] + (u[s] - 323)*hat_FH[s] <= 0)

    @constraint(IOP, [s = 1:NoS], a[s] == sum(p[1, j]*u[s]^(j-1) for j = 1:4))
    @constraint(IOP, [s = 1:NoS], b[s] == sum(p[2, j]*u[s]^(j-1) for j = 1:4))

    @constraint(IOP, [s = 1:NoS], w1[s] == a[s]*λ[s, 2])
    @constraint(IOP, [s = 1:NoS], w2[s] == b[s]*λ[s, 2])

    @constraint(IOP, [s = 1:NoS], (hat_Qc[s]/2 + 553 - T3)*λ[s, 1] <= dual2[1, s])
    @constraint(IOP, [s = 1:NoS], -((hat_Qc[s]/2 + 553 - T3)*λ[s, 1]) <= dual2[1, s])

    @constraint(IOP, [s = 1:NoS], (-λ[s, 2]*10 - λ[s, 2]*hat_Qc[s] + (u[s]- T3 - 170)*λ[s, 2]*hat_FH[s] + w1[s]*hat_Qc[s] + w2[s]*hat_FH[s]) <= dual2[2, s])
    @constraint(IOP, [s = 1:NoS], -(-λ[s, 2]*10 - λ[s, 2]*hat_Qc[s] + (u[s]- T3 - 170)*λ[s, 2]*hat_FH[s] + w1[s]*hat_Qc[s] + w2[s]*hat_FH[s]) <= dual2[2, s])

    @constraint(IOP, [s = 1:NoS], (2*T3 - 786 - hat_Qc[s] + (u[s] - 393)*hat_FH[s])*λ[s, 3] <= dual2[3, s])
    @constraint(IOP, [s = 1:NoS], -((2*T3 - 786 - hat_Qc[s] + (u[s] - 393)*hat_FH[s])*λ[s, 3]) <= dual2[3, s])

    @constraint(IOP, [s = 1:NoS], ((2*T3 - 1026 - hat_Qc[s] + (u[s] - 313)*hat_FH[s])*λ[s, 4]) <= dual2[4, s])
    @constraint(IOP, [s = 1:NoS], -((2*T3 - 1026 - hat_Qc[s] + (u[s] - 313)*hat_FH[s])*λ[s, 4]) <= dual2[4, s])

    @constraint(IOP, [s = 1:NoS], ((2*T3 - 1026 - hat_Qc[s] + (u[s] - 323)*hat_FH[s])*λ[s, 5]) <= dual2[5, s])
    @constraint(IOP, [s = 1:NoS], -((2*T3 - 1026 - hat_Qc[s] + (u[s] - 323)*hat_FH[s])*λ[s, 5]) <= dual2[5, s])
    if count >= 2
        @constraint(IOP, [s = NoS-count+2:NoS], -10 - hat_Qc[s] + (u[s] - T3 - 170 + 0.5*hat_Qc[s])*hat_FH[s] >= dual4[s])
    end
    @constraint(IOP, [i = 1:2, j = 1:4], p[i, j] <= t[i, j])
    @constraint(IOP, [i = 1:2, j = 1:4], -p[i, j] <= t[i, j])
    @constraint(IOP, Q[3] == 1)

    @objective(IOP, Min, 1000*sum((hat_Qc[s] - Qc[s])^2 for s=1:NoS) + sum((hat_FH[s] - FH[s])^2 for s=1:NoS)
    + (sum(dual1) + sum(dual2) - sum(dual3) - sum(dual4)) + 1e-5*sum(t))
    optimize!(IOP)
    # Qc = value.(hat_Qc)
    # FH = value.(hat_FH)
    # init_p = value.(p)
    # init_a = value.(a)
    # init_b = value.(b)
    # init_Q = value.(Q)
    # init_λ = value.(λ)

    # println(Qc)
    # println(FH)
    # set_optimizer(IOP, Ipopt.Optimizer)

    # @JuMP.variables(IOP, begin
    # λ[i = 1:NoS, j = 1:5] >= 0, (start = init_λ[i, j])
    # hat_Qc[i = 1:NoS], (start = Qc[i])
    # hat_FH[i = 1:NoS], (start = FH[i])
    # a[i = 1:NoS], (start = init_a[i])
    # b[i = 1:NoS], (start = init_b[i])
    # p[i = 1:2, j = 1:4], (start = init_p[i, j])
    # Q[i = 1:4], (start = init_Q[i])
    # dual1[1:2, 1:NoS] >= 0
    # dual2[1:5, 1:NoS] >= 0
    # dual3[1:NoS] <= 0
    # dual4 <= 0
    # w1[1:NoS]
    # w2[1:NoS]
    # end)
    # optimize!(IOP)
    return value.(Q), value.(p), value.(dual1), value.(dual2), value.(dual3)
end

I = 20
rng = MersenneTwister(1234)
Temp_range = 586:0.2:593 # Finer the grid, better the feasible region plot. 
T3 = 388

# T3 = rand(rng, -100:100, I)/10 .+ 388
# Training dataset
T5 = rand(rng, 10:40, I)/10 .+ 586
# T5[end] = 585
# println(T5)
global p, Q
# Test dataset
T5_Validation = deepcopy(Temp_range)
true_Qc = zeros(Float64, length(T5_Validation))
true_FH = zeros(Float64, length(T5_Validation))

for i in 1:length(T5_Validation)
    true_Qc[i], true_FH[i] = ForwardProblem(T3, T5_Validation[i])
end
Pred_Error = []
Infeasibility_Measure = []
Inf_Temp = []
for iter in 1:1
    global p, Q
    # Plot the feasible region in the background
    init = 600
    min = zeros(Float64, length(Temp_range))
    max = zeros(Float64, length(Temp_range))
    i = 0
    for T in Temp_range
        i += 1
        min[i], max[i] = BoundingProblem(T3, T, init)
    end
    figure()
    PyPlot.fill_between(Temp_range, min, max, facecolor = "mistyrose")

    init = 0
    min = zeros(Float64, length(Temp_range))
    max = zeros(Float64, length(Temp_range))
    i = 0
    for T in Temp_range
        i += 1
        min[i], max[i] = BoundingProblem(T3, T, init)
    end
    PyPlot.fill_between(Temp_range, min, max, facecolor = "mistyrose")

    Opt_Qc = zeros(Float64, length(T5))
    Opt_FH = zeros(Float64, length(T5))

    for i in 1:length(T5)
        Opt_Qc[i], Opt_FH[i] = ForwardProblem(T3, T5[i])
    end

    # q = InitializeQ(T5, Opt_Qc, Opt_FH, p)
    # # Initialize Parameters
    if iter == 1
        p = InitializeP(T5, Opt_Qc, Opt_FH)
        Q = [1, 1, 1, 1]
    end
    Q, p, d1, d2, d3 = InverseProblem(T5, Opt_Qc, Opt_FH, p, iter, Q)
    println(Q)
    println(p)

    min = zeros(Float64, length(Temp_range))
    max = zeros(Float64, length(Temp_range))
    i = 0
    for T5 in Temp_range
        i += 1
        min[i], max[i] = SurrogateBoundingProblem(Q, p, T5)
    end
    PyPlot.fill_between(Temp_range, max, min, facecolor = "blue", alpha = 0.2)
    xlim([minimum(Temp_range), maximum(Temp_range)])
    ylim([0, 700])
    labels = minimum(Temp_range):4:maximum(Temp_range)
    xticks(labels, string.(labels))
    gcf()
    # legend(fancybox = true, framealpha = 0.5)
    # savefig(string("Qc_",iter,".png"))
    savefig(string("FeasibleRegion_",iter,".png"))

    # figure()
    # Perform prediciton error analysis
    bar_Qc = zeros(Float64, length(T5))
    bar_FH = zeros(Float64, length(T5))

    # Check accuracy
    for i in 1:length(T5)
        bar_Qc[i], bar_FH[i] = SurrogateFP(Q, p, T5[i])
    end

    push!(Pred_Error, (1/length(T5))*norm((Opt_Qc - bar_Qc), 1))

    PyPlot.scatter(T5, bar_Qc, color = "blue", marker = "s", label = "Surrogate Model")
    PyPlot.scatter(T5, Opt_Qc, color = "firebrick", marker = "o", label = "True Model")
    # PyPlot.scatter(T5, 0, color = "black")
    println(T5)
    xlim([minimum(Temp_range), maximum(Temp_range)])
    ylim([0, 700])
    labels = minimum(Temp_range):4:maximum(Temp_range)
    xticks(labels, string.(labels), fontsize = 14)
    yticks(fontsize = 14)
    gcf()
    legend(fancybox = true, framealpha = 0.5, fontsize = 14)
    xlabel(L"T_5 \, (K)", fontsize = 18)
    ylabel(L"Q_c \, (kW)", fontsize = 18)
    savefig(string("Test_Qc_",iter,".png"), bbox_inches="tight")
    # savefig(string("FeasibleRegion_",iter,".png"))

    temp, infeasibility = SeparationProblem(Q, p)
    push!(Inf_Temp, temp)
    push!(Infeasibility_Measure, infeasibility)
    println("Infeasibility is", infeasibility)
    if infeasibility >= -1
        break
    end
    push!(T5, temp)
end
