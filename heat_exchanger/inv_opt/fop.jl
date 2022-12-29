function forward_problem(T5)
    T3 = 388
    FP = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0,
     "NonConvex" => 2))
    @JuMP.variables(FP, begin
    Qc >= 0
    FH >= 0
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
