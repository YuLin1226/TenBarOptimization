function [c , ceq] = NonlinearConstrains(r)
    [stress , Q] = FEM(r);
    yield_stress = 250e6;
    stress_con = abs(stress) - yield_stress;
    dis_con = sqrt(Q(3)^2 + Q(4)^2) - 0.02;
    c = [stress_con ; dis_con];
    ceq = [];
end