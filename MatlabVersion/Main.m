clear;
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0.1,0.1];
ub = [0.75,0.75];
nonlcon = @NonlinearConstrains;
fun = @objective;
x0 = [0.7 , 0.7];
options = optimset('display', 'iter');
[x, fval, exitflag, output] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, nonlcon, options)