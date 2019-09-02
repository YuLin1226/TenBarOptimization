clear;
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0,0];
ub = [];
nonlcon = @NonlinearConstrains;
fun = @objective;
x0 = [0.5 , 0.5];
[x, fval, exitflag, output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon)