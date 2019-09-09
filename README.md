# TenBarOptimization
## Matlab Version
Introduction
Files in the folder :
1. FEM.m : Calculate the stress and displacement of this project by Finite Element Method.
2. Main.m : Optimize the project.
3. NonlinearConstrains.m : Constraints function for the Main.m file.
4. objective.m : Target function for the Main.m file.
5. StiffnessMatrix.m : Matrix forming function for the FEM.m file.

How to USE ?
>>
Open the Main.m file, and press the run button.

The optimization result will be like :
```
                                            First-order      Norm of
 Iter F-count            f(x)  Feasibility   optimality         step
    0       3    1.289128e+06    0.000e+00    1.882e+06
    1      11    1.112332e+06    0.000e+00    1.693e+06    7.036e-02
    2      15    6.798781e+05    0.000e+00    1.670e+06    3.395e-01
    3      23    6.328793e+05    0.000e+00    1.303e+06    9.307e-02
    4      31    5.899769e+05    0.000e+00    1.004e+06    6.194e-02
    5      40    2.855083e+05    0.000e+00    9.348e+04    2.032e-01
    6      55    2.749923e+05    0.000e+00    5.703e+04    1.330e-02
    7      59    2.231628e+05    0.000e+00    4.707e+05    4.573e-02
    8      80    2.226433e+05    0.000e+00    1.661e+05    2.236e-03
    9      86    2.224451e+05    0.000e+00    1.644e+05    3.332e-04
   10      97    2.219872e+05    0.000e+00    1.427e+05    2.214e-03
   11     110    2.213107e+05    0.000e+00    1.075e+05    3.825e-03
   12     115    2.210286e+05    0.000e+00    1.060e+05    4.190e-04
   13     126    2.205210e+05    0.000e+00    7.158e+04    3.700e-03
   14     130    2.125854e+05    0.000e+00    4.801e+04    1.080e-02
   15     133    2.124152e+05    0.000e+00    2.589e+03    6.789e-04
   16     136    2.124059e+05    8.858e-01    3.409e+02    1.860e-04
   17     139    2.124061e+05    0.000e+00    1.991e+00    8.942e-06
   18     142    2.124061e+05    0.000e+00    2.745e-02    1.217e-08

Local minimum found that satisfies the constraints.

Optimization completed because the objective function is non-decreasing in 
feasible directions, to within the default value of the optimality tolerance,
and constraints are satisfied to within the default value of the constraint tolerance.

<stopping criteria details>

x =

    0.3000    0.2663

fval =

   2.1241e+05

exitflag =

     1

output = 

  struct with fields:

         iterations: 18
          funcCount: 142
    constrviolation: 0
           stepsize: 1.2166e-08
          algorithm: 'interior-point'
      firstorderopt: 0.0275
       cgiterations: 28
            message: 'Local minimum found that satisfies the constraints.…'
```      
## Python Version
Introduction :
The structure is the same as the Matlab version, but all the functions are written in one file, Main.m.

How to USE ?
>>
Open the Main.py file, and press the run button.

The optimization result will be like :
```
Optimization terminated successfully.    (Exit mode 0)
            Current function value: 212406.03352685508
            Iterations: 20
            Function evaluations: 121
            Gradient evaluations: 18
     fun: 212406.03352685508
     jac: array([812557.50976562, 679880.109375  ])
 message: 'Optimization terminated successfully.'
    nfev: 121
     nit: 20
    njev: 18
  status: 0
 success: True
       x: array([0.30002294, 0.26626172])
```
