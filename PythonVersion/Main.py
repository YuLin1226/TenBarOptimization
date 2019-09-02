import math
import numpy as np
import scipy.optimize as sco

def StiffnessMatrix(angle):
    C = math.cos(angle)
    S = math.sin(angle)
    SM = np.array([
            [C*C , C*S , -C*C , -C*S],
            [C*S , S*S , -C*S , -S*S],
            [-C*C , -C*S , C*C , C*S],
            [-C*S , -S*S , C*S , S*S]
            ])
    return SM

def FEM(r):
    Beam = np.array([1,2,3,4,5,6,7,8,9,10])
    Element = np.array([[3,5] , [1,3] , [4,6] , [2,4] , [3,4] , [1,2] , [4,5] , [3,6] , [2,3] , [1,4]] )
    Node = np.array([1,2,3,4,5,6])
    Angle = np.array([math.radians(0),math.radians(0),math.radians(0),
                      math.radians(0),math.radians(90),math.radians(90),
                      math.radians(135),math.radians(45),math.radians(135),math.radians(45),])
    Force = np.array([0,0,0,-1e7,0,0,0,-1e7,0,0,0,0])
    Force = Force.transpose()
    L1 = np.ones(6)
    L2 = np.ones(4)*math.sqrt(2)
    Length = np.hstack((L1,L2))*9.14  #9.14 = element unit length
    E = 200e9
    A1 = np.ones(6)*r[0]**2
    A2 = np.ones(4)*r[1]**2
    Area = np.hstack((A1,A2))*math.pi
    
    K = np.zeros((2*len(Node),2*len(Node)))
    index = np.zeros(4)
    for n in range(len(Beam)):
        index[0] = 2*Element[n][0]-2
        index[1] = 2*Element[n][0]-1
        index[2] = 2*Element[n][1]-2
        index[3] = 2*Element[n][1]-1
        SM = StiffnessMatrix(Angle[n])*E*Area[n]/Length[n]
        for i in range(4):
            for j in range(4):
                Matrix = np.zeros( ( 2*len(Node),2*len(Node) ) )
                row = int(index[i])
                column = int(index[j])
                Matrix[row][column] = SM[i][j] #array index number must be "int"
                K += Matrix
    BC_node = np.array([5,6])
    
    Ks = K
    
#    To delete the elements in arrays, it can be done in another way like K = K[0:8 , 0:8]
#    ,which means we choose the 0~8 rows in K and 0~8 columns in K.
    K = np.delete(K , [ 2*BC_node[1]-2 , 2*BC_node[1]-1 ] , axis = 0) #axis = 0 : delete the row
    K = np.delete(K , [ 2*BC_node[0]-2 , 2*BC_node[0]-1 ] , axis = 0)
    K = np.delete(K , [ 2*BC_node[1]-2 , 2*BC_node[1]-1 ] , axis = 1) #axis = 1 : delete the column
    K = np.delete(K , [ 2*BC_node[0]-2 , 2*BC_node[0]-1 ] , axis = 1)
    
    Force = np.delete(Force , [ 2*BC_node[1]-2 , 2*BC_node[1]-1 ] , axis = 0)
    Force = np.delete(Force , [ 2*BC_node[0]-2 , 2*BC_node[0]-1 ] , axis = 0)
    iK = np.linalg.inv(K)
    Q = iK.dot(Force)
    Q = np.hstack(( Q , np.zeros(4) ))
    stress = np.zeros((len(Beam) , 1))
    for i in range(len(Beam)):
        S1 = np.array([ -math.cos(Angle[i]) , -math.sin(Angle[i]) , math.cos(Angle[i]) , math.sin(Angle[i]) ])
        D = np.array([[ Q[2*Element[i][0]-2] ],
                      [ Q[2*Element[i][0]-1] ],
                      [ Q[2*Element[i][1]-2] ],
                      [ Q[2*Element[i][1]-1] ] ])
        stress[i] = E/Length[i]*S1.dot(D)
#    method 1 to delete elements in 2D arrays        
#    for i in range(2*BC_node[0] - 2):
#        Ks = np.delete(Ks , 2*BC_node[0] - 3 - i ,axis = 0)
    
#    method 2 to slice the arrays
    Ks = Ks[8:,0:]
    print(Ks.shape)
    R = Ks.dot(Q)
    return Q,stress

def cons_d2(r):
    Q,_ = FEM(r)
    con_d2 = 0.02 - math.sqrt( Q[2]**2 + Q[3]**2 ) 
    return con_d2

def cons_stress(r):
    _,stress = FEM(r)
    con_stress = (250e6 - np.absolute(stress))/1e+8
    con_stress = con_stress.flatten().tolist()
    return con_stress

def fun(r):
    weight = (6*r[0]**2 + 4*r[1]**2*math.sqrt(2) )*math.pi*7860*9.14
    return weight

if __name__ == "__main__":
    r0 = [0.1,0.1]
    cons = [{'type': 'ineq','fun': cons_d2},
            {'type': 'ineq','fun': cons_stress}]
    
#    type >> decides types of constraints
#    fun  >> calls constraints functions, hence constraints must be the callable functions.
    
    bnds = ((0.001 , None),(0.001 , None))
    res = sco.minimize(fun , r0 , method = 'SLSQP' ,
                       constraints = cons,
                       bounds = bnds,
                       options={'disp': True})
    
#    The optimization sometimes fails due to the scales of constraints.
#    e.g. 
#    In this case, the scale of D2 constraints equals to 0.001,
#    but the one of stress constraints equals to 1e+6.
#    Fortunately, it can be solved by using multipication to change the scales.
    
    
    result = print(res)


        
        
        
        
        
        
        
        
        
