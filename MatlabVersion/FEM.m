function [stress , Q] = FEM(r)
%% Parameter 
Beam = [1,2,3,4,5,6,7,8,9,10];
Element = [3,5 ; 1,3 ; 4,6 ; 2,4 ; 3,4 ; 1,2 ; 4,5 ; 3,6 ; 2,3 ; 1,4];
Node = [1,2,3,4,5,6];

Force = zeros( 2*length(Node) , 1);         % Unit: N
Force(4) = -1e7;
Force(8) = -1e7;

angle_rad = zeros( length(Beam) , 1);       % Unit: radius
angle_rad(1) = deg2rad( 0 );
angle_rad(2) = deg2rad( 0 );
angle_rad(3) = deg2rad( 0 );
angle_rad(4) = deg2rad( 0 );
angle_rad(5) = deg2rad( 90 );
angle_rad(6) = deg2rad( 90 );
angle_rad(7) = deg2rad( 135 );
angle_rad(8) = deg2rad( 45 );
angle_rad(9) = deg2rad( 135 );
angle_rad(10) = deg2rad( 45 );

L = [ ones(6,1) ; ones(4,1)*sqrt(2) ]*9.14; % Unit: meter

E = 200e9;      % Unit: Pa
Area = pi*[ones(6,1)*r(1)^2 ; ones(4,1)*r(2)^2];    % Unit: m^2

%% Stiffness Matrix
K = zeros( 2*length(Node) , 2*length(Node) );
for n = 1 : length(Beam)
    row = [2*Element(n,1)-1 ; 2*Element(n,1) ; 2*Element(n,2)-1 ; 2*Element(n,2)];
    column = [2*Element(n,1)-1 ; 2*Element(n,1) ; 2*Element(n,2)-1 ; 2*Element(n,2)];
    S = StiffnessMatrix( angle_rad(n) ) * E * Area(n) / L(n);
    for i = 1 : 4
        for j = 1:4
            Matrix = zeros( 2*length(Node) , 2*length(Node) );
            Matrix( row(i) , column(j) ) = S(i,j);
            K = K + Matrix;
        end    
    end
end
%% Boundary Conditions Applied
BC_Node = [5,6];
Ks = K;
K( 2*BC_Node(2)-1 : 2*BC_Node(2) , : ) = [];
K( 2*BC_Node(1)-1 : 2*BC_Node(1) , : ) = [];
K( : , 2*BC_Node(2)-1 : 2*BC_Node(2) ) = [];
K( : , 2*BC_Node(1)-1 : 2*BC_Node(1) ) = [];

Force( 2*BC_Node(2)-1 : 2*BC_Node(2) , : ) = [];
Force( 2*BC_Node(1)-1 : 2*BC_Node(1) , : ) = [];
%% Displacement Calculation
Q = inv(K)*Force;
Q = [Q ; zeros(4,1)];
%% Stress Calculation
stress = zeros(length(Beam),1);
for n = 1 : length(Beam)
    stress(n) = E/L(n) * ...
        [-cos(angle_rad(n)) , -sin(angle_rad(n)) , cos(angle_rad(n)) , sin(angle_rad(n)) ] * ...
        [ Q(2*Element(n,1)-1) ; Q(2*Element(n,1)) ; Q(2*Element(n,2)-1) ; Q(2*Element(n,2)) ];
end

%% Reaction Calculation
Ks( 1 : BC_Node(1)*2-2 , : ) = [];
R = Ks * Q;
end