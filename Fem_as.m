clc;clear all;
%% Loading Raw data from individual data files

E= load('modulus.txt');
A= load('test.txt');
alpha= load('alp.txt');
del= load('delta.txt');
Node= load('nod.txt');
elcon= load('connection.txt'); % elements connection
DBC = load('DispBC.txt'); % Displacement Boundary condition
FBC = load('ForBC.txt');  % Force Boundary Condition
numnode= size(Node,1); % number of nodes
numelm= size(elcon,1); % number of elements
Kg= zeros(2*numnode); % Global stiffness matrix
Ug= zeros(2*numnode,1); % Global Disp vector
R= zeros(2*numnode,1); % Reaction Force

for i=1:numelm
    n1= elcon(i,1);
    n2= elcon(i,2);
    x1= Node(n1,1);
    x2= Node(n2,1);
    y1= Node(n1,2);
    y2= Node(n2,2);
    rot= atan2(y2-y1,x2-x1);
    L= sqrt((x2-x1)^2+(y2-y1)^2);
    S= sin(rot);
    C= cos(rot);
    
    kl= (E(i)*A(i)/L)* [C^2 C*S -C^2 -C*S;
                        C*S S^2 -C*S -S^2;
                        -C^2 -C*S C^2 C*S;
                        -C*S -S^2 C*S S^2];
     
     k1= 2*n1-1; k2= 2*n1;  % defining DOF
     k3= 2*n2-1; k4= 2*n2;
     
     Kg(k1:k2,k1:k2)= Kg(k1:k2,k1:k2) +kl(1:2,1:2);
     Kg(k1:k2,k3:k4)= Kg(k1:k2,k3:k4) +kl(1:2,3:4);
     Kg(k3:k4,k1:k2)= Kg(k3:k4,k1:k2) +kl(3:4,1:2);
     Kg(k3:k4,k3:k4)= Kg(k3:k4,k3:k4) +kl(3:4,3:4);    
end
  Kg  % Global Stiffness Matrix

c1= A(1)*E(1)*alpha(1)*del*[-1;1]; %Thermal force vector for ele1
c2= A(2)*E(2)*alpha(2)*del*[-1;1]; %Thermal force vector for ele2
c3= A(3)*E(3)*alpha(3)*del*[-1;1]; %Thermal force vector for ele3
fT= [0+c1(1);0+c2(1);0+c3(2);0+c1(2);5+c3(1);0+c2(2)] %Global force vector

%% Applying elimination approach

f= [fT(3);fT(5)];  
k= [Kg(3,3) Kg(3,5);Kg(5,3) Kg(5,5)];
u= k\f

Ug = [0;0;u(1);0;u(2);0] % Global Displacement vector
Fg = Kg*Ug;
R= Fg-fT %reaction forces

u1= [Ug(1);Ug(2);Ug(3);Ug(4)];%displacement of ele1 
u2= [Ug(1);Ug(2);Ug(5);Ug(6)]; %displacement of ele2
u3= [Ug(3);Ug(4);Ug(5);Ug(6)];  %displacement of ele3

%% stress in each element

for i=1
    n1= elcon(i,1);
    n2= elcon(i,2);
    x1= Node(n1,1);
    x2= Node(n2,1);
    y1= Node(n1,2);
    y2= Node(n2,2);
    rot= atan2(y2-y1,x2-x1);
    L= sqrt((x2-x1)^2+(y2-y1)^2);
    S= sin(rot);
    C= cos(rot);
    sig1= E(i)/L*[-C -S C S]*u1;
    m1= E(i)*alpha(i)*del;
    sigma1= sig1-m1;
end

for i=2
    n1= elcon(i,1);
    n3= elcon(i,2);
    x1= Node(n1,1);
    x3= Node(n3,1);
    y1= Node(n1,2);
    y3= Node(n3,2);
    rot= atan2(y3-y1,x3-x1);
    L= sqrt((x3-x1)^2+(y3-y1)^2);
    S= sin(rot);
    C= cos(rot);
    sig2 = E(i)/L*[-C -S C S]*u2;
    m2 = E(i)*alpha(i)*del;
    sigma2 = sig2-m2;
end

for i=3
    n3= elcon(i,1);
    n2= elcon(i,2);
    x3= Node(n3,1);
    x2= Node(n2,1);
    y3= Node(n3,2);
    y2= Node(n2,2);
    rot= atan2(y2-y3,x2-x3);
    L= sqrt((x2-x3)^2+(y2-y3)^2);
    S= sin(rot);
    C= cos(rot);
    sig3 = E(i)/L*[-C -S C S]*u3;
    m3 = E(i)*alpha(i)*del;
    sigma3= sig3-m3;
end
sigma1 % Display Stress in Element 1
sigma2 % Display Stress in Element 2
sigma3 % Display Stress in Element 3