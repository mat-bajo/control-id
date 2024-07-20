function sysc = SCL_cont(sys,n0,n,eigctrb,eigobsv)
% Input : 
% sys system open loop
% n0 number of unstable poles
% n size of the controller
% eigctrb n0 poles to be reached with the control gain
% eigobsv n0 poles to be reached with the observer gain
% Output : 
% sysc system controller

%% Spliting
% original system
A = sys.A;
B = sys.B;
C = sys.C;
% active part of size n0
A0 = A(1:n0,1:n0);
B0 = B(1:n0,1);
C0 = C(1,1:n0);
% passive part of size n-n0
A1 = A(n0+1:n,n0+1:n);
B1 = B(n0+1:n,1);
C1 = C(1,n0+1:n);

%% Controller
% Pole placement
K0 = -acker(A0,B0,eigctrb); % control gain
G0 = -acker(A0',C0',eigobsv)'; % observer gain
% State, input and output matrices
F = [A0+G0*C0 G0*C1; zeros(n-n0,n0) A1]; 
G = [-G0; zeros(n-n0,1)];
H = [B0;B1];
K = [K0 zeros(1,n-n0)];
sysc = ss(F+H*K,G,K,0);
end

