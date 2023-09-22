function sysc = control_id(sys,n0,n1,eigctrb,eigobsv)
% Input : sysol system open loop
% Output : sysc controller

A = sys.A;
B = sys.B;
C = sys.C;

%% Spliting
A0 = A(1:n0,1:n0);
A1 = A(n0+1:n1,n0+1:n1);
B0 = B(1:n0,1);
B1 = B(n0+1:n1,1);
C0 = C(1,1:n0);
C1 = C(1,n0+1:n1);

%% Controller
K0 = -acker(A0,B0,eigctrb);
G0 = -acker(A0,C0',eigobsv)';
F = [A0+G0*C0 G0*C1; zeros(n1-n0,n0) A1];
G = [-G0; zeros(n1-n0,1)];
H = [B0;B1];
K = [K0 zeros(1,n1-n0)];
sysc = ss(F+H*K,G,K,0);


end

