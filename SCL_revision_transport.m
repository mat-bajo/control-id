clear
%close

% --------------------------------------------------------------------------
% ODE-Transport system
% --------------------------------------------------------------------------

% Parameters
%ODE part
Aode = [-1 2; -1 0];
Bu = [0;1];
Cy = [0 1];
nx = size(Aode,1);
sysol = ss(Aode,Bu,Cy,0); % Open-loop ODE part
% PDE part
Bode = [-2; 1]; % Retarded delay part
Eode = 1/10*[-0.2; 0.1]; % Neutral delay part
Code = [0 1];

T = 10; % Simulation time
h = 1; % Delay

% Eigenvalues approximation (Pade)
syms x; symH = pade(exp(-h*x),x,0,'Order',[9 10]);  
[symNum,symDen] = numden(symH);
TFnum = sym2poly(symNum);
TFden = sym2poly(symDen);
Hrde = tf(TFnum,TFden);
Hnde = tf([TFnum 0],TFden);
sysrde = ss(Hrde); sysnde = ss(Hnde);
Arde = sysrde.A; Brde = sysrde.B; Crde = sysrde.C; 
coef = sysnde.B(1)/Brde(1); Cnde = coef*sysnde.C; Dnde = coef*sysnde.D;
Apade = [Aode-Eode*Dnde*Code Bode*Crde-Eode*Cnde; Brde*Code Arde];
Bpade = [Bu; zeros(10,1)];
Cpade = [Cy zeros(1,10)];
[Vpade,Jpade] = jordan(Apade);
syspade = ss(Jpade,Vpade\Bpade,Cpade*Vpade,0); % both eigenvalues and eigenvectors

% Control design
n0 = 1; n = 3;
eigctrb = -1; % Eigenvalues of the controller part
eigobsv = -1; % Eigenvalues of the observer part
sysc = SCL_cont(syspade,n0,n,eigctrb,eigobsv); % Synthesis
syscl = lft(sysol,sysc); % Closed-loop ODE part
syscl.B = zeros(nx+n); syscl.B(1:nx,1:nx) = Bode*Code; % Retarded delay part
syscl.C = zeros(nx+n); syscl.C(1:nx,1:nx) = Eode*Code; % Neutral delay part

% --------------------------------------------------------------------------
% Simulation
% --------------------------------------------------------------------------

xinit = [0; 1];
outol = ddensd(@(t,y,ydel,ypdel) ddetransport(t,y,ydel,ypdel,Aode,Bode*Code,Eode*Code),@(t,y) dely(t,y,h), @(t,y) delyp(t,y,h),@(t) ddeinittransport(t,xinit),[0 T]);
outcl0 = ddensd(@(t,y,ydel,ypdel) ddetransport(t,y,ydel,ypdel,syscl.A(1:nx+n0,1:nx+n0),syscl.B(1:nx+n0,1:nx+n0),syscl.C(1:nx+n0,1:nx+n0)),@(t,y) dely(t,y,h), @(t,y) delyp(t,y,h),@(t) ddeinittransport(t,[xinit;zeros(n0,1)]),[0 T]); % ones(n0,1)
outcl = ddensd(@(t,y,ydel,ypdel) ddetransport(t,y,ydel,ypdel,syscl.A,syscl.B,syscl.C),@(t,y) dely(t,y,h), @(t,y) delyp(t,y,h),@(t) ddeinittransport(t,[xinit;zeros(n,1)]),[0 T]); % 1/n*ones(n,1)
figure(1)
plot(outol.x,sqrt(outol.y(1,:).^2+outol.y(2,:).^2),'Linewidth',4); hold on
plot(outcl0.x,sqrt(outcl0.y(1,:).^2+outcl0.y(2,:).^2),'-o','Linewidth',2,'MarkerIndices',1:5:length(outcl0.y(1,:))); hold on
plot(outcl.x,sqrt(outcl.y(1,:).^2+outcl.y(2,:).^2),'-x','Linewidth',2,'MarkerIndices',1:15:length(outcl.y(1,:))); hold on
set(gca,'Fontsize',18)
xlabel('Time $t$','Interpreter','Latex'); ylabel('Norm of the state','Interpreter','Latex')
grid on; box on;
legend('Open-loop','Closed-loop $n=1$','Closed-loop $n=3$','Interpreter','Latex')
ylim([0 3])

% --------------------------------------------------------------------------
% Validation of theorem
% --------------------------------------------------------------------------
K0 = sysc.B(1:n0,:); G0 = sysc.C(:,1:n0); P0 = eye(2*n0);
bk = zeros(n0); cg = zeros(n0);
for k = n+1:10
    eigenvalue = Jpade(k,k);
    eigenvector = Vpade(k,1:2)';%inv(eigenvalue*eye(nx)-Aode)*(Bode-eigenvalue*Eode);
    bk = bk + K0'*Bu'*eigenvector*eigenvector'*Bu*K0;
    cg = cg + eigenvector'*Cy'*[-G0;G0]'*P0*[-G0;G0]* Cy*eigenvector;
end
lambda = abs(real(Jpade(n,n)));
rho = 16*max(eig(bk))*max(eig(cg))/((lambda-1)*lambda);
if rho < 1
    disp('Closed-loop system is guaranteed to be exponentially stable')
end


function dydt = ddetransport(t,y,ydel,ypdel,A,B,E)
% Differential equations function 
    dydt = A*y+B*ydel-E*ypdel;
end
%-------------------------------------------
function dy = dely(t,y,h) % delay for y
    dy = t-h;
end
%-------------------------------------------
function dyp = delyp(t,y,h) % delay for y'
    dyp = t-h;
end
%-------------------------------------------
function s = ddeinittransport(t,xinit)
    s = xinit;%*(t==0);
end
%-------------------------------------------


