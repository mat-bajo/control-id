clear
close

% --------------------------------------------------------------------------
% ODE-Diffusion system
% --------------------------------------------------------------------------

% Parameters
% ODE part
Aode = [0 1; -4 -4];
Bu = [0; 1];
Cy = [0 1];
nx = size(Aode,1);
sysol = ss(Aode,Bu,Cy,0); % Open-loop ODE part
% PDE part
Bode = [0; -1];
Code = [0 1];
nu = 1/9;
lambda = 1.5;

th = linspace(0,1,25); dth = th(2)-th(1);
t = linspace(0,10,50);
tho = 24;

% Eigenvalues approximation (Pade)
syms x; symH = pade(-x/sqrt(nu)/sinh(x/sqrt(nu)),x,0,'Order',[20 20]);  % Pade at order 20 for the PDE part
%cosh(x/2/sqrt(nu))*
[symNum,symDen] = numden(symH);
TFnum = sym2poly(symNum);
TFnum(TFnum==0) = [];
TFden = sym2poly(symDen);
TFden(TFden==0) = []; % (pi/3)^2 puis (2*pi/3)^2
nmax = 10; % 20/2 to take into account the square root of lambda
TFnum2 = zeros(1,nmax+1);
TFden2 = zeros(1,nmax+1);
for ni = 0:nmax
    for k = 0:ni
        TFnum2(end-k) = TFnum2(end-k) + TFnum(end-ni)*nchoosek(ni,k)*(-lambda)^(ni-k);
        TFden2(end-k) = TFden2(end-k) + TFden(end-ni)*nchoosek(ni,k)*(-lambda)^(ni-k);
    end
end
Hpde = tf(TFnum2,TFden2);
syspde = ss(Hpde); Apde = syspde.A; Bpde = syspde.B; Cpde = syspde.C; Dpde = syspde.D;
Apade = [Aode+Bode*Dpde*Code Bode*Cpde; Bpde*Code Apde];
Bpade = [Bu; zeros(nmax,1)];
Cpade = [Cy zeros(1,nmax)];
[Vpade,Jpade] = jordan(Apade);
syspade = ss(Jpade,Vpade\Bpade,Cpade*Vpade,0); % both eigenvalues and eigenvectors

% Control design
n0 = 2; n = 4;
eigctrb = [-1 -1.5];%[-1+i -1-i]; % Eigenvalues of the controller part
eigobsv = eigctrb; % Eigenvalues of the observer part
sysc = SCL_cont(syspade,n0,n,eigctrb,eigobsv); % Synthesis
syscl = lft(sysol,sysc); % Closed-loop ODE part
syscl.B = zeros(nx+n,1); syscl.B(1:nx,1) = Bode; % Closed-loop PDE output
syscl.C = zeros(1,nx+n); syscl.C(1,1:nx) = Code;  % Closed-loop PDE intput
Aclpade = [syscl.A+syscl.B*Dpde*syscl.C syscl.B*Cpde; Bpde*syscl.C Apde];

% --------------------------------------------------------------------------
% Simulation
% --------------------------------------------------------------------------

% Open-loop
xinit = [0;1];
xode(:,1) = xinit;
sysode = ss(Aode,Bode,Code,0); 
pinit = [-sysode.C*xinit sysode.C*xinit];
for ind = 1:length(t)-1
    ti = [t(ind) (t(ind)+t(ind+1))/2 t(ind+1)];
    solpde = pdepe(0,@(x,t,u,dudx) pde(x,t,u,dudx,nu,lambda),@(x) init(x,pinit),@(xl,ul,xr,ur,t) ode(xl,ul,xr,ur,t,sysode.C*xode(:,ind)),th,ti);
    pinit = polyfit(th,solpde(end,:),5);
    ypde = (solpde(end,tho+1)-solpde(end,tho-1))/(2*dth);
    solode = ode45(@(t,y) sysode.A*y+sysode.B*ypde, ti, xode(:,ind));
    xode(:,ind+1) = solode.y(:,end);
end
figure(1)
plot(t,sqrt(xode(1,:).^2+xode(2,:).^2),'Linewidth',4); hold on
clear xode

% Closed-loop n0
xinit = [0;1;zeros(n0,1)];
xode(:,1) = xinit;
sysn0 = ss(syscl.A(1:nx+n0,1:nx+n0),syscl.B(1:nx+n0,1),syscl.C(1,1:nx+n0),0); 
pinit = [-sysn0.C*xinit sysn0.C*xinit];
for ind = 1:length(t)-1
    ti = [t(ind) (t(ind)+t(ind+1))/2 t(ind+1)];
    solpde = pdepe(0,@(x,t,u,dudx) pde(x,t,u,dudx,nu,lambda),@(x) init(x,pinit),@(xl,ul,xr,ur,t) ode(xl,ul,xr,ur,t,sysn0.C*xode(:,ind)),th,ti);
    pinit = polyfit(th,solpde(end,:),5);
    ypde = (solpde(end,tho+1)-solpde(end,tho-1))/(2*dth);
    solode = ode45(@(t,y) sysn0.A*y+sysn0.B*ypde, ti, xode(:,ind));
    xode(:,ind+1) = solode.y(:,end);
end

plot(t,sqrt(xode(1,:).^2+xode(2,:).^2),'-o','Linewidth',2);
clear xode

% Closed-loop n
xinit = [0;1;zeros(n,1)];
xode(:,1) = xinit;
sysn = ss(syscl.A,syscl.B,syscl.C,0); 
pinit = [-sysn.C*xinit sysn.C*xinit];
for ind = 1:length(t)-1
    ti = [t(ind) (t(ind)+t(ind+1))/2 t(ind+1)];
    solpde = pdepe(0,@(x,t,u,dudx) pde(x,t,u,dudx,nu,lambda),@(x) init(x,pinit),@(xl,ul,xr,ur,t) ode(xl,ul,xr,ur,t,sysn.C*xode(:,ind)),th,ti);
    pinit = polyfit(th,solpde(end,:),5);
    ypde = (solpde(end,tho+1)-solpde(end,tho-1))/(2*dth);
    solode = ode45(@(t,y) sysn.A*y+sysn.B*ypde, ti, xode(:,ind));
    xode(:,ind+1) = solode.y(:,end);
end
plot(t,sqrt(xode(1,:).^2+xode(2,:).^2),'-x','Linewidth',2); 

set(gca,'Fontsize',18)
xlabel('Time $t$','Interpreter','Latex'); ylabel('Norm of the state','Interpreter','Latex')
grid on; box on;
legend('Open-loop','Closed-loop $n=2$','Closed-loop $n=4$','Interpreter','Latex')
ylim([0 4])


function [c,f,s] = pde(x,t,u,dudx,nu,lambda)
c = 1;
f = nu*dudx;
s = lambda*u;
end
%----------------------------------------------
function u0 = init(x,p)
   u0 = polyval(p,x);
end
%----------------------------------------------
function [pl,ql,pr,qr] = ode(xl,ul,xr,ur,t,bc)
pl = ul-bc;
ql = 0;
pr = ur;
qr = 0; 
end