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
coef = sysnde.B(1)/Brde(1); Cnde = 1/coef*sysnde.C; Dnde = 1/coef*sysnde.D;
Apade = [Aode-Eode*Dnde*Code Bode*Crde-Eode*Cnde; Brde*Code Arde];
Bpade = [Bu; zeros(10,1)];
Cpade = [Cy zeros(1,10)];
[Vpade,Jpade] = jordan(Apade);
syspade = ss(Jpade,Vpade\Bpade,Cpade*Vpade,0); % both eigenvalues and eigenvectors

% Control design
n0 = 1; n = 4;
eigctrb = -1; % Eigenvalues of the controller part
eigobsv = -1; % Eigenvalues of the observer part
sysc = SCL_cont(syspade,n0,n,eigctrb,eigobsv); % Synthesis
syscl = lft(sysol,sysc); % Closed-loop ODE part
syscl.B = zeros(nx+n,1); syscl.B(1:nx,1) = Bode; % Retarded delay part
sysclE = zeros(nx+n,1); sysclE(1:nx,1) = Eode; % Neutral delay part
syscl.C = zeros(1,nx+n); syscl.C(1,1:nx) = Code; % Output part
Acl0pade = [syscl.A(1:nx+n0,1:nx+n0)-sysclE(1:nx+n0,1)*Dnde*syscl.C(1,1:nx+n0) syscl.B(1:nx+n0,1)*Crde-sysclE(1:nx+n0,1)*Cnde; Brde*syscl.C(1,1:nx+n0) Arde];
Aclpade = [syscl.A-sysclE*Dnde*syscl.C syscl.B*Crde-sysclE*Cnde; Brde*syscl.C Arde];

% --------------------------------------------------------------------------
% Simulation
% --------------------------------------------------------------------------

xinit = [0; 1];
outol = ddensd(@(t,y,ydel,ypdel) ddetransport(t,y,ydel,ypdel,Aode,Bode*Code,Eode*Code),@(t,y) dely(t,y,h), @(t,y) delyp(t,y,h),@(t) ddeinittransport(t,xinit),[0 T]);
outcl0 = ddensd(@(t,y,ydel,ypdel) ddetransport(t,y,ydel,ypdel,syscl.A(1:nx+n0,1:nx+n0),syscl.B(1:nx+n0,1)*syscl.C(1,1:nx+n0),sysclE(1:nx+n0,1)*syscl.C(1,1:nx+n0)),@(t,y) dely(t,y,h), @(t,y) delyp(t,y,h),@(t) ddeinittransport(t,[xinit;zeros(n0,1)]),[0 T]); % ones(n0,1)
outcl = ddensd(@(t,y,ydel,ypdel) ddetransport(t,y,ydel,ypdel,syscl.A,syscl.B*syscl.C,sysclE*syscl.C),@(t,y) dely(t,y,h), @(t,y) delyp(t,y,h),@(t) ddeinittransport(t,[xinit;zeros(n,1)]),[0 T]); % 1/n*ones(n,1)
figure(1)
plot(outol.x,sqrt(real(outol.y(1,:).^2+outol.y(2,:).^2)),'Linewidth',4); hold on
plot(outcl0.x,sqrt(real(outcl0.y(1,:).^2+outcl0.y(2,:).^2)),'-o','Linewidth',2,'MarkerIndices',1:5:length(outcl0.y(1,:))); hold on
plot(outcl.x,sqrt(real(outcl.y(1,:).^2+outcl.y(2,:).^2)),'-x','Linewidth',2,'MarkerIndices',1:15:length(outcl.y(1,:))); hold on
set(gca,'Fontsize',18)
xlabel('Time $t$','Interpreter','Latex'); ylabel('Norm of the state','Interpreter','Latex')
grid on; box on;
legend('Open-loop','Closed-loop $n=1$','Closed-loop $n=4$','Interpreter','Latex')
ylim([0 3])

% --------------------------------------------------------------------------
% Visualisation des spectres
% --------------------------------------------------------------------------

figure(2)
scatter(real(eig(syspade.A)), imag(eig(syspade.A)),100,'filled'); hold on
scatter(real(eig(Acl0pade)), imag(eig(Acl0pade)),50,'o','Linewidth',2); hold on
scatter(real(eig(Aclpade)), imag(eig(Aclpade)),50,'x','Linewidth',2); 
set(gca,'Fontsize',18)
xlabel('Real part','Interpreter','Latex'); ylabel('Imag part','Interpreter','Latex')
grid on; box on;
legend('Open-loop','Closed-loop $n=1$','Closed-loop $n=4$','Interpreter','Latex')
xlim([-2.5 1.5]); ylim([-12 12])

% --------------------------------------------------------------------------
% Visualisation des spectres with QPmR
% --------------------------------------------------------------------------

% Characteristic equation 
box = [-4 1 -75 75];
%syms s; polyol = simplify(det(s*(eye(2)+Eode*Code*exp(-h*s))-Aode-Bode*Code*exp(-h*s))); % open loop
Pol = [1 1 2;1/100 -97/100 -3]; D = [0;1]; Rol = QPmR(box,Pol,D); % open loop
%syms s; polycl0 = simplify(det(s*(eye(2+n0)+sysclE(1:nx+n0,1)*syscl.C(1,1:nx+n0)*exp(-h*s))-syscl.A(1:nx+n0,1:nx+n0)-syscl.B(1:nx+n0,1)*syscl.C(1,1:nx+n0)*exp(-h*s))); % closed loop n0
Pcl0 = [1 7492675639051575/2251799813685248 5537578738330717038242197086341/633825300114114700748351602688 5745103538990914570830178361477/633825300114114700748351602688; 
    1/100 -213183706102102729/225179981368524800 -1183904899166108119/225179981368524800 -15722627476098981/2251799813685248 ];
Rcl0= QPmR(box,Pcl0,D); % closed loop n0
%syms s; polycl = simplify(det(s*(eye(2+n)+sysclE*syscl.C*exp(-h*s))-syscl.A-syscl.B*syscl.C*exp(-h*s))); % closed loop n0
p15 = 1486658069589317470032358459338142224730610758620508028524780370157948108800000;
p10 = 54177705732422385738152543826991788747779031631584040878140401354622768463186125 + 3553046236796577379018784033427509078649334547040311935028428800i;
p12 = 58749537359073947131522669199908511913919984131999751078902846375555569706598400 + 873610775149894237707403305227332422465227190745769000828928000i;
p13 = 27842018694557730242470545613640556144989284141222124889155338917033814261760000 - 549061934503927367620943229095615484183245778965091232147046400i;
p14 = 8433293101151868996761424156074727374217012674309497702806918855883600848486400 - 353850408744467084181045482390806632676194303809295734538240000i;
p21 = - (53835157952928246664536014604011116121305990168717744070224067682565646756872192 + 1137592847567959178504417016097915348623772829847164795662368768i);
p16 = 180925139433306555349329664076074856020734351040063381311652475012364265062400;
p11 = (90741930142049463633552401745221621409032183185575004268418810324139411245598925 + 1386100066879368122629596523304136836310564323208564972296601600i);
p22 = -(31926499725622548096619019057334885337687433372585366288256406093758005814231040 - 1268671637404060614374931181643222029383669515064210395970404352i);
p23 = -(9560467921089717157601463834048554652707045655872891363262865763865597789077504 - 341282781224538468821215140452034345180837959943458907526004736i);
p24 = - (1749656870421702671492272555132796414862992444890854047055103007844653938507776 + 3538504087444670841810454823908066326761943038092957345382400i) ;
p25 = - 162440055948747249542019486201171936653013556433057033400171621810537498673152;
p26 = 1809251394333065553493296640760748560207343510400633813116524750123642650624;
p20 = - (31369281555460640415609667247573212650992362457608745779002230835629806898380800 + 5329569355194866068528176050141263617974001820560467902542643200i);
Pcl = [p16 p15 p14 p13 p12 p11 p10; p26 p25 p24 p23 p22 p21 p20]/180925139433306555349329664076074856020734351040063381311652475012364265062400;
Rcl = QPmR(box,Pcl,D);

figure(3)
scatter(real(Rol),imag(Rol),100,'filled'); hold on;
scatter(real(Rcl0),imag(Rcl0),50,'o','Linewidth',2);  hold on;
scatter(real(Rcl),imag(Rcl),50,'x','Linewidth',2); 
set(gca,'Fontsize',18)
xlabel('Real part','Interpreter','Latex'); ylabel('Imag part','Interpreter','Latex')
grid on; box on;
legend('Open-loop','Closed-loop $n=1$','Closed-loop $n=4$','Interpreter','Latex')

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


