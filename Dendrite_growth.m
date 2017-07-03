clear
clc
clf

setstep =250;

PI = 3.14159;
RR = 8.3145;

%*********************** Calculation Domain *****************************
nd = 750;  
dx = 30.0e-9;
dy = 30.0e-9; 
al = dx*(nd-1);  


%************************ material parameters  ***************************

Tm = 1728.0;                    
cndct = 84.01;   
speht = 5.42e+06;   
rlate  = 2.350e+09; 

%********************* input factors ************************************

Tini = 1337.8; % 0.9                   % modify here before running
% Tini = 1424.5; % 0.7
% Tini = 1511.2; % 0.5
% Tini = 1597.9; % 0.3

%astre = 0.05;                          % modify here before running
%astre = 0.015;
%astre = 0.005;

%********************** initial configuration ****************************

phi = zeros(nd,nd);
T = zeros(nd,nd);

delta = 3.0*dx; 

for i = 1:nd
for j = 1:nd
    if ((i-375)*(i-375) + (j-375)*(j-375) < 20.0) phi(i,j) = 0.9; end
end
end

for i = 1:nd
for j = 1: nd
T(i,j) = Tini + phi(i,j)*(Tm - Tini);
end
end

%*********************** independent varibales ***************************

gamma = 0.37;  % interfacial energy
ram = 0.1;    
bbb = 2.0*log((1.0 + (1.0 - 2.0*ram))/(1.0 - (1.0 - 2.0*ram)))/2.0;
j_fold = 4.0;  % modes of anisotropy
skine = 2.0;  % interface kinetics [m/Ks]
aaa = sqrt(3.0*delta*gamma/bbb);  % gradient coefficient
www = 6.0*gamma*bbb/delta;   % energy barrier
pmobi = bbb*Tm*skine/(3.0*delta*rlate); % mobility of phase field
th0 = 0.0;    % angle of prioritized growth direction
anois = 0.1;   % amplitude of noise

%*********************** time interval and steps **************************

dtp = dx*dx/(5.0*pmobi*aaa*aaa);   
dtt = dx*dx/(5.0*cndct/speht);
if(dtp>dtt) delt=dtt;
else delt = dtp; end

outstep = setstep;
instep = setstep;
  
%***********************  implementation euqations  ***********************

for m = 1: outstep
for n = 1: instep

for i = 2:nd-1
for j = 2:nd-1
    
if (i == nd-1)phi(i+1,j) = phi(i-1,j); end    %Neumann boundary conditions
if (i == 2) phi(i-1,j) = phi(i+1,j); end
if (j == nd-1)phi(i,j+1) = phi(i,j-1); end
if (j == 2) phi(i,j+1) = phi(i,j-1); end
if (i == nd-1 && j == nd-1) phi(i+1,j+1) = phi(i-1,j-1); end
if (i == 2  && j == 2) phi(i-1,j-1) = phi(i+1,j+1); end

if (i == nd-1)T(i+1,j) = T(i-1,j); end
if (i == 2) T(i-1,j) = T(i+1,j); end
if (j == nd-1)T(i,j+1) = T(i,j-1); end
if (j == 2) T(i,j+1) = T(i,j-1); end
if (i == nd-1 && j == nd-1) T(i+1,j+1) = T(i-1,j-1); end
if (i == 2  && j == 2) T(i-1,j-1) = T(i+1,j+1); end

%first order derivatives
dxphi = (phi(i+1,j) - phi(i-1,j))/2.0/dx;  
dyphi = (phi(i,j+1) - phi(i,j-1))/2.0/dy; 
%second order derivatives
dxxphi = (phi(i+1,j) + phi(i-1,j) - 2.0*phi(i,j))/dx/dx;  
dyyphi = (phi(i,j+1) + phi(i,j-1) - 2.0*phi(i,j))/dy/dy;
dxyphi = (phi(i+1,j+1) + phi(i-1,j-1) - phi(i-1,j+1) - phi(i+1,j-1))/4.0/dx/dy;

% angle of normal line for interface
th = atan(dyphi/(dxphi + 1.0e-20)); 
% gradient coefficient
ep = aaa*(1.0 + astre*cos(j_fold*(th-th0))); 
%  1st order derivative to angle
ep1p = -aaa*astre*j_fold*sin(j_fold*(th-th0));
% 2nd order derivative to angle
ep2p = -aaa*astre*j_fold*j_fold*cos(j_fold*(th-th0));

phikais = - ep*ep*(dxxphi + dyyphi) - ep*ep1p*((dyyphi - dxxphi)*sin(2.0*th) + 2.0*dxyphi*cos(2.0*th))+ 0.5*(ep1p*ep1p + ep*ep2p)*(2.0*dxyphi*sin(2.0*th) - (dxxphi + dyyphi) - (dyyphi - dxxphi)*cos(2.0*th));
dF = 15.0*rlate*(T(i,j) - Tm)*phi(i,j)*(1.0-phi(i,j))/2.0/www/Tm;
phikai = 4.0*www*phi(i,j)*(1.0 - phi(i,j))*(0.5 - phi(i,j) + dF + anois*(rand(1) - 0.5));
% kinectic equation
%phi(i,j) = phi(i,j) + phiddtt*delt;    %explicit method
phiddtt = -pmobi*(phikai + phikais); 
% addtion of thermal flunctuation
phi(i,j) = phi(i,j) + phiddtt*delt + anois*(rand(1) - 0.5)*phi(i,j)*(1.0 - phi(i,j)); 
% correction of boundary value
if (phi(i,j) >= 1.0) phi(i,j) = 1.0;end    
if (phi(i,j) <= 0.0) phi(i,j) = 0.0;end    

Tddtt = (cndct*((T(i+1,j) + T(i-1,j) - 2.0*T(i,j))/dx/dx + (T(i,j+1) + T(i,j-1) - 2.0*T(i,j))/dx/dy) + 30.0*phi(i,j)*phi(i,j)*(1.0 - phi(i,j))*(1.0 - phi(i,j))*rlate*phiddtt)/speht;

T(i,j) = T(i,j) + Tddtt*delt; %explicit method
end
end
mesh(phi);
view(2);
end
filename = sprintf('input_factors_%d.jpg',m);          %modify here before running
saveas(gcf, filename,'png');
end
