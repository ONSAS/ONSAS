%%%%%% Dinámica No-Lineal: Péndulo con Barra elástica (Green)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%%% Parametros Estructura
L = 2; %m
A = .7854/100^2; %m2 (fi 1cm)
E = 210000e6 %Pa (Acero)
m = 214; %kg
c = 0; %kg/s (amortiguamiento)
g = 9.81; %m/s2

%% Defino Vector de fuerzas Internas: fint(u)
%%    u = [u1 , u2]^T

Fint = @(u) E*A*L*(u(1)^2-2*L*u(2)+u(2)^2)/2/L^4*[u(1);-(L-u(2))];

%% Defino Matriz Tangente: Kt(u) = d Fint / du

KT = @(u) E*A/2/L^3 * ([u(1),-(L-u(2))]'*[2*u(1), -2*L+2*u(2)] + (u(1)^2-2*L*u(2)+u(2)^2)*eye(2,2));

%% Defino Vector de fuerzas Externas: gravedad

ft = @(t) [0;-m*g]; %N

%% Defino Matriz de Masa Concentrada

M = [m 0 ; 0 m];

%% Defino Matriz de Amortiguamiento

C = [c 0 ; 0 c];

%% Defino Condiciones Iniciales y Tiempo Final

t0 = 0;
tf = 20 %sec
%u0 = [0;-m*g/(E*A/L)];
%v0 = [2;.010];
u0 = [L;L];
v0 = [0;0];
ac0 = M\(ft(t0)-C*v0-Fint(u0)); % de ec de movimiento Mu.. + Fint(u) = ft


% Inicializacion HHT

dt = .20;
alfa = 0.175; delta = 1/2+alfa; beta = (1+alfa)^2/4;  %delta = 1/2; beta = 1/4; alfa = 0;

u(:,1) = u0; % u(0)
v(:,1) = v0; %v(0)
a(:,1) = ac0; %a(0)

% Comienza Marcha en el Tiempo usando HHT
t(1) = t0;

k=1;

epsg(k) = (u(1,k)^2-2*L*u(2,k)+u(2,k)^2)/2/L^2;
ftol = 1e-6;
Maxiter = 15;

while t<tf %% Paso temporal de HHT
 j = 1;
 ferr = +inf;
 ut = u(:,k); vt = v(:,k); at = a(:,k);
 aktdt = a(:,k);
 
 while and(j<Maxiter,ferr>ftol) %% Iteración tipo N-R para Equilibrio
    uktdt = ut + dt*vt + dt^2*((1/2-beta)*at+beta*aktdt);
    vktdt = vt + dt*((1-delta)*at+delta*aktdt);
    Keff = (1-alfa)*beta*dt^2*KT(uktdt) + M;
    feff = -( M*aktdt + (1-alfa)*Fint(uktdt) + alfa*Fint(ut) - (1-alfa)*ft(t(k)+dt) - alfa*ft(t(k)) );
    da = Keff\feff;
    aktdt = aktdt + da;
    ferr = norm(feff)/norm(ft(t(k)+dt));
    j = j+1;
  end
  j
  u(:,k+1) = ut + dt*vt + dt^2*((1/2-beta)*at+beta*aktdt);
  a(:,k+1) = aktdt;
  v(:,k+1) = vt + dt*((1-delta)*at+delta*aktdt);
  epsg(k+1) = (u(1,k+1)^2-2*L*u(2,k+1)+u(2,k+1)^2)/2/L^2;
  t(k+1)=t(k)+dt;
  
  k=k+1;
end


subplot(1,2,1)
plot(u(1,:),u(2,:),'-r',0,2,'*k')
axis([-1.1*L 1.1*L -0.1*L 2.1*L],'square')
xlabel('x_1')
ylabel('x_2')
title('Solución Con HHT')

subplot(1,2,2)
plot(t(2:end),E*A*epsg(2:end))
xlabel('t (sec)')
ylabel('Normal (N)')

%print -dpng pendulo.png
