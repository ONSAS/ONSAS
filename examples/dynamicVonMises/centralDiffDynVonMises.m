
function [u, normalForce, t] = centralDiffDynVonMises(rho, Lx, L0, Lc, Ic, Ac, E, kc, m, c, g, tf, dt)

  mb = L0*Ac*rho ;
  Lz = sqrt( L0^2 - Lx^2 ); %m

%% Defino Vector de fuerzas Internas: fint(u)  - u = [u1 , u2]^T
Fint = @(u) E*Ac*L0*(u(1)^2+u(2)^2-2*Lx*u(1)+2*Lz*u(2))/2/L0^4*[-Lx+u(1);Lz+u(2)]+[kc*u(1);0];

%% Defino Vector de fuerzas Externas: gravedad
ft = @(t) [0;-(m+mb)/2*g]; %N

%% Defino Matriz de Masa Concentrada
M = [mb 0 ; 0 (mb+m)/2] ;

%% Defino Matriz de Amortiguamiento
C = [c/10 0 ; 0 c];

%% Defino Condiciones Iniciales
t0 = 0;
u0 = [0;0];
v0 = [0;0];
ac0 = M\(ft(t0)-C*v0-Fint(u0)); % de ec de movimiento Mu.. + Fint(u) = ft

% Inicializacion Difrerencia Centrada

a0 = 1/dt^2; a1=1/2/dt; a2=2*a0; a3=1/a2;

nTimes = tf/dt

Meff = a0*M+a1*C;
M2 = a0*M-a1*C;

% Comienza Marcha en el Tiempo usando Diferencia Centrada

k = 2 ; % index at which equilibrium is done. in this case t=0

epsg = @(u) (u(1)^2+2*Lz*u(2)-2*Lx*u(1)+u(2)^2)/2/L0^2;

u      = zeros(2, nTimes+2) ;
acc    = zeros(2, nTimes+2) ;
vel    = zeros(2, nTimes+2) ;
normalForce = zeros(1, nTimes+2) ;
t      = -dt:dt:tf ;

u(:,k-1) = u0 - dt*v0 + a3*ac0; % solution u(-dt) at time -dt
u(:,k  ) = u0; % u(0) %

while k<=(nTimes+1)
    feff      = ft(t(k)) - Fint(u(:,k)) +a2*M*u(:,k) - M2*u(:,k-1) ;
    u(:,k+1) = Meff\feff ; % sets solution at t+dt
    acc(:,k) = a0*(u(:,k+1)-2*u(:,k)+u(:,k-1)); % computers acc and vel at t
    vel(:,k) = a1*(u(:,k+1)-u(:,k-1));
    normalForce(k+1) = E * Ac * epsg( u(:,k+1 ) ) ;
    k=k+1;
end
