
% numerical example
% cantilever beam loaded with a vertical force at the free end

function [ Mn, Fn, matdes ] = softHinge1DOF_numericSol(l, A, E, Inertia, Mc, My, Mu, kh1, kh2, Ks)
        
freedofs = [2 4 6]; % u2 v2 theta2

% at the beginning..., there was no softening hinge
soft_hinge_boolean = false ;

% Gauss-Lobatto Quadrature with 3 integration points [a (a+b)/2 b]
npi = 3 ;
xpi = [0 l/2 l] ;
wpi = [1/3 4/3 1/3] * l * 0.5 ;

nu   = 0.3 ;
tol1 = 1e-8;
tol2 = 1e-8 ;
tolk = 15 ;

% initial values
dn   = [0 0 0 0 0 0]' ;
Fint = [0 0 0 0 0 0]' ;
tM   = 0 ;

kpn  = zeros(npi,1) ;
xin1 = zeros(npi,1) ;
xin2 = 0 ;

kpn1  = zeros(npi,1) ;
xin11 = zeros(npi,1) ;
xin21 = 0 ;

khat1 = zeros(npi,1) ;

M1 = zeros(npi,1) ;
Fn = zeros(npi,1) ;

alfan = 0 ;

xd = 0 ;
xdi = 1 ;

Final_force = 1000 ; % value of the final force
%Final_force = 4 ; % value of the final force

load_case = [0 0 0 1 0 0]' ; % load applied in vertical direction (Y)
load_factors = 0:Final_force ;

% --- element params ---
elemParams = [l A Inertia] ;

% --- elastoplastic params ---
elastoplasticParams = [E Mc My Mu kh1 kh2 Ks] ;

matdes = zeros (6, Final_force+1) ;

matdes(:,1) = dn ;

gxin = zeros(Final_force, 1) ;
gxin2 = zeros(Final_force, 1) ;
gkpn = zeros(Final_force, 1) ;

Mn = zeros(Final_force, 1) ;
TM = zeros(Final_force, 1) ;

Alf = zeros(Final_force, 1) ;

for ind = 2:length(load_factors)

    curr_load_factor = load_factors(ind) ;

    Fext = load_case * curr_load_factor ;

    dnk = matdes(:,ind-1) ;

    % iteration vars
    converged_boolean = false ;
    k = 0 ; % set iterations zero

    gxin(ind-1,1)   = xin1(1)   ;
    gxin2(ind-1,1)  = xin2      ;
    gkpn(ind-1,1)   = kpn(1)    ;
    Mn(ind-1,1)     = M1(1)     ;
    TM(ind-1,1)     = tM        ;
    Fn(ind-1,1)     = Fint(4)   ;
    Alf(ind-1,1)    = alfan     ;

    while converged_boolean == false && k < tolk && tM >= 0

        k = k + 1 ;

        [soft_hinge_boolean, Fint, M1, Kelement, kpn1, xin11, xin21, alfan1, xd, xdi, tM] = framePlastic(soft_hinge_boolean, dnk, kpn, xin1, xin2, alfan, xd, xdi, tM, elemParams, elastoplasticParams) ;

        residualForce = Fext - Fint ;

        Krelement = Kelement(freedofs,freedofs) ;
 
        residualForceRed = residualForce(freedofs) ;
      
        % system of equilibrium equations

        deltadred = Krelement\residualForceRed ;

        % /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\
        
        deltad = zeros(6,1) ;      
        deltad(freedofs) = deltadred ;

        dnk1 = dnk + deltad ;

        dnk     = dnk1   ;

        kpn     = kpn1   ;
        xin1    = xin11  ;

        xin2    = xin21  ;
        alfan   = alfan1 ;

        norm1 = norm(deltadred) ;
        norm2 = norm(residualForceRed) ;

        converged_boolean = norm1 < tol1 || norm2 < tol2 ;

    end

    matdes(:,ind) = dnk1 ;

end

# lw = 2.5 ; ms = 0.2 ; plotfontsize = 14 ;

# figure('Name','Cantilever Beam / Plasticity','NumberTitle','off');
# hold on, grid on
# plot(abs(matdes(6,1:length(load_factors)-1)), Mn,'b-x' , 'linewidth', lw, 'markersize', ms, "Color", "#EDB120") ;
# plot(abs(matdes(4,1:length(load_factors)-1)), Mn, 'k-o' , 'linewidth', lw, 'markersize', ms, "Color", "#0072BD") ;
# labx = xlabel('Generalized displacements in free node (m, rad)');   laby = ylabel('Moment in plastic hinge (KN.m)') ;
# legend('Degree of Freedom y','Degree of Freedom \theta','location','Southeast') ;
# set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
# set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
# title('Cantilever Beam / Plasticity') ;

# figure('Name','Cantilever Beam / Plasticity','NumberTitle','off');
# hold on, grid on
# plot(abs(matdes(6,1:length(load_factors)-1)), Fn,'b-x' , 'linewidth', lw, 'markersize', ms, "Color", "#EDB120") ;
# plot(abs(matdes(4,1:length(load_factors)-1)), Fn,'b-x' , 'linewidth', lw, 'markersize', ms, "Color", "#A2142F") ;
# labx = xlabel('Generalized displacements in free node (m, rad)');   laby = ylabel('Load Applied (KN)') ;
# legend('Degree of Freedom y','Degree of Freedom \theta','location','Southeast') ;
# set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
# set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
# title('Cantilever Beam / Plasticity') ;