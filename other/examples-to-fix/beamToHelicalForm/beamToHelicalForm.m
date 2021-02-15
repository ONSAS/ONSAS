% ------------------------------------------------------------------------------ 
% ------               ONSAS example file: beam bent to a helical form    ------
% ------------------------------------------------------------------------------

clear all, close all
dirOnsas = [ pwd '/..' ] ;
addpath( dirOnsas ) ;
problemName = 'beam2helicalform' ;

l = 10   ;    ty = .05 ;  tz = .05 ;  Nelem = 10     ;
%Ibrahimbegovic poroperties
EA = 10^4 ; nu = -0.5 ;  rho = 0 ; % E = G

% values consistent with given geometrical properties of the problem
A   =  1    ; I   =  1e-2 ; J   = I ; 
E = EA/A ;
G = E/(2*(1+nu));

%Verification OK

Nodes = [ (0:(Nelem))'*l/Nelem zeros(Nelem+1,2) ] ;

auxconec = [ (ones(Nelem,1)*[ 1 2 0 1 0]) (1:(Nelem))' (2:(Nelem+1))' ] ;

Conec = cell(2+Nelem,1) ;

Conec{1, 1} = [ 0 1 0 0 1                     1       ] ; % fixed node
Conec{2, 1} = [ 0 1 1 0 0                     Nelem+1 ] ; % loaded node
for i=1:Nelem
  Conec{2+i, 1} =  auxconec(i,:) ;
end

% ======================================================================
%% --- MELCS parameters ---

 
materialsParams = { [ rho 1 E nu] } ;

elementsParams  = { 1; 3} ;

loadsParams   = {[ 1 1   0 0 0 0 1 0 ]} ;

% --- cross section ---

crossSecsParams = {[ 1 A J I I ]} ;  % problem data given

springsParams    = {[ inf  inf  inf  inf  inf  inf ]} ;

% ======================================================================
%% Loads 
controlDofs     = [ Nelem+1  3  1 ]   ;
nLoadSteps      = 30                    ;
targetLoadFactrForce = 0.1             ;
targetLoadFactrMoment= 2.5*pi*10     ;
targetLoadFactr = targetLoadFactrForce ;  

global nLoadSteps targetLoadFactrForce targetLoadFactrMoment Nelem ;
userLoadsFilename = 'myBeam2HelicalLoadFunc';

%% Params
% ======================================================================
stopTolIts      = 50      ;
stopTolDeltau   = 1e-8       ;
stopTolForces   = 1e-8       ;
numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts targetLoadFactr nLoadSteps ] ; 
analyticFunc = @(w) 2*l/( l / ( E * I ))./w .*(sin((w * l / (2* E * I ))).^2) ;

%% Booleans
% ======================================================================
angExponUpdate      = 1    ;
storeBoolean        = 1     ;
plotParamsVector    = [ 3 ] ; 
printFlag           = 2     ;

ONSAS ;

%% Ploteo 
close all
%Plot visual params
lw = 3.5; ms = 5; plotfontsize = 22 ;

figure
loadFactors = 0:targetLoadFactrMoment/nLoadSteps:targetLoadFactrMoment;
plot( loadFactors, analyticFunc(loadFactors) ,'b-' , 'linewidth', lw,'markersize',ms )
hold on, grid on
 controlDisps(end)=[]
 plot( loadFactors, controlDisps ,'rx' , 'linewidth', lw,'markersize',ms )

legend ('AnalyticSol','NumericalSol')
laby = ylabel('Disp_z (m)');   labx = xlabel('Moment_y (N.m)') ;
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;


figure
plot( matUs(end-1,:) ,'k' , 'linewidth', lw,'markersize',ms )

