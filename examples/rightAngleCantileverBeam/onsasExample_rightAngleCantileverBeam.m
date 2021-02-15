% ------------------------------------------------------------------------------
% example Right-angle cantilever
%
% ------------------------------------------------------------------------------

clear all, close all

dirOnsas = [ pwd '/../../src' ] ; % set ONSAS.m directory
addpath( dirOnsas ); % add ONSAS directory to path

acdir    = pwd ;
problemName = 'rightAngleCantileverBeam' ;

% --- Scalar parameters --------------
E   =  1e6  ; nu  = -0.5  ; L   =  10   ; rho =  1    ;

% values consistent with given geometrical properties of the problem
A   =  1    ; I   =  1e-3 ; J   = I ; 
% ------------------------------------


% --- MELCS params ----

% materials
materialsParams = {[ rho 1 E nu ]} ;

% elements
elementsParams  = { 1; 3} ; %Nodo=1 Truss=1 Beam=3

% loads
loadsParams     = {[ 1 1   0 0 0  0 1 0 ]} ;
loadFactorsFunc = @(t) 50*t*(t<1) + (100-50*t)*(t>=1)*(t<2) + 0 ;

% crosssections
crossSecsParams = {[ 1 A J I I 20 10 10 ]} ;  % problem data given

% springs
springsParams    = {[ inf  inf  inf  inf  inf  inf ]} ;

% --------------------


% --- Nodes Matrix ---
nElemsPerBeam = 10 ;
auxCoords     = linspace( 0, L, nElemsPerBeam+1 )' ;

Nodes = [ zeros(nElemsPerBeam+1,1)  auxCoords                zeros(nElemsPerBeam+1,1) ; ...
          -auxCoords(2:end)         ones(nElemsPerBeam,1)*L  zeros(nElemsPerBeam  ,1) ] ;

% ----------------------
          
% --- Conec Cell -------

% add node elements
Conec = { [ 0 1 0 0 1	  1               ]   ; ... % fixed node
          [ 0 1 1 0 0   nElemsPerBeam+1 ] } ;     % loaded node

aux = (1:(2*nElemsPerBeam+1))' ;
% add frame elements
for i=1:(length(aux)-1)
  Conec{i+2,1} = [ 1 2 0 1 0  aux(i) aux(i+1) ] ;
end
% ----------------------


% --- Method and tolerances	----------
% tolerances
stopTolDeltau = 0     ;  stopTolForces = 1e-7 ;  stopTolIts = 30 ;

 timeIncr      =  0.25 ;  finalTime     = 20   ;    
% timeIncr      =  0.25 ;  finalTime     = 2   ;    

%~ alphaHHT = 0 ; % Newmark case
alphaHHT = -0.05 ; % HHT

numericalMethodParams = [ 4 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts alphaHHT ] ;
% ------------------------------

% --- Plot parameters ---
controlDofs      = [ nElemsPerBeam+1 3  1 ] ;
plotParamsVector = [ 3 ];

ONSAS ;

lw  = 2   ; ms  = 5.5 ;
lw2 = 3.2 ; ms2 = 23 ;
plotfontsize = 22 ;

figure
%~ plot(controlDisps,'b--o','linewidth',lw,'markersize',ms);
plot(timesVec, controlDisps,'b--o','linewidth',lw,'markersize',ms);
grid on
labx = xlabel('Time (s)');   laby = ylabel(sprintf('Displacement node: %2i dof %1i', controlDofs(1), controlDofs(2) ) ) ;
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;

figure
grid on
plot(timesVec, loadFactors,'r','linewidth',lw,'markersize',ms);
