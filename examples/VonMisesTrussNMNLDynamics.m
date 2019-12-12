%% Example VonMisesNMNLDynamics
% An example for nonlinear dynamic analysis of the VonMises Truss.
% Parameters as in example 4.3.1 from Bazzano and Perez Zerpa book: https://www.colibri.udelar.edu.uy/jspui/bitstream/20.500.12008/22106/1/Bazzano_P%C3%A9rezZerpa_Introducci%C3%B3n_al_An%C3%A1lisis_No_Lineal_de_Estructuras_2017.pdf
%%

% uncomment to delete variables and close windows
clear all, close all

%% General data
dirOnsas = [ pwd '/..' ] ;
problemName = 'VonMisesTrussNMNLDynamics' ;

%% Structural properties

Es  = 200e9 ;
nu  = 0.3   ;
rho = 7850  ;
hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 Es nu rho ] ;

a = .0032 ;
b = .0254 ;
A = a*b ;  %
secGeomProps = [ A 0 0 0 ] ;

Lx = .187 ; 
Lz = .084 ;
Lc = .240 ;

Ic = b*a^3 / 12 ;
kc = 3*Es*Ic/(Lc^3) ;

Nodes = [ 0  0 0  ; ...
          Lx 0 Lz ] ;

Conec = [ 1 2 0 0 1 1 1 ] ;

% in global system of coordinates
nodalSprings = [ 1  kc   0  inf  0  inf 0 ; ...
                 2  inf  0  inf  0    0 0 ] ;

%~ ceda    =  0.02 ;

%~ Nodes
%~ Conec
%~ nodalSprings
%~ stop

%~ %Parameters
%~ targetLoadFactr=1.5;

%~ nodalDamping = 0 ;

% method
timeIncr   =  0.00001    ;
finalTime  =  0.5             ;
nLoadSteps = finalTime/timeIncr ;
DeltaNW    =  0.5               ;
AlphaNW    =  0.25              ;

% tolerances
stopTolDeltau = 1e-8           ; 
stopTolForces = 1e-8           ;
stopTolIts    = 100            ;
% ------------------------------------

p0       = 0.1 ;
loadFactorsFunc = @(t) p0 * t ; 

nodalVariableLoads   = [ 2  0  0  0  0  -1  0 ];

controlDofInfo = [ 2 5 -1 ] ;
% ------------------------------

% ------------------------------
% analysis parameters
dynamicAnalysisBoolean   = 1 ; 

% initial conditions
%~ u0   = 0.2;
%~ udot0= 0;

%~ nonHomogeneousInitialCondU0 = [ 2 1 u0 ] ;
%~ nonHomogeneousInitialCondUdot0 =[ 2 1 udot0 ] ;

numericalMethodParams = [ 4 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts DeltaNW AlphaNW] ;

%~ plotParamsVector = [2 5 ];
plotParamsVector = [3 ]; sectPar = [ 12 a b ];
printflag = 0 ;
% ------------------------------------------------

acdir = pwd ;
cd(dirOnsas);
ONSAS
cd(acdir) ;
