% ------------------------------------
% example simple plane chain
% ------------------------------------

clear all, close all

dirOnsas = [ pwd '/../..' ] ; % set ONSAS.m directory
addpath( dirOnsas ); % add ONSAS directory to path

problemName       = 'chain' ;

% -- scalar params -----
Es  = 210e9 ; % Pa
nu  = 0     ; %
b   = 0.05  ; % m - width of square section
l   = 2     ; % m
h   = 1     ;
rho = 8050  ; % kg/m^3

nodalDispDamping = 0 ;

% ----- geometry -----------------
Nelem = 10 ;

coordsXs = (0:(Nelem))'*l/Nelem ;

Nodes = [ coordsXs zeros(size(coordsXs)) h*abs( coordsXs-l/2 )-l/2 ] ;

auxconec = [ (ones(Nelem,1)*[ 1 2 0 1 0]) (1:(Nelem))' (2:(Nelem+1))' ] ;

Conec = cell(2+Nelem,1) ;

Conec{1, 1} = [ 0 1 0 0 1                     1       ] ; % fixed node
Conec{2, 1} = [ 0 1 0 0 1                     Nelem+1 ] ; % fixed node
for i=1:Nelem
  Conec{2+i, 1} =  auxconec(i,:) ;
end

%  --- MELCS params ---
materialsParams = { [ rho 2 Es nu ] } ; %truss
elementsParams  = { 1 ; [ 2 0] } ;

% loadsParams
booleanSelfWeightZ = 1 ;

crossSecsParams = { [2 b b ] } ;

springsParams = {[  inf  0  inf  0  inf 0 ] } ;
% -------------------

timeIncr      =  0.01              ;
finalTime     =  1.5                 ;
nLoadSteps    = finalTime/timeIncr ;
stopTolDeltau = 0                  ; 
stopTolForces = 1e-8               ;
stopTolIts    = 30                 ;

alphaHHT = -0.05 ;

numericalMethodParams = [ 4 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts alphaHHT ] ;

% ------------------------------------

controlDofs = [ Nelem/2 5 1 ] ;
plotParamsVector = [ 3 150 ];
printFlag = 0 ;
storeBoolean = 1;

ONSAS

disp('verificacion: ')
sumaFuerzaPesoNum = -sum(BCsData.constantFext)
sumaFuerzaPesoAna = sqrt( (l/2)^2 + h^2 )*2 * b*b * rho * 9.81
