% ------------------------------------
% Test problem: continous beam with internal releases.
% ------------------------------------
inputONSASversion = '0.1.8' ;
problemName = 'Pendulum' ;

% Constitutive properties
E  = 210e9 ; nu = 0.3 ; rho= 4200;
hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 E nu rho] ;

% Geometrical properties
l = 300e-2; 
EA=10e9; A=EA/E; 
d=sqrt(4*A/pi);
Iy = pi*d^4/64; ; Iz =Iy ; It = 2*Iy ;
secGeomProps = [ A Iy Iz It ] ;

M=rho*A*l;
g=9.81;

%delta mass mat and nelem-------------------------
if exist('previouslyDefinedSelectedFileVar')~=0
  % ver como se hace voler atras
 
  Casoparticular=load ('../MIOPendulum/caso.txt');
  
  deltamassMat = Casoparticular(1);
  Nelem        = Casoparticular(2);

else
  deltamassMat=1e-3;
  Nelem =20;
end
%-----------------------------------
% Nodes and Conec 
Nodes       = [ (0:(Nelem))'*l/Nelem zeros(Nelem+1,2) ] ;
nnodes      = size(Nodes,1);
auxelemType = 2;
Conec       = [ (1:(Nelem))' (2:(Nelem+1))'  zeros(Nelem,2) (ones(Nelem,1)*[ 1 1 auxelemType]) ] ;
Nnodes=size(Conec,1)+1;

% Boundary conditions

nodalSprings = [   1  inf  0    inf   0    inf       0  ];    

% Nodal loads
ElemWheight               =  g * rho * A * l/Nelem ; 
nodalConstantLoads        = [ (1:Nnodes)'  zeros(Nnodes,4) -ElemWheight*ones(Nnodes,1) zeros(Nnodes,1) ] ;
nodalConstantLoads(1,6)   = nodalConstantLoads(1,6)*1/2;    %reduce whight at the end
nodalConstantLoads(end,6) = nodalConstantLoads(end,6)*1/2;
nodalConstantLoads(end,6) = nodalConstantLoads(end,6) ;

% Analysis parameters
nonLinearAnalysisBoolean = 1 ; 
dynamicAnalysisBoolean   = 1 ; 

% dynamic method params
timeIncr  =  0.05               ;
finalTime =  4             ;
nLoadSteps=finalTime/timeIncr  ;
DeltaNW   =  0.5               ;
AlphaNW   =  0.25              ;


% tolerances
stopTolDeltau = 1e-12           ; 
stopTolForces = 1e-12           ;
stopTolIts    = 1000            ; 

% initial conditions
u0   = 0
udot0= 0;
nonHomogeneousInitialCondU0 = [ 2 5 u0 ] ;
nonHomogeneousInitialCondUdot0 =[ 2 5 udot0 ] ;

if exist('previouslyDefinedSelectedFileVar')~=0
  plotParamsVector = [ 0 ] ;
 else
  plotParamsVector = [ 3 ] ;
end




printflag = 1 ;


controlDofInfo = [ nnodes 5 -1 ]     ;

numericalMethodParams = [ 4 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts DeltaNW AlphaNW] ;


