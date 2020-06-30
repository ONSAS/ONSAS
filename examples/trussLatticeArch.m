%% Truss Lattice Arch example
% In this example a planar arch (in plane x-z) is considered, using 
% truss elements in a modular or lattice configuration.
% At this version the rotation of each section is not done TO DO.
%%

clear all, close all

%% General data
dirOnsas = [ pwd '/..' ] ;
problemName = 'trussLatticeArch' ;

%% Structural properties

Es = 210e9 ;
nu = 0 ;
hyperElasParams = cell(1,1) ;
hyperElasParams{1} = [1 Es nu] ;

% each row shows the properties of each section: A, Iy Iz and J
A = 0.02^2 ; 
secGeomProps = [ A 0 0 0  ] ;

numModules =  20 ;
lengthX    = 10 ;
angle = 10 * (pi/180) ;
radius = lengthX*.5 / sin(angle);

angles = linspace(-angle,+angle,numModules+1);
zs = radius * ( cos(angles) - cos(angle) ) ;
factors = ( cos(angles) - cos(angle) ) ;

lxmod = lengthX/numModules ;
height     = .5 ;
hsup = height*2/3;
hinf = height*1/3;
ladotriang = height*2/sqrt(3);

Nodes = [] ;


for i=0:numModules
  Nodes = [ Nodes; ...
            (radius+hsup)*sin(angles(i+1))               0  (radius+hsup)*cos(angles(i+1))   ; ...
            (radius-hinf)*sin(angles(i+1))  -ladotriang*.5  (radius-hinf)*cos(angles(i+1))   ; ...
            (radius-hinf)*sin(angles(i+1))  +ladotriang*.5  (radius-hinf)*cos(angles(i+1))   ] ;
end

numnodes = size(Nodes,1);

% in global system of coordinates
nodalSprings = [ 1  inf  0  inf  0  inf 0 ; ...
                 2  inf  0  inf  0  inf 0 ; ...
                 3  inf  0  inf  0  inf 0 ; ...
                 numnodes-2  inf  0  inf  0  inf 0 ; ...
                 numnodes-1  inf  0  inf  0  inf 0 ; ...
                 numnodes    inf  0  inf  0  inf 0 ];

Conec = [ 1 2 0 0 1 1 1 ; ...
          2 3 0 0 1 1 1 ; ...
          3 1 0 0 1 1 1 ] ;

for i=0:(numModules-1)
  Conec = [ Conec ; ...
            i*3+1  i*3+4 0 0 1 1 1 ; ...
            i*3+2  i*3+5 0 0 1 1 1 ; ...
            i*3+3  i*3+6 0 0 1 1 1 ; ...
            %
            i*3+4  i*3+5 0 0 1 1 1 ; ...
            i*3+5  i*3+6 0 0 1 1 1 ; ...
            i*3+6  i*3+4 0 0 1 1 1 ; ...
            %
            i*3+1  i*3+6 0 0 1 1 1 ; ...
            i*3+2  i*3+4 0 0 1 1 1 ; ...
            i*3+3  i*3+5 0 0 1 1 1 ; ...
            %
            i*3+1  i*3+5 0 0 1 1 1 ; ...
            i*3+2  i*3+6 0 0 1 1 1 ; ...
            i*3+3  i*3+4 0 0 1 1 1 ; ...
            %
            ] ;      
end

%% Loading parameters

nodalVariableLoads   = [ round(numnodes*.5)-1  0  0  0  0 -1  0 ];

%% Analysis parameters

% [ node nodaldof scalefactor(positive or negative) ]
controlDofInfo = [ round(numnodes*.5)-1 5 -1 ] ;

% analysis parameters
nonLinearAnalysisBoolean = 1 ;  dynamicAnalysisBoolean   = 0 ; 
LBAAnalyFlag             = 0 ; 

targetLoadFactr  = 1e8    ;

stopTolIts       = 30     ;
stopTolDeltau    = 1.0e-10 ;
stopTolForces    = 1.0e-6 ;
nLoadSteps       = 80    ;
incremArcLen     = .05    ;

numericalMethodParams = [ 2 stopTolDeltau stopTolForces stopTolIts ...
                            targetLoadFactr nLoadSteps incremArcLen ] ; 
%~ numericalMethodParams = [ 1 stopTolDeltau stopTolForces stopTolIts ...
                            %~ targetLoadFactr nLoadSteps ] ; 

stabilityAnalysisBoolean = 0 ;

%% Output parameters
printflag = 0 ;
%~ plotParamsVector = [ 2 4 ];
plotParamsVector = [ 3 30 ];
sectPar = [12 sqrt(A) sqrt(A) ] ;

reportBoolean = 1 ;

%% ONSAS execution
% move to onsas directory and ONSAS execution

acdir = pwd ;
cd(dirOnsas);
ONSAS
cd(acdir) ;

