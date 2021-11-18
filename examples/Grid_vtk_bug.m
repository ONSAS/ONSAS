% Orthogonal grid

close all, clear all, close all
% add path
addpath( genpath( [ pwd '/../../../../repos/ONSAS.m/src'] ) );

% material scalar parameters
E = 30e6 ; % kN/m2
nu = 0.2 ;
% mesh
l = 5 ; %m
width = 1 ;
thk = 0.1 ;
%~ nrows = 20 ;
nrows = 4 ;
ncolumns = 3 ;
%~ ncolumns = 5 ;

ty = thk  ;
tz = 2*thk ; 



% MEBI parameters

% Materials
% ----------------------------------------------------------------------
materials.hyperElasModel  = '1DrotEngStrain' ;
materials.hyperElasParams = [ E nu ] ;
% Elements
% ----------------------------------------------------------------------
% Types
elements(1).elemType = 'node'  ;

% Columns
elements(2).elemType = 'frame' ;

% rows
elements(3).elemType = 'frame' ;

% Sections
% Columns
elements(2).elemTypeGeometry = [2 ty tz ] ;
elements(2).elemTypeParams   = 1          ;

% Rows
elements(3).elemTypeGeometry = [2 ty tz*2 ] ;
elements(3).elemTypeParams   = 1          ;


% Boundary conditions
% ----------------------------------------------------------------------
% Pinned support
boundaryConds(1).imposDispDofs = [ 1 3 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 ] ;
% Roller support
boundaryConds(2).imposDispDofs = [ 1 3 6 ] ;
boundaryConds(2).imposDispVals = [ 0 0 0 ] ;
% Load
P = -1 ;
imp = P/500 ;

boundaryConds(2).loadsCoordSys = 'global'        ;
boundaryConds(2).loadsBaseVals = [ 0 imp 0 0 P 0 ] ;

% Initial conditions
% ----------------------------------------------------------------------
initialConds                = struct() ;


% Mesh

% Nodes coords
% ----------------------------------------------------------------------
mesh.nodesCoords = [ ] ;


for i = 1:nrows

	zCoord = l * (i-1) / (nrows-1) ;
	%~ yCoord = sin( pi * zCoord / l ) ;
	yCoord = 0 ;
	
	for j = 1:ncolumns	
		xCoord = width * (j-1) / (ncolumns-1) ;
		vec = [ xCoord yCoord zCoord ] ;
		mesh.nodesCoords = [ mesh.nodesCoords ; vec ] ;
	end	
end

% Conec cell
% ----------------------------------------------------------------------
mesh.conecCell = {  } ;

% First BC
for i = 1:ncolumns
	%~ mesh.conecCell{ (1:ncolumns), 1 } = 																	[ zeros(ncolumns,1) ones(ncolumns,1) 	 ones(ncolumns,1) zeros(ncolumns,1) (1:ncolumns)' ] ;
	mesh.conecCell{ i, 1 } = [ 0 1 1 0 i ] ;
end

% Second BC
for i = 1:ncolumns
	mesh.conecCell{ ncolumns+i, 1 } = [ 0 1 2 0 ncolumns*(nrows-1)+i ] ;
end

% Row elements
for i=1:nrows
	for j=1:(ncolumns-1)

		mesh.conecCell{ ncolumns*2 + (ncolumns-1)*(i-1) + j, 1 } = [ 1 3 0 0 ncolumns*(i-1)+(j) ncolumns*(i-1)+(j+1) ] ;
	end
end

% Column elements
for j=1:ncolumns
	for i=1:(nrows-1)
	
		mesh.conecCell{ ncolumns*2 + nrows*(ncolumns-1) + (nrows-1)*(j-1)+i, 1 } = [ 1 2 0 0 ncolumns*(i-1)+(j) ncolumns*i+(j) ] ;
	end
end

% Analysis settings

% Parameters
% ----------------------------------------------------------------------
analysisSettings.methodName    = 'newtonRaphson' ;
boundaryConds(2).loadsTimeFact = @(t) 10*t     ;

analysisSettings.deltaT        =   .1  ;
analysisSettings.finalTime     =   0.5   ;
analysisSettings.incremArcLen  =   0.005   ;
analysisSettings.stopTolDeltau =   1e-6 ;
analysisSettings.stopTolForces =   1e-6 ;
analysisSettings.stopTolIts    =   10   ;
analysisSettings.posVariableLoadBC = 2 ;
analysisSettings.iniDeltaLamb = 0.1 ;


otherParams.problemName = 'GridVTK';
%~ otherParams.controlDofs = [ numElements+1  5 ] ;
otherParams.plotsFormat = 'vtk' ;

A = ty*tz ;
I = width*ty^3/12 ;

Pcrit = pi()^2*E*I/l^2 

[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;


centralNode = ceil((ncolumns/2+1))*(nrows/2+1) ;
Dof      		 = centralNode*6 - 3 	;
controlDisps =  matUs(Dof, :) 						;
loadFactors  =  loadFactorsMat(:, 2) 			;

lw = 2.0 ; ms = 11 ; plotfontsize = 22 ;
figure
grid on
plot( controlDisps, loadFactors, 'k-o' , 'linewidth', lw,'markersize',ms )
labx = xlabel('Displacement');   laby = ylabel('$\lambda$') ;
