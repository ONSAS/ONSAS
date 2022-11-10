% Cantilever with elastoplastic perfect constitutive model
% ----------------------------------------------------------------------
%
% Previous definitions

close all, clear all

onsasPath = '/../../src' ;
addpath( genpath( [ pwd onsasPath ] ) ) ; % add ONSAS directory to path

% MEBI parameters: Material-Element-BoundaryConditions-InitialConditions
% ----------------------------------------------------------------------

% Materials
% Constitutive model parameters
E  = 210e6 ; nu = 0.3 ; % 
sigmaY = 250e3 ; 
K = E/3 ; % Linear hardening parameter 

materials(1).hyperElasModel = 'biLinear' ;
materials(1).hyperElasParams = [ E, nu, sigmaY, K ] ;

% Elements
elements(1).elemType  = 'node'  ;
elements(2).elemType  = 'frame' ;

% Section properties
ty = 0.1 ; 
tz = 0.1 ; % cross-section widths
elements(2).elemCrossSecParams = { 'rectangle'; [ ty tz ] };
elements(2).massMatType     =  'consistent' ; 

% BoundaryConditions
% Supports
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;

% Loads
% Applied nodal loads
P = 5 ; % applied nodal load
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 -P 0 ] ;

% InitialConditions
% empty struct
initialConds = struct() ;

% Nodes
% ----------------------------------------------------------------------
L = 1 ; %
nnodesMesh = 11 ;
xcoords = linspace(0,L,nnodesMesh)' ;
ycoords = zeros(length(xcoords),1) ;
zcoords = zeros(length(xcoords),1) ;
 
mesh.nodesCoords = [ xcoords ycoords zcoords ] ;

% Conec Cell
% ----------------------------------------------------------------------
mesh.conecCell = { } ;

% Auxiliar
vec1 = (1:1:(length(xcoords)-1))' ;
vec2 = (2:1:(length(xcoords)))' ;
vec3 = ones(length(vec1),1) ;
loadedNode = length(xcoords) ;

% Node elements
mesh.conecCell{1, 1 } = [ 0 1 1 0  1 ] ;
mesh.conecCell{2, 1 } = [ 0 1 2 0  loadedNode ] ;

% Frame elements
for i=1:(nnodesMesh-1) 
	mesh.conecCell{i+2,1} = [ 1 2 0 0 i i+1 ] ;
end

% Call ONSAS Function
Conec = myCell2Mat( mesh.conecCell ) ;
Conec(find(Conec(:,1)==0),:) = [] ;

% Analysis settings
% ----------------------------------------------------------------------
% NR settings
analysisSettings = struct() ;
analysisSettings.finalTime = 20 ;
analysisSettings.stopTolIts = 15 ;

% Other parameters
% ----------------------------------------------------------------------
otherParams.problemName = 'cantileverLinearHardening' ;
otherParams.plots_format = 'vtk' ;

% Auxiliar variables
% ----------------------------------------------------------------------

nelems 				= nnodesMesh-1 	;
ndofs 				= 6 						;
numNodes 			= 2 						;
LocBendXZdofs = [ 5 4 11 10 ] ;

% Global variables
% ----------------------------------------------------------------------

global ne
	ne = 10 ;

% ONSAS Solve
% ----------------------------------------------------------------------
[matUs, loadFactorsMat, matFint ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

% Post process
% ----------------------------------------------------------------------
kappaHistElem = zeros(nelems,analysisSettings.finalTime+1) ;
RXYXZ = eye(4) ; RXYXZ(2,2) = -1; RXYXZ(4,4) = -1;
for j=1:nelems
	nodeselem   = Conec( j, (4+1):(4+numNodes) )' ;
  dofselem    = nodes2dofs( nodeselem , ndofs )        ;
	elemCoords  = reshape( mesh.nodesCoords( Conec( j, (4+1):(4+numNodes) )' , : )', 1, 3*numNodes) ;
	
	[ loc, l ] = beamParameters( elemCoords ) ;
	R = RotationMatrix(ndofs, loc) ;
	R = R(LocBendXZdofs,LocBendXZdofs) ;
	Uke = R'*matUs(dofselem(LocBendXZdofs),1:end) ;  
	Be = bendingInterFuns (0 , l, 2 ) * RXYXZ ;
	kappae = Be*Uke ;
	kappaHistElem(j,:) = kappae ;
end

% M-k auxiliar 
elem 	= 1 ;

kappaVec = abs(kappaHistElem(elem,1:end)) ;
aux  = cell2mat(matFint) ;
mVec = abs(aux(elem,:)) ;
mVec = mVec(4:24:end) ;

% Analytical solution
% ----------------------------------------------------------------------
I = ty*tz^3/12 ;
My = sigmaY * ty*tz^2/6 ;
kappa_e = 2*sigmaY/(E*tz) ;

Ma = zeros(1,length(kappaVec)) ;
C = E*K / (E+K) ;
for i=1:length(kappaVec)
	if kappaVec(i) < kappa_e
		Ma(i) = E*I*kappaVec(i) ;
	else	
		Ma(i) = sigmaY*ty*tz^2/12 * (3 - kappa_e^2/kappaVec(i)^2 + kappaVec(i)/kappa_e*C/E*( 2 - 3*kappa_e/kappaVec(i) + kappa_e^3/kappaVec(i)^3 ) ) ;
	end
end

% Plots
% ----------------------------------------------------------------------

% Sets
lw1 = 2.0 ; ms1 = 11 ; plotFontSize = 22 ;

plotLoadVec = [1:(analysisSettings.finalTime+1)] ;
scaleFactor = 1 ;

Nodes = mesh.nodesCoords ;

M_kappa = figure ;
hold on, grid on

% NR
plot(kappaVec/kappa_e, mVec/My, 'b-o', 'linewidth', lw1, 'markersize', ms1)
% Ana
plot(kappaVec/kappa_e, Ma/My, 'g-x', 'linewidth', lw1, 'markersize', ms1)

legend ({'NR', 'Analytic'}, 'location', 'east');
labx = xlabel('\kappa / \kappa_e'); laby = ylabel('M / M_y') ;

tit = title('M-\kappa');
set(labx, 'fontsize', plotFontSize*.8);
set(laby, 'fontsize', plotFontSize*.8);
set(tit, 'fontsize', plotFontSize);

verifBoolean = max ((mVec(2:end)-Ma(2:end))./Ma(2:end) ) < 0.03 ;
