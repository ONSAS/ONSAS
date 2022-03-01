
close all, clear all ; addpath( genpath( [ pwd '/../../src'] ) );

% structural parameters
E  = 210e9 ;		% Young's modulus (Pa)
nu = 0.3 ;		% Poisson ratio
a  = .05; 		% Rectangular section width (mm)
b  = .01; 		% Rectangular section height (mm)
A  = a*b; 		% Rectangular section Area (mm2)
R  = 1 ;
fiApoyo = 30 ;

nEM = 10; 			% numero elementos finitos mitad de arco

%Material Definitions
materials(1).hyperElasModel  = '1DrotEngStrain' ;
materials(1).hyperElasParams = [ E nu ] ;

%Element Definitions
elements(1).elemType = 'node' ;
elements(2).elemType = 'frame';
elements(2).elemCrossSecParams{1,1} = 'rectangle' ;
elements(2).elemCrossSecParams{2,1} = [a b ]      ;

%Boundary Conditions Definitions
boundaryConds(1).imposDispDofs = [ 1 3 5 ] ;	% pinned node
boundaryConds(1).imposDispVals = [ 0 0 0 ] ;
boundaryConds(2).imposDispDofs =  3 ;			% node restrained in X-Y plane
boundaryConds(2).imposDispVals =  0 ;
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsTimeFact = @(t) t ;

boundaryConds(2).loadsBaseVals = 1e4*[ .001 0 0 0 -1 0 ] ;	% Fx = -12.8e3 (N), My = 12.8e3*ecc (Nmm)
																% Note that LBA Load is: Pcr = 12.8e3 (N)

%Nodal Coordinates Definitions. X(mm)	Y(mm)	Z(mm)

anguloTotalArco = (180-fiApoyo) - fiApoyo
deltaTheta = anguloTotalArco / ( 2*nEM )

for i=1:(2*nEM+1)
	anguloi = (180-fiApoyo) - (i-1)*deltaTheta ;
	mesh.nodesCoords(i,:) = [  R*cos( pi / 180 * anguloi )  0 R*sin( pi / 180 * anguloi )   ] ; 	% Coords of Nodes in vertical member, incl corner node #(m+1)
end

largo = 0 ;
for j=1:2*nEM
  largo = largo + norm( mesh.nodesCoords(j+1,:) - mesh.nodesCoords(j,:) );
end
perimArcoTeo = R*anguloTotalArco*pi/180 ;

%Conectivity Definitions
mesh.conecCell = cell(5,1) ;
         					   %M E B I / Node
mesh.conecCell{ 1 } = [ 0 1 1 0  1   ] ; 	% Node at coord (0,0,0)
mesh.conecCell{ 2 } = [ 0 1 1 0  2*nEM+1 ] ;	% Node at coord (L,0,L)
mesh.conecCell{ 3 } = [ 0 1 2 0  nEM+1   ] ;	%

for j=1:2*nEM					   %M E B I / Nodes
	mesh.conecCell{ 4+j-1 } = [ 1 2 0 0  j j+1 ] ; % frame finite elements
end

initialConds = [];	% no initial conditions for static analysis

%Static Analysis Parameter Definitions
analysisSettings.methodName    = 'arcLength' ;
analysisSettings.deltaT        =   0.02  ;
analysisSettings.finalTime      =   1   ;
analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   30   ;

analysisSettings.incremArcLen = .02 ;
analysisSettings.iniDeltaLamb = boundaryConds(2).loadsTimeFact(1)/100 ;
analysisSettings.posVariableLoadBC = 2 ;

otherParams.problemName = 'circularArch_AL';
otherParams.plotsFormat = 'vtk' ;

%Analysis case 1: Solution of Stable Branch with ArcLength
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
controlDispsALstab =  -matUs(6*nEM+5,:) ; % rotation wrt Y-axis at node m+1
loadFactorsALstab  =  loadFactorsMat(:,2) ;

return
%Plot Load Displacement Curves
lw = 1.0 ; ms = 5 ; plotfontsize = 12 ;
figure
plot( controlDispsALstab, loadFactorsALstab, 'k-' , 'linewidth', lw,'markersize',ms )
hold on
plot( controlDispsALunstab, loadFactorsALunstab, 'r-' , 'linewidth', lw,'markersize',ms )
hold off

labx = xlabel('rotation @Corner node (rad)');   laby = ylabel('\lambda(t)') ;
legend( 'stable branch', 'unstable branch','location','southeast') ;
set(gca, 'linewidth', 1.0, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Roorda Frame / Load-Displacement curves / With imperfectons') ;
grid on ;
print('output/RoordaFrame.png','-dpng')
