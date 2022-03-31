% Static Arclength Analysis of Roorda Frame
% 2021.10.06
close all, clear all ; addpath( genpath( [ pwd '/../../src'] ) );

% structural parameters
E = 210e3 ;	% Young's modulus (MPa)
nu = 0.3  ;	% Poisson ratio
a = 50 ; 		% Rectangular section width (mm)
b = 10 ; 		% Rectangular section height (mm)
A = a*b; 		% Rectangular section Area (mm2)
L = 1000 ;	% Frame members length (mm)
m = 4; 			% number of finite elements per member

vececc = [ 1 .25 ] ;
vecpAL = [ 2 .5 ] ;

figure, grid on, hold on

for indecc = 1:length(vececc)
  ecc = vececc(indecc);

%Material Definitions
materials(1).hyperElasModel  = '1DrotEngStrain' ;
materials(1).hyperElasParams = [ E nu ] ;

%Element Definitions
elements(1).elemType = 'node' ;
elements(2).elemType = 'frame';
elements(2).elemCrossSecParams = {'rectangle'; [2 b b ] } ;


%Boundary Conditions Definitions
boundaryConds(1).imposDispDofs = [ 1 3 5 ] ;	% pinned node
boundaryConds(1).imposDispVals = [ 0 0 0 ] ;
boundaryConds(2).imposDispDofs =  3 ;			% node restrained in X-Y plane
boundaryConds(2).imposDispVals =  0 ;
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsTimeFact = @(t) t ;

boundaryConds(2).loadsBaseVals = 12.8e3*[ 0 0 0 ecc*1 -1 0 ] ;	% Fx = -12.8e3 (N), My = 12.8e3*ecc (Nmm)
																% Note that LBA Load is: Pcr = 12.8e3 (N)

%Nodal Coordinates Definitions. X(mm)	Y(mm)	Z(mm)
dl = L/m;												% length of finite elements (mm)

for i=1:m+1
	mesh.nodesCoords(i,:) = [   0 	 0  (i-1)*dl ] ; 	% Coords of Nodes in vertical member, incl corner node #(m+1)
end

for i=m+2:2*m+1
	mesh.nodesCoords(i,:) = [   (i-m-1)*dl 0  L  ] ;  	% Coords of Nodes in horizontal member, excl corner node #(m+1)
end

%Conectivity Definitions
mesh.conecCell = cell(5,1) ;
					   %M E B I / Node
mesh.conecCell{ 1 } = [ 0 1 1 0  1   ] ; 	% Node at coord (0,0,0)
mesh.conecCell{ 2 } = [ 0 1 1 0  2*m+1 ] ;	% Node at coord (L,0,L)
mesh.conecCell{ 3 } = [ 0 1 2 0  m+1   ] ;	% Corner Node   (0,0,L)

for j=1:2*m					   %M E B I / Nodes
	mesh.conecCell{ 4+j-1 } = [ 1 2 0 0  j j+1 ] ; % frame finite elements
end

initialConds = [];	% no initial conditions for static analysis

%Static Analysis Parameter Definitions
analysisSettings.methodName    = 'arcLength' ;
analysisSettings.deltaT        =   0.02  ;
analysisSettings.finalTime     =   1.1   ;
analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   30   ;

analysisSettings.finalTime     = 1.1 ;
analysisSettings.incremArcLen = vecpAL(indecc) ;
analysisSettings.iniDeltaLamb = boundaryConds(2).loadsTimeFact(1)/100 ;
analysisSettings.posVariableLoadBC = 2 ;

otherParams.problemName = 'staticRoordaFrame_Stable_AL';
otherParams.plotsFormat = 'vtk' ;

%Analysis case 1: Solution of Stable Branch with ArcLength
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
controlDispsALstab =  -matUs(6*m+4,:) ; % rotation wrt Y-axis at node m+1
loadFactorsALstab  =  loadFactorsMat(:,2) ;


%Analysis case 2: Solution of Unstable Branch with ArcLength
otherParams.problemName       = 'staticRoordaFrame_Unstable_AL' ;

boundaryConds(2).loadsBaseVals = 12.8e3*[ 0 0 0 -ecc -1 0 ] ; % choose opposite eccentricity

[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
controlDispsALunstab =  -matUs(6*m+4,:) ; % rotation wrt Y-axis at node m+1
loadFactorsALunstab  =  loadFactorsMat(:,2) ;

%Plot Load Displacement Curves
lw = 1.0 ; ms = 5 ; plotfontsize = 12 ;

plot( controlDispsALstab, loadFactorsALstab, 'k-x' , 'linewidth', lw,'markersize',ms )
plot( controlDispsALunstab, loadFactorsALunstab, 'r-o' , 'linewidth', lw,'markersize',ms )

labx = xlabel('rotation @Corner node (rad)');   laby = ylabel('\lambda(t)') ;
legend( 'stable branch', 'unstable branch','location','southeast') ;
set(gca, 'linewidth', 1.0, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Roorda Frame / Load-Displacement curves / With imperfectons') ;
grid on ;

end

print('output/RoordaFrame.png','-dpng')
