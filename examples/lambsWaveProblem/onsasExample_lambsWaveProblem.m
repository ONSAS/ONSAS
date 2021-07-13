%md## Lamb's wave problem
%md

clear all, close all
% add path
addpath( genpath( [ pwd '/../../src'] ) ) ;
% scalar parameters
E = 1.8773e10 ; nu = 0.25 ; thickness = 1 ; rho = 2200 ;

global spitMatrices
spitMatrices = true


%md### MEBI parameters

%md#### materials

%md since only one material is considered, the structs defined for the materials contain only one entr
materials.hyperElasModel  = {'linearElastic'} ;
materials.hyperElasParams = { [ E nu ] }      ;
materials.density         = { rho }      ;

%md#### elements
%md In this model two kinds of elements are used: tetrahedrons for the solid and triangles for introducing the external loads. Since two kinds of elements are used, the structs have length 2:
elements.elemType = { 'node', 'edge', 'triangle' } ;
%md since triangle and tetrahedron elements dont have specific parameters the struct entries contain empty vectors
elements.elemTypeParams = { []; [] ; 2  } ;
elements.elemTypeGeometry = { []; thickness ; thickness } ;
%md
%md#### boundaryConds
%md in this case four BCs are considered, one corresponding to a load and three to displacements.
boundaryConds.loadsCoordSys = {[]     ; []  ; 'global'  } ;
boundaryConds.loadsTimeFact = { []    ; []  ; @(t) 1e6*( 1*( t <= .15 ) - 3*( t <= .1 ) + 3*( t <= .05 ) ) } ;
boundaryConds.loadsBaseVals = { []    ; []  ; [ 0 0  1 0  0 0 ]  } ;
boundaryConds.imposDispDofs = { [1 3] ; [1] ; []  } ;
boundaryConds.imposDispVals = { [0 0] ; [0] ; []  } ;
%md
%md#### initialConds
%md since no initial non-homogeneous initial conditions are used, an empty struct is used .
initialConds = struct();
%md

%md### Mesh
%md An 8-node mesh is considered with its connectivity matrix
%md
%md```@raw html
%md<img src="https://raw.githubusercontent.com/ONSAS/ONSAS_docs/master/docs/src/solidCubeMeshHTML.svg" alt="structure diagram" width="500"/>
%md```
%md```@raw latex
%md\begin{center}
%md\def\svgwidth{0.6\textwidth}
%md\input{solidCubeMeshPDF.pdf_tex}
%md\end{center}
%md```

[ mesh.nodesCoords, mesh.conecCell ] = meshFileReader( 'lambsProblemDomain.msh' ) ;

%md
%md### Analysis parameters
%md
nAproxElem = length( mesh.conecCell )

CFL = 0.125 ;
cP  = 3200  ;
h   = sqrt( 2*3200^2 / nAproxElem ) * .5
dt  = CFL * h / cP

analysisSettings.methodName    = 'newmark' ;
analysisSettings.stopTolIts    = 30      ;
analysisSettings.stopTolDeltau = 1.0e-12 ;
analysisSettings.stopTolForces = 1.0e-12 ;
analysisSettings.finalTime      = 1    ;
analysisSettings.alphaNM       = 0.25    ;
analysisSettings.deltaNM       = 0.5     ;
analysisSettings.deltaT        = dt      ;
%md
%md
%md### Output parameters
otherParams.problemName = 'lambsProblem' ;
otherParams.plotsFormat = 'vtk' ;
%md
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md


return

% Analytic solution
p = abs(p) ; Rext = 0.15 ; Rint = 0.1 ;
a = ( (1+nu)*(1-2*nu)*Rint^2 ) / ( E*(Rext^2-Rint^2) ) ;
b = ( (1+nu)*Rint^2*Rext^2 )   / ( E*(Rext^2-Rint^2) ) ;
analyticCheckTolerance = 1e-3 ;
analyticValRInt = a*Rint + b/Rint ;
analyticValRExt = a*Rext + b/Rext ;
verifBoolean = ( ( matUs(         1, 2) - analyticValRInt ) < analyticCheckTolerance ) && ...
               ( ( matUs( (8-1)*6+3, 2) - analyticValRExt ) < analyticCheckTolerance )

%md
return
%md### Results
%md
%md```math
%md\lambda(t) = \frac{1}{p} \frac{E}{2}  \left( \left( 1+\frac{u}{Lx} \right)^3 - \left( 1+ \frac{u}{Lx} \right) \right)
%md```

analyticCheckTolerance = 1e-6 ;
analyticFunc           = @(w) 1/p * E * 0.5 * ( (1 + w/Lx).^3 - (1+w/Lx) ) ;
disps = matUs(6*6+1,:) ;
analyticVals = analyticFunc(disps) ;
%
verifBoolean = ( norm( analyticVals - loadFactorsMat') / norm( analyticVals) ) < analyticCheckTolerance

%md
%md### plot
%md
lw = 2.0 ; ms = 11 ; plotfontsize = 22 ;
figure, hold on, grid on
plot( disps, loadFactorsMat, 'k-o' , 'linewidth', lw,'markersize',ms )
plot( disps, analyticVals, 'b-x' , 'linewidth', lw,'markersize',ms )
labx = xlabel('Displacement');   laby = ylabel('$\lambda$') ;
legend('Numeric','Analytic','location','North')
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
print('verifUniaxial.png','-dpng')
%md

return

% --------------------------------------------------------

clear iniMatUs


controlDispsValsCase1         = controlDisps  ;
loadFactorAnalyticalValsCase1 = analyticVals  ;
loadFactorNumericalValsCase1  = numericalVals ;

close all

% --------------------------------------------------------
% solid model using gmsh mesh, local tension load and complex step
% --------------------------------------------------------

problemName = 'uniaxialExtension_GMSH_ComplexStep' ;

[ Nodes, Conec ] = meshFileReader( 'geometry_uniaxialExtension.msh' ) ;

loadsParams{1,1}    = [ 0 1  0 0 0 0 p 0 ] ; % local coords appliend tension

elementsParams{2,1} = [ 4 1 ] ; % complex step constitutive tensor

plotParamsVector = [ 0 ] ;
analyticSolFlag        = 0 ;

% run ONSAS
ONSAS

controlDispsValsCase2         = controlDisps  ;
loadFactorNumericalValsCase2  = numericalVals ;


% --------------------------------------------------------
% truss element model
% --------------------------------------------------------

problemName = 'uniaxialExtension_truss' ;

Nodes = [ 0    0    0 ; ...
          Lx   0    0   ...
        ] ;

Conec = {[ 0 1 0 0 1   1   ] ; ... % fixed node
         [ 0 1 1 0 2   2   ] ; ... % loaded node
         [ 1 2 0 1 0   1 2 ]   ... % truss element
        } ;

% ======================================================================
% --- MELCS parameters ---

materialsParams = cell(1,1) ; % M
elementsParams  = cell(1,1) ; % E
loadsParams     = cell(1,1) ; % L
crossSecsParams = cell(1,1) ; % C
springsParams   = cell(1,1) ; % S

% --- Material parameters ---
E = 1 ; nu = 0.3 ;
materialsParams{1} = [ 0 2 E nu ] ;

% --- Element parameters ---
elementsParams = { 1  ; [ 2 0 ]} ;

% --- Load parameters ---
loadsParams{1,1} = [ 1 1  p 0 0 0 0 0 ] ;

% --- CrossSection parameters ---
crossSecsParams = { [ 2 Ly Lz] } ; %

% ----------------------------------------------------------------------
% --- springsAndSupports parameters ---
springsParams{1, 1} = [ inf 0  inf 0   inf 0 ] ;
springsParams{2, 1} = [ 0   0  inf 0   inf 0 ] ;

% ======================================================================

plotParamsVector       = [ 0 ] ;


controlDofs = [ 2 1 1 ] ;

%% run ONSAS
ONSAS

controlDispsValsCase3         = controlDisps  ;
%~ loadFactorNumericalValsCase3  = numericalVals .* (1+controlDisps) / Lx ;
loadFactorNumericalValsCase3  = numericalVals ;


% ----------------------------------------------------------------------
% --- plots ---
lw = 2.0 ; ms = 10 ; plotfontsize = 22 ;

figure, grid on, hold on

plot( controlDispsValsCase1, ...
      loadFactorAnalyticalValsCase1 ,'b-o' , 'linewidth', lw,'markersize',ms )

plot( controlDispsValsCase1, ...
      loadFactorNumericalValsCase1  ,'k-s' , 'linewidth', lw,'markersize',ms)

plot( controlDispsValsCase2, ...
      loadFactorNumericalValsCase2  ,'r-x' , 'linewidth', lw,'markersize',ms)

plot( controlDispsValsCase3, ...
      loadFactorNumericalValsCase3  ,'g--' , 'linewidth', lw,'markersize',ms)

labx = xlabel('Displacement');   laby = ylabel('$\lambda$') ;
legend('analytic Sol','numerical Sol 1','numerical Sol 2','numerical Sol 3','location','North')
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
%~ print( [ 'plotsExtensionSVK' ] ,'-dpdflatex','-tight') ;
print( [ 'plotsExtensionSVK.png' ] ,'-dpng') ;

% ----------------------------------------------------------------------
