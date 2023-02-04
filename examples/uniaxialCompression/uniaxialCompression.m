%md# Uniaxial Compression example
%md
%md[![Octave script](https://img.shields.io/badge/script-url-blue)](https://github.com/ONSAS/ONSAS.m/blob/master/examples/uniaxialCompression/uniaxialCompression.m)
%md
%mdIn this tutorial example an hyper-elastic solid is submitted to a uniaxial compression test. The geometry and tension applied are shown in the figure, where the $Lx$, $Ly$ and $Lz$ are the dimensions, and the compression $p$ is applied on the face $x=Lx$, as nominal tension.
%md
%md```@raw html
%md<img src="../../assets/uniaxialCompression/diagramSolidUniaxialCompression.svg" alt="structure diagram" width="500"/>
%md```
%md---
%md
%md
%md ## Analytic solution
%md
%mdLet us consider that a uniform deformation is produced, with a nonzero axial stretch $\alpha$ and nonzero transversal stretch $\beta$. The corresponding deformation gradient and Green-Lagrange strain tensor are given by:
%md
%md```math
%md\textbf{F} = \left[ \begin{matrix} \alpha & 0 & 0 \\ 0 & \beta & 0 \\ 0 & 0 & \beta \end{matrix} \right]
%md\qquad
%md\textbf{E} = \frac{1}{2}(\textbf{C} - \textbf{I}) = \left[  \begin{matrix} \frac{1}{2} \left(\alpha^2 -1\right) & 0 & 0 \\ 0 &  \frac{1}{2} \left(\beta^2 -1\right) & 0 \\ 0 & 0 &  \frac{1}{2} \left(\beta^2 -1\right) \end{matrix}
%md\right]
%md```
%mdThe neo-Hookean elastic strain energy potential $\Psi$ is given by:
%md```math
%md\Psi(\mathbf{C=\textbf{F}^T\textbf{F}}) = \frac{\mu}{2}(I_1 -3) + K ( J^2 -1 )
%md```
%mdwhere $I_1 = \mathrm{tr}(\mathbf{C})$ is the first invariant, $J = \sqrt{\det(\mathbf{C})}$ and $K$ and $\mu$ are the bulk and shear material parameters, respectively.
%md
%mdThe second Piola-Kirchhoff tensor is given by:
%md```math
%md\textbf{S}( \mathbf{C} ) = \mu (\textbf{I} - \mathbf{C}^{-1}) + K (J(J-1)\mathbf{C}^{-1}) 
%md```
%mdthen, using the relation $\textbf{P}=\textbf{F}\textbf{S}$, the $P_{yy}=P_{zz}$ component is computed: and set to zero (using the boundary conditions)
%md```math
%mdP_{yy}( \mu,K ) =
%md2\beta \left( 
%md    \frac{\mu}{2} -  \frac{\mu}{\beta^2} + K (\alpha^2\beta^2 -\alpha) 
%md \right) = 0
%md```math
%mdthus, considering the axial displacement $u_x$ and $u_y$ at $\mathbf{X}=(L_x,L_y,L_z)$ and using the stretch definitions $\alpha = (1+u_x/Lx)$ and $\beta = (1+u_y/Ly)$, we obtain the axial component of the nominal stress:
%md```math
%mdP_{xx}( \mu,K ) =
%md2\alpha \left( 
%md    \frac{\mu}{2} -  \frac{\mu}{\alpha^2} + \frac{K\beta^2}{\alpha} (\beta^2 \alpha -1) 
%md \right) = p
%md```
%md
%md## Numerical solution: case 1
%md---
%mdBefore defining the structs, the workspace is cleaned, the ONSAS directory is added to the path and scalar geometry and material parameters are defined.
clear all, close all
% add path
addpath( genpath( [ pwd '/../../src'] ) ) ;
% scalar parameters
E = 1 ; nu = 0.3 ; p = -5 ; Lx = 2 ; Ly = 1 ; Lz = 1 ;
%md
%md
%md### MEBI parameters
%md
%md#### materials
%md The material of the solid considered is a neo-Hookean model with $\lambda$, $\mu$ and bulk($K$) parameters:
lambda = E*nu/((1+nu)*(1-2*nu)) ; mu = E/(2*(1+nu)) ; bulk = E / ( 3*(1-2*nu) ) ;
%md since only one material is considered, a scalar struct is defined as follows
materials.hyperElasModel = 'NHC' ;
materials.hyperElasParams = [ mu bulk ] ;
%md
%md#### elements
%md In this model two kinds of elements are used: `tetrahedron` for the solid and `triangle` for introducing the external loads. Since two kinds of elements are used, the struct have length 2:
elements(1).elemType = 'triangle' ;
elements(2).elemType = 'tetrahedron' ;
elements(2).elemTypeParams = [ 2 ] ;
%md
%md#### boundaryConds
%md in this case four BCs are considered, one corresponding to a load and three to displacements.
%md the first BC introduced is a load, then the coordinate system, loadfactor time function and base load vector are defined
boundaryConds(1).loadsCoordSys = 'global';
boundaryConds(1).loadsTimeFact = @(t) p*t ;
boundaryConds(1).loadsBaseVals = [ 1 0 0 0 0 0 ] ;
%md the other BCs have imposed displacements
boundaryConds(2).imposDispDofs = [1] ;
boundaryConds(2).imposDispVals =  0  ;
%
boundaryConds(3).imposDispDofs = [3] ;
boundaryConds(3).imposDispVals =  0  ;
%
boundaryConds(4).imposDispDofs = [5] ;
boundaryConds(4).imposDispVals =  0  ;
%
%md
%md#### initialConds
%md since no initial non-homogeneous initial conditions are used, an empty struct is used .
initialConds = struct();
%md
%md### Mesh
%md A simple hand-made 8-node mesh, with 6 tetrahedrons is considered
%md
%md```@raw html
%md<img src="https://raw.githubusercontent.com/ONSAS/ONSAS_docs/master/docs/src/solidCubeMeshHTML.svg" alt="mesh diagram" width="500"/>
%md```
%md```@raw latex
%md\begin{center}
%md\def\svgwidth{0.6\textwidth}
%md\input{solidCubeMeshPDF.pdf_tex}
%md\end{center}
%md```
%md The node coordinates matrix is given by the following
mesh.nodesCoords = [ 0    0    0 ; ...
                     0    0   Lz ; ...
                     0   Ly   Lz ; ...
                     0   Ly    0 ; ...
                     Lx   0    0 ; ...
                     Lx   0   Lz ; ...
                     Lx  Ly   Lz ; ...
                     Lx  Ly    0 ] ;
%md and the connectivity cell is defined as follows with the four MEBI parameters for each element followed by the indexes of the nodes of each element. All the eight triangle elements are considered with no material (since they are used only to include load) and the following six elements are solid neo-Hookean material tetrahedrons.
mesh.conecCell = {[ 0 1 1 0    5 8 6   ]; ... % loaded face
                  [ 0 1 1 0    6 8 7   ]; ... % loaded face
                  [ 0 1 2 0    4 1 2   ]; ... % x=0 supp face
                  [ 0 1 2 0    4 2 3   ]; ... % x=0 supp face
                  [ 0 1 3 0    6 2 1   ]; ... % y=0 supp face
                  [ 0 1 3 0    6 1 5   ]; ... % y=0 supp face
                  [ 0 1 4 0    1 4 5   ]; ... % z=0 supp face
                  [ 0 1 4 0    4 8 5   ]; ... % z=0 supp face
                  [ 1 2 0 0    1 4 2 6 ]; ... % tetrahedron
                  [ 1 2 0 0    6 2 3 4 ]; ... % tetrahedron
                  [ 1 2 0 0    4 3 6 7 ]; ... % tetrahedron
                  [ 1 2 0 0    4 1 5 6 ]; ... % tetrahedron
                  [ 1 2 0 0    4 6 5 8 ]; ... % tetrahedron
                  [ 1 2 0 0    4 7 6 8 ]  ... % tetrahedron
                } ;
%md
%md### Analysis parameters
%md
analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.stopTolIts    = 30     ;
analysisSettings.stopTolDeltau = 1.0e-8 ;
analysisSettings.stopTolForces = 1.0e-8 ;
analysisSettings.finalTime      = 1      ;
analysisSettings.deltaT        = .1     ;
%md
%md### Output parameters
otherParams.plotsFormat = 'vtk' ;
otherParams.problemName = 'uniaxialCompression_HandMadeMesh' ;
%md
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md

controlDispsValsCase1 = matUs(6*6+1,:) ;
loadFactorsCase1 = loadFactorsMat ;
%md
%md## Numerical solution: case 2
%mdIn this analysis case, the mesh information is read from a gmsh-generated
%mdmesh file, the pressure is applied using local coordinates and the stiffness
%md matrix is computed using the complex-step method.
%md
otherParams.problemName = 'uniaxialCompression_GMSH_ComplexStep' ;
[ mesh.nodesCoords, mesh.conecCell ] = meshFileReader( 'geometry_uniaxialExtension.msh' ) ;
boundaryConds(1).loadsCoordSys = 'local';
boundaryConds(1).loadsBaseVals = [0 0 0 0 1 0 ] ;

[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

%md
%md### Analytic solution computation
%md
%md Displacements of node 7 are extracted:
%md
controlDispsValsCase2         = matUs(6*6+1,:) ;  ;
loadFactorsCase2  = loadFactorsMat ;
%md
%md The values of $\beta$ and $\alpha$ for each load step are compted:
%md
alphas         = (Lx + matUs(6*6+1,:)) / Lx ;
betas          = (Ly + matUs(6*6+3,:)) / Ly ;
analyticFunc      = @(alphas,betas) mu * alphas - mu*1./alphas + bulk * betas.^2 .* ( alphas .* betas.^2 -1) ;
analyticVals = analyticFunc( alphas, betas ) ;
%md## Plot
%mdThe numerical and analytic solutions are plotted.
lw = 2.0 ; ms = 11 ; plotfontsize = 18 ;
figure, hold on, grid on
plot( controlDispsValsCase1, loadFactorsCase1, 'r-x' , 'linewidth', lw,'markersize',ms )
plot( controlDispsValsCase2, loadFactorsCase2,  'k-o' , 'linewidth', lw,'markersize',ms )
plot( controlDispsValsCase2, analyticVals,  'g-s' , 'linewidth', lw,'markersize',ms )
labx = xlabel('Displacement');   laby = ylabel('\lambda(t)') ;
legend( 'Numeric-1', 'Numeric-2', 'Analytic-NHC' , 'location', 'SouthEast' )
set(gca, 'linewidth', 1.0, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
%md
%md```@raw html
%md<img src="https://raw.githubusercontent.com/ONSAS/ONSAS.docs/master/docs/src/verifCompression.png" alt="plot check" width="500"/>
%md```
