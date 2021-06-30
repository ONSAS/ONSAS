%md## Example Linear Plane Strain
%md
%mdIn this example a hollow cylinder submitted to an internal pressure $p_i=0.1$ GPa as shown in the Figure is considered.
%\begin{figure}[htb]
% 	\centering
% 	\includegraphics[width=350px]{./images/tikzCylinder.png}
% 	\caption{Diagram of cylinder submitted to internal pressure.}
% 	\label{fig:infcylinder}
% \end{figure}
% %
% The material is considered isotropic and homogeneous with elasticity modulus $E=210$ GPa and Poisson's coefficient $\nu=0.3$. The length of the cylinder is $L_z = 0.75$ m and the internal and external radii are $R_i=0.2$ m and $R_e=0.24$ m, respectively. %
%
% The pressure is radial and applied on the internal surface, given by $\bfp =  - p_i \, \bfe_r$. The displacement vector $\bfu$ is $\bfu=u_r\bfe_r + u_z\bfe_z $. The boundary conditions correspond to a plane strain state:
%
% \begin{eqnarray}
% &u_z(x,y,z=0)=0\;\;\quad\forall\; x,y \\
% &u_z(x,y,z=L_z)=0\;\;\quad\forall\; x,y
% \end{eqnarray}
%
%
% %\paragraph{The analytical solution}
% The analytical solution is obtained considering the Navier equation:
%
% \begin{equation}\label{eq:navier}
% (\lambda+\mu)\triangledown(\triangledown\cdot \bfu ) + \mu\triangle \bfu +\bfb = (3\lambda + 2\mu)\alpha\triangledown\theta
% \end{equation}
%
% Due to the symmetry of the problem, the following identities can be considered:
%
% \begin{equation}
% \begin{array}{l}\label{eq:cylinderrel}
% u_r = u_r(r) \\
% u_z = u_z(z) \\
% \triangle \bfu = \triangledown(\triangledown\cdot\bfu)
% \end{array}
% \end{equation}
%
% Given the identities in \autoref{eq:cylinderrel} the Navier equation is reduced to:
%
% \begin{equation}
% (\lambda+2\mu)\triangledown(\triangledown\cdot\bfu)+\bfb = (3\lambda+2\mu)\alpha\triangledown\theta
% \end{equation}
%
% Solving the differential equation, applying the displacement and force conditions and using the constitutive relationship, the solution displacement vector $\bfu$ can be obtained. The components $u_r$ and $u_z$ are given by:
%
% \begin{equation}\label{eq:infcylinder}
% u_r(r) = Ar + \dfrac{B}{r} \quad \text{and} \quad u_z(z) = 0,
% \end{equation}
% where:
% \begin{equation}
% A = \dfrac{(1+\nu)(1-2\nu)R_i^2p_i}{E(R_e^2-R_i^2)}, \quad
% B = \dfrac{(1+\nu)R_i^2R_e^2p_i}{E(R_e^2-R_i^2)}
% \end{equation}
%
% %\paragraph{The numerical results}
% The numerical results were obtained using \href{https://github.com/onsas/onsas/blob/development/input/onsas_input_TEST_Cylinder.m}{this} input file.
% %
% Due to the symmetry of the problem, only a quarter of cylinder is modeled. %
% The Finite Elements mesh is formed by tetrahedron elements and it was generated using the open-source software gmsh \cite{Geuzaine2009a}. %
%
% The output vtk files can be visualized using Paraview, showing for instance displacements in \autoref{fig:cylinderOnsasOutput} or stresses \autoref{fig:cylinderOnsasOutputVM}. %
%
% \begin{figure}[htb]
% 	\centering
% 	\includegraphics[width=500px]{./images/cylinderSolid.png}
% 	\caption{Undeformed and deformed configuration of the cylinder with a scale factor.}
% 	\label{fig:cylinderOnsasOutput}
% \end{figure}
%
% \begin{figure}[htb]
% 	\centering
% 	\includegraphics[width=500px]{./images/cylinderVM.png}
% 	\caption{Von Mises stress in reference configuration.}
% 	\label{fig:cylinderOnsasOutputVM}
% \end{figure}
%
% The magnitude of the displacements are represented by the deformed color scale and the reference configuration is shown in gray with blue (triangle edges). %
%
% The displacements provided by ONSAS match the analytical solution.
%



% Analytic solution
% analyticSolFlag = 3 ; p = abs(p) ; r = Rext ;
% a = ( (1+nu)*(1-2*nu)*Rint^2*p ) / ( E*(Rext^2-Rint^2) ) ;
% b = ( (1+nu)*Rint^2*Rext^2*p )   / ( E*(Rext^2-Rint^2) ) ;
% analytSol = a*r + b/r ; analyticSolDofs = [ 6*(6-1)+1 ] ;
% analyticCheckTolerance = 1e-3 ;



%md### Numerical solution
%mdBefore defining the structs, the workspace is cleaned, the ONSAS directory is added to the path and scalar geometry and material parameters are defined.
clear all, close all
% add path
addpath( genpath( [ pwd '/../../src'] ) ) ;
% scalar parameters
E = 1 ; nu = 0.3 ; p = 3 ; Lx = 2 ; Ly = 3 ; thickness = 0.5 ;

%md
%md
%md### MEBI parameters
%md
%md#### materials
%md The material of the solid considered is the Saint-Venant-Kirchhoff with Lamé parameters computed as

%md since only one material is considered, the structs defined for the materials contain only one entr
materials.hyperElasModel = {'linearElastic'} ;
materials.hyperElasParams = { [ E nu ] } ;
%md
%md#### elements
%md In this model two kinds of elements are used: tetrahedrons for the solid and triangles for introducing the external loads. Since two kinds of elements are used, the structs have length 2:
elements.elemType = { 'edge', 'triangle' } ;
%md since triangle and tetrahedron elements dont have specific parameters the struct entries contain empty vectors
elements.elemTypeParams = { [];[] } ;
elements.elemTypeGeometry = { thickness ; thickness } ;
%md
%md#### boundaryConds
%md in this case four BCs are considered, one corresponding to a load and three to displacements.
boundaryConds.loadsCoordSys = {'global'        ; [ ] ; [ ]  } ;
boundaryConds.loadsTimeFact = { @(t) t         ; [ ] ; [ ]  } ;
boundaryConds.loadsBaseVals = { [p 0 0 0 0 0 ] ; [ ] ; [ ]  } ;
boundaryConds.imposDispDofs = { []             ; [1] ; [3]  } ;
boundaryConds.imposDispVals = { []             ; [0] ; [0]  } ;
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
%md The connectivity matrix is given by the following matrix
mesh.nodesCoords = [ 0    0    0 ; ...
                     Lx   0    0 ; ...
                     0   Ly    0 ; ...
                     Lx  Ly    0 ] ;
%md and the connectivity cell is defined as follows with the MEBI integer parameters for each element. All the eight triangle elements are considered with no material (since they are used only to include load) and the following six elements are solid SVK material tetrahedrons.
%md
mesh.conecCell = {[ 0 1 3 0    1 2   ]; ... % constrained edge
                  [ 0 1 2 0    1 3   ]; ... % constrained edge
                  [ 0 1 1 0    2 4   ]; ... % loaded edge
                  [ 1 2 0 0    2 3 1 ]; ... % triangle
                  [ 1 2 0 0    2 4 3 ]  ... % triangle
                } ;
%md
%md### Analysis parameters
%md
analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.stopTolIts    = 30      ;
analysisSettings.stopTolDeltau = 1.0e-12 ;
analysisSettings.stopTolForces = 1.0e-12 ;
analysisSettings.finalTime      = 2       ;
analysisSettings.deltaT        = 1      ;
%md
%md
%md### Output parameters
otherParams.plotParamsVector = [ 3 ] ;
otherParams.problemName = 'uniaxialExtension_Manual' ;
%~ printflag = 2 ;
%md
[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md
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
