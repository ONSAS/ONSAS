%md# Elastplastic perfect cylinder plane strain  example  

%md### Numerical solution
%md Before defining the structs, the workspace is cleaned, the ONSAS directory is added to the path and scalar geometry and material parameters are defined:
close all, if ~strcmp( getenv('TESTS_RUN'), 'yes'), clear all, end
% add path
addpath( genpath( [ pwd filesep '..' filesep '..' filesep 'src' ] ) ) ;
% scalar parameters
global a
global b

b = 200 ; 
a = 100 ;
L = .75 ;

p = 0.01 ;
E = 210 ; nu = 0.3 ; H = 0 ; sigmaY0 = 0.24 ;
%md
%md
%md### MEB parameters
%md
%md#### materials
%md The constitutive behavior of the material considered is isotropic hardening.
%md Since only one material is considered, the structs defined for the materials contain only one entry:
materials.hyperElasModel  = 'isotropicHardening' ;
materials.hyperElasParams =  [ E nu H sigmaY0 ] ;
%md
%md#### elements
%md 
%md In this plane model, three kinds of elements are used: `triangle` for the solid, `edges` to add pressure loads and `nodes` to set additional boundary conditions for the numerical resolution. Since three kinds of elements are used, the struct has length 3: 
elements(1).elemType           = 'node'    ;
elements(2).elemType           = 'edge'    ;
elements(2).elemCrossSecParams = L         ;
elements(3).elemType           = 'triangle';
elements(3).elemTypeParams     = 2         ;
elements(3).elemCrossSecParams = L         ;
%md where `elemCrossSecParams` field sets the thickness of the edge and `elemTypeParams` sets the plane strain triangle element.  
%md
%md#### boundaryConds
%md Three BCs are considered, one corresponding to a load and two for displacements.
%md The first two BCs constrain displacements in $x$ and $y$ global directions respectively:
boundaryConds(1).imposDispDofs = [1] ;
boundaryConds(1).imposDispVals = [0] ;
boundaryConds(2).imposDispDofs = [3] ;
boundaryConds(2).imposDispVals = [0] ;
%md then the third BC corresponds to the pressure. It is introduced in `local` coordinates so the first entry is along the edge (tangent) and the second towards the normal vector obtained by rotating the tangent vector 90 degrees in global axis $z$:
boundaryConds(3).loadsCoordSys = 'local' ;
boundaryConds(3).loadsTimeFact = @(t) t*(t<=19) + (t-(t-19)*2)*(t>19)  ;
boundaryConds(3).loadsBaseVals = [ 0 p ]  ;
%md
%md#### initialConds
%md Any non-homogeneous initial conditions are considered, thereafter an empty struct is set:
initialConds = struct();
%md
%md### Mesh
%md The mesh can be read from the msh file. However, if any changes to the mesh are desired, the .geo file can be edited and the msh file can be re-generated using GMSH.
%md
%md```@raw html
%md<img src="../../assets/linearCylinderPlaneStrain/meshCylinderPlaneStrain.png" alt="mesh plot" width="500"/>
%md```
%md The element properties are set using labels into GMSH follwing the MEBI nomenclature. First `triangle` elements have linear elastic material so entry $1$ of the _materialṣ_ struct is assigned. Then for both `node` and `edge` elements any material is set. 
%md Next displacement boundary conditions are assigned to the element, since the problem is modeled into $x-y$ plane, a constrain to avoid rotation along $z$ is necessary. This is done fixing $y$ and $x$ displacements (using `boundaryConds(1)` and `boundaryConds(2)` as labels) on points 2 3 4 5.
%md Finally the internal pressure is applied on the `edge` elements linked with curves from one to four (Circles 1-4 in Figure). In accordance with the orientation of the curve set in GMSH, the normal vector obtained in local coordinates is $e_r$ so the internal pressure is assigned using `boundaryConds(3)`. Once the mesh is created is read using:
base_msh='';
if strcmp( getenv('TESTS_RUN'),'yes') && isfolder('examples'),
  base_msh=['.' filesep 'examples' filesep 'EPPCylinderPlaneStrain' filesep];
end
mesh = struct();
[ mesh.nodesCoords, mesh.conecCell ] = meshFileReader( [ base_msh 'ringEPP_2.msh'] ) ;

Conec = myCell2Mat( mesh.conecCell ) ;
elems = size(Conec,1) 
%md
%md### Analysis parameters
%md
%md The Newton-Raphson method is employed to solve 19 load steps. The ratio between `finalTime` and `deltaT` sets the number of load steps used to evaluate `boundaryConds(3).loadsTimeFact` function:  
analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.stopTolIts    = 30      ;
analysisSettings.stopTolDeltau = 1.0e-8 ;
analysisSettings.stopTolForces = 1.0e-6 ;
analysisSettings.finalTime     = 19       ;
analysisSettings.deltaT        = 1      ;
%md
%md### Output parameters
%md
otherParams.problemName = 'EPPPlaneStrain' ;
otherParams.plots_format = 'vtk' ;
%md The ONSAS software is executed for the parameters defined above and the displacement solution of each load(time) step is saved in `matUs`matrix:
%md
[matUs, loadFactorsMat, ~, cellStress ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%md
%md## Verification
%mdThe numerical and analytic solutions are compared for the external surface (since all the elements on the same surface have the same analytic solution):

global Y

Y = 2*sigmaY0 / sqrt(3) ;

p0 = Y/2 * (1-a^2/b^2) ; % Yielding pressure

% Plastic front value
for i = 1:length(cvals)
	p = pressure_vals(i) ;
	if i == 1  
		val = fsolve(@(c)c_val(c,p,Y,a,b), a) ;
	else
		val = fsolve(@(c)c_val(c,p,Y,a,b), cvals(i-1)) ;
	end
	cvals(i) = val ;
end

pressure_vals = loadFactorsMat(:,3)*p ;

cvals = zeros(length(pressure_vals),1) ;
ubAna = zeros(length(pressure_vals),1) ;

% Analytic radial displacement at outer surface
for i = 1:length(cvals)
	p = pressure_vals(i) ;
	if p < p0
		ubAna(i) = 2*p*b / ( E*( b^2/a^2-1 ) ) * (1-nu^2) ;
	else
		c = cvals(i) ;
		ubAna(i) = Y*c^2/(E*b) * (1-nu^2) ;
	end	
end

% Plot parameters
lw = 2.0 ; ms = 11 ; plotFontSize = 10 ;
fig = figure;
hold on, grid on

node = 5 ;
dofX = node * 6 - 5 ;
ubNum = matUs(dofX, :) ; 

plot(ubNum, pressure_vals, 'b-o', 'linewidth', lw,'markersize', ms)
plot(ubAna, pressure_vals, 'g-x', 'linewidth', lw,'markersize', ms)

legend ({'FEM', 'Analytic',}, 'location', 'east');
labx = xlabel('u_b'); laby = ylabel('p') ;
tit = title('p-u_b');
set(labx, 'fontsize', plotFontSize*.8);
set(laby, 'fontsize', plotFontSize*.8);
set(tit, 'fontsize', plotFontSize);


% Plastic front
fig2 = figure;
hold on, grid on

cvals(1:11) = a ;
plot(cvals,pressure_vals, 'r-s', 'linewidth', lw,'markersize', ms)

legend ({ 'Plastic front'}, 'location', 'east');


% Plastic front - cylinder
fig3 = figure;
hold on, grid on

theta = [0:0.01:2*pi] ;
aux = zeros(length(cvals)) ;
plot(a*cos(theta), a*sin(theta) , 'b', 'linewidth', lw,'markersize', ms)
plot(cvals, aux, 'r-s', 'linewidth', lw,'markersize', ms)
plot(b*cos(theta), b*sin(theta) , 'b', 'linewidth', lw,'markersize', ms)

for i = 12:length( cvals )
	plot(cvals(i)*cos(theta), cvals(i)*sin(theta), 'g--', 'linewidth', lw,'markersize', ms)
end

legend ({ 'Cylinder' , 'Plastic front'}, 'location', 'northeast');
axis("equal")


% Check solution
analyticCheckTolerance = 1e-2 ;
verifBoolean = ( ( ubNum(end) - ubAna(end) ) < analyticCheckTolerance ) ;

%~ % 2nd run
%~ [ mesh.nodesCoords, mesh.conecCell ] = meshFileReader( 'ringEPP_2.msh' ) ;

%~ Conec = myCell2Mat( mesh.conecCell ) ;
%~ elems = size(Conec,1) 

%~ [matUs2, loadFactorsMat, ~, cellStress ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

%~ ubNum_2 = matUs2(dofX, :) ; 

%~ fig2 = figure, hold on, grid on
%~ plot(ubNum_2, pressure_vals, 'b-o', 'linewidth', lw,'markersize', ms)
%~ plot(ubAna, pressure_vals, 'g-x', 'linewidth', lw,'markersize', ms)

%~ legend ({'FEM', 'Analytic'}, 'location', 'east');

%~ labx = xlabel('u_b'); laby = ylabel('p') ;

%~ tit = title('p-u_b');
%~ set(labx, 'fontsize', plotFontSize*.8);
%~ set(laby, 'fontsize', plotFontSize*.8);
%~ set(tit, 'fontsize', plotFontSize);


%~ ubNumVec = [ubNum(end), ubNum_2(end)]

%~ errVec = (ubNumVec-ubAna(end)) / ubAna(end) *100 



%%
%%

% Implicit function
function val = c_val(c,p,Y,a,b)
	val = p/Y-( log(c/a)+1/2*(1 - c^2/b^2) ) ;
end
