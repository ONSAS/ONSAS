%md# Beam truss joint example
close all, clear all
problemName = 'beamTrussJoint' ;
addpath( genpath( [ pwd '/../../src'] ) );

%mdThe goal of this example is to provide a minimal validation of the integration between truss and frame elements in the same model. The structure considered is formed by two elements (one truss (t) and one beam (b)) and small displacements are considered.

%md Truss geometrical and material properties are:
Et = 1e9 ; nu = 3; dt = .05; At = pi*dt^2/4 ;  lt = 1 ; nut = 0.3 ;  
%md and frame geometrical and material properties are:
Eb = Et/3 ;db = 5*dt ; Ab = pi*db^2/4 ;  lb = .5 ; nub = 0.3 ; Ib = pi*db^4/64 ;

%md##Numerical solution
%md### MEBI parameters
%md### materials
%mdSince the example contains two different type of materials the fields of the `materials` struct will have two entries. Although the structure develops small displacements a Rotated Engineering strain material constitutive behavior is considered.
materials(1).hyperElasModel  = '1DrotEngStrain' ;
materials(1).hyperElasParams = [ Et nu ] ;
%
materials(2).hyperElasModel  = '1DrotEngStrain' ;
materials(2).hyperElasParams = [ Eb nu ] ;

%md### elements
%md
%mdThree different types of elements are considered: node, frame and truss, defined as follows:
elements(1).elemType = 'node'  ;
elements(2).elemType = 'truss' ;
elements(3).elemType = 'frame' ;
%mdIn order to add the struct of geometry the assign to the node is an empty input (because it has not geometrical properties), and truss and beam elements will be set as:
elements(2).elemTypeGeometry = [3 dt] ;
elements(3).elemTypeGeometry = [3 db] ;

%mdTruss number of elements
%
numElemT  = 1            ;
numNodesT = numElemT + 1 ;
%md and beam:
numElemB  = 10           ;
numNodesB = numElemB + 1 ;

%md
%md### boundaryConds
%md
%mdThe fixed frame BC:
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
%mdloaded BC:
boundaryConds(2).imposDispDofs = [ 2 3 6 ] ;
boundaryConds(2).imposDispVals = [ 0 0 0 ] ;
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsTimeFact = @(t) 1e4 * t ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 1 0 ] ;
%mdbase truss suppor BC:
boundaryConds(3).imposDispDofs = [ 1 2 3 4 5 ] ;
boundaryConds(3).imposDispVals = [ 0 0 0 0 0 ] ;

%md#### initial Conditions
%md homogeneous initial conditions are considered, then an empty struct is set:
initialConds                = struct() ;

%md
%md### mesh parameters
%md
%mdThe coordinates of the nodes of the mesh are given by the matrix:
mesh.nodesCoords = [(0:(numElemB))'*lb/numElemB zeros(numElemB+1,1) zeros(numElemB+1,1) 	       ;	 
		            lb*ones(numElemT,1) 	    zeros(numElemT,1)   -(1:(numElemT))'*lt/numElemT ] ;
%md where the columns 1,2 and 3 correspond to $x$, $y$ and $z$ coordinates, respectively, and the row $i$-th corresponds to the coordinates of node $i$.

%md The conectivity struct using MEBI nomenclature is defined by using the following auxiliar Element and Nodes matrix:
%mdAuxliar conecNodes:
auxConecNodes = [0 1 1 0  1                     ;
                 0 1 2 0  numNodesB             ;
                 0 1 3 0  numNodesB + numElemT  ] ;  

%mdand the auxiliar conecNodes is:
auxConecElem  = [   %MEBI frame elements
                    [ (ones(numElemB,1)*2 )  (ones(numElemB,1)*3)   (zeros(numElemB,1))  (zeros(numElemB,1)) ...
                    %ElemNodes..
                    (1:(numElemB))'         (2:numElemB + 1)' ] ;
                    %MEBI truss elements
                    [ (ones(numElemT,1)*1 )  (ones(numElemT,1)*2)   (zeros(numElemT,1))  (zeros(numElemT,1)) ... 
                     %ElemNodes...
                    (numElemB + 1: numElemB + numElemT)'    (numElemB + 2:numElemB + numElemT + 1)'   ] ... 
                ] ;  


%Build conectivity cell:
mesh.conecCell = cell( size(auxConecElem,1) + size(auxConecNodes,1) ) ;
for i = 1:size(auxConecNodes,1)
    mesh.conecCell{i,1}= auxConecNodes(i,:) ;
end
for i =  1:numElemB + numElemT
    mesh.conecCell{size(auxConecNodes,1) + i,1} = auxConecElem(i,:);
end                                          
%md### analysisSettings
%md The method used in the analysis is the Newton-Raphson, then the field `methodName` must be introduced as:
analysisSettings.methodName    = 'newtonRaphson' ;
%md and the following parameters correspond to the iterative numerical analysis settings
analysisSettings.deltaT        =   0.1  ;
analysisSettings.finalTime     =   1    ;
analysisSettings.stopTolDeltau =   1e-6 ;
analysisSettings.stopTolForces =   1e-6 ;
analysisSettings.stopTolIts    =   10   ;

%md### otherParams
otherParams.problemName = problemName   ;
otherParams.plotsFormat = ''            ;

[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;


%md the analytic solution is computed at the joint where uz is given by:
analyticFunc            = @(w)( Et*At/lt + 3*Eb*Ib/lb^3 )*w ;
beamTruss_stiffRatio    = Et*At/lt/(3*Eb*Ib/lb^3)           ;
%mdand the numerical: 
%md 
controlDof  = (numNodesB)*6 - 1     ;
dispZnum    = matUs(controlDof,:)   ;
#mdverification boolean is computed as follows:
%md
difLoadEngRot   = loadFactorsMat(:,2)' - analyticFunc(dispZnum) ;
verifBoolean    = ( ( norm( difLoadEngRot    ) / norm( loadFactorsMat(:,2)  ) ) <  1e-4 ) ;
%md
%md### Plots
%mdOutput diplacments and load factor function are plotted:
figure
hold on, grid on
lw = 2.0 ; lw2 = 1.0 ; ms = 11 ; plotfontsize = 18 ;
plot( dispZnum, analyticFunc(dispZnum), 'b-x' , 'linewidth', lw,'markersize',ms)    ;
plot( dispZnum, loadFactorsMat(:,2),    'r-s' , 'linewidth', lw,'markersize',ms )   ;
labx = xlabel('Displacement u_z(t)');   laby = ylabel('\lambda(t)')                 ;
set(gca, 'linewidth', lw2, 'fontsize', plotfontsize )                               ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize)            ;

