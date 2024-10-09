
close all, clear all
addpath( genpath( [ pwd '/../../src'] ) );

% scalar parameters
E = 210e9 ;  A = 2.5e-3 ;
Kplas = E*.1 ;
sigma_Y_0 = 25e6 ;

Lx = 2 ; Ly = 2 ;
n_cells_x = 8 ;

% 
n_cells_y = round(n_cells_x*Ly/Lx) ;

materials = struct();
materials.modelName  = 'plastic-rotEngStr' ;
materials.modelParams = [ E Kplas sigma_Y_0 ] ;

elements = struct();
elements(1).elemType = 'node' ;
elements(2).elemType = 'truss';
elements(2).elemCrossSecParams = { 'circle' , sqrt(A*4/pi) } ;

boundaryConds = struct();
boundaryConds(1).imposDispDofs = [ 1 3 5 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 ] ;

boundaryConds(2).imposDispDofs = [ 5 ] ;
boundaryConds(2).imposDispVals = [ 0 ] ;

boundaryConds(3).imposDispDofs = [ 5 ] ;
boundaryConds(3).imposDispVals = [ 0 ] ;
boundaryConds(3).loadsCoordSys = 'global'         ;
boundaryConds(3).loadsTimeFact = @(t) 1*t     ;
boundaryConds(3).loadsBaseVals = [ 0 0 -1 0 0 0 ] ;

mesh = struct();
n_nodes_x = n_cells_x +1;
n_nodes_y = n_cells_y +1;

xs = linspace(0,Lx,n_nodes_x)' ;
ys = linspace(0,Ly,n_nodes_y)' ;

aux = zeros(n_nodes_x*n_nodes_y,3) ;
for j=1:n_nodes_y
    aux( (n_nodes_x*(j-1)+1):n_nodes_x*j,:) = [ xs ones(n_nodes_x,1)*[ys(j) 0] ]; 
end

mesh.nodesCoords = aux;

aux={} ;
for i=1:n_nodes_x
    aux{ end+1 } = [ 0 1 1    i   ] ;
end

for i=(n_nodes_x+1):(n_nodes_x*n_nodes_y)
    aux{ end+1 } = [ 0 1 2    i   ] ;
end

assert(mod(n_cells_x,2)==0)

loadednode = n_nodes_x*n_nodes_y - n_cells_x/2 ;

aux{ loadednode } = [ 0 1 3    loadednode ] ;
aux{ loadednode-1 } = [ 0 1 3    loadednode-1 ] ;
aux{ loadednode+1 } = [ 0 1 3    loadednode+1 ] ;

for i=1:n_cells_x
    for j=1:n_cells_y
        n1 = n_nodes_x*(j-1)+i+0 ;
        n2 = n_nodes_x*(j-1)+i+1 ;
        n3 = n_nodes_x*(j  )+i+0 ;
        n4 = n_nodes_x*(j  )+i+1 ;

        aux{ end+1 } = [ 1 2 0  n1 n2 ] ;
        aux{ end+1 } = [ 1 2 0  n1 n3 ] ;
        aux{ end+1 } = [ 1 2 0  n2 n3 ] ;
        aux{ end+1 } = [ 1 2 0  n1 n4 ] ;
    end
end

j=n_cells_y;
for i=1:n_cells_x
    aux{ end+1 } = [ 1 2 0  n_nodes_x*(j)+i  n_nodes_x*(j)+i+1 ] ;
end

for j=1:n_cells_y
    aux{ end+1 } = [ 1 2 0  n_nodes_x*(j-1)+(n_nodes_x) n_nodes_x*(j)+(n_nodes_x) ] ;
end

mesh.conecCell = aux; 

initialConds                = struct() ;

analysisSettings = struct();
analysisSettings.methodName    = 'arcLength' ;
analysisSettings.deltaT        =   1  ;
analysisSettings.finalTime     =   300    ;

analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   15   ;
analysisSettings.incremArcLen  =   [1e-4*ones(1,20) 1e-3*ones(1,280) ]   ;

analysisSettings.posVariableLoadBC = 3 ;
analysisSettings.ALdominantDOF = [ (loadednode-1)*6+3 -1];

otherParams = struct();
otherParams.problemName = 'plasticLattice';
otherParams.plots_format = 'vtk' ;
# otherParams.plots_deltaTs_separation = 2 ;

[ modelCurrSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;


%
%mdAfter that the structs are used to perform the numerical time analysis
[matUs, loadFactorsMat, modelSolutions ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData ) ;

deltas = -matUs(6+5,:)' ;
eles = sqrt( x2^2 + (z2-deltas).^2 ) ;

valsLin = -2* ( E * (eles - L) ./ L ) * A .* (z2-deltas ) ./ eles  ;

Etan = E*Kplas / ( E + Kplas) ;
sigmas_hard = sigma_Y_0 + Etan * ( abs( (eles - L)) ./ L - sigma_Y_0/E ) ;

valsFlu = 2*( sigmas_hard * A ) .* ( ( z2 - deltas ) ./ eles )  ;
valsP = min( valsLin, valsFlu) ;


% softening
Kplas = -E*.05 ;
materials.modelParams = [ E Kplas sigma_Y_0 ] ;

analysisSettings.methodName    = 'arcLength' ;
analysisSettings.iniDeltaLamb = boundaryConds(2).loadsTimeFact(.2)/100 ;
analysisSettings.incremArcLen = [ 2e-4 4e-5*ones(1,100)];

[ modelCurrSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
%
%mdAfter that the structs are used to perform the numerical time analysis
[matUsB, loadFactorsMatB, modelSolutions ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData ) ;

deltasB = -matUsB(6+5,:)' ;

eles = sqrt( x2^2 + (z2-deltasB).^2 ) ;

valsLin = -2*E*A * (z2-deltasB ) ./ eles .* (eles - L) ./ L ;

Etan = E*Kplas / ( E + Kplas) ;
sigmas_hard = sigma_Y_0 + Etan * ( abs( (eles - L)) ./ L - sigma_Y_0/E ) ;

valsFlu = 2*( sigmas_hard * A ) .* ( ( z2 - deltasB ) ./ eles )  ;
valsPB = min( valsLin, valsFlu) ;

verifBoolean = ( norm( valsP  - loadFactorsMat(:,2) ) < 1e-4*norm(valsP ) ) && ...
               ( norm( valsPB - loadFactorsMatB(:,2)) < 1e-4*norm(valsPB) ) ;

figure
plot( deltas , valsP, 'b-x' )
grid on, hold on
plot( deltas , loadFactorsMat(:,2), 'r-o' )
plot( deltasB , loadFactorsMatB(:,2), 'g-o' )
plot( deltasB , valsPB, 'k-*' )
legend('analytic-hard','numeric-hard', 'analytic-soft','numeric-soft')
