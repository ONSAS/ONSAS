%  A Plastic Frame Analysis / Softening Hinges
close all, if ~strcmp( getenv('TESTS_RUN'), 'yes'), clear, end
addpath( genpath( [ pwd '/../../src' ] ) ) ; % add ONSAS directory to path

% assumed XY plane

% /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\
% material
EI  = 77650 ;       % KN.m^2
kh1 = 29400 ;       % KN.m^2
kh2 = 2730 ;
Ks  = -kh1 ;        % KN.m

nu = 0.3 ;          % Poisson's ratio

% geometry
L1 = 3 ;              % m
L2 = 6 ;
L3 = 3 ;
ty = 0.3 ;              % width cross section
tz = 0.3 ;              % height cross section

Inertia = tz*ty^3/12 ;  % m^4

E = EI/Inertia ;        % KN/m^2 [KPa]

A  = ty*tz ;            % m^2
Mc = 37.9 ;             % KN.m
My = 268 ;
Mu = 374 ;

% /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\

materials = struct();
materials(1).modelName = 'plastic-2Dframe'  ;
materials.modelParams = [ E Mc My Mu kh1 kh2 Ks nu ] ;
% elements
elements = struct();
elements(1).elemType  = 'node'  ;
elements(2).elemType  = 'frame' ;
elements(2).elemCrossSecParams = {'generic' ; [A 1 Inertia Inertia] } ;
% Boundary Conditions
% Supports
boundaryConds = struct();
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;

boundaryConds(3).imposDispDofs = [ 2 4 5 ] ;
boundaryConds(3).imposDispVals = [ 0 0 0 ] ;

% Loads
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsBaseVals = [ 1 0 0 0 0 0 ] ;
boundaryConds(2).loadsTimeFact = @(t) t ;

boundaryConds(2).imposDispDofs = [ 2 4 5 ] ;
boundaryConds(2).imposDispVals = [ 0 0 0 ] ;

% Mesh
% Mesh nodes
mesh = struct();
mesh.nodesCoords = [ 0      0       0 ; ...
                     0      L1/2    0 ; ...
                     0      L1      0 ; ...
                     0      L2      0 ; ...
					 L3     L2      0 ; ...
                     L3     L1      0 ; ...
                     L3     0       0 ] ;
% Conec Cell
mesh.conecCell = { } ;
% nodes
mesh.conecCell{1, 1 } = [ 0 1 1   1 ] ;
mesh.conecCell{7, 1 } = [ 0 1 1   7 ] ;

mesh.conecCell{2, 1 } = [ 0 1 3   2 ] ;
mesh.conecCell{3, 1 } = [ 0 1 3   3 ] ;
mesh.conecCell{4, 1 } = [ 0 1 2   4 ] ;
mesh.conecCell{5, 1 } = [ 0 1 3   5 ] ;
mesh.conecCell{6, 1 } = [ 0 1 3   6 ] ;

% and frame elements
mesh.conecCell{8, 1 }  = [ 1 2 0   1 2 ] ;
mesh.conecCell{9, 1 }  = [ 1 2 0   2 3 ] ;
mesh.conecCell{10, 1 } = [ 1 2 0   3 4 ] ;
mesh.conecCell{11, 1 } = [ 1 2 0   3 6 ] ;
mesh.conecCell{12, 1 } = [ 1 2 0   4 5 ] ;
mesh.conecCell{13, 1 } = [ 1 2 0   5 6 ] ;
mesh.conecCell{14, 1 } = [ 1 2 0   6 7 ] ;


% InitialConditions
% empty struct
initialConds = struct() ;

% Analysis settings
analysisSettings                    = {} ;
analysisSettings.methodName         = 'arcLength' ;
analysisSettings.deltaT             = 1 ;
analysisSettings.incremArcLen       = [1e-3*ones(1,160) 1e-4*ones(1,250)] ;
analysisSettings.finalTime          = length(analysisSettings.incremArcLen) ;
analysisSettings.iniDeltaLamb       = 1 ;
analysisSettings.posVariableLoadBC  = 2 ;
analysisSettings.stopTolDeltau      = 1e-10 ;
analysisSettings.stopTolForces      = 1e-10 ;
analysisSettings.stopTolIts         = 200 ;
analysisSettings.ALdominantDOF      = [4*6+1 1] ;

%
otherParams = struct() ;
otherParams.problemName = 'plastic_2dframe' ;
% otherParams.plots_format = 'vtk' ;

[ modelCurrSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

[matUs, loadFactorsMat, modelSolutions ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData ) ;

rotations = matUs((3)*6+6,:) ;
displacements = matUs((3)*6+1,:) ; % node with horizontal load applied
loadfactors = loadFactorsMat(:,2) ;

Hinges = zeros(7,3) ;

moments_hist = zeros(4,length(modelSolutions)) ;
for i =1:length(modelSolutions)
    for jj = 1:7
    
        aux = modelSolutions{i}.localInternalForces(jj) ;
        moments_hist(:,i) = [ aux.Mz; aux.Mz2; aux.Mz3; aux.tM ] ;

        if abs(moments_hist(1,i)) > Mu && any(Hinges(jj,:)) == false
        
            Hinges(jj,1) = true ;

        elseif abs(moments_hist(2,i)) > Mu && any(Hinges(jj,:)) == false
        
            Hinges(jj,2) = true ;

        elseif abs(moments_hist(3,i)) > Mu && any(Hinges(jj,:)) == false
        
            Hinges(jj,3) = true ;

        end
    end
end

disp(Hinges) ;

moments_hist = zeros(4,length(modelSolutions)) ;
for i =1:length(modelSolutions)
    aux = modelSolutions{i}.localInternalForces(4) ;
    moments_hist(:,i) = [ aux.Mz; aux.Mz2; aux.Mz3; aux.tM ] ;
end
Mn1_numericONSAS = moments_hist(1,:) ;
Mn2_numericONSAS = moments_hist(2,:) ;
Mn3_numericONSAS = moments_hist(3,:) ;
tMn_numericONSAS = moments_hist(4,:) ;

% Plots

lw = 2 ; ms = 1 ; plotfontsize = 14 ; step = 10 ;

figure('Name','Frame / Plasticity (load factors)','NumberTitle','off') ;
hold on, grid on

plot(abs(displacements), loadfactors, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#0072BD") ;

labx = xlabel('Displacements (m)') ;
laby = ylabel('\lambdaF') ;

legend('ONSAS \lambdaF[y]', 'location', 'Southeast') ;

set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Frame / Plasticity (load factors)') ;

figure('Name','Frame / Plasticity (Moments)','NumberTitle','off') ;
hold on, grid on

plot(abs(displacements(1:step:length(tMn_numericONSAS))), tMn_numericONSAS(1:step:length(tMn_numericONSAS)), '-x', 'linewidth', lw, 'markersize', ms*6, "Color", "#0072BD") ;
plot(abs(displacements), Mn1_numericONSAS, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#EDB120") ;
plot(abs(displacements), Mn2_numericONSAS, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#77AC30") ;
plot(abs(displacements), Mn3_numericONSAS, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#A2142F") ;

labx = xlabel('Displacements (m)') ;
laby = ylabel('Moments tMn') ;

legend('ONSAS tMn [y]', 'ONSAS Mn1 [y]', 'ONSAS Mn2 [y]', 'ONSAS Mn3 [y]', 'location', 'Southeast') ;

set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Frame / Plasticity (load factors)') ;

print('-f1','../../../Tesis/tex/imagenes/plasticFrameVLoadFactors.png','-dpng') ;
print('-f2','../../../Tesis/tex/imagenes/plasticFrameVMomentstMn.png','-dpng') ;