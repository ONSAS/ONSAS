%  A Plastic Frame Analysis / Softening Hinges
close all, if ~strcmp( getenv('TESTS_RUN'), 'yes'), clear, end
addpath( genpath( [ pwd '/../../src' ] ) ) ; % add ONSAS directory to path

% assumed XY plane

% /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\
% material
EI  = 77650 ;           % KN.m^2
kh1 = 29400 ;           % KN.m^2
kh2 = 273 ;
Ks  = -18000 ;            % KN.m

nu = 0.3 ;              % Poisson's ratio

% geometry
L1 = 2.5 ;                % m
ty = 0.3 ;              % width cross section
tz = 0.4 ;              % height cross section
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
boundaryConds(2).loadsBaseVals = [ 0 0 -1 0 0 0 ] ;
boundaryConds(2).loadsTimeFact = @(t) t ;

boundaryConds(2).imposDispDofs = [ 2 4 5 ] ;
boundaryConds(2).imposDispVals = [ 0 0 0 ] ;

% Mesh
% Mesh nodes
mesh = struct();
mesh.nodesCoords = [ 0      0           0 ; ...
                     L1/4   0           0 ; ...
                     L1/2   0           0 ; ...
                     L1*3/4 0           0 ; ...
                     L1     0           0 ];

% Conec Cell
mesh.conecCell = { } ;
% nodes
mesh.conecCell{ 1, 1 } = [ 0 1 1   1 ] ;
mesh.conecCell{ 5, 1 } = [ 0 1 1   5 ] ;

mesh.conecCell{2, 1 }  = [ 0 1 3   2 ] ;
mesh.conecCell{3, 1 }  = [ 0 1 2   3 ] ; % node with vertical load applied
mesh.conecCell{4, 1 }  = [ 0 1 3   4 ] ;

% and frame elements
mesh.conecCell{6, 1 } = [ 1 2 0   1  2 ] ;
mesh.conecCell{7, 1 } = [ 1 2 0   2  3 ] ; 
mesh.conecCell{8, 1 } = [ 1 2 0   3  4 ] ;
mesh.conecCell{9, 1 } = [ 1 2 0   4  5 ] ;

% InitialConditions
% empty struct
initialConds = struct() ;

% Analysis settings
analysisSettings                    = {} ;
analysisSettings.methodName         = 'arcLength' ;
analysisSettings.deltaT             = 1 ;
analysisSettings.incremArcLen       = [1e-4*ones(1,1000)] ;
analysisSettings.finalTime          = length(analysisSettings.incremArcLen) ;
analysisSettings.iniDeltaLamb       = 1 ;
analysisSettings.posVariableLoadBC  = 2 ;
analysisSettings.stopTolDeltau      = 1e-14 ;
analysisSettings.stopTolForces      = 1e-8 ;
analysisSettings.stopTolIts         = 30 ;
analysisSettings.ALdominantDOF      = [1*6+3 -1] ;

%
otherParams = struct() ;
otherParams.problemName = 'plastic_2dframe' ;
% otherParams.plots_format = 'vtk' ;

[ modelCurrSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

[ matUs, loadFactorsMat, modelSolutions ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData ) ;

rotations = matUs((2)*6+6,:) ;
displacements = matUs((2)*6+3,:) ; % node with horizontal load applied
loadfactors = loadFactorsMat(:,2) ;

Hinges = zeros(4,3) ;

moments_hist = zeros(4,length(modelSolutions)) ;
for i =1:length(modelSolutions)
    for jj = 1:4
    
        aux = modelSolutions{i}.localInternalForces(jj) ;
        moments_hist(:,i) = [ aux.Mz; aux.Mz2; aux.Mz3; aux.tM ] ;

        if abs(moments_hist(1,i)) >= Mu && Hinges(jj,1) == false
        
            Hinges(jj,1) = true ;

        elseif abs(moments_hist(2,i)) >= Mu && Hinges(jj,2) == false
        
            Hinges(jj,2) = true ;

        elseif abs(moments_hist(3,i)) >= Mu && Hinges(jj,3) == false
        
            Hinges(jj,3) = true ;

        end
    end
end

disp(Hinges) ;

moments_hist = zeros(4,length(modelSolutions)) ;
for i =1:length(modelSolutions)
    aux = modelSolutions{i}.localInternalForces(1) ;
    moments_hist(:,i) = [ aux.Mz; aux.Mz2; aux.Mz3; aux.tM ] ;
end
Mn1_numericONSAS = moments_hist(1,:) ;
Mn2_numericONSAS = moments_hist(2,:) ;
Mn3_numericONSAS = moments_hist(3,:) ;
tMn_numericONSAS = moments_hist(4,:) ;

% Plots

lw = 2 ; ms = 1 ; plotfontsize = 14 ;

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

plot(abs(displacements), abs(tMn_numericONSAS), '-x', 'linewidth', lw, 'markersize', ms, "Color", "#0072BD") ;

labx = xlabel('Displacements (m)') ;
laby = ylabel('Moments (KN.m)') ;

legend('ONSAS tMn [y]', 'location', 'Southeast') ;

set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Frame / Plasticity (load factors)') ;

print('-f1','../../../Tesis/tex/imagenes/plasticFrameSimpleIILoadFactors.pdf','-dpdf') ;
print('-f2','../../../Tesis/tex/imagenes/plasticFrameSimpleIIMomentstMn.pdf','-dpdf') ;