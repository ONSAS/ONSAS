%  A Plastic Frame Analysis / Softening Hinges
close all, if ~strcmp( getenv('TESTS_RUN'), 'yes'), clear, end
addpath( genpath( [ pwd '/../../src' ] ) ) ; % add ONSAS directory to path

% assumed XY plane

% /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\
% material
EI  = 10000 ;           % KN.m^2
kh1 = 29400 ;           % KN.m^2
kh2 = 2730 ;
Ks  = -kh1 ;            % KN.m

nu = 0.3 ;              % Poisson's ratio

% geometry
L1 = 3 ;                % m
L2 = 3 ;
L3 = 3 ;
ty = 0.1 ;              % width cross section
tz = 0.1 ;              % height cross section
Inertia = tz*ty^3/12 ;  % m^4

E = EI/Inertia ;        % KN/m^2 [KPa]

A  = ty*tz ;            % m^2
Mc = 37.9 ;             % KN.m
My = 50 ;
Mu = 70 ;

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
mesh.nodesCoords = [ 0      0           0 ; ...
                     0      L1/2        0 ; ...
                     0      L1          0 ; ...
                     0      L1 + L2/2   0 ; ...
                     0      L1 + L2     0 ; ...
                     L3/2   L1 + L2     0 ; ...
                     L3/2   L1          0 ; ...
                     L3     L1 + L2     0 ; ...
                     L3     L1 + L2/2   0 ; ...
                     L3     L1          0 ; ...
                     L3     L1/2        0 ; ...
                     L3     0           0 ] ;

% Conec Cell
mesh.conecCell = { } ;
% nodes
mesh.conecCell{ 1, 1 } = [ 0 1 1   1 ] ;
mesh.conecCell{12, 1 } = [ 0 1 1  12 ] ;

mesh.conecCell{2, 1 }  = [ 0 1 3   2 ] ;
mesh.conecCell{3, 1 }  = [ 0 1 3   3 ] ;
mesh.conecCell{4, 1 }  = [ 0 1 3   4 ] ;
mesh.conecCell{5, 1 }  = [ 0 1 2   5 ] ; % node with horizontal load applied
mesh.conecCell{6, 1 }  = [ 0 1 3   6 ] ;
mesh.conecCell{7, 1 }  = [ 0 1 3   7 ] ;
mesh.conecCell{8, 1 }  = [ 0 1 3   8 ] ;
mesh.conecCell{9, 1 }  = [ 0 1 3   9 ] ;
mesh.conecCell{10, 1 } = [ 0 1 3  10 ] ;
mesh.conecCell{11, 1 } = [ 0 1 3  11 ] ;

% and frame elements
mesh.conecCell{13, 1 } = [ 1 2 0   1  2 ] ;
mesh.conecCell{14, 1 } = [ 1 2 0   2  3 ] ;
mesh.conecCell{15, 1 } = [ 1 2 0   3  4 ] ;
mesh.conecCell{16, 1 } = [ 1 2 0   4  5 ] ;
mesh.conecCell{17, 1 } = [ 1 2 0   5  6 ] ;
mesh.conecCell{18, 1 } = [ 1 2 0   6  8 ] ;
mesh.conecCell{19, 1 } = [ 1 2 0   3  7 ] ;
mesh.conecCell{20, 1 } = [ 1 2 0   7 10 ] ;
mesh.conecCell{21, 1 } = [ 1 2 0   8  9 ] ;
mesh.conecCell{22, 1 } = [ 1 2 0   9 10 ] ;
mesh.conecCell{23, 1 } = [ 1 2 0  10 11 ] ;
mesh.conecCell{24, 1 } = [ 1 2 0  11 12 ] ;


% InitialConditions
% empty struct
initialConds = struct() ;

% Analysis settings
analysisSettings                    = {} ;
analysisSettings.methodName         = 'arcLength' ;
analysisSettings.deltaT             = 1 ;
analysisSettings.incremArcLen       = [1e-4*ones(1,600) 1e-5*ones(1,400)] ;
analysisSettings.finalTime          = length(analysisSettings.incremArcLen) ;
analysisSettings.iniDeltaLamb       = 1 ;
analysisSettings.posVariableLoadBC  = 2 ;
analysisSettings.stopTolDeltau      = 1e-14 ;
analysisSettings.stopTolForces      = 1e-8 ;
analysisSettings.stopTolIts         = 50 ;
analysisSettings.ALdominantDOF      = [4*6+1 1] ;

%
otherParams = struct() ;
otherParams.problemName = 'plastic_2dframe' ;
% otherParams.plots_format = 'vtk' ;

[ modelCurrSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

[matUs, loadFactorsMat, modelSolutions ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData ) ;

rotations = matUs((4)*6+6,:) ;
displacements = matUs((4)*6+1,:) ; % node with horizontal load applied
loadfactors = loadFactorsMat(:,2) ;

Hinges = zeros(12,3) ;

moments_hist = zeros(4,length(modelSolutions)) ;
for i =1:length(modelSolutions)
    for jj = 1:12
    
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
    aux = modelSolutions{i}.localInternalForces(7) ;
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

plot(abs(displacements), Mn1_numericONSAS, '-x', 'linewidth', lw, 'markersize', ms, "Color", "#0072BD") ;

labx = xlabel('Displacements (m)') ;
laby = ylabel('Moments (KN.m)') ;

legend('ONSAS tMn [y]', 'location', 'Southeast') ;

set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize ) ;
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;
title('Frame / Plasticity (load factors)') ;

print('-f1','../../../Tesis/tex/imagenes/plasticFrameLoadFactors.pdf','-dpdf') ;
print('-f2','../../../Tesis/tex/imagenes/plasticFrameMomentstMn.pdf','-dpdf') ;