
clear all, close all
addpath( genpath( [ pwd '/../../src'] ) ) ;

% scalar parameters
E = 1e5 ; nu = 0.45 ;
% tension applied
p = 3 ;

lambda = 2*E*nu/((1+nu)*(1-2*nu)) ; mu = E/(2*(1+nu)) ;
materials(1).hyperElasModel = 'SVK' ;
materials(1).hyperElasParams = [ lambda mu ] ;

lambda = E*nu/((1+nu)*(1-2*nu)) ; mu = E/(2*(1+nu)) ;
materials(2).hyperElasModel = 'SVK' ;
materials(2).hyperElasParams = [ lambda mu ] ;

elements(1).elemType = 'triangle' ;
elements(2).elemType = 'tetrahedron' ;



boundaryConds(1).loadsCoordSys = 'global';
boundaryConds(1).loadsTimeFact = @(t) t ;
boundaryConds(1).loadsBaseVals = [ p 0 0 0 0 0 ] ;
%md the other BCs have imposed displacements
boundaryConds(2).imposDispDofs = [1 3 5] ;
boundaryConds(2).imposDispVals =  [0 0 0 ]  ;

initialConds = struct();

[ mesh.nodesCoords, mesh.conecCell ] = meshFileReader( 'geometry_semiSphereWithInclusion.msh' ) ;

otherParams.problemName = 'semiSphereWithInclusion' ;
otherParams.plotsFormat = 'vtk' ;


analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.stopTolIts    = 30     ;
analysisSettings.stopTolDeltau = 1.0e-8 ;
analysisSettings.stopTolForces = 1.0e-8 ;
analysisSettings.finalTime      = 1      ;
analysisSettings.deltaT        = .125   ;

[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
