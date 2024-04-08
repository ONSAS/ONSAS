%md# Rectangular Plate test
%md
close all, if ~strcmp( getenv('TESTS_RUN'), 'yes'), clear all, end
addpath( genpath( [ pwd '/../../src'] ) );
%md Scalar parameters
E = 1e3 ;
nu = 0.3 ;
tz = 1 ;

materials = struct() ;
materials.modelName  = 'elastic-linear' ;
materials.modelParams =  [ E nu ]       ;

elements = struct() ;
elements(1).elemType           = 'edge'    ;
elements(2).elemCrossSecParams = tz         ;
elements(3).elemType           = 'triangle-plate';
elements(2).elemCrossSecParams = {'thickness', tz } ;

boundaryConds                  = struct() ;

boundaryConds(1).imposDispDofs =  [ 1 3 5 ]                  ; % fixed nodes: 1 2 4
boundaryConds(1).imposDispVals =  [ 0 0 0 ] ;

boundaryConds(2).loadsCoordSys = 'global'                   ;
boundaryConds(2).loadsTimeFact = @(t) t  ;
boundaryConds(2).loadsBaseVals = [ 0 0 0 0 -1 0 ]            ; % forces in node 1

base_msh='';
if strcmp( getenv('TESTS_RUN'),'yes') && isfolder('examples'),
  base_msh=['.' filesep 'examples' filesep 'rectangularPlate' filesep];
end
mesh = struct();
[ mesh.nodesCoords, mesh.conecCell ] = meshFileReader( [ base_msh 'rectangularPlate.msh'] ) ;
