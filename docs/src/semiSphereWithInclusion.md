# Semi-sphere with inclusion
In this example a semi-spherical solid is considered with a spherical inclusion inside. The constitutive behavior of the solids is assumed given by a Neo-Hookean law, with the inclusion being stiffer than the matrix solid.

```@raw html
<img src="https://raw.githubusercontent.com/ONSAS/ONSAS.m/master/docs/src/assets/semiSphere.png" alt="plot" width="500"/>
```

```
clear all, close all
addpath( genpath( [ pwd '/../../src'] ) ) ;

% scalar parameters
E = 500 ; nu = 0.45 ;
% tension applied

lambda = 2*E*nu/((1+nu)*(1-2*nu)) ; mu = E/(2*(1+nu)) ;
materials(1).hyperElasModel = 'NHC' ;
materials(1).hyperElasParams = [ lambda mu ] ;

lambda = E*nu/((1+nu)*(1-2*nu)) ; mu = E/(2*(1+nu)) ;
materials(2).hyperElasModel = 'NHC' ;
materials(2).hyperElasParams = [ lambda mu ] ;

elements(1).elemType = 'triangle' ;
elements(2).elemType = 'tetrahedron' ;



boundaryConds(1).loadsCoordSys = 'global';
boundaryConds(1).loadsTimeFact = @(t) t*200 ;
boundaryConds(1).loadsBaseVals = [ 0 0 -1 0 0 0 ] ;
```
 the other BCs have imposed displacements
```
boundaryConds(2).imposDispDofs = [1 3 5] ;
boundaryConds(2).imposDispVals = [0 0 0 ]  ;

initialConds = struct();

[ mesh.nodesCoords, mesh.conecCell ] = meshFileReader( 'geometry_semiSphereWithInclusion.msh' ) ;

otherParams.problemName = 'semiSphereWithInclusion' ;
otherParams.plotsFormat = 'vtk' ;

analysisSettings.methodName    = 'newtonRaphson' ;
analysisSettings.stopTolIts    = 30     ;
analysisSettings.stopTolDeltau = 1.0e-10 ;
analysisSettings.stopTolForces = 1.0e-10 ;
analysisSettings.finalTime      = 1      ;
analysisSettings.deltaT        = .1   ;

[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

maxy = max( mesh.nodesCoords(:,2) )

inds = find( mesh.nodesCoords(:,2)==maxy )

controlNode = inds(1);

disps = -matUs(6*(controlNode-1)+3,:) ;

figure
plot( disps, loadFactorsMat, 'b-x' )
