
function vtkMainWriter
[ vtkNodes, vtkDispMat, vtkNormalForces, vtkStress ] = vtkGeometry( ...
   modelProperties.Conec, crossSecsParamsMat, modelNextSol.U, normalForces, modelProperties.Nodes, elementsParamsMat, modelNextSol.Stress  ) ;

% data
[ cellPointData, cellCellData, filename ] = vtkData( outputDir, problemName, indplot, vtkNormalForces, vtkStress, vtkDispMat ) ;

% writes file
vtkWriter( filename, vtkNodes, vtkConec , cellPointData, cellCellData ) ;
