
function plotTemp(NodesCoord, Conec, timeIndex, temps)

vtkWriter( sprintf('tempPlot_%04i.vtk',timeIndex), NodesCoord, [ 4*ones(size(Conec,1),1) Conec], {'SCALARS','temp',temps }, {} )
