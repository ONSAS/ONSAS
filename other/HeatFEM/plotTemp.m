
function plotTemp(NodesCoord, Conec, timeIndex, temps, geometryType )

if geometryType == 1
  % plot settings
  MS = 10 ; LW = 1.5 ;
  plot( NodesCoord, temps,'b-o', 'markersize', MS,'linewidth',LW )
  
elseif geometryType == 2
  vtkWriter( sprintf('tempPlot_%04i.vtk',timeIndex), NodesCoord, [ 4*ones(size(Conec,1),1) Conec], {'SCALARS','temp',temps }, {} )
end
