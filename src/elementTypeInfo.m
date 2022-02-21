function [numNodes, dofsStep] = elementTypeInfo ( elemType )

%~ if strcmp( elemType, 'node') %non physical real element
  %~ numNodes = 1 ;
  %~ dofsStep = 1 ;

if strcmp( elemType, 'truss') || strcmp( elemType, 'edge')
  numNodes = 2 ;
  dofsStep = 2 ;

elseif strcmp( elemType, 'frame')
  numNodes = 2 ;
  dofsStep = 1 ;

elseif strcmp( elemType, 'tetrahedron')
  numNodes = 4 ;
  dofsStep = 2 ;

elseif strcmp( elemType, 'triangle')
  numNodes = 3 ;
  dofsStep = 2 ;

end
