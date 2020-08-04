function [numNodes, dofsStep] = elementTypeInfo ( elemType )

switch elemType
case 1 % node
  numNodes = 1 ;
  dofsStep = 1 ;

case 2 % truss
  numNodes = 2 ;
  dofsStep = 2 ;

case 3
  numNodes = 2 ;
  dofsStep = 1 ;

case 4
  numNodes = 4 ;
  dofsStep = 2 ;

end
