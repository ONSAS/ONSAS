function [numNodes, dofsStep] = elementTypeInfo ( elemType )

    switch elemType
    case 1
      numNodes = 2 ;
      dofsStep = 2 ;
    case 2
      numNodes = 2 ;
      dofsStep = 1 ;
    case 3
      numNodes = 4 ;
      dofsStep = 2 ;
    end
