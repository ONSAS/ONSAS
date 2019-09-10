// Test for nodes2dofs C++ conversion with mkoctfile
// mkoctfile nodes2dofsCpp.cpp

#include <stdio.h>
#include<iostream>

#include <octave/oct.h>

using namespace std;

DEFUN_DLD (nodes2dofsCpp, args, , "Nodes to DOFs conversion in C++.")
{
  
  int nodeNumber  = int( args(0).array_value()(0) ) ;
  int dofsPerNode = int( args(1).array_value()(0) ) ;

  NDArray dofsVec(dofsPerNode) ;
  
  for (int i=0; i < dofsPerNode; i++){
    dofsVec(i) = (double) ( (nodeNumber-1)*dofsPerNode + (i+1) ) ;
    cout << dofsVec(i) << "\n";
  }
  
  return octave_value( dofsVec );
}
