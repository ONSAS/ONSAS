// Test for nodes2dofs C++ conversion with mkoctfile
// mkoctfile nodes2dofsCpp.cpp

#include <octave/oct.h>
#include <array>
#include <iostream>

using namespace std;

DEFUN_DLD (nodes2dofsCpp, args, nargout, "Nodes to DOFs conversion in C++.")
{

  octave_stdout << "Hello World has "
                << args.length () << " input arguments and "
                << nargout << " output arguments.\n";
	
	// Array with nodes									
	NDArray nodesNumber = args(0).array_value ();
	// Length of array
	octave_idx_type n = nodesNumber.numel() ;
	// Dofs per node
	int dofsPerNode = int( args(1).array_value()(0) ) ;
	NDArray dofsVec(dim_vector(n*dofsPerNode,1)) ;
	// Aux vec
	NDArray dofs(dim_vector(dofsPerNode,1)) ;
	for (int i=0; i < dofsPerNode; i++){
		dofs(i,1) = (double) (i+1) ;
	}

  for (int i=0; i < n; i++){
    for (int j=0; j < dofsPerNode; j++){
			dofsVec( (i)*dofsPerNode + j ) = (double) (dofsPerNode*(nodesNumber(i)-1)+j+1);
			std::cout << dofsVec((i)*dofsPerNode + j) << " DOF \n";
		}
  }
	
  return octave_value( dofsVec );
  return octave_value_list ();
}


