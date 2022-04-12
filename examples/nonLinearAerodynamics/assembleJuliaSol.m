%This function is used to assemble julia solution of the non linear cantilever Aero differential equation. The solver used is DiffentialEquations.jl withh DiffEq.jl file. 

function [dSolJulia xdefJulia] = assembleJuliaSol(elements,mesh)
  % Build conec matrix
  conecElemMatrix =  [] ; 
  for indexElem = 2:length(mesh.conecCell)
    conecElemMatrix = [ conecElemMatrix; mesh.conecCell{indexElem,:} ];
  end
  % number of elements is:
  numElements = size(conecElemMatrix,2) - 1 ; 
  % span length
  l = mesh.nodesCoords(end,1) ;
  % Read DiffEq.jl output
  solCell     = dlmread ('output/solJDiffEq.csv', ',', 1, 0) ;
  xJulia      = solCell(:,1) ;
  thetaZJulia = solCell(:,4) ;
  uYJulia     = solCell(:,5) ;

  % Fill vector solution 
  dSolJulia          = zeros( 6*size( xJulia, 1 ), 1 ) ;
  dSolJulia(3:6:end) = -1 * uYJulia'                   ;
  dSolJulia(6:6:end) = -1 * thetaZJulia'               ;

  % Create and fill solution with the same size of the vector locating their respective coordinates
  dSol = zeros( (numElements + 1) * 6 ,1 );
  numElemJulia = size( dSolJulia(1:6:end) ) - 1;

  % Find the element displacements solution
  for elem = 1:size(conecElemMatrix,1)
    % locate dofs and element coordinates
    % element coords
    x = mesh.nodesCoords(conecElemMatrix(elem,5:6)',1) ;
    y = mesh.nodesCoords(conecElemMatrix(elem,5:6)',2) ;	
    z = mesh.nodesCoords(conecElemMatrix(elem,5:6)',3) ;
    elemCoords = [x(1);y(1);z(1);x(2);y(2);z(2)]; 

    % element dofs
    dofElemVec        = 6*conecElemMatrix( elem,5:6 )-5                             ;
    dofElem           = [ dofElemVec(1) : dofElemVec(2) + 5 ]'                      ; 
    xelemCoords       = elemCoords(1:3:end)                                         ;
    solJuliaIndex     = round (xelemCoords * numElemJulia / l)                      ;
    dofNode1ElemJulia = ( solJuliaIndex(1)*6 + 1 : 1:( solJuliaIndex(1)+1 )*6)'     ;
    dofNode2ElemJulia = ( solJuliaIndex(2)*6 + 1 : 1 :( solJuliaIndex(2) + 1 )*6 )' ;
    dofElemJulia      = [ dofNode1ElemJulia; dofNode2ElemJulia ]                    ;

    % Read displacments from julia solution
    UeSol             = dSolJulia( dofElemJulia ) ;
    dSol(dofElem)     = UeSol ;
  end
  % xdefJulia
  xdefJulia = xJulia +  dSolJulia(1:6:end) ;
end
