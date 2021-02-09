
function [ fun , vol ] = DerivFun( tetcoordmat , varargin )

  if nargin == 1
    derivOrder = 0 ;
  else
    aux4Matlab =  varargin(1)   ;
    derivOrder =  aux4Matlab{1} ;
  end

  A        = zeros(4,4)   ;
  A(:,1)   = 1.0          ;
  A(:,2:4) = tetcoordmat' ;

  invA = inv(A) ;

  vol = det(A) / 6.0 ;

  fun = invA(2:4,:) ;
