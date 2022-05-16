
function [ funder, jacobianmat, vol, tetCoordMat ] = computeFuncDerivVolTetraSolid( elemCoords )

tetCoordMat        = reshape( elemCoords', 3, 4 ) ;

eleCoordMat = tetCoordMat               ;

xi = 0.25 ;  wi = 1/6  ;

% matriz de derivadas de fun forma respecto a coordenadas isoparametricas
deriv = shapeFuns( xi, xi , xi , 1 ) ;

% jacobiano que relaciona coordenadas materiales con isoparametricas
jacobianmat = eleCoordMat * deriv'  ;

vol = analyDet( jacobianmat ) / 6.0 ;

if vol<0,  vol, error('Element with negative volume, check connectivity.'), end

funder = inv(jacobianmat)' * deriv ;
