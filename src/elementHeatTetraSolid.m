function [ KTe, Me, Qhe ] = elementHeatTetraSolid( elemCoords, materialParams )

  tetCoordMat        = reshape( elemCoords', 3, 4 ) ;

  %~ [ deriv , vol ] =  DerivFun( tetCoordMat ) ;

  %~ tetVol = vol ;
  %~ BMat = BMats ( deriv ) ;

  % material params
  rho   = materialParams(1) ;  
  cSpHe = materialParams(2) ;
  kCond = materialParams(3) ;


  eleCoordMat = tetCoordMat               ; 

  xi = 0.25 ;  wi = 1/6  ;

  % matriz de derivadas de fun forma respecto a coordenadas isoparametricas
  deriv = shapeFuns( xi, xi , xi , 1 ) ;

  % jacobiano que relaciona coordenadas materiales con isoparametricas
  jacobianmat = eleCoordMat * deriv'  ;

  vol = analyDet( jacobianmat ) / 6.0 ;

  if vol<0,  vol, error('Element with negative volume, check connectivity.'), end

  funder = inv(jacobianmat)' * deriv ;
  
  KTe = kCond * funder' * funder * vol ;
  Me  = rho * cSpHe * vol/4 * eye(4) ;
  Qhe = vol / 4 * ones(4,1);
