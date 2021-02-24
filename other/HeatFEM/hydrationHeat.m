function Qint = hydrationHeat( t )
  global rho
  
  Q0 = 330.0 ;
  a = 0.69 / 24.0 ;
  b = 0.56 ;
  
  Qint = Q0 * ( -1 * (-a*b*t.^( b-1 ) * exp( -a * t.^b ) ) ) * rho ;
end
