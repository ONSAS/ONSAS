% Copyright (C) 2021, Mauricio Vanzulli, Jorge M. Perez Zerpa.
%
% This file is part of the paper a consitent corrotational aerodynamic element.
%
% NREL’s S809 Airfoil (s809-nr).http://airfoiltools.com/airfoil/details?
% airfoil=s809-nr, 2017. Accesse: 2017-07-12. 
% Mark Drela. Xfoil: An analysis and design system for low reynolds number airfoils.
% In Low Reynolds number aerodynamics, pages 1–12. Springer, 1989.
function c_m = momentCoefS809(betaRel)
  if -13 <= rad2deg( betaRel ) && rad2deg( betaRel ) <= -8 
      c_m = 0.401070 * betaRel + 0.026000 ;
  elseif -8 < rad2deg( betaRel ) && rad2deg( betaRel ) <= 8
      c_m = -0.089525 * betaRel  -0.042500;
  elseif 8 <rad2deg( betaRel ) && rad2deg( betaRel ) < 20
      c_m = -1.276618 * betaRel ^ 4 + -1.878836 * betaRel ^ 3 + 0.697222 * betaRel ^ 2 + 0.137372 * betaRel  - 0.054041 ;
  elseif rad2deg( betaRel ) >= 20
      c_m = -0.02 ;
  elseif rad2deg( betaRel ) < -13
      c_m = -0.06 ;
  end
    c_m = 0 ;
end

function computeLiftPolyFitCoefs()
    % Code to compute polyfit coeficients
    % Angle between -13 and -8
    momentPoints = [ -0.065 -0.03 ]             ;
    anglePoints = deg2rad([-13 -8  ])           ;
    p1 = polyfit(anglePoints, momentPoints, 1 ) ;
    % Angle between -8 and 8
    momentPoints2 = [ -0.03 -0.055 ]            ;
    anglePoints2 = deg2rad([-8 8 ])             ;
    p2 = polyfit(anglePoints2, momentPoints2, 1) ;
    % Angle between 8 and 20
    momentPoints3 = [ -0.055 -0.02 -0.01 -0.02 ] ;
    anglePoints3 = deg2rad([-8 10 15 20 ])       ;
    p3 = polyfit(anglePoints3, momentPoints3, 4 );
end