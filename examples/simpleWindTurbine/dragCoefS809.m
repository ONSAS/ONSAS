% Copyright (C) 2021, Mauricio Vanzulli, Jorge M. Perez Zerpa.
%
% This file is part of the paper a consitent corrotational aerodynamic element.
%
% NREL’s S809 Airfoil (s809-nr).http://airfoiltools.com/airfoil/details?
% airfoil=s809-nr, 2017. Accesse: 2017-07-12. 
% Mark Drela. Xfoil: An analysis and design system for low reynolds number airfoils.
% In Low Reynolds number aerodynamics, pages 1–12. Springer, 1989.
function c_d = dragCoefS809(betaRel)
  if abs(rad2deg( betaRel ) ) <= 8 
      c_d = 0.01 ;
  elseif 8 < rad2deg( betaRel ) && rad2deg( betaRel ) <= 20
      c_d = 0.8207016 * betaRel ^ 2 + 0.8207016 * betaRel + -0.01 ;
  elseif -13 <=rad2deg( betaRel ) && rad2deg( betaRel ) < -8
      c_d = -0.458366 * betaRel + -0.054000 ;
  elseif rad2deg( betaRel ) > 20
      c_d = 0.1 ;
  elseif rad2deg( betaRel ) < -13
      c_d = 0.05 ;
  else 
  end
  c_d = 0 ;
end

function computeDragPolyFitCoefs()
    % Code to compute polyfit coeficients
    % Angle between 8 and 20
    dragPoints = [ .01 .02 .1 ]               ;
    anglePoints = deg2rad([8 10 20])          ;
    p1 = polyfit(anglePoints, dragPoints, 2 ) ;

    % Angle between -13 and -8
    dragPoints2 = [ 0.05 .01 ]                  ;
    anglePoints2 = deg2rad([-13 -8])            ;
    p2 = polyfit(anglePoints2, dragPoints2, 1 ) ;
end