% Copyright (C) 2021, Mauricio Vanzulli, Jorge M. Perez Zerpa.
%
% This file is part of the paper a consitent corrotational aerodynamic element.
%
% NREL’s S809 Airfoil (s809-nr).http://airfoiltools.com/airfoil/details?
% airfoil=s809-nr, 2017. Accesse: 2017-07-12. 
% Mark Drela. Xfoil: An analysis and design system for low reynolds number airfoils.
% In Low Reynolds number aerodynamics, pages 1–12. Springer, 1989.
function c_l = liftCoefS809(betaRel)
  if abs(rad2deg( betaRel ) ) <= 8 
      c_l = 5.01338 * betaRel + 0.21667 ;
  elseif 8 < rad2deg( betaRel ) && rad2deg( betaRel ) <= 20
      c_l = -6.83918 * betaRel ^ 2 + 5.01338 * betaRel + 0.3333 ;
  elseif -15 <=rad2deg( betaRel ) && rad2deg( betaRel ) < -8
      c_l = 9.37945 * betaRel ^ 2 + 5.81143 * betaRel + 0.12857 ;
  elseif rad2deg( betaRel ) > 20
      c_l = 1.25 ;
  elseif rad2deg( betaRel ) < -15
      c_l = -0.75 ;
  end
  c_l = 0.2 ;
end

function computeLiftPolyFitCoefs()
    % Code to compute polyfit coeficients
    % Angle between -8 and 8
    liftPoints = [ -0.5 0.25 0.9]              ;
    anglePoints = deg2rad([-8 0 8 ])           ;
    p1 = polyfit(anglePoints, liftPoints, 1 )  ;
    % Angle between 8 and 20
    liftPoints2 = [ 0.9 1 1.25]                ;
    anglePoints2 = deg2rad([8 10 20 ])         ;
    p2 = polyfit(anglePoints2, liftPoints2, 2) ;
    % Angle between -8 and -15
    liftPoints3 = [ -0.5 -0.6 -0.75 ]          ;
    anglePoints3 = deg2rad([-8 -10 -15 ])      ;
    p3 = polyfit(anglePoints3, liftPoints3, 2 );
end