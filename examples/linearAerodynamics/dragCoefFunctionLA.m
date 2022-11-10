function C_d = dragCoefFunctionLA( betaRel, Re )
    % Aeroelastic stability of a 3DOF system based on quasi-steady theory with reference to inertial coupling
    % Emulation of drag extrated from the reference above
    if 0<betaRel< deg2rad(210)
      omega_drag = 0.8 * 2 * pi / deg2rad(180);
      C_d = +1.2 - 0.6* cos( omega_drag * betaRel);
    else 
      betaRel
      error ('The angle must be beteween 0 and 210 to use this aerodynamic coefficents')
    end
end