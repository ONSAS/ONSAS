function C_l = liftCoefFunctionLA( betaRel, Re )
    % Aeroelastic stability of a 3DOF system based on quasi-steady theory with reference to inertial coupling
    % Emulation of lift extracted from the reference above
    if max(betaRel) < deg2rad(210)
      omega_drag = 2*pi / deg2rad(180);
      C_l = 0.5 + 1.5 * sin( omega_drag * betaRel);
    else 
      betaRel
      error ("The angle must be beteween 0º and 210 to use this aerodynamic coefficents")
    end
end
