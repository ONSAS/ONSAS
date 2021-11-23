function C_l = liftCoefFunction (betaRel)
    % Aeroelastic stability of a 3DOF system based on quasi-steady theory with reference to inertial coupling
    % Emulation of lift extracted from the reference above
    if max(betaRel) < deg2rad(210)
        omega_drag = 2*pi / deg2rad(180);
        C_l = 0.5 + 1.5 * sin( omega_drag * betaRel);
    end
end
