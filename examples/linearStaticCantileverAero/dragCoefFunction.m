function C_d = dragCoefFunction (betaRel)
    % Aeroelastic stability of a 3DOF system based on quasi-steady theory with reference to inertial coupling
    % Emulation of drag extrated from the reference above
    if max(betaRel)< deg2rad(210)
        omega_drag = 0.8 * 2 * pi / deg2rad(180);
        C_d = +1.2 - 0.6* cos( omega_drag * betaRel);
        else
    end
    % C_d = 0;
end