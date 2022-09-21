function C_m = momentCoefFunctionLA( betaRel, Re )
    % Aeroelastic stability of a 3DOF system based on quasi-steady theory with reference to inertial coupling
    % Emulation of aerodinamic moment extracted from the reference above
    if max(betaRel) < deg2rad(210)
      betaRel = -rad2deg(betaRel);
      C_m = betaRel / 25  .*(betaRel<25) +  (-(betaRel-25)/100 + 1 ).*(betaRel>25) + 0.15;
    end
end