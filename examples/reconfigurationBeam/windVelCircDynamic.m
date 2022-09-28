function windVel = windVelCircDynamic(x,t)

[~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, uy_vec ]  = loadParametersCirc() ;


% compute the wind vel vector according to the cy vector
windVel = [0 uy_vec(end -3) 0 ]'; 

end