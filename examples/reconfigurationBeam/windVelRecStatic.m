function windVel = windVelRecStatic(x,t)

% load parameters
 [~, ~, ~, ~, ~, ~, ~, uy_vec  ] = loadParamteters() ;
% compute the wind vel vector according to the cy vector
windVel = [0 uy_vec(t) 0 ]'; 

end