function ys = test_fun_to_integrate( xinput )

ys = 20 * xinput .^ 4 ;


% ys = zeros( size(xinput) ) ;
% for j=1:length( xinput)
%   x = xinput(j);

%   if x > 0.5
%     y = 100 ;
%     %y = 0 ;
%   elseif x>0,
%     y = 200 * x ;
%   elseif x>-1,
%     y = 50 * x ;
%   else
%     y = -50 ;
%   end

%   ys(j) = y ;
% end