% ========================================================================================
% ========================================================================================
function resultBoolean = gaussIntegrationTest()

  plots_boolean = false ;

  point_nums_to_test = [ 1:10 12 14 16   ] ;
%  point_nums_to_test = [ 1:10    ]
  int_values = zeros( size(point_nums_to_test) ) ;

  a = -1.5 ;  b =  1 ;
  %analyInt = -50*.5 + -50*1*.5 + 100*.5*.5 + 100*.5  ;
  analyInt = 20 * ( b^5/5 - a^5/5 )  ;

  xs = -1.5:.02:1;
  ys = 0;
  for i=1:length(xs)
    ys(i) = test_fun_to_integrate(xs(i));
  end

  if plots_boolean
    figure
    plot( xs, ys, 'g-x' )
    title('test_fun_to_integrate')
  end

  for j=1:length( point_nums_to_test)
    [xIntPoints, wIntPoints] = GaussPointsAndWeights ( point_nums_to_test(j) ) ;

    for k=1:length(xIntPoints)
      int_values(j) = int_values(j) + ...
                      wIntPoints(k) * test_fun_to_integrate( (b-a)/2 * xIntPoints(k) + (a+b)/2 ) * ( (b-a) /2 ) ;
    end
  end
  
  if plots_boolean
    numericalInt = quadl( 'test_fun_to_integrate', a, b )
    figure
    plot(point_nums_to_test, int_values,'b-x')
    hold on, grid on
    plot(point_nums_to_test, numericalInt*ones(size(point_nums_to_test)),'r-o')
  end

  resultBoolean = max( abs( int_values(3:end) - analyInt ) ) / abs( analyInt ) < 1e-8 ;

% ========================================================================================
% ========================================================================================
function ys = test_fun_to_integrate( xinput )
  ys = 20 * xinput .^ 4 ;
