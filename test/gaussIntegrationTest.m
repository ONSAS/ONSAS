
function resultBoolean = gaussIntegrationTest()

  point_nums_to_test = [ 1:10 ]
  int_values = zeros( size(point_nums_to_test) ) ;

  a = -1.5 ;  b =  1 ;
  analyInt = -50*.5 + -50*1*-5 + 100*.5*.5 

  xs = -1.5:.02:1;
  ys = 0;

  for i=1:length(xs)
    ys(i) = fun_to_integrate(xs(i));
  end

  plot( xs, ys )

  for j=1:length( point_nums_to_test)
    [xIntPoints, wIntPoints] = GaussPointsAndWeights ( point_nums_to_test(j) ) ;

    for k=1:length(xIntPoints)
      int_values(j) = int_values(j) + ...
                      wIntPoints(k) * fun_to_integrate( (b-a)/2 * xIntPoints(k) + (a+b)/2 ) * ( (b-a) /2 ) ;
    end
  end

  figure
  plot(point_nums_to_test, int_values)

resultBoolean = true
function y = fun_to_integrate( x )

  if x > 0.5
    y = 0 ;
  elseif x>0,
    y = 200 * x ;
  elseif x>-1,
    y = 50 * x ;
  else
    y = -50 ;
  end