% Implicit function
function val = c_val(c,p,Y,a,b)
	val = p/Y-( log(c/a)+1/2*(1 - c^2/b^2) ) ;
end
