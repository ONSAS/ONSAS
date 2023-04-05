% Implicit function
function val = c_val(c,p,Y,Ri,Re)
	val = p/Y-( log(c/Ri)+1/2*(1 - c^2/Re^2) ) ;
end
