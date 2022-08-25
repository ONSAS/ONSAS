function epsk = epsVal(z, B, R, Ut)

	epsk = -z'.* B * R * Ut ;

end
