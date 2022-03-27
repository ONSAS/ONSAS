
function [Area, J, Iyy, Izz, Jrho ] = crossSectionProps ( elemCrossSecParams, rho )
% --- cross section ---
if strcmp( elemCrossSecParams{1}, 'generic' ) %general section
    elemCrossSecParamsVec = elemCrossSecParams{2} ;
    Area = elemCrossSecParamsVec( 1 ) ;
    J    = elemCrossSecParamsVec( 2 ) ;
    Iyy  = elemCrossSecParamsVec( 3 ) ;
    Izz  = elemCrossSecParamsVec( 4 ) ;
    %
    if length( elemCrossSecParamsVec ) > 5
        Jrho =  diag( elemCrossSecParamsVec( 5:7 ) ) ;
    else
        Jrho = rho * diag( [ J Iyy Izz ] ) ;
    end

elseif strcmp( elemCrossSecParams{1}, 'rectangle' )
    elemCrossSecParamsVec = elemCrossSecParams{2} ;
    Area = elemCrossSecParamsVec( 1 ) *elemCrossSecParamsVec( 2 )        ;
    Iyy  = elemCrossSecParamsVec( 1 ) *elemCrossSecParamsVec( 2 ) ^ 3/12 ;
    Izz  = elemCrossSecParamsVec( 1 ) *elemCrossSecParamsVec( 2 ) ^ 3/12 ;

    if elemCrossSecParamsVec( 1 )== elemCrossSecParamsVec( 2 )
        J    = 1/3*0.40147*elemCrossSecParamsVec(2)^4 ;
    else
        h = max( elemCrossSecParamsVec(1:2) ) ;
        b = min( elemCrossSecParamsVec(1:2) ) ;
        r = h/b ;
        beta = -1e-5*r^6 + 4e-4*r^5 - 5.8e-3*r^4 + 4.1e-2*r^3 - 0.1625*r^2 + 0.3628*r - 0.0949 ;
        J = beta * b^3*h ;
    end
    Jrho = rho * diag( [ J Iyy Izz ] ) ;

elseif strcmp( elemCrossSecParams{1}, 'circle' )
    diameter = elemCrossSecParams{2} ;
    Area = pi*diameter^2/4             ;
    Iyy  = pi*diameter^4/64            ;
    Izz  = Iyy                         ;
    J    = Iyy + Izz                   ;
    Jrho = rho * diag( [ J Iyy Izz ] ) ;
else
  error(' section type not implemented yet, please create an issue')
end
