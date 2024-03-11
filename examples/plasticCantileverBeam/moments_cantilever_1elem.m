
% los deplazamientos que recibe sería los dnp1 v1, v2, theta1, theta2 , alpha xd

% lo viejo es kappas plas n

function Ms = moments_cantilever_1elem( v1, v2, theta1, theta2 , xd, alpha, L, kappas_plas_n )

    x = [0;L/2;L];

    Bv1 = -6/L^2*(1-2*x/L) ;
    Bv2 =  6/L^2*(1-2*x/L) ;

    Bt1 = -2/L*(2-3*x/L) ;
    Bt2 = -2/L*(1-3*x/L) ;

    %    Gtecho = -1 ;  xd 

%%%% CREO QUE NO    while % no convergio

        kappas_techo = Bv1*v1 + Bv2*v2 + Bt1*theta1 + Bt2*theta2  % +Gtecho * alpha

    kappas_plas_n1 = kappas_plas_n ; 

    Mnp1_test = E*I*( kappas_techo - kappas_plas_n1 )

    phi_test = 

    if phi_test <0

        jfñka

    end

%%%s    end
    moments
    Ms = zeros(3,1);

    


