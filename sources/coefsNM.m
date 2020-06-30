function [a0,a1,a2,a3,a4,a5,a6,a7] = coefsNM( AlphaNW, deltaNW, deltaT)

    a0 =   1.0     / ( AlphaNW * (deltaT)^2 )                ;
    a1 =   deltaNW / ( AlphaNW *  deltaT    )                ;
    a2 =   1.0     / ( AlphaNW *  deltaT    )                ;
    a3 =   1.0     / ( AlphaNW * 2          ) - 1            ;
    a4 =   deltaNW /   AlphaNW                - 1            ;
    a5 = ( deltaNW / ( AlphaNW * 2 )          - 1 ) * deltaT ;
    a6 = ( 1 - deltaNW                            ) * deltaT ;
    a7 = deltaNW                                    * deltaT ;
  
