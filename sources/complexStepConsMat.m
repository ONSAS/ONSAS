function ConsMat = complexStepConsMat( stressFun, consParams, Egreen )
  
  compStep = 1e-8 ;

  for i=1:6

    direction = zeros( 3,3) ;
    
    if i<4
      direction(i,i) = 1 ;
    elseif i==4
      direction(2,3) = 1 ;
      direction(3,2) = 1 ;
    elseif i==5
      direction(1,3) = 1 ;
      direction(3,1) = 1 ;
    elseif i==6
      direction(1,2) = 1 ;
      direction(2,1) = 1 ;
    end

    EgreenComp = Egreen  + compStep * direction * j ;
  
    Scomp = feval( stressFun, consParams, EgreenComp, 0 ) ;
    
    dsde = mat2voigt( ( imag( Scomp ) / compStep ) , 1 ) ;
  
    ConsMat(1:3,i) = dsde(1:3) ;
    ConsMat(4:6,i) = dsde(4:6)*0.5 ;
    %~ ConsMat(:,i) = dsde ;
  end
