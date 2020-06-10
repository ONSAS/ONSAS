function Fext = computeFext( constantFext, variableFext, loadFactor, userLoadsFilename )

  if strcmp( userLoadsFilename , '')
    FextUser = zeros(size(constantFext)) ;
  else
    FextUser = feval( userLoadsFilename, nextLoadFactor)  ;
  end

  Fext  = variableFext * loadFactor + constantFext  + FextUser  ;
