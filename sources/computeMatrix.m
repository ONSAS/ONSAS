% ======================================================================
function systemDeltauMatrix = computeMatrix( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk, neumdofs, solutionMethod , bendStiff)

  % computes static tangent matrix
  [~, KT ] = assemblyFintVecTangMat( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk, bendStiff, 2 ) ;

  % performs one iteration
  if solutionMethod == 1 || solutionMethod == 2
    systemDeltauMatrix = KT ( neumdofs, neumdofs ) ;
  end
