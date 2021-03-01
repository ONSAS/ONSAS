% function for testing ONSAS using moxunit
% ----------------------------------------

function test_suite=runTestProblems_moxunit
  % initialize tests
  try
    test_functions=localfunctions();
  catch
  end
  initTestSuite;
  
function staticVonMisesTruss
  onsasExample_staticVonMisesTruss
  assertEqual( verifBoolean, 1 );

%~ function uniaxialExtension_test
  %~ onsasExample_uniaxialExtension_test
  %~ assertEqual( verifBoolean, 1 );

%~ function uniformCurvatureCantilever_test
  %~ onsasExample_uniformCurvatureCantilever_test
  %~ assertEqual( verifBoolean, 1 );
