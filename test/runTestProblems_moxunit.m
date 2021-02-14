% function for testing ONSAS using moxunit
% ----------------------------------------

function test_suite=runTestProblems_moxunit
  % initialize tests
  try
    test_functions=localfunctions();
  catch
  end
  initTestSuite;
  
function staticVonMisesTruss_test
  onsasExample_staticVonMisesTruss_test
  assertEqual( verifBoolean, 1 );

function uniaxialExtension_test
  onsasExample_uniaxialExtension_test
  assertEqual( verifBoolean, 1 );

function uniformeCurvature_test
  onsasExample_uniformCurvature_test
  assertEqual( verifBoolean, 1 );
