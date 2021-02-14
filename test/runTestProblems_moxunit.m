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
