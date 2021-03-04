% function for testing ONSAS using moxunit
% ----------------------------------------

function test_suite=runTestProblems_moxunit
  % initialize tests
  try
    test_functions=localfunctions()
  catch
  end
  initTestSuite;
  
function test_1
  onsasExample_staticVonMisesTruss
  verifBoolean  
  assertEqual( verifBoolean, 1 );

function test_2
  example_HeatAnalytic
  assertEqual( verifBoolean, 1 );
