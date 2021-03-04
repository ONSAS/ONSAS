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
  verifBoolean  
  assertEqual( verifBoolean, 1 );

function heatAnalytic
  example_HeatAnalytic
  assertEqual( verifBoolean, 1 );
