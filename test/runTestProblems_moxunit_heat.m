% function for testing ONSAS using moxunit
% ----------------------------------------

function test_suite=runTestProblems_moxunit_heat
  % initialize tests
  try
    test_functions=localfunctions()
  catch
  end
  initTestSuite;
  
function test_1
  example_HeatAnalytic
  verifBoolean  
  assertEqual( verifBoolean, true );

function test_2
  example_HeatRobiAndNeum
  verifBoolean  
  assertEqual( verifBoolean, true );
