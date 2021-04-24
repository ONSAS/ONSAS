% function for testing ONSAS using moxunit
% ----------------------------------------

function test_suite=runTestProblems_moxunit_disp
  % initialize tests
  try
    test_functions=localfunctions()
  catch
  end
  initTestSuite;
  
function test_1
  onsasExample_staticVonMisesTruss
  assertEqual( verifBoolean, true );

function test_2
  onsasExample_uniformCurvatureCantilever
  assertEqual( verifBoolean, true );
