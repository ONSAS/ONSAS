  handlesFig.anaPanel = uipanel ('parent', handlesFig.inputPanel, "position", vecAnaPanel) ;
  
  % From bottom geometry
  set(handlesFig.anaPanel, 'units', 'pixels')
  vec = get(handlesFig.anaPanel, 'position') ;
  PPanaFromBottom = vec(4)-ref ;
  PPanaFromLeft = 0.11 * anchoPantalla ;
  
  
  ##  , 'bordertype', 'none'
  handlesFig.anaTitle = uicontrol('parent', handlesFig.anaPanel, 'style', 'text', 'string', 'Analysis options', 'HorizontalAlignment', 'left', 'fontsize', 11, 'position', [PPFromLeft PPanaFromBottom buttonWidth buttonHeight]) ; 
  % Analysis radiobuttons 
  handlesFig.linAna = uicontrol('parent', handlesFig.anaPanel, 'style', 'radiobutton', 'string', 'Linear Analysis', 'position', [PPFromLeft PPanaFromBottom-buttonHeight buttonWidth buttonHeight ], 'callback', {@anaChoice, fig, buttonHeight, buttonWidth} ) ;
  handlesFig.dynNonLinAna = uicontrol('parent', handlesFig.anaPanel, 'style', 'radiobutton', 'string', 'Dynamic NonLinear Analysis', 'position', [PPFromLeft PPanaFromBottom-2*buttonHeight buttonWidth buttonHeight ], 'callback', {@anaChoice, fig, buttonHeight, buttonWidth} ) ;
  handlesFig.nonLinAna = uicontrol('parent', handlesFig.anaPanel, 'style', 'radiobutton', 'string', 'NonLinear Analysis', 'position', [PPFromLeft PPanaFromBottom-3*buttonHeight buttonWidth buttonHeight ], 'callback', {@anaChoice, fig, buttonHeight, buttonWidth} ) ;
  %~ handlesFig.dynNonLinAna = uicontrol('parent', handlesFig.anaPanel, 'style', 'radiobutton', 'string', 'Dynamic NonLinear Analysis', 'position', [PPFromLeft PPanaFromBottom-4*buttonHeight buttonWidth buttonHeight ], 'callback', {@anaChoice, fig, buttonHeight, buttonWidth} ) ;
  % ------------------------------------------------------------------------------------------------------------------------
  % auxiliar
  handlesFig.visible = 'off' ;
  handlesFig.visibleLin = 'off' ;

  % Numerical Method options
  handlesFig.numericalTitle = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible, 'style', 'text', 'string', 'Numerical Method parameters','HorizontalAlignment', 'left', 'fontsize', 11, 'position', [PPFromLeft PPanaFromBottom-5*buttonHeight 1.15*buttonWidth buttonHeight]) ;
  handlesFig.numMethod = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible, 'units', 'pixels', 'style', 'popupmenu', 'position', [PPFromLeft PPanaFromBottom-6*buttonHeight buttonWidth buttonHeight], 'callback', {@numMethodChoice, fig}  );
  
  
  % --------------------- Linear Parameters ------------------------------
  % delta T
  handlesFig.deltaTText = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visibleLin,'style', 'text', 'string', 'Delta t:', 'horizontalalignment', 'left', 'visible', 'off', 'position', [PPFromLeft PPanaFromBottom-7*buttonHeight textW textH]) ;
  handlesFig.deltaT = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visibleLin,'style', 'edit', 'string', handlesFig.prevValdeltaTLin, 'position', [PPanaFromLeft PPanaFromBottom-7*buttonHeight editW editH], 'tag' , 'delta T', 'callback', {@editBoxes, fig}) ;
  % final Time
  handlesFig.finalTimeLinText = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible,'style', 'text', 'string', 'Final time:', 'horizontalalignment', 'left', 'visible', 'off','position', [PPFromLeft PPanaFromBottom-8*buttonHeight textW textH] ) ;
  handlesFig.finalTimeLin = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible,'style', 'edit', 'string', handlesFig.prevValfinalTimeLin, 'position', [PPanaFromLeft PPanaFromBottom-8*buttonHeight editW editH], 'tag' , 'finalTime', 'callback', {@editBoxes, fig}) ;
  % indexTime Sol
  handlesFig.indexTimeSolText = uicontrol('parent', handlesFig.anaPanel, 'style', 'text', 'string', 'Index Time Sol', 'horizontalalignment', 'left', 'visible', handlesFig.visible, 'position', [PPFromLeft PPanaFromBottom-9*buttonHeight textW textH] );  
  handlesFig.indexTimeSol = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible, "style", "edit", 'string', handlesFig.prevValIndexTimeSol, "position",[PPanaFromLeft PPanaFromBottom-9*buttonHeight editW editH], 'tag', 'indexTimeSol', 'callback', {@editBoxes, fig}) ;
  
  % --------------------- Non Linear Parameters --------------------------
  % DeltaU 
  handlesFig.tolDeltaUText = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible,'style', 'text', 'string', 'Tolerance Delta u:', 'horizontalalignment', 'left', 'visible', 'off', 'position', [PPFromLeft PPanaFromBottom-7*buttonHeight textW textH]) ;
  handlesFig.tolDeltaU = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible,'style', 'edit', 'string', handlesFig.prevValTolDeltaU, 'position', [PPanaFromLeft PPanaFromBottom-7*buttonHeight editW editH], 'tag' , 'tolDeltaU', 'callback', {@editBoxes, fig}) ;
  % tolForces    
  handlesFig.tolForcesText = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible,'style', 'text', 'string', 'Tolerance Forces:', 'horizontalalignment', 'left', 'visible', 'off','position', [PPFromLeft PPanaFromBottom-8*buttonHeight textW textH] ) ;
  handlesFig.tolForces = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible,'style', 'edit', 'string', handlesFig.prevValTolForces, 'position', [PPanaFromLeft PPanaFromBottom-8*buttonHeight editW editH], 'tag' , 'tolForces', 'callback', {@editBoxes, fig}) ;
  % tolIter    
  handlesFig.tolIterText = uicontrol('parent', handlesFig.anaPanel, 'style', 'text', 'string', 'Tolerance Iterations:', 'horizontalalignment', 'left', 'visible', handlesFig.visible, 'position', [PPFromLeft PPanaFromBottom-9*buttonHeight textW textH] );  
  handlesFig.tolIter = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible, "style", "edit", 'string', handlesFig.prevValTolIter, "position",[PPanaFromLeft PPanaFromBottom-9*buttonHeight editW editH], 'tag', 'tolIter', 'callback', {@editBoxes, fig}) ;
  % target Load Fac
  handlesFig.targetLoadFacText = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible, 'style', 'text', 'string', 'Target load factor:', 'horizontalalignment', 'left', 'visible', 'off', 'position', [PPFromLeft PPanaFromBottom-10*buttonHeight textW textH] ) ;
  handlesFig.targetLoadFac = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible, "style", "edit", 'string', handlesFig.prevValLoadFac, "position", [PPanaFromLeft PPanaFromBottom-10*buttonHeight editW editH], 'tag', 'loadFac', 'callback', {@editBoxes, fig}) ;
  % load Steps  
  handlesFig.loadStepsText = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible, 'style', 'text', 'string', 'Load steps:', 'horizontalalignment', 'left', 'visible', 'off', 'position', [PPFromLeft PPanaFromBottom-11*buttonHeight textW textH] ) ;
  handlesFig.loadSteps = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible, "style", "edit", 'string', handlesFig.prevValLoadSteps, "position", [PPanaFromLeft PPanaFromBottom-11*buttonHeight editW editH], 'tag', 'loadSteps', 'callback', {@editBoxes, fig}) ;
  % increment Arc Length
  handlesFig.incremALText = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible, 'style', 'text', 'string', 'Increment AL:', 'horizontalalignment', 'left', 'visible', 'off', 'position', [PPFromLeft PPanaFromBottom-12*buttonHeight textW textH] ) ;
  handlesFig.incremAL = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible, "style", "edit", 'string', handlesFig.prevValIncremAL, "position", [PPanaFromLeft PPanaFromBottom-12*buttonHeight editW editH], 'tag', 'incremAL', 'callback', {@editBoxes, fig}) ;
  
  % --------------------- Dynamic Parameters --------------------------
  % time Increment 
  handlesFig.timeIncrText = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible,'style', 'text', 'string', 'Time increment:', 'horizontalalignment', 'left', 'visible', 'off', 'position', [PPFromLeft PPanaFromBottom-7*buttonHeight textW textH]) ;
  handlesFig.timeIncr = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible,'style', 'edit', 'string', handlesFig.prevValTimeIncr, 'position', [PPanaFromLeft PPanaFromBottom-7*buttonHeight editW editH], 'tag' , 'timeIncr', 'callback', {@editBoxes, fig}) ;
  % final time    
  handlesFig.finalTimeText = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible,'style', 'text', 'string', 'Final time:', 'horizontalalignment', 'left', 'visible', 'off','position', [PPFromLeft PPanaFromBottom-8*buttonHeight textW textH] ) ;
  handlesFig.finalTime = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible,'style', 'edit', 'string', handlesFig.prevValFinalTimeDyn, 'position', [PPanaFromLeft PPanaFromBottom-8*buttonHeight editW editH], 'tag' , 'finalTimeDyn', 'callback', {@editBoxes, fig}) ;
  % stopTol DeltaU    
  handlesFig.tolDeltaUDynText = uicontrol('parent', handlesFig.anaPanel, 'style', 'text', 'string', 'Tolerance Delta U:', 'horizontalalignment', 'left', 'visible', handlesFig.visible, 'position', [PPFromLeft PPanaFromBottom-9*buttonHeight textW textH] );  
  handlesFig.tolDeltaUDyn = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible, "style", "edit", 'string', handlesFig.prevValTolDeltaUDyn, "position",[PPanaFromLeft PPanaFromBottom-9*buttonHeight editW editH], 'tag', 'tolDeltaUDyn', 'callback', {@editBoxes, fig}) ;
  % stopTol forces
  handlesFig.tolForcesDynText = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible, 'style', 'text', 'string', 'Tolerance Forces:', 'horizontalalignment', 'left', 'visible', 'off', 'position', [PPFromLeft PPanaFromBottom-10*buttonHeight textW textH] ) ;
  handlesFig.tolForcesDyn = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible, "style", "edit", 'string', handlesFig.prevValTolForcesDyn, "position", [PPanaFromLeft PPanaFromBottom-10*buttonHeight editW editH], 'tag', 'tolForcesDyn', 'callback', {@editBoxes, fig}) ;
  % stopTol Iterations  
  handlesFig.tolIterDynText = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible, 'style', 'text', 'string', 'Tolerance Iterations:', 'horizontalalignment', 'left', 'visible', 'off', 'position', [PPFromLeft PPanaFromBottom-11*buttonHeight textW textH] ) ;
  handlesFig.tolIterDyn = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible, "style", "edit", 'string', handlesFig.prevValTolIterDyn, "position", [PPanaFromLeft PPanaFromBottom-11*buttonHeight editW editH], 'tag', 'tolIterDyn', 'callback', {@editBoxes, fig}) ;
  % DeltaNW
  handlesFig.deltaNWText = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible, 'style', 'text', 'string', 'Delta NW:', 'horizontalalignment', 'left', 'visible', 'off', 'position', [PPFromLeft PPanaFromBottom-12*buttonHeight textW textH] ) ;
  handlesFig.deltaNW = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible, "style", "edit", 'string', handlesFig.prevValDeltaNW, "position", [PPanaFromLeft PPanaFromBottom-12*buttonHeight editW editH], 'tag', 'deltaNW', 'callback', {@editBoxes, fig}) ;
  % AlphaNW
  handlesFig.alphaNWText = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible, 'style', 'text', 'string', 'Alpha NW:', 'horizontalalignment', 'left', 'visible', 'off', 'position', [PPFromLeft PPanaFromBottom-13*buttonHeight textW textH] ) ;
  handlesFig.alphaNW = uicontrol('parent', handlesFig.anaPanel, 'visible', handlesFig.visible, "style", "edit", 'string', handlesFig.prevValAlphaNW, "position", [PPanaFromLeft PPanaFromBottom-13*buttonHeight editW editH], 'tag', 'alphaNW', 'callback', {@editBoxes, fig}) ;
  % ------------------------------------------------------------------------------------------------------------------------
