
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.FlowSystem.PumpRates;

import com.Sts.Actions.Wizards.SensorLoad.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

public class StsPumpRatesWizard extends StsSensorLoadWizard // implements Runnable
{
    public StsSensorSelect sensorSelect = null;
    public StsSensorTime sensorTime = null;
    public StsSensorStatic sensorStatic = null;
    public StsSensorRelative sensorRelative = null;
    public StsPumpRatesLoad pumpLoad = null;

    byte vUnits = StsParameters.DIST_NONE;
    byte hUnits = StsParameters.DIST_NONE;
    byte binaryHorzUnits = StsParameters.DIST_NONE;
    byte binaryVertUnits = StsParameters.DIST_NONE;
    
    private StsSensorFile[] sensorFiles = null;

    private StsWizardStep[] mySteps =
    {
        sensorSelect = new StsSensorSelect(this),
        sensorTime = new StsSensorTime(this),
        sensorStatic = new StsSensorStatic(this),
        sensorRelative = new StsSensorRelative(this),
        pumpLoad = new StsPumpRatesLoad(this)
    };

    public StsPumpRatesWizard(StsActionManager actionManager)
    {
        super(actionManager);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("Pump Rates Wizard");
        super.initialize();
        disableFinish();
        if(!super.start()) return false;
        return true;
    }
    
    public void gotoLoad()
    {
        super.prepareFileSets();
        pumpLoad.constructPanel();
        gotoStep(pumpLoad);
    }
    
    public void initialize()
    {
        pumpLoad.setPumpRatesFactory(new StsPumpRatesFactory());
        hUnits = model.getProject().getXyUnits();
        vUnits = model.getProject().getDepthUnits();
    }

}
