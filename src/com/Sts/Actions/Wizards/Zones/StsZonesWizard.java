
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Zones;


import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.Model.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;


public class StsZonesWizard extends StsWizard
{
	private StsDefineZones selectZones = new StsDefineZones(this);
	private StsEditZones editZones = new StsEditZones(this);
    private StsWizardStep[] mySteps =
    {
    	selectZones,
        editZones
    };

	public StsZonesWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
    	dialog.setTitle("Define Zones and Zone Properties");
    	return super.start();
    }

    public void previous()
    {
        gotoPreviousStep();
    }

    public void next()
    {
        gotoNextStep();
    }

    public boolean end()
    {
        if(success) model.setActionStatus(StsZonesWizard.class.getName(), StsModel.STARTED);
        boolean ok = super.end();
		model.getDatabase().saveTransactions();
		return ok;
    }

    public StsZone[] getSelectedZones() { return selectZones.getSelectedZones(); }
}