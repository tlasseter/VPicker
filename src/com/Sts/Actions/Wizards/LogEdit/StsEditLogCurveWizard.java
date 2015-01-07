//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.LogEdit;

import com.Sts.Actions.Wizards.EditTd.StsEditTdCurve;
import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.DBTypes.StsWell;
import com.Sts.DBTypes.StsWellClass;
import com.Sts.MVC.StsActionManager;
import com.Sts.MVC.StsModel;

public class StsEditLogCurveWizard extends StsWizard
{
    private StsEditLogCurve editLogCurve = new StsEditLogCurve(this);
    private StsWizardStep[] wizardSteps = {editLogCurve};
    boolean wellNamesWereOn = false;
    StsWellClass wellClass;

	public StsEditLogCurveWizard(StsActionManager actionManager)
	{
		super(actionManager);
        addSteps(wizardSteps);
        wellClass = (StsWellClass)model.getStsClass(StsWell.class);
        wellNamesWereOn = wellClass.getDisplayNames();
        if(!wellNamesWereOn) wellClass.setDisplayNames(true);
	}

    public void previous()
    {
            gotoPreviousStep();
    }

    public void next()
    {
            gotoNextStep();
    }

    public boolean start()
    {
        if (model.getCreateStsClass(StsWell.class).getSize() == 0)
        {
            logMessage("No wells, found. Terminating action.");
            return false;
        }

        System.runFinalization();
        System.gc();
        dialog.setTitle("Edit Well Log Curve");
//        dialog.getContentPane().setSize(500, 600);
        return super.start();
    }

    public boolean end()
    {
        if(!wellNamesWereOn) wellClass.setDisplayNames(false);
        if (success)
        {
            model.setActionStatus(getClass().getName(), StsModel.STARTED);
        }
        return super.end();
    }
}
