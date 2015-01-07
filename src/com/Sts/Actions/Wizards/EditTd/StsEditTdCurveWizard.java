//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.EditTd;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;

public class StsEditTdCurveWizard extends StsWizard
{
    private StsEditTdCurve editTdCurve = new StsEditTdCurve(this);
    private StsWizardStep[] wizardSteps = { editTdCurve };
    boolean wellNamesWereOn = false;
    StsWellClass wellClass;

	public StsEditTdCurveWizard(StsActionManager actionManager)
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
		/* jbw remove limitataion
        StsObject seismicVolume = model.getCurrentObject(StsSeismicVolume.class);
        if(seismicVolume == null)
        {
            logMessage("No seismic volume is loaded.  Terminating action.");
            return false;
        }
		*/
        if (model.getCreateStsClass(StsWell.class).getSize() == 0)
        {
            logMessage("No wells, found. Terminating action.");
            return false;
        }

        System.runFinalization();
        System.gc();
        dialog.setTitle("Edit Well TD Curve");
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
