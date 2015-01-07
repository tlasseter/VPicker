package com.Sts.Actions.Wizards.SimulationFile;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Export.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.Beans.*;


public class StsTranslateEclipseWizard extends StsWizard
{
	private StsTranslateEclipse translateEclipse = new StsTranslateEclipse(this);
    private StsWizardStep[] mySteps = { translateEclipse };

	public StsTranslateEclipseWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
        disableFinish();
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        StsBuiltModel builtModel = (StsBuiltModel)model.getCurrentObject(StsBuiltModel.class);
        if(builtModel == null)
        {
            //TODO replace with simple OK dialog
             StsYesNoDialog.questionValue(model.win3d,"Model has not been constructed yet. Run Model Construction wizard first.");
             return false;
        }
        dialog.setTitle("Export Simulation File");
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
}