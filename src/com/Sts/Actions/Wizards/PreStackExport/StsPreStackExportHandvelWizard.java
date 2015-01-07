package com.Sts.Actions.Wizards.PreStackExport;

import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;

public abstract class StsPreStackExportHandvelWizard extends StsWizard 
{
	private StsExportHandvel[] handVelSteps;
	
    public StsPreStackExportHandvelWizard(StsActionManager actionManager)
    {
        super(actionManager);
        handVelSteps = createHandVelSteps();
        addSteps(handVelSteps);
    }

    protected abstract StsExportHandvel[] createHandVelSteps();

	public boolean start()
    {
        System.runFinalization();
        System.gc();

        dialog.setTitle("Export Pre-Stack Data");
        dialog.getContentPane().setSize(400, 300);
        if(!super.start()) return false;
        
        return true;
    }

    public boolean end()
    {
        return super.end();
    }

    public void previous()
    {
        gotoPreviousStep();
    }

    public void next()
    {
    	gotoNextStep();
    }
    
    @Override
    public void finish() 
    {
    	handVelSteps[0].exportHandVels();
    	super.finish();
    }
}

