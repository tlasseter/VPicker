
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SubVolumeMgmt;

import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;

public class StsSubVolumeMgmtWizard extends StsWizard
{
    public StsSubVolumeMgmtActivate svActivate = null;
    public StsSubVolumeMgmtApply svApply = null;

    private StsWizardStep[] mySteps =
    {
        svActivate = new StsSubVolumeMgmtActivate(this),
        svApply = new StsSubVolumeMgmtApply(this)
    };

    public StsSubVolumeMgmtWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("Activate and Apply Sub-Volumes");
        dialog.getContentPane().setSize(300, 300);

        if(!super.start())
            return false;

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
        if(currentStep == svActivate)
            enableFinish();
        gotoNextStep();
    }

    public void finish()
    {
        super.finish();
    }
}
