package com.Sts.Actions.Wizards.InteractiveFilter;

import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;

public class StsInteractiveFilterWizard extends StsWizard
{
    private StsFilterVolume filterVelocityModel;

    private StsWizardStep[] mySteps =
    {
        filterVelocityModel = new StsFilterVolume(this)
    };

    public StsInteractiveFilterWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();

        dialog.setTitle("Interactive Filter");
        dialog.getContentPane().setSize(400, 200);
        if(!super.start()) return false;

        checkAddToolbar();

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
}
