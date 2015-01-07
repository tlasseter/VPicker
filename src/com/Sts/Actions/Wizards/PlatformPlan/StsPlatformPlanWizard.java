package com.Sts.Actions.Wizards.PlatformPlan;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;

public class StsPlatformPlanWizard extends StsWizard
{
    private StsSelectPlatform selectPlatform;
    private StsConfigurePlatform configurePlatform;
    private StsSlotConfigure configureSlots;
    private StsWellAssignment assignWells;

    private StsPlatform currentPlatform = null;

    private StsWizardStep[] mySteps =
    {
        selectPlatform = new StsSelectPlatform(this),
        configurePlatform = new StsConfigurePlatform(this),
        configureSlots = new StsSlotConfigure(this),
        assignWells = new StsWellAssignment(this)
   };

    public StsPlatformPlanWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("Platform Definition");
        dialog.getContentPane().setSize(500, 700);
        return super.start();
    }

    private void initialize()
    {
    }

    public void next()
    {
        if(currentStep == selectPlatform)
        {
            currentPlatform = selectPlatform.getSelectedPlatform();
            currentPlatform.configurePlatform();
        }
        gotoNextStep();
    }

    public void createPlatform()
    {
        currentPlatform = new StsPlatform(false);
        gotoNextStep();
    }
    public StsPlatform getPlatform() { return currentPlatform; }

    public void cancel()
    {
        super.cancel();
    }

    public boolean end()
    {
        if (!super.end())return false;
        model.setActionStatus(getClass().getName(), StsModel.STARTED);
        return true;
    }

    public void previous()
    {
        gotoPreviousStep();
    }

    public void finish()
    {
        currentPlatform.addToModel();
        super.finish();
        return;
    }

    static void main(String[] args)
    {
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
        StsPlatformPlanWizard platformWizard = new StsPlatformPlanWizard(actionManager);
        platformWizard.start();
    }
}
