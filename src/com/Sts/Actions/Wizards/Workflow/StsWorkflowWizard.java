package com.Sts.Actions.Wizards.Workflow;

import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;

public class StsWorkflowWizard extends StsWizard
{
    public StsDefineWorkflow defineWorkflow;

    private StsWizardStep[] mySteps =
    {
        defineWorkflow = new StsDefineWorkflow(this)
    };

    public StsWorkflowWizard(StsActionManager actionManager)
    {
        super(actionManager);
        if(Main.singleWorkflow)
        {
            new StsMessage(this.frame, StsMessage.INFO,
                           "While many capabilites exist beyond those currently licensed, \n" +
                           "your license is limited to one workflow. If you desire access \n" +
                           "to expanded capabilites and workflows, please contact the us \n" +
                           "to purchase the additional capabilites.");
            return;
        }
        addSteps(mySteps);
    }

    public boolean start()
    {

        System.runFinalization();
        System.gc();

        dialog.setTitle("Select Workflow Screen");
        dialog.getContentPane().setSize(300, 300);
        if(!super.start()) return false;
        enableFinish();
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

    static void main(String[] args)
    {
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
        StsWorkflowWizard workflowWizard = new StsWorkflowWizard(actionManager);
        workflowWizard.start();
    }
}
