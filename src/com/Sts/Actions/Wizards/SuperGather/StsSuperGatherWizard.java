package com.Sts.Actions.Wizards.SuperGather;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;

public class StsSuperGatherWizard extends StsWizard implements Runnable
{
	protected StsStackParameter stackParameter;

    protected StsPreStackLineSet currentLineSet = null;
	protected StsPreStackVelocityModel velocityModel = null;

    private StsWizardStep[] mySteps =
	{
            stackParameter = new StsStackParameter(this)
    };

	public StsSuperGatherWizard(StsActionManager actionManager)
	{
		super(actionManager, 500, 600);
		addSteps(mySteps);
 	}

	public boolean start()
	{
		System.runFinalization();
		System.gc();

		dialog.setTitle("Define Super Gather");
		dialog.getContentPane().setSize(500, 600);

		if (!super.start())
			return false;

        return true;
	}

	public boolean checkStartAction()
    {
        currentLineSet = StsPreStackLineSetClass.currentProjectPreStackLineSet;
        velocityModel = currentLineSet.getVelocityModel();
        return (velocityModel != null);
    }

	public boolean end()
	{
        model.enableDisplay();
        model.setActionStatus( getClass().getName(), StsModel.STARTED);
        return super.end();
	}

	public void next()
	{
            gotoNextStep();
	}

	public void previous()
	{
		gotoPreviousStep();
	}

    public void cancel()
    {
        super.cancel();
    }

    public void finish()
    {
        super.finish();
        model.viewObjectChangedAndRepaint(this, currentLineSet.superGatherProperties);
    }

    public StsPreStackVelocityModel getVelocityModel(){return velocityModel;}
    public void setVelocityModel(StsPreStackVelocityModel model) {velocityModel = model;}
    public StsPreStackLineSet getPreStackVolume() {return currentLineSet;}
    public void setPreStackVolume(StsPreStackLineSet volume) {currentLineSet = volume;}

    static void main(String[] args)
	{
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
		StsSuperGatherWizard superGatherWizard = new StsSuperGatherWizard(actionManager);
		superGatherWizard.start();
	}
}
