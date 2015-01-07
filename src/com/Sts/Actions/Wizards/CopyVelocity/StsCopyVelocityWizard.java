package com.Sts.Actions.Wizards.CopyVelocity;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Toolbars.*;

public class StsCopyVelocityWizard extends StsWizard implements Runnable
{
	protected StsCopyVelocity copyVelocity;

    protected StsPreStackLineSet currentLineSet = null;
    protected StsPreStackVelocityModel velocityModel = null;
    protected StsPreStackVelocityModel velocityModelCopy = null;
    protected StsPreStackLineSet lineSetCopy = null;
    private StsWizardStep[] mySteps =
	{
            copyVelocity = new StsCopyVelocity(this)
    };

	public StsCopyVelocityWizard(StsActionManager actionManager)
	{
		super(actionManager, 500, 600);
		addSteps(mySteps);
 	}

	public boolean start()
	{
		System.runFinalization();
		System.gc();

		dialog.setTitle("Copy Velocity Model");
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
        createVelocityModel();

        if( lineSetCopy != null)
        {
            setCursor();
            addToolbar();
            lineSetCopy.lineSetClass.setShowVelStat( true);
            model.viewObjectChangedAndRepaint(this, lineSetCopy);
        }
    }

    private void addToolbar()
    {
        StsVelocityAnalysisToolbar.checkAddToolbar(model, lineSetCopy, true);
    }

    private void setCursor()
    {
        lineSetCopy.setCurrentDataRowCol(model.win3d, 0, 0);
        lineSetCopy.nextProfile();
    }

    private void createVelocityModel()
    {
        // Assign the current velocity model as inputVelocityModel and create an outputVelocity;
        // clone all profiles from input velocityModel
        velocityModelCopy = (StsPreStackVelocityModel )velocityModel.clone();

        lineSetCopy = (StsPreStackLineSet )currentLineSet.clone();
        lineSetCopy.initializeSuperGathers();

        lineSetCopy.setName( copyVelocity.panel.getOutputVelocityName() + lineSetCopy.getIndex());
        lineSetCopy.setVelocityModel( velocityModelCopy);

        velocityModelCopy.setPreStackLineSet( lineSetCopy);
        velocityModelCopy.setName( copyVelocity.panel.getOutputVelocityName() + "-OutputVelocity" + velocityModel.nameExtension);

        StsObjectList velocityProfileList = new StsObjectList(0);
        StsObjectRefList inProfiles  = velocityModel.getVelocityProfiles();
        for (int n = 0; n < inProfiles.getSize(); n++)
        {
            StsVelocityProfile inProfile = (StsVelocityProfile )inProfiles.getElement(n);
            StsVelocityProfile velocityProfile = new StsVelocityProfile( inProfile.row, inProfile.col, inProfile);
            velocityProfileList.add( velocityProfile);
        }

        velocityModelCopy.setVelocityProfiles( velocityProfileList.convertListToRefList(model, "velocityProfiles", velocityModelCopy));
        velocityModelCopy.initializeInterpolation();
        velocityModelCopy.checkInitializeInterpolation();
        velocityModelCopy.initializeInterpolateVelocityProfiles();
   }

    public StsPreStackVelocityModel getVelocityModel(){return velocityModel;}
    public void setVelocityModel(StsPreStackVelocityModel model) {velocityModel = model;}
    public StsPreStackLineSet getPreStackVolume() {return currentLineSet;}
    public void setPreStackVolume(StsPreStackLineSet volume) {currentLineSet = volume;}

    static void main(String[] args)
	{
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
		StsCopyVelocityWizard copyVelocityWizard = new StsCopyVelocityWizard(actionManager);
		copyVelocityWizard.start();
	}
}
