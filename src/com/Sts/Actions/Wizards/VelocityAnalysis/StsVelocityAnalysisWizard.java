package com.Sts.Actions.Wizards.VelocityAnalysis;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;

public class StsVelocityAnalysisWizard extends StsWizard implements Runnable
{
    protected StsDefineVelocityAnalysis defineVelocityAnalysis;
    protected StsIdentifyVelocityModels identifyVelocityModels;
    protected StsExtractVelocityProfile extractProfile;
    protected StsPreStackLineSet currentLineSet = null;
	protected StsPreStackVelocityModel velocityModel = null;
    protected StsPreStackVelocityModel inputVelocityModel = null;

    private Object inputVelocityModelOrVolume = StsPreStackLineSet.NO_MODEL;

    private StsWizardStep[] mySteps =
	{
			defineVelocityAnalysis = new StsDefineVelocityAnalysis(this),
			identifyVelocityModels = new StsIdentifyVelocityModels(this),
            extractProfile = new StsExtractVelocityProfile(this)
     };

	public StsVelocityAnalysisWizard(StsActionManager actionManager)
	{
		super(actionManager, 500, 500);
		addSteps(mySteps);
 	}

	public boolean start()
	{
		System.runFinalization();
		System.gc();

		dialog.setTitle("Define Velocity Analysis");
		dialog.getContentPane().setSize(500, 600);

		if (!super.start())
			return false;

        return true;
	}

	/** If we already have a velocity model check with user about creating a new one;
	 *  if the user wants a new one, or one doesn't exist: start action.
    */
	public boolean checkStartAction()
	{
		currentLineSet = StsPreStackLineSetClass.currentProjectPreStackLineSet;
		if (currentLineSet == null)
			return false;

        if( currentLineSet.getVelocityModel() != null)
        {
            if( currentLineSet.getVelocityModel().getVelocityProfiles().getSize() > 0)
                return StsYesNoDialog.questionValue(this.frame, "Replace existing velocity model?\n\n Yes - New model\n No - Cancel request");
        }
        return true;
    }

	public boolean end()
	{
        model.enableDisplay();
        model.setActionStatus( getClass().getName(), StsModel.STARTED);
        return super.end();
	}

	public void next()
	{
        if( currentStep == identifyVelocityModels)
        {
            createVelocityModel();
            if( getCreateProfiles() ||
                inputVelocityModelOrVolume == StsPreStackLineSet.HAND_VEL)
            {
                extractProfile.constructPanel();
                gotoStep( extractProfile);
            }
            else
                finish();
        }
        else
            gotoNextStep();
	}

	public void previous()
	{
		gotoPreviousStep();
	}

    public void finish()
    {
        super.finish();

        if (success)
        {
            setCursor();
            addToolbar();
            currentLineSet.lineSetClass.setShowVelStat(true);
        }
        else
        {
            if( velocityModel != null)
                velocityModel.delete();
            if( inputVelocityModel != null)
                inputVelocityModel.delete();
        }
        model.win3dDisplayAll();
    }

    private void createVelocityModel()
    {
        // Use case 1 -- Dataset with no predecessor and no velocity model
        if( inputVelocityModelOrVolume == StsPreStackLineSet.NO_MODEL)
        {
            // This is the first velocity model
            velocityModel = currentLineSet.constructVelocityModel();
            velocityModel.initializeVelocityProfiles();
            currentLineSet.setVelocityModel(velocityModel);
            return;
        }

        // Use case 2 -- We have an input velocityVolume with or without handVels OR we have a 2d associated line set
        // Create an input and output velocity model;
        if( inputVelocityModelOrVolume instanceof StsPreStackVelocityModel ||inputVelocityModelOrVolume instanceof StsSeismic || inputVelocityModelOrVolume == StsPreStackLineSet.HAND_VEL)
        {
            inputVelocityModel = currentLineSet.constructVelocityModel();
            inputVelocityModel.initializeVelocityProfiles();
            inputVelocityModel.setName(currentLineSet.getName() + "-InputVelocity-" + currentLineSet.getIndex() + velocityModel.nameExtension);
            if( inputVelocityModelOrVolume instanceof StsSeismicVolume)
                inputVelocityModel.setSeismicVelocityVolume((StsSeismicVolume)inputVelocityModelOrVolume);
            currentLineSet.setInputVelocityModel(inputVelocityModel);

            velocityModel = currentLineSet.constructVelocityModel();
            velocityModel.initializeVelocityProfiles();
            velocityModel.setName(currentLineSet.getName() + "-OutputVelocity-" + currentLineSet.getIndex() + velocityModel.nameExtension);
            currentLineSet.setVelocityModel(velocityModel);
        }
    }

	protected void addToolbar()
	{
		StsVelocityAnalysisToolbar.checkAddToolbar(model, currentLineSet, true);
	}

	protected void setCursor()
	{
        currentLineSet.setCurrentDataRowCol(model.win3d, 0, 0);
        model.getCursor3d().setCol(0);  //current velocity profile now stored in row/col of cursor3d
        model.getCursor3d().setRow(0);
        currentLineSet.nextProfile(); 
        boolean is3d = StsPreStackLineSetClass.currentProjectPreStackLineSet instanceof StsPreStackLineSet3d;
        if (!is3d)
            currentLineSet.previousProfile(); //a bit of a hack - for some reason "nextProfile()" actually takes you to 2nd profile 1st time through on 2d
	}

    public double getVelScaleMultiplier() { return identifyVelocityModels.panel.getVelScaleMultiplier(); }
    public StsPreStackVelocityModel getVelocityModel(){return velocityModel;}
    public void setVelocityModel(StsPreStackVelocityModel model) {velocityModel = model;}
    public StsPreStackLineSet getPreStackVolume() {return currentLineSet;}
    public void setPreStackVolume(StsPreStackLineSet volume) {currentLineSet = volume;}
    public void setInputVelocityModelOrVolume(Object velocityModel) {inputVelocityModelOrVolume = velocityModel;}
    public Object getInputVelocityModelOrVolume() {return inputVelocityModelOrVolume;}
    public StsPreStackVelocityModel getInputVelocityModel(){return inputVelocityModel;}
    public void setInputVelocityModel(StsPreStackVelocityModel model) {inputVelocityModel = model;}
    public StsPreStackLineSet getCurrentLineSet() {return currentLineSet;}
    public boolean getCreateProfiles() {return identifyVelocityModels.panel.getCreateProfiles();}
    public int getNumVelPts() {return identifyVelocityModels.panel.getNumVelPts();}    

    static void main(String[] args)
	{
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
		StsVelocityAnalysisWizard velAnalysisWizard = new StsVelocityAnalysisWizard(actionManager);
		velAnalysisWizard.start();
	}
}