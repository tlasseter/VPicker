package com.Sts.Actions.Wizards.WellPlan;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

public class StsLoadWellPlanWizard extends StsWizard
{
    byte vUnits = StsParameters.DIST_NONE;
    byte hUnits = StsParameters.DIST_NONE;
    private StsWellPlanSelect wellSelect;
    private StsWellPlanLoad loadWells;

    private StsWellPlanSet planSet = null;
    private StsWellPlan wellPlan = null;

    private StsWizardStep[] mySteps =
    {
        wellSelect = new StsWellPlanSelect(this),
        loadWells = new StsWellPlanLoad(this),
//        definePath = new StsDefinePath(this),
//        relateWell = new StsWellRelate(this),
//        defineVertices = new StsDefineVertices(this)
   };

    public StsLoadWellPlanWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);

    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();

        dialog.setTitle("Define Well Plan");
        dialog.getContentPane().setSize(500, 500);
        initialize();
        if(!super.start()) return false;

        return true;
    }

    public boolean checkPrerequisites()
    {
    	// Check that a velocity model exists.
    	if(model.getProject().velocityModel != null)
    		return true;
    	else
    	{
    		reasonForFailure = "Velocity model must exist in project prior to loading a well plan.";
    		return false;
    	}
    }
    
    private void initialize()
    {
        hUnits = model.getProject().getXyUnits();
        vUnits = model.getProject().getDepthUnits();
    }

    public boolean end()
    {
        if (success)
        {
            model.setActionStatus(getClass().getName(), StsModel.STARTED);
        }
        return super.end();
    }

    public void previous()
    {
    /*
        if(currentStep == relateWell)
        {
            if(!selectPlan.isNew())
            {
                new StsMessage(getModel().win3d, StsMessage.WARNING, "Unable to select another plan, when defining plan. Finish and re-enter wizard.");
            }
            else
                gotoPreviousStep();
        }
        else if(currentStep == definePath)
        {
            new StsMessage(getModel().win3d, StsMessage.WARNING, "Unable to select existing plan, when defining new plan. Finish and re-enter wizard.");
        }
        else
    */
            gotoPreviousStep();
    }

    public void next()
    {
        if(currentStep == wellSelect)
        {
           StsFile[] selectedFiles = wellSelect.panel.getSelectedFiles();
           loadWells.setSelectedFiles(selectedFiles);
           gotoStep(loadWells);
        }
        else if(currentStep == loadWells)
        {
//            finish();
            gotoNextStep();
        }
    }
/*
    public StsWellPlan getWellPlan()
    {
        return (StsWellPlan)planSet.getPlan();
    }
*/
    public String getPathName()
    {
        if(planSet != null)
            return planSet.getName();
        else
            return null;
    }

    public StsColor getPathColor()
    {
        if(planSet != null)
            return planSet.getStsColor();
        else
            return null;
    }

    public StsWellPlanSet getWellPlanSet()
    {
        return planSet;
    }

    public void setWellPlanSet(StsWellPlanSet planSet)
    {
        this.planSet = planSet;
    }

	public StsWellPlan getWellPlan()
	{
		return this.wellPlan;
	}

	public void setWellPlan(StsWellPlan wellPlan)
	{
		this.wellPlan = wellPlan;
    }

    public Object[] getWellPlans()
    {
        if(planSet == null)
            return null;
        else
            return planSet.getWellPlans();
    }

    public void finish()
    {
        super.finish();
        return;
    }
/*
    public void constructAvailableFiles()
    {
        StsWellImport.initializeWellFilenameSets();
        StsWellImport.constructWellFilenameSets(StsWellImport.LASFILES);
        StsWellImport.constructWellFilenameSets(StsWellImport.ASCIIFILES);
       StsWellImport.constructWellFilenameSets(StsWellImport.UTFILES);
        StsWellImport.compressWellSet();
   }
*/
}
