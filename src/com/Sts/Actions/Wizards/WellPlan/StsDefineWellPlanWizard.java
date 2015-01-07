package com.Sts.Actions.Wizards.WellPlan;

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

public class StsDefineWellPlanWizard extends StsWizard
{
    byte vUnits = StsParameters.DIST_NONE;
    byte hUnits = StsParameters.DIST_NONE;
    private StsPlatformOrWell platformOrWell;
    private StsDefinePlatform definePlatform;
    private StsSelectPlatformSlot selectPlatformSlot;
    private StsDefineKbAndKickoff defineKBandKickoff;
    private StsDefineTarget defineTarget;
    private StsAddMidPoint addMidPoint;
	private StsDefineActualWellStart defineActualWellStart;

    private String planName = null;
    private boolean isCanceled = false;
    private StsWellPlanSet wellPlanSet = null;
	private boolean isNewPlanSet;
	private StsWellPlan wellPlan = null;
    private StsWell drillingWell = null;
//    private String wellType = StsPlatformOrWellPanel.DEFINE_PLATFORM;
    private StsPlatform platform = null;

    private StsWizardStep[] mySteps =
    {
        platformOrWell = new StsPlatformOrWell(this),
        definePlatform = new StsDefinePlatform(this),
        selectPlatformSlot = new StsSelectPlatformSlot(this),
        defineKBandKickoff = new StsDefineKbAndKickoff(this),
		defineActualWellStart = new StsDefineActualWellStart(this),
        defineTarget = new StsDefineTarget(this),
        addMidPoint = new StsAddMidPoint(this)
   };

   StsPlatformOrWellPanel panel = platformOrWell.panel;

    public StsDefineWellPlanWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();

        dialog.setTitle("Define Well Plan");
        initialize();
        if(!super.start()) return false;

        return true;
    }

    private void initialize()
    {
        StsWellImport.initialize(model);
    }

    public void next()
    {
        model.getProject().setProjectTimeToCurrentTime();
        if(currentStep == platformOrWell)
        {
			wellPlanSet = panel.getWellPlanSet();
            wellPlanSet.setDrawCurtain(false);
			isNewPlanSet = panel.getIsNewPlan();
			if(isNewPlanSet)
			{
				wellPlan = wellPlanSet.getCurrentWellPlan();
				platform = panel.getCreatePlatform();
				wellPlan.setPlatform(platform);
			}
			else
			{
				wellPlan = wellPlanSet.copyLastWellPlan();
			}
			wellPlan.initializePlan();
            wellPlan.addTimeStampName();
			drillingWell = panel.getDrillingWell();
			if(drillingWell != null)
			{
				wellPlan.setDrillingWell(drillingWell);
				gotoStep(defineActualWellStart);
			}
			else
			{
				if (platform != null)
				{
					if (panel.getIsNewPlatform())
						gotoStep(definePlatform);
					else
						gotoStep(selectPlatformSlot);
				}
				else
					gotoStep(defineKBandKickoff);
			}
        }
		else if(currentStep == defineActualWellStart && !isNewPlanSet)
		     finish();
		else if(currentStep == defineKBandKickoff)
		{
			gotoStep(defineTarget);
		}
        else if(currentStep == defineTarget && wellPlan.getTrajectoryType() == StsWellPlan.STRAIGHT_STRING)
            finish();
        else
            gotoNextStep();
        model.getProject().setProjectTimeToCurrentTime();
    }

	public boolean getIsNewPlanSet() { return isNewPlanSet; }
/*
    public void cancel()
    {
        platform.deleteWellfromPlatform(wellPlan.getName());
        wellPlan.delete();
        if(platform != null && platform.getNWells() == 0)
        {
            platform.delete();
            platform = null;
        }
        wellPlan = null;
        isCanceled = true;
        super.cancel();
        model.refreshObjectPanel();
    }
*/
    public boolean end()
    {
        if (success)
        {
            model.setActionStatus(getClass().getName(), StsModel.STARTED);
            wellPlan.wizardFinish();
			// the first wellPlan constructed for this set is just made persistent above, so it can be added to the persistent wellPlans object
			// Otherwise it has already been added when it was copied
            if(!wellPlan.isPersistent())  wellPlan.addToModel();
            if(isNewPlanSet)
				wellPlanSet.addWellPlan(wellPlan);
			else
				wellPlanSet.setPlan(wellPlan);

            wellPlanSet.setDrawCurtain(true);
            model.refreshObjectPanel();
        }
        return super.end();
    }

    public void previous()
    {
        currentStep.cancel();

		if(currentStep == defineTarget)
		{
			if (drillingWell != null)
				gotoStep(defineActualWellStart);
			else
				gotoStep(defineKBandKickoff);
		}
        else if(currentStep == selectPlatformSlot)
        {
			if (panel.getWellType() == StsPlatformOrWellPanel.PLATFORM_WELL)
			{
				gotoStep(platformOrWell);
			}
			else
				gotoPreviousStep();
		}
		else if(currentStep == defineKBandKickoff)
		{
			if(platform == null)
				gotoStep(platformOrWell);
			else
				gotoStep(selectPlatformSlot);
		}
    }

    public void setPlatform(Object StsPlatform) { this.platform = platform; }
    public StsPlatform getPlatform() { return platform; }


    public String getPathName()
    {
        if(wellPlanSet != null)
            return wellPlanSet.getName();
        else
            return null;
    }

    public boolean isCanceled()
    {
        return isCanceled;
    }

    public StsColor getPathColor()
    {
        if(wellPlanSet != null)
            return wellPlanSet.getStsColor();
        else
            return null;
    }

    public StsWellPlanSet getWellPlanSet()
    {
        return wellPlanSet;
    }

    public void setWellPlanSet(StsWellPlanSet planSet)
    {
        this.wellPlanSet = planSet;
    }

	public StsWellPlan getWellPlan()
	{
		return wellPlan;
	}

	public void setWellPlan(StsWellPlan wellPlan)
	{
		this.wellPlan = wellPlan;
    }

    public Object[] getWellPlans()
    {
        if(wellPlanSet == null)
            return null;
        else
            return wellPlanSet.getWellPlans();
    }

    public boolean addWellToPlatform()
    {
        return platform.addWellToCurrentSlot(wellPlanSet.getName());
    }
    public void setWellXY(double[] xy)
    {
        wellPlan.setXOrigin(xy[0]);
        wellPlan.setYOrigin(xy[1]);
    }

	public StsWell getDrillingWell() { return drillingWell; }

    public void setZKB(float z)
    {
        wellPlan.setZKB(z);
    }

    public void finish()
    {
        success = true; // TODO probably should do some checking to make sure wellPlan construction is successful
        super.finish();
        return;
    }
    public void constructAvailableFiles()
    {
        StsWellImport.initializeWellFilenameSets();
        StsWellImport.constructWellFilenameSets(StsWellImport.UTFILES);
        StsWellImport.compressWellSet();
   }

    static void main(String[] args)
    {
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
        StsDefineWellPlanWizard drillWizard = new StsDefineWellPlanWizard(actionManager);
        drillWizard.start();
    }
}
