package com.Sts.Actions.Wizards.WellPlan;

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.Well.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

public class StsLoadUtWellPlanWizard extends StsWizard
{
    byte vUnits = StsParameters.DIST_NONE;
    byte hUnits = StsParameters.DIST_NONE;
    private StsUtWellSelect wellSelect;
    private StsDefineUtWells defineUtWells;
    private StsUtWellPlanLoad loadWells;

    private StsWellPlanSet planSet = null;
    private StsWellPlan wellPlan = null;

    private StsWizardStep[] mySteps =
    {
 //       selectPlan = new StsSelectPlan(this),
 //       initialPlan = new StsInitialPlan(this),
        wellSelect = new StsUtWellSelect(this),
        defineUtWells = new StsDefineUtWells(this),
        loadWells = new StsUtWellPlanLoad(this),
//        definePath = new StsDefinePath(this),
//        relateWell = new StsWellRelate(this),
//        defineVertices = new StsDefineVertices(this)
   };

    public StsLoadUtWellPlanWizard(StsActionManager actionManager)
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

    private void initialize()
    {
        StsWellImport.initialize(model);
        loadWells.setWellFactory(new StsWellPlanFactory());
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
    /*
        if(currentStep == selectPlan)
        {
            if(selectPlan.isNew())
            {
                planSet = new StsWellPlanSet("WellPlanSet");
//                setWellPlan((StsWellPlan)planSet.getPlan());
                gotoStep(initialPlan);
            }
            else
            {
                enableFinish();
                gotoStep(relateWell);
            }

        }
        else if(currentStep == initialPlan)
        {
            if(initialPlan.isImportPlan())
            {
                gotoStep(wellSelect);
            }
            else
            {
                gotoStep(definePath);
            }
        }
        else if(currentStep == wellSelect)
    */
        if(currentStep == wellSelect)
        {
           String[] selectedWellnames = wellSelect.getSelectedWellnames();
           String[] selectedWellFilenames = StsWellImport.getSelectedWellFilenames(selectedWellnames, StsWellImport.UTFILES, true);
           loadWells.setSelectedWellnames(selectedWellnames);
           if(selectedWellFilenames != null)
           {
               defineUtWells.setSelectedWellnames(selectedWellnames);
               gotoStep(defineUtWells);
           }
           else
               gotoStep(loadWells);
        }
        else if(currentStep == defineUtWells)
        {
            StsPoint[] topPoints = defineUtWells.getTopHoleLocations();
            byte[] types = defineUtWells.getTypes();
            StsWellImport.setTopHoleLocations(topPoints, types);
 //           enableFinish();
            gotoStep(loadWells);
        }
        else if(currentStep == loadWells)
        {
//            finish();
            gotoNextStep();
        }
    /*
        else if(currentStep == definePath)
        {
            gotoStep(relateWell);
        }
        else if(currentStep == relateWell)
        {
//            enableFinish();
             gotoStep(defineVertices);
        }
        else if(currentStep == defineVertices)
         {
             enableFinish();
             gotoNextStep();
        }

        else
            gotoNextStep();
     */
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
    public void constructAvailableFiles()
    {
        StsWellImport.initializeWellFilenameSets();
        StsWellImport.constructWellFilenameSets(StsWellImport.LASFILES);
        StsWellImport.constructWellFilenameSets(StsWellImport.ASCIIFILES);
        StsWellImport.constructWellFilenameSets(StsWellImport.UTFILES);
        StsWellImport.compressWellSet();
   }
}
