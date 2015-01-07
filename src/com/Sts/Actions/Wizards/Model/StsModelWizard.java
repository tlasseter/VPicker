
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Model;


import com.Sts.Actions.Build.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;


public class StsModelWizard extends StsWizard
{
	private StsSelectZones selectZones = new StsSelectZones(this);
	private StsEditZones editZones = new StsEditZones(this);
	private StsEditFaultGaps editFaultGaps = new StsEditFaultGaps(this);
	private StsSurfaceAutopick surfaceAutopick = new StsSurfaceAutopick(this);
	private StsBuild3dGrid build3dGrid = new StsBuild3dGrid(this);
    private StsWizardStep[] mySteps =
    {
    	selectZones,
        editZones,
        editFaultGaps,
		surfaceAutopick,
        build3dGrid
    };

    StsBuiltModel builtModel = null;
    private boolean skipEditFaultGaps = false;
	private boolean skipSurfaceAutopick = false;

	public StsModelWizard(StsActionManager actionManager)
    {
        super(actionManager);

        // see if we have fault sections
        try
        {
            int nFaultSections = model.getNObjectsOfType(StsSection.class, StsSection.FAULT);
            if (nFaultSections == 0) skipEditFaultGaps = true;

        }
        catch (NullPointerException e) { skipEditFaultGaps = true; }

        addSteps(mySteps);
    }

	public boolean checkStartAction()
	{
        StsBuiltModel builtModel = (StsBuiltModel)model.getCurrentObject(StsBuiltModel.class);
        if(builtModel == null) return true;
		boolean yes = StsYesNoDialog.questionValue(model.win3d,"Model already built.\n  Do you wish to delete it and rebuild?.");
		if(yes)
		{
			builtModel.delete();
            StsEclipseModel eclipseModel = (StsEclipseModel)model.getCurrentObject(StsEclipseModel.class);
            if(eclipseModel != null) eclipseModel.deleteModel(model);
        }
		return yes;
    }

    public boolean start()
    {
		// model.getDatabase().blockTransactions();
        model.win3d.removeToolbar(StsBuildFrameToolbar.NAME);
        model.setActionStatus(StsBuildFrame.class.getName(), StsModel.ENDED);
        System.runFinalization();
        System.gc();
    	dialog.setTitle("Build Model");
    	return super.start();
    }

    public void previous()
    {
		if(currentStep == build3dGrid)
		{
			if(skipSurfaceAutopick)
			{
				if(skipEditFaultGaps)
					gotoStep(editZones);
				else
					gotoStep(editFaultGaps);
			}
			else
				gotoStep(surfaceAutopick);
		}
		else if(currentStep == surfaceAutopick)
		{
			if(skipEditFaultGaps)
				gotoStep(editZones);
			else
				gotoStep(editFaultGaps);
		}
		else
			gotoPreviousStep();
    }

    public void next()
    {
        if(currentStep == selectZones)
        {
            StsZone[] zones = selectZones.getSelectedZones();
            if(zones != null)
            {
                model.deleteStsClass(StsBuiltModel.class);
                if(!StsBlock.constructBlocks(model)) return;
                StsBlock[] blocks = (StsBlock[])model.getCastObjectList(StsBlock.class);
                builtModel = new StsBuiltModel(zones, blocks);
                gotoStep(editZones);
            }
        }
        else if(currentStep == editZones)
        {
            if(skipEditFaultGaps)
            {
				if(skipSurfaceAutopick)
				{
					gotoStep(build3dGrid);
					enableFinish();
				}
				else
					gotoStep(surfaceAutopick);
            }
            else
                gotoStep(editFaultGaps);
        }
        else if(currentStep == editFaultGaps)
        {
			if(skipSurfaceAutopick)
			{
				gotoStep(build3dGrid);
				enableFinish();
			}
			else
				gotoStep(surfaceAutopick);

        }
		else if(currentStep == surfaceAutopick)
		{
			gotoStep(build3dGrid);
			enableFinish();
		}
    }

    public void finish()
    {
        actionManager.endCurrentAction();
    }

    public boolean end()
    {
        success = super.end();
        if(success)
        {
            model.setActionStatus(StsModelWizard.class.getName(), StsModel.STARTED);
            model.commit();
        }
        return success;
    }
    
    public StsZone[] getSelectedZones() { return selectZones.getSelectedZones(); }
    public StsBuiltModel getBuiltModel() { return builtModel; }
}
