
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Horizons;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.*;

public class StsHorizonsWizard extends StsWizard
{
	private static final int SELECT_SURFACES = 0;
	private static final int CORRELATE_MARKERS = 1;
	private static final int DEFINE_GRID = 2;
	private static final int BUILD_GRIDS = 3;
//	private static final int BUILD_BOUNDARY = 4;

	private StsSelectSurfaces selectSurfaces;
	private StsDefineGrid defineGrid;
	private StsBuildHorizonGrids buildGrids;
	private StsCorrelateMarkers correlateMarkers;
//	private StsBuildBoundary buildBoundary = new StsBuildBoundary(this);
    private StsWizardStep[] mySteps =
    {
    	selectSurfaces = new StsSelectSurfaces(this),
        defineGrid = new StsDefineGrid(this),
        buildGrids = new StsBuildHorizonGrids(this),
		correlateMarkers = new StsCorrelateMarkers(this)
    };
    private boolean skipCorrelateMarkers = false;
	private StsSurface[] surfaces;
	private StsModelSurface[] horizons;
	public StsHorizonsWizard(StsActionManager actionManager)
    {
        super(actionManager, 500, 400);
        defineGrid.setModel(model);
        addSteps(mySteps);

        // see if we have well markers
        try
        {
            if (model.getNObjects(StsMarker.class) < 1)
                skipCorrelateMarkers = true;
        }
        catch (NullPointerException e) { skipCorrelateMarkers = true; }

        // check that we have a zInc defined for project

        model.getProject().checkCreateZInc();
    }
/*
    public boolean checkStartAction()
    {
        StsMessage.questionValue(model.win3d,"NOT AVAILABLE YET.\n\nNeed to delete model, faults and boundary to define new horizons.\n\n       Proceed?\n");
        return true;
//        return false;
    }
 */
    public boolean start()
    {
    	dialog.setTitle("Build Horizons");
        surfaces = (StsSurface[])model.getCastObjectList(StsSurface.class);
        if(surfaces == null || surfaces.length == 0)
        {

            new StsMessage(model.win3d, StsMessage.WARNING, "No imported surfaces available.");
            cancel();
            return false;
        }

    	return super.start();
    }

	public boolean checkPrerequisites()
	{
		StsSurfaceClass surfaceClass = (StsSurfaceClass) model.getStsClass("com.Sts.DBTypes.StsSurface");
		StsObject[] surfaces = surfaceClass.getSurfaces();
		if (surfaces.length == 0) {
			reasonForFailure = "No surfaces available for conversion to horizons. Must have loaded or created at least one.";
			return false;
		}
        return true;
    }

    public void next()
    {
        if(isTest) gotoNextStep();

        if(currentStep == selectSurfaces)
        {
		    surfaces = selectSurfaces.getSelectedSurfaces();
            if(surfaces != null)
                gotoNextStep();
        }
        else if(currentStep == defineGrid)
        {
            StsGridDefinition gridDef = defineGrid.getGridDefinition();
            buildGrids.setGridDef(gridDef);
            model.getProject().setGridDefinition(gridDef);
        	gotoNextStep();
        }
        else if(currentStep == buildGrids)
		{
			if(skipCorrelateMarkers) finish();
			else
			{
				StsModelSurface[] horizons = buildGrids.getHorizons();
				correlateMarkers.setModelSurfaces(horizons);
				enableFinish();
				gotoNextStep();
			}
		}
		else if(currentStep == correlateMarkers)
			finish();
    }

    public void previous()
    {
        gotoPreviousStep();
    }
/*
    public void finish()
    {
    	next();
		super.finish();
//		if(currentStep == buildGrids) super.finish();
    }
*/
    public boolean end()
    {
		if(!super.end()) return false;
        model.setActionStatus(getClass().getName(), StsModel.STARTED);
		return true;
    }

	public StsSurface[] getSurfaces() { return surfaces; }
}

