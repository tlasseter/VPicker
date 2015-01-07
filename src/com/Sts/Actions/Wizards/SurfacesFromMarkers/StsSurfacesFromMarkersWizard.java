
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SurfacesFromMarkers;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsSurfacesFromMarkersWizard extends StsWizard
{
	private StsSelectMarkers selectMarkers;
    private StsBuildSurfaces buildSurfaces;

	public StsRotatedGridBoundingBox boundingBox;


	public StsSurfacesFromMarkersWizard(StsActionManager actionManager)
    {
        super(actionManager);
        setupBoundingBox();
        selectMarkers = new StsSelectMarkers(this);
		buildSurfaces = new StsBuildSurfaces(this);
		addSteps(new StsWizardStep[] { selectMarkers, buildSurfaces } );
        try
        {
            if (model.getNObjects(StsMarker.class) < 1)
            {
                new StsMessage(frame, StsMessage.ERROR, "No wells with markers exist.");
                return;
            }
        }
        catch (NullPointerException e)
        {
            StsException.outputException("Problem accessing markers in project", e, StsException.WARNING);
            return;
        }
    }

    private void setupBoundingBox()
    {
        boundingBox = model.getProject().getRotatedBoundingBox().getGridBoundingBoxClone();
        boundingBox.adjustBoundingBox(200);      
    }

    public boolean start()
    {
        byte zDomain = model.getProject().getZDomain();
        boolean ok;
        if(zDomain == StsProject.TD_DEPTH)
            ok = StsYesNoDialog.questionValue(model.win3d, "Surfaces will be built from markers in depth. Is this ok?.  If not switch to time.");
        else
            ok = StsYesNoDialog.questionValue(model.win3d, "Surfaces will be built from markers in time. Is this ok?.  If not switch to depth.");
        if(!ok) return false;

        dialog.setTitle("Build Surfaces from Markers");
    	return super.start();
    }

    public void next()
    {
        gotoNextStep();
    }

    public void previous()
    {
        gotoPreviousStep();
    }

    public boolean end()
    {
		if(!super.end()) return false;
        model.setActionStatus(getClass().getName(), StsModel.STARTED);
		return true;
    }

    public StsMarker[] getSelectedMarkers() { return selectMarkers.getSelectedMarkers(); }

	public double getXOrigin() { return boundingBox.xOrigin; }
	public double getYOrigin() { return boundingBox.yOrigin; }
	public float getXInc() { return boundingBox.xInc; }
	public float getYInc() { return boundingBox.yInc; }
	public int getNRows() { return boundingBox.nRows; }
	public int getNCols() { return boundingBox.nCols; }

	public void setXOrigin(double value) { boundingBox.xOrigin = value; }
	public void setYOrigin(double value) { boundingBox.yOrigin = value; }
	public void setXInc(float value) { boundingBox.xInc = value; }
	public void setYInc(float value) { boundingBox.yInc = value; }
	public void setNRows(int value) { boundingBox.nRows = value; }
	public void setNCols(int value) { boundingBox.nCols = value; }
}

