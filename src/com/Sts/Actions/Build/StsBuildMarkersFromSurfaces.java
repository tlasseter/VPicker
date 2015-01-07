
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Build;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.RFUI.*;
import com.Sts.Utilities.*;

public class StsBuildMarkersFromSurfaces extends StsAction implements Runnable
{
    private int nMarkersBuilt = 0;
    private boolean success = false;

 	public StsBuildMarkersFromSurfaces(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public void run()
    {
        try
        {
            // get wells
            StsClass wells = model.getCreateStsClass(StsWell.class);
            int nWells = wells.getSize();

            if (nWells == 0)
            {
                logMessage("No wells found to intersect with.");
    	        actionManager.endCurrentAction();
                return;
            }

            // get surfaces
            StsSelectStsObjects selector = StsSelectStsObjects.constructor(model,
                    StsModelSurface.class, model.getName(), "Select surfaces:", false, false);
			if(selector == null) return;

            StsObject[] surfaces = selector.selectObjects();
            if (surfaces==null)
            {
                logMessage("No surfaces selected to intersect with.");
    	        actionManager.endCurrentAction();
                return;
            }

            // build the markers
            statusArea.setMaximum(surfaces.length);
            statusArea.setProgress(0);
            statusArea.addProgress();
            statusArea.sleep(10);
            for (int i=0; i<surfaces.length; i++)
            {
                StsModelSurface s = (StsModelSurface)surfaces[i];
                for (int n = 0; n < nWells; n++)
                {
                    StsWell well = (StsWell)wells.getElement(n);
                    if (doIntersection(s, well)) nMarkersBuilt++;
                }
                statusArea.setProgress(i+1); // update progress
            }

            success = true;
			actionManager.endCurrentAction();
        }
        catch(Exception e)
        {
            statusArea.removeProgress();
            StsException.outputException("StsBuildMarkersFromSurfaces failed.",
                    e, StsException.WARNING);
        }
    }

    public boolean end()
    {
        if (success)
        {
            logMessage("Built " + nMarkersBuilt +
                    " new well markers by surface/well intersections.");
        }
        else
        {
            logMessage("Unable to build markers from surfaces.");
        }
        return success;
    }

    /** intersect a well & surface */
    public boolean doIntersection(StsModelSurface modelSurface, StsWell well)
    {
        if (modelSurface == null || well == null) return false;
        StsWellMarker wellMarker = well.getMarker(modelSurface.getName());
        if (wellMarker== null)
            wellMarker = StsWellMarker.constructor(modelSurface, well); // adds itself to instance list
        if(wellMarker == null) return false;
//        wellMarker.setSurfaceLocation();
        return true;
    }
}
