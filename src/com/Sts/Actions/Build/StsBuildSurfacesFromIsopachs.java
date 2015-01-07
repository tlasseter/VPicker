
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

public class StsBuildSurfacesFromIsopachs extends StsAction implements Runnable
{
    private String surfaceName = null;
    private boolean success = false;
    private int nSurfacesBuilt = 0;

 	public StsBuildSurfacesFromIsopachs(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public void run()
    {
        try
        {
            statusArea.setTitle("Build surfaces:");
            logMessage("Choose surface to isopach up/down from.");

            StsObject[] surfaces = model.getObjectListOfType(StsModelSurface.class, StsModelSurface.IMPORTED);

            StsSelectStsObjects selector = StsSelectStsObjects.constructor(model, surfaces, model.getName(),
                                                "Select a surface grid to isopach from:", true);
			if(selector == null) return;
            StsObject surfaceObj = selector.selectObject();
            if (surfaceObj == null)
            {
                statusArea.textOnly();
                logMessage("No surface selected.");
    	        actionManager.endCurrentAction();
                return;
            }
            StsModelSurface surface = (StsModelSurface)surfaceObj;
			surfaceName = surface.getName();
            StsMarker surfaceMarker = surface.getMarker();
            if (surfaceMarker == null)
            {
                statusArea.textOnly();
                logMessage("Surface has no marker correlated with it.");
    	        actionManager.endCurrentAction();
                return;
            }

            // get markers
            logMessage("Choose markers to derive isopach surfaces from.");
            selector = StsSelectStsObjects.constructor(model,
                    StsMarker.class, model.getName(), "Select markers:", false, false);
			if(selector == null) return;
            StsObject[] markers = selector.selectObjects();
            if (markers==null)
            {
                statusArea.textOnly();
                logMessage("No markers selected.");
    	        actionManager.endCurrentAction();
                return;
            }

            // build the surfaces
            statusArea.setMaximum(markers.length);
            statusArea.setProgress(0);
            statusArea.addProgress();
            statusArea.sleep(10);
            StsSurfaceClass surfaceClass = (StsSurfaceClass)model.getCreateStsClass(StsModelSurface.class);
            for (int i=0; i<markers.length; i++)
            {
                StsMarker marker = (StsMarker)markers[i];
                String name = marker.getName();
                if (surfaceClass.getObjectWithName(name) != null || marker == surfaceMarker)
                {
                    logMessage("Already have surface: " + name + ".  Skipping to next marker.");
                    statusArea.setProgress(i+1); // update progress
                    continue;
                }
                StsModelSurface newSurface = StsModelSurface.constructModelSurface(name, marker.getStsColor(), StsModelSurface.IMPORTED);
                newSurface.copyModelSurface(surface);
                newSurface.setPrefixString("Surface: " + name);
                if (!newSurface.setGridPoints(surface, marker, StsParameters.largeFloat, 4.0f))
                {
                    logMessage("Grid creation for surface:  " + name + " failed.");
                    statusArea.setProgress(i+1); // update progress
                    model.delete(newSurface);
                }
//                newSurface.persistPoints();
                newSurface.constructGrid();
                //newSurface.saveValuesToDataStore();
                newSurface.setMarker(marker);
                marker.setModelSurface(newSurface);
                nSurfacesBuilt++;
                logMessage("Built new surface: " + name);
                statusArea.setProgress(i+1); // update progress
            }

            success = true;
        	actionManager.endCurrentAction();
        }
        catch(Exception e)
        {
            statusArea.removeProgress();
            StsException.outputException("StsBuildSurfacesFromIsopachs failed.",
                    e, StsException.WARNING);
        }
    }

    public boolean end()
    {
        statusArea.textOnly();
        if (success)
        {
            logMessage("Built " + nSurfacesBuilt +
                    " new surface(s) from marker isopachs and surface"
                    + surfaceName + ".");
        }
        else
        {
            logMessage("Unable to build surfaces from marker isopachs");
        }
        return success;
    }

}
