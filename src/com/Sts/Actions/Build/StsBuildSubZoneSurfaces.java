
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
import com.Sts.Utilities.*;

public class StsBuildSubZoneSurfaces extends StsAction implements Runnable
{
    private boolean success = false;
    private int nSurfacesBuilt = 0;

 	public StsBuildSubZoneSurfaces(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public void run()
    {
        try
        {
            StsClass zones = model.getCreateStsClass(StsZone.class);
            int nZones = zones.getSize();

            if(nZones <= 0)
            {
                logMessage("No zones have been built.");
    	        actionManager.endCurrentAction();
            }
        /*
            zones.forEach("ConstructSubZoneSurfaces", null);

            status.setMaximum(nZones);
            StsStatusArea.staticSetProgress(0);
            status.addProgress();
            status.sleep(10);
            for (int i=0; i<nZones; i++)
            {

                StsSurface newSurface = new StsSurface(model, name,
                        marker.getStsColor(), StsSurface.IMPORTED);
                model.refreshModelState();
                StsGrid newGrid = surfaceGrid.cloneParameters(newSurface);
                newSurface.setGrid(newGrid);
                newGrid.setPrefixString("Surface: " + name);
                if (!newGrid.setGridPoints(surface, marker,
                        StsParameters.largeFloat, 4.0f))
                {
                    logMessage("Grid creation for surface:  " +
                            name + " failed.");
                    StsStatusArea.staticSetProgress(i+1); // update progress
                    model.delete(newSurface);
                }
                newGrid.persistPoints();
                newGrid.constructGrid();
                //newGrid.saveValuesToDataStore();
                newSurface.setMarker(marker);
                marker.setParentSurface(newSurface);
                nSurfacesBuilt++;
                logMessage("Built new surface: " + name);
                StsStatusArea.staticSetProgress(i+1); // update progress
            }
        */
            success = true;
        	actionManager.endCurrentAction();
        }
        catch(Exception e)
        {
            statusArea.removeProgress();
            StsException.outputException("StsBuildSubZoneSurfaces failed.",
                    e, StsException.WARNING);
        }
    }

    public boolean end()
    {
        statusArea.textOnly();
        if (success)
        {
            logMessage("Built " + nSurfacesBuilt +
                    " new subZone surface(s).");
        }
        else
        {
            logMessage("Unable to build subZone surfaces.");
        }
        return success;
    }

}
