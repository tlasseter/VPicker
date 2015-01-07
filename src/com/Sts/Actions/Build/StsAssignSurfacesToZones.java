
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

public class StsAssignSurfacesToZones extends StsAction
{
 	public StsAssignSurfacesToZones(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public boolean start()
    {
        try
        {
            // retrieve zones
            StsClass zones = model.getCreateStsClass(StsZone.class);
            if (zones==null || zones.getSize()==0)
            {
                logMessage("No zones found.");
                return true;
            }
            StsObject[] zoneObjects = zones.getElements();

            // retrieve surfaces
            StsSurfaceClass surfaceClass = (StsSurfaceClass)model.getCreateStsClass(StsModelSurface.class);
            StsObject[] surfaceObjects = surfaceClass.getObjectList();
            if (surfaceObjects.length == 0)
            {
                logMessage("No surfaces found.");
                return true;
            }
            // reorder the surfaces
            StsSurfaceAssignmentDialog dialog = new StsSurfaceAssignmentDialog(model.win3d,
                    surfaceObjects, zoneObjects);
            dialog.setVisible(true);
            int[] selectedIndices = dialog.getSelectedIndices();
            if (selectedIndices == null)
            {
                logMessage("No surfaces found to assign to zones.");
                return true;
            }

            int nIndices = selectedIndices.length;
            if (nIndices != zoneObjects.length+1)
            {
                logMessage("Invalid number of surfaces found to assign to zones.");
                return false;
            }

            // add surfaces to zones
            StsModelSurface[] zoneSurfaces = new StsModelSurface[nIndices];
            for (int i=0; i<nIndices; i++)
                zoneSurfaces[i] = (StsModelSurface)surfaceObjects[selectedIndices[i]];

            for (int i=0; i<zoneObjects.length; i++)
            {
                StsZone zone = (StsZone)zoneObjects[i];
                StsModelSurface top = zoneSurfaces[i];
                if (top != null)
                {
                    zone.setTopSurface(top);
                    top.setStsColor(zone.getStsColor());
                }
                StsModelSurface base = zoneSurfaces[i+1];
                if (base != null) zone.setBotSurface(base);
            }
        }
        catch(Exception e)
        {
            System.out.println("Exception in StsAssignSurfacesToZones.start()\n" + e);
            return false;
        }
        return true;
    }
}
