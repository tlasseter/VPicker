
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

public class StsBuildZonesFromSurfaces extends StsAction
{
    StsModelSurface[] surfaces;
    String[] surfaceNames;

 	public StsBuildZonesFromSurfaces(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public boolean start()
    {
        try
        {
			surfaces = (StsModelSurface[])model.getCastObjectList(StsModelSurface.class);
            int nSurfaces = surfaces.length;
            if (nSurfaces < 2)
            {
                logMessage("At least two surfaces are required to define a zone.");
    	        actionManager.endCurrentAction();
                return true;
            }
            surfaceNames = StsMainObject.getNamesFromObjects(surfaces);
            while (true)
            {
    		    StsZoneBuilderDialog dialog = new StsZoneBuilderDialog(model.win3d, surfaceNames);
                dialog.setLocationRelativeTo(model.win3d);
                dialog.setVisible(true);
                surfaceNames = dialog.getSelectedItems();
                if (surfaceNames==null)  // user cancelled
                {
                    StsMessageFiles.logMessage("Surface selection cancelled.");
			        actionManager.endCurrentAction();
                    return true;
                }
                if (surfaceNames.length==1)
                {
                    StsMessageFiles.logMessage("Only 1 surface chosen.  Please reselect.");
                    dialog.dispose();
                    continue;
                }
                break;
            }
            surfaceNames = StsOrderedListDialog.parseNames(surfaceNames);
            int nZones = buildZones();
            if (nZones<1)
            {
                StsMessageFiles.logMessage("Error! Unable to define zones.");
    	        actionManager.endCurrentAction();
                return false;
            }

            //StsZone.buildModelSurfaces(model);

            if (nZones==1)
            {
                StsMessageFiles.logMessage("Successfully defined 1 zone.");
            }
            else
            {
                StsMessageFiles.logMessage("Successfully defined " + nZones
                        + " zones.");
            }
			actionManager.endCurrentAction();
        }
        catch(Exception e)
        {
            StsException.outputException("StsBuildZonesFromSurfaces failed.", e, StsException.WARNING);
            return false;
        }
        return true;
    }

    /** build zones from surfaces (bottom to top)  top to bottom??  */
    private int buildZones() throws StsException
    {
        if (surfaceNames==null || surfaceNames.length < 2) return 0;
        int nZonesBuilt = 0;

        StsModelSurface topModelSurface = null;
        StsModelSurface baseModelSurface = null;
        baseModelSurface = getSurfaceWithName(surfaceNames[0]);
        StsZoneClass zoneClass = (StsZoneClass)model.getCreateStsClass(StsZone.class);
        for (int i=1; i<surfaces.length; i++)
        {
            topModelSurface = baseModelSurface;
            baseModelSurface = getSurfaceWithName(surfaceNames[i]);
            if(baseModelSurface == null) return 0;
            if (zoneClass.getExistingZone(topModelSurface, baseModelSurface) != null) continue;
            StsZone zone = new StsZone(topModelSurface.getName(), topModelSurface, baseModelSurface);
            StsMarker.tryToLinkIntervalBelow(topModelSurface.getMarker(), topModelSurface);
            nZonesBuilt++;

            zone.setTopModelSurface(topModelSurface);
            topModelSurface.setZoneBelow(zone);
            zone.setBaseModelSurface(baseModelSurface);
            baseModelSurface.setZoneAbove(zone);
        }
        return nZonesBuilt;
    }

    private StsModelSurface getSurfaceWithName(String name)
    {
        return (StsModelSurface)StsMainObject.getListObjectWithName(surfaces, name);
    }
}
