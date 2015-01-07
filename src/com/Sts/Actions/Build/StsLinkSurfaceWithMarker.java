
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

public class StsLinkSurfaceWithMarker extends StsAction
{
 	public StsLinkSurfaceWithMarker(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public boolean start()
    {
        try
        {
            statusArea.setTitle("Link surfaces & markers:");

            // check special cases
			StsModelSurface[] surfaces = (StsModelSurface[])model.getCastObjectList(StsModelSurface.class);
            int nSurfaces = surfaces.length;
            if (nSurfaces < 1)
            {
                statusArea.textOnly();
                logMessage("No surface found.");
              	actionManager.endCurrentAction();
                return true;
            }
            int nMarkers = model.getNObjects(StsMarker.class);
            if (nMarkers < 1)
            {
                statusArea.textOnly();
                logMessage("No markers found.");
              	actionManager.endCurrentAction();
                return true;
            }

            String[] surfaceNames = StsMainObject.getNamesFromObjects(surfaces);
            while (true)  // keep asking until user cancels
            {
                // retrieve surfaces and select one
                StsSurfaceLinkDialog dialog = new StsSurfaceLinkDialog(model.win3d, surfaceNames);
                dialog.setLocationRelativeTo(model.win3d);
                dialog.setVisible(true);
                StsModelSurface surface = dialog.getSelectedSurface();
                if (surface==null)
                {
                    statusArea.textOnly();
                    logMessage("Surface - marker linking finished.");
                   	actionManager.endCurrentAction();
                    return true;  // user cancelled
                }

                logMessage("Selected surface:  " + surface.getName());
                StsMarker currentMarker = surface.getMarker();

                // retrieve marker and select one
                StsMarkerLinkDialog dialog2 = new StsMarkerLinkDialog(model.win3d,
                        "Surface: " + surface.getName());
                if (currentMarker != null)
                {
                    dialog2.setSelectedItem(currentMarker.getName());
                }
                dialog2.setLocationRelativeTo(model.win3d);
                dialog2.setVisible(true);
                StsMarker newMarker = dialog2.getSelectedMarker();
                if (newMarker==null) continue;

                logMessage("Selected marker:  " + newMarker.getName());
                //status.sleep(1000);
                if (newMarker != currentMarker)
                {
                    surface.setMarker(newMarker);
                    newMarker.setModelSurface(surface);
                    logMessage("Linked surface: " + surface.getName() +
                            " with marker:  " + newMarker.getName());
                    StsMarker.tryToBuildNewIntervals(model, newMarker, surface);
                }
                else
                {
                    logMessage("Previously linked surface: " + surface.getName() +
                            " with marker:  " + newMarker.getName());
                }
                if (nSurfaces==1 || nMarkers==1)  // user can't cancel!
                {
                   	actionManager.endCurrentAction();
                    return true;
                }
            }
        }
        catch(Exception e)
        {
            System.out.println("Exception in StsLinkSurfaceWithMarker.start()\n" + e);
			actionManager.endCurrentAction();
            return false;
        }
    }
}
