
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

public class StsReorderSurfaces extends StsAction
{
 	public StsReorderSurfaces(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public boolean start()
    {
        try
        {
            statusArea.setTitle("Reorder surfaces:");
            logMessage("Move surfaces into stratigraphic order, then press Okay.");

			StsModelSurface[] surfaces = (StsModelSurface[])model.getCastObjectList(StsModelSurface.class);
            int nSurfaces = surfaces.length;
            if (nSurfaces == 0)
            {
                statusArea.textOnly();
                logMessage("No surfaces found to reorder.");
                return true;
            }
            else if (nSurfaces < 2)
            {
                statusArea.textOnly();
                logMessage("At least two surfaces are required to do reordering.");
                return true;
            }

            // build the dialog
            String[] surfaceNames = StsMainObject.getNamesFromObjects(surfaces);
			StsOrderedSurfacesDialog dialog = new StsOrderedSurfacesDialog(model.win3d, surfaceNames);
            dialog.setLocationRelativeTo(model.win3d);
            dialog.setVisible(true);
			actionManager.endCurrentAction();
        }
        catch(Exception e)
        {
            System.out.println("Exception in StsReorderSurface.start()\n" + e);
            return false;
        }
        statusArea.textOnly();
        logMessage("Done reordering surfaces.");
        return true;
    }
}
