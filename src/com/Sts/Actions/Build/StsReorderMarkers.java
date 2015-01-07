
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

public class StsReorderMarkers extends StsAction
{
 	public StsReorderMarkers(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public boolean start()
    {
        try
        {
            statusArea.setTitle("Reorder markers:");
            logMessage("Move markers into stratigraphic order, then press Okay.");
            StsClass markers = model.getCreateStsClass(StsMarker.class);
            if (markers == null)
            {
                statusArea.textOnly();
                logMessage("No markers found to reorder.");
                return true;
            }
            else if (markers.getSize() < 2)
            {
                statusArea.textOnly();
                logMessage("At least two markers are required to do reordering.");
                return true;
            }

            // build the dialog
			StsOrderedMarkersDialog dialog = new StsOrderedMarkersDialog(model.win3d);
            dialog.setLocationRelativeTo(model.win3d);
            dialog.setVisible(true);
			actionManager.endCurrentAction();
        }
        catch(Exception e)
        {
            System.out.println("Exception in StsReorderMarker.start()\n" + e);
            return false;
        }
        statusArea.textOnly();
        logMessage("Done reordering markers.");
        return true;
    }
}
