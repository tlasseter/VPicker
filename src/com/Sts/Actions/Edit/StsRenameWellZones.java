
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Edit;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.RFUI.*;

public class StsRenameWellZones extends StsAction
{
 	public StsRenameWellZones(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public boolean start()
    {
        try
        {
            StsRenameDialog d = new StsRenameDialog(model,
                    StsWellZoneSet.class, model.getName(), "Rename Well Zones:", true);
            d.setVisible(true);
        }
        catch(Exception e)
        {
            System.out.println("Exception in StsRenameWellZones.start()\n" + e);
            return false;
        }

	    actionManager.endCurrentAction();
        return true;
    }
}
