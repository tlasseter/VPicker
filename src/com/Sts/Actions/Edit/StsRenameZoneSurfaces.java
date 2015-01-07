
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

public class StsRenameZoneSurfaces extends StsAction
{
 	public StsRenameZoneSurfaces(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public boolean start()
    {
        try
        {
            StsObject[] surfaces = model.getObjectListOfType(StsModelSurface.class, StsModelSurface.MODEL);
            StsRenameDialog d = new StsRenameDialog(model, surfaces,
                    model.getName(), "Rename Zone Top/Base Surfaces:", false);
            d.setVisible(true);
        }
        catch(Exception e)
        {
            System.out.println("Exception in StsRenameZoneSurfaces.start()\n" + e);
            return false;
        }

	    actionManager.endCurrentAction();
        return true;
    }
}
