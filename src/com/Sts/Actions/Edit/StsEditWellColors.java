
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

public class StsEditWellColors extends StsAction
{
 	public StsEditWellColors(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public boolean start()
    {
        try
        {
            StsColorSelectDialog d = new StsColorSelectDialog(model,
                    StsWell.class, model.getName(), "Change Well Colors:", true);
            d.setVisible(true);
        }
        catch(Exception e)
        {
            System.out.println("Exception in StsEditWellColors.start()\n" + e);
            return false;
        }

	    actionManager.endCurrentAction();
        return true;
    }
}
