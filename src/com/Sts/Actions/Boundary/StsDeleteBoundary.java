
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Boundary;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;

// Here we just want to install a toolbar and then terminate the action
public class StsDeleteBoundary extends StsAction
{

    public StsDeleteBoundary(StsActionManager actionManager)
    {
        super(actionManager);
    }

    public boolean start()
    {

        StsBuiltModel.checkDeleteBoundary(model);
        actionManager.endCurrentAction();
        return true;
    }

    public boolean end()
    {
        model.setActionStatus(StsBuildBoundary.class.getName(), StsModel.CAN_START);
        return true;
    }
}
