
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
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;

// Here we just want to install a toolbar and then terminate the action
public class StsBuildFrame extends StsAction
{
	StsBuildFrameToolbar buildFrameToolbar;

 	public StsBuildFrame(StsActionManager actionManager)
    {
        super(actionManager);

        buildFrameToolbar = (StsBuildFrameToolbar)model.win3d.getToolbarNamed(StsBuildFrameToolbar.NAME);
        if(buildFrameToolbar != null)
        {
            buildFrameToolbar.setVisible(true);
            return;
        }
        buildFrameToolbar = new StsBuildFrameToolbar(model.win3d);
        model.win3d.addToolbar(buildFrameToolbar);
    }

    public boolean checkStartAction()
    {
        StsBuiltModel builtModel = (StsBuiltModel)model.getCurrentObject(StsBuiltModel.class);
        if(builtModel == null) return true;
        boolean yes = StsYesNoDialog.questionValue(model.win3d,"Model already built.\n  Do you wish to delete it so additional faults can be defined?.");
        if(yes)
        {
            builtModel.delete();
            model.setActionStatus(StsBuildFrame.class.getName(), StsModel.CAN_START);
        }
        return yes;
    }

	public boolean start()
    {
        return false; 
    }

	public void checkAddToolbar()
	{
		model.win3d.addToolbar(buildFrameToolbar);
	}
}
