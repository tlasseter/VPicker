
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001i
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Sections;

import com.Sts.Actions.*;
import com.Sts.MVC.*;
import com.Sts.UI.Toolbars.*;

// Here we just want to install a toolbar and then terminate the action
public class StsBuildVerticalSections extends StsAction
{
    StsBuildVerticalSectionsToolbar buildFrameToolbar;

     public StsBuildVerticalSections(StsActionManager actionManager)
    {
        super(actionManager);
        buildFrameToolbar = new StsBuildVerticalSectionsToolbar(model.win3d);
    }

    public boolean start() { return false; }

    public void checkAddToolbar()
    {
        model.win3d.addToolbar(buildFrameToolbar);
    }
}

