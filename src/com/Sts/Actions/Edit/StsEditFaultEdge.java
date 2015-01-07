
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Edit;

import com.Sts.Actions.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;

public class StsEditFaultEdge extends StsAction
{
    /** constructor */
    public StsEditFaultEdge(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public boolean start()
    {
        /*
        StsStatusArea.staticSetTitle("Edit fault cuts:");
        logMessage("Select a fault cut to edit.");
        addEndButton();
        return true;
        */

        logMessage("Fault cut editing not yet implemented...");
        return false;
    }

    public boolean end()
    {
        return true;
    }

    /** mouse action for 3d window */
   	public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
    {
        return true;
    }
}
