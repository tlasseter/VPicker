
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions;

import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;

public class StsSurfacePointDebug extends StsSurfaceAction
{
    public StsSurfacePointDebug(StsActionManager actionManager)
    {
        super(actionManager);
        title = new String("Debug point reporting: ");
    }

    public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
    {
        StsMousePoint mousePoint;

    	int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
       	if (leftButtonState == StsMouse.PRESSED || leftButtonState == StsMouse.DRAGGED)
        {
            if(surface == null) return false;
            mousePoint = mouse.getMousePoint();
            surface.getSurfacePosition(mouse, true, (StsGLPanel3d)glPanel);
        }
        else if(StsMouse.isButtonStateReleasedOrClicked(leftButtonState))
        {
            mousePoint = mouse.getMousePoint();
            StsGridPoint gridPoint = surface.getSurfacePosition(mouse, true, (StsGLPanel3d)glPanel);
            if(gridPoint != null) surface.debugSurfacePoint(gridPoint);
            actionManager.endCurrentAction();
        }
        return true;
    }
}
