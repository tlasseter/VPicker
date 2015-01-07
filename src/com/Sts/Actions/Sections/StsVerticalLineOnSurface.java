
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Sections;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

public class StsVerticalLineOnSurface extends StsAction
{
    StsLine line;
    StsCursor3d cursor3d;

    /** construct a action handler for building a StsLine on a cursor. */
    public StsVerticalLineOnSurface(StsActionManager actionManager)
    {
        super(actionManager);
        addEndAllButton();
    }

    public boolean start()
    {
        try
        {
            logMessage("Pick point on top visible surface to create vertical line.");
        }
        catch(Exception e)
        {
            if (line != null) line.setDrawVertices(false);
             StsException.outputException(e, StsException.WARNING);
            return false;
        }
        return true;
    }


    /** mouse action for 3d window */
       public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
    {
        StsSurfaceClass surfaceClass = (StsSurfaceClass) model.getStsClass(StsSurface.class);
        StsSurface surface = surfaceClass.getTopVisibleSurface();
        if (surface == null) return false;

        int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

        if(leftButtonState == StsMouse.PRESSED || leftButtonState == StsMouse.DRAGGED)
        {
            surface.getSurfacePosition(mouse, true, (StsGLPanel3d)glPanel);
        }
        else if(StsMouse.isButtonStateReleasedOrClicked(leftButtonState))  // add a point
        {
            StsGridPoint gridPoint = surface.getSurfacePosition(mouse, true, (StsGLPanel3d)glPanel);
            if (gridPoint == null)return true;
            StsLine line = StsLine.buildVertical(gridPoint, StsParameters.FAULT);
            if(line != null)
            {
                logMessage(line.getLabel() + " built successfully.");
                model.win3d.win3dDisplay();
                return true;
            }
            else
            {
                logMessage("Failed to build a vertical line.");
                return true;
            }
        }
        return true;
    }

    /** add a line when the mouse action finishes */
    public boolean end()
    {
        statusArea.textOnly();

        try { Thread.sleep(1000); }
        catch(Exception e) { }
        return true;
    }
}
