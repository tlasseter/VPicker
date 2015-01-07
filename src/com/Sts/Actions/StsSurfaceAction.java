
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;

abstract public class StsSurfaceAction extends StsAction
{
    protected StsSurface surface;
    protected String title = new String("Action label not implemented: ");
    protected String label = null;

    public StsSurfaceAction(StsActionManager actionManager, boolean canInterrupt)
    {
        super(actionManager, canInterrupt);
    }

    public StsSurfaceAction(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public boolean start()
    {
        statusArea.setTitle(title);

        try
        {
            getVisibleSurface();
            if(surface == null) return false;

            label = "surface: " + surface.getName();

            if(!surface.checkIsLoaded())
            {
                statusArea.textOnly();
                logMessage("Unable to load " + label + ".");
                return false;
            }

            infoMessage("Hold down left mouse button and move around surface." +
                    "  Release button to terminate.");
        }
        catch(Exception e) { return false;}

        return true;
    }

    private void getVisibleSurface()
    {
        // We might have horizons but not zones yet, so check for horizons
        StsModelSurfaceClass modelSurfaceClass = (StsModelSurfaceClass)model.getCreateStsClass(StsModelSurface.class);
        surface = modelSurfaceClass.selectSurface(true, StsModelSurface.class);
        if (surface != null) return;

		StsSurfaceClass surfaceClass = (StsSurfaceClass)model.getCreateStsClass(StsSurface.class);
        surface = surfaceClass.selectSurface(true, StsSurface.class);
        if (surface != null) return;

        statusArea.textOnly();
        logMessage("No surface available or selected.");
        return;
    }

    /** mouse action for 3d window */
   	abstract public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel);

/*
   	public boolean performMouseAction(StsMouse mouse)
	{
        StsMousePoint mousePoint;

    	int leftButtonState = mouse.getButtonState(StsMouse.LEFT);
       	if (leftButtonState == StsMouse.PRESSED || leftButtonState == StsMouse.DRAGGED)
        {
            mousePoint = mouse.getMousePoint();
            performMousePressedOrDraggedAction(mousePoint);
        }
        else if (StsMouse.isButtonStateReleasedOrClicked(leftButtonState))
        {
            mousePoint = mouse.getMousePoint();
            performMouseReleasedAction(mousePoint);
            glPanel.getActionManager().endCurrentAction();
        }
        return true;
 	}

    protected void performMousePressedOrDraggedAction(StsMousePoint mousePoint)
    {
        StsException.systemError("StsSurfaceAction.performMousePressedOrDraggedAction() failed." +
            " not implemented by subclass: " + this.getClass().toString());
    }

    protected void performMouseReleasedAction(StsMousePoint mousePoint)
    {
        StsException.systemError("StsSurfaceAction.performMouseReleasedAction() failed." +
            " not implemented by subclass: " + this.getClass().toString());
    }
*/
  	public boolean end()
    {
        statusArea.textOnly();
    	logMessage(" ");
        return true;
    }
}
