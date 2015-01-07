
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

public class StsSurfaceReporting extends StsSurfaceAction
{
    public StsSurfaceReporting(StsActionManager actionManager)
    {
        super(actionManager, true);
        title = new String("Surface reporting: ");
        if(surface == null) return;
    }

    public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
    {
        StsGLPanel3d glPanel3d = (StsGLPanel3d)glPanel;
        StsView view = glPanel3d.getView();
        if(!(view instanceof StsView3d)) return true;
        ((StsView3d)view).surfaceReadout(surface, mouse);
        return true;
    }

  	public boolean end()
    {
        StsGLPanel3d glPanel3d = model.getGlPanel3d();
        StsView view = glPanel3d.getView();
        if(view == null) return false;
        view.cursorButtonState = StsMouse.CLEARED;
        statusArea.textOnly();
    	logMessage(" ");
        return true;
    }
}
