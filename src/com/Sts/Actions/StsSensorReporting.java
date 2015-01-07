
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions;

import com.Sts.DBTypes.StsSensorClass;
import com.Sts.MVC.StsActionManager;
import com.Sts.MVC.StsGLPanel;
import com.Sts.MVC.View3d.StsGLPanel3d;
import com.Sts.MVC.View3d.StsView;
import com.Sts.MVC.View3d.StsView3d;
import com.Sts.Types.StsMouse;

public class StsSensorReporting extends StsAction
{
    public StsSensorReporting(StsActionManager actionManager)
    {
        super(actionManager, true);
    }

    public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
    {
        StsGLPanel3d glPanel3d = (StsGLPanel3d)glPanel;
        StsView view = glPanel3d.getView();
        if(!(view instanceof StsView3d)) return true;
        ((StsView3d)view).sensorReadout(mouse);
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
