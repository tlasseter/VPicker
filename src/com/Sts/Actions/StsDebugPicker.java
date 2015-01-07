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

public class StsDebugPicker extends StsAction
{
    StsGLPanel3d glPanel3d;

    public StsDebugPicker(StsActionManager actionManager)
    {
        super(actionManager);
//        glPanel3d = model.win3d.glPanel3d;
    }

    public boolean start()
    {
//        glPanel3d.debugPicker = true;
//        glPanel3d.debugPicking = true;
        return true;
    }

    public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
    {
        if(glPanel3d != null && glPanel3d != glPanel)
            glPanel3d.debugPicker = false;
        glPanel3d = (StsGLPanel3d)glPanel;
        glPanel3d.debugPicker = true;
        StsView view = glPanel3d.getView();
        if(!(view instanceof StsView3d)) return true;
        int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
        if (StsMouse.isButtonStateReleasedOrClicked(leftButtonState))
            glPanel3d.debugPicking = true;
        else
            glPanel3d.debugPicking = false;
        System.out.println("Left mouse released.");
        glPanel3d.projectionChanged = true;
        glPanel3d.repaint();
        return true;
    }

    public boolean end()
    {
        if(glPanel3d != null)
        {
            glPanel3d.debugPicker = false;
            glPanel3d.projectionChanged = true;
        }
        return true;
    }
}