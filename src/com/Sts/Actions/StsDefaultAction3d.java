package com.Sts.Actions;

import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Toolbars.StsMouseActionToolbar;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsDefaultAction3d extends StsAction
{
	StsView currentView = null;
    StsGLPanel3d glPanel3d;

    public StsDefaultAction3d(StsGLPanel3d glPanel3d, StsActionManager actionManager)
	{
		super(actionManager, true);
        this.glPanel3d = glPanel3d;
        currentView = glPanel3d.getView();
	}

	public boolean start()
	{
		return true;
	}

	public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
	{
		StsView currentView = glPanel3d.getView();

 //       if(!mouse.isButtonStateReleased(StsMouse.LEFT)) return true;
        if(!((StsMouseActionToolbar)getModel().win3d.getToolbarNamed(StsMouseActionToolbar.NAME)).isReadoutSelected())
        {
            if(currentView instanceof StsView3d)
			    return currentView.moveCursor3d(mouse, (StsGLPanel3d)glPanel);

            if(currentView instanceof StsView2d)
                return currentView.moveCursor3d(mouse, (StsGLPanel3d)glPanel);
        }
        return true;
	}

	public boolean end()
	{
		StsView view = glPanel3d.getView();
		if(view == null) return false;
		view.cursorButtonState = StsMouse.CLEARED;
		statusArea.textOnly();
		logMessage(" ");
		return true;
	}
}
