package com.Sts.Actions;

import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsSetCursor3dAction extends StsAction
{
	StsView currentView = null;

	public StsSetCursor3dAction(StsActionManager actionManager)
	{
		super(actionManager, true);
//		currentView = glPanel3d.getCurrentView();
	}

	public boolean start()
	{
		return true;
	}

	public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
	{
//		StsView currentView = glPanel3d.getCurrentView();

 //       if(!mouse.isButtonStateReleased(StsMouse.LEFT)) return true;

		if(currentView instanceof StsView3d)
			return currentView.moveCursor3d(mouse, (StsGLPanel3d)glPanel);

        if(currentView instanceof StsView2d)
            return currentView.moveCursor3d(mouse, (StsGLPanel3d)glPanel);

		return true;
	}

	public boolean end()
	{
//		StsView view = glPanel3d.getCurrentView();
//		if(view == null) return false;
//		view.cursorButtonState = StsMouse.CLEARED;
		statusArea.textOnly();
		logMessage(" ");
		return true;
	}
}
