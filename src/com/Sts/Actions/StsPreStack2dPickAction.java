package com.Sts.Actions;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */
public class StsPreStack2dPickAction extends StsAction
{
	StsPreStackLineSet2d prestackVolume2d;

	public StsPreStack2dPickAction(StsActionManager actionManager, StsPreStackLineSet2d prestackVolume2d)
	{
		super(actionManager);
		this.prestackVolume2d = prestackVolume2d;
	}

	public void setVolume(StsPreStackLineSet2d prestackVolume2d)
	{
		this.prestackVolume2d = prestackVolume2d;
	}

	public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
	{
		if(prestackVolume2d == null) return true;
		if(!mouse.isButtonStateReleased(StsMouse.LEFT)) return true;
		mouse.clearButtonState(StsMouse.LEFT);
		return prestackVolume2d.performMouseAction(mouse, glPanel);
	}
}
