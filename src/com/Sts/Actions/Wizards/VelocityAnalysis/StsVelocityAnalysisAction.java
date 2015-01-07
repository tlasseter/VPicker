package com.Sts.Actions.Wizards.VelocityAnalysis;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;

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
public class StsVelocityAnalysisAction extends StsAction implements Runnable
{
	StsPreStackLineSet currentVolume = null;
    StsPreStackVelocityModel velocityVolume = null;
	StsVelocityAnalysisToolbar velocityAnalysisToolbar = null;

	public StsVelocityAnalysisAction(StsActionManager actionManager)
	{
		super(actionManager);
	}

	public void run()
	{
		start();
	}

	public boolean start()
	{
		currentVolume = StsPreStackLineSetClass.currentProjectPreStackLineSet;
        if(currentVolume == null) return false;
        velocityVolume = currentVolume.getVelocityModel();
		if(velocityVolume == null)
		{
			actionManager.endCurrentAction();
			actionManager.startAction(StsVelocityAnalysisWizard.class);
			return false;
		}
		new StsMessage(model.win3d, StsMessage.INFO, "Velocity model exists.  Edit properties on object panel.", true);
		setCursor();
        currentVolume.lineSetClass.setShowVelStat(true);
        StsVelocityAnalysisToolbar.checkAddToolbar(model, currentVolume, true);
		actionManager.endCurrentAction();
		return true;
	}

	private void setCursor()
	{
        currentVolume.setCurrentDataRowCol(model.win3d, 0, 0);
        currentVolume.nextProfile();
    }

	public boolean end()
	{
		return true;
	}
}
