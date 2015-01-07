package com.Sts.Actions.Wizards.Velocity;

import java.awt.Dimension;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.Actions.Wizards.WizardHeaders.StsHeaderPanel;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.StsException;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author T.J.Lasseter
 * @version c51c
 */
public class StsCreateVelocityVolume extends StsWizardStep implements Runnable
{
	StsProgressPanel panel;
	StsHeaderPanel header, info;

	public StsCreateVelocityVolume(StsWizard wizard)
	{
		super(wizard);
		panel = StsProgressPanel.constructor(5, 50);
		header = new StsHeaderPanel();
		setPanels(panel, header);
                panel.setPreferredSize(new Dimension(300, 300));
		header.setTitle("Create Velocity Volume");
		header.setSubtitle("Using Velocity Functions");
		//header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/WellLoad.html");
        header.setInfoText(wizardDialog,"Progress bar indicates status.");
	}

	public boolean start()
	{
        disablePrevious();
        disableNext();
		run();
		return true;
	}

	public void run()
	{
		try
        {
			success = ((StsSeisVelWizard)wizard).constructVelocityVolume(panel);
            if(success)
            {
                panel.appendLine("Complete. Press Finish> Button");
                wizard.enableFinish();
                panel.finished();
            }
		}
		catch(Exception e)
		{
			StsException.outputException("StsCreateVelocityVolume.run() failed.", e, StsException.WARNING);
		}
    }

	public boolean end()
	{
        return success;
	}
}
