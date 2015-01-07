package com.Sts.Actions.Wizards.Velocity;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;

import java.awt.*;

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
public class StsCreateVelocityModel extends StsWizardStep implements Runnable
{
	StsProgressPanel panel;
	StsHeaderPanel header, info;

	public StsCreateVelocityModel(StsWizard wizard)
	{
		super(wizard);
		panel = StsProgressPanel.constructor(5, 50);
		header = new StsHeaderPanel();
		setPanels(panel, header);
        panel.setPreferredSize(new Dimension(300, 300));
		header.setTitle("Create Velocity Model");
		header.setSubtitle("Using Surface-Marker Ties");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Velocity");
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
			success = ((StsVelocityWizard)wizard).constructVelocityModel(panel);
            if(success)
            {
                panel.appendLine("Complete. Press Finish> Button");
                wizard.enableFinish();
                panel.finished();
            }
		}
		catch(Exception e)
		{
			StsException.outputException("StsCreateVelocityModel.run() failed.", e, StsException.WARNING);
		}
    }

	public boolean end()
	{
        model.getProject().checkSetDepthRangeForVelocityModel();
        return success;
	}
}
