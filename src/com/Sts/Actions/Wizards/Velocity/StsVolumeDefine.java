package com.Sts.Actions.Wizards.Velocity;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

import java.awt.*;

/**
 * <p>Title: jS2S development</p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author T.J.Lasseter
 * @version c51c
 */
public class StsVolumeDefine extends StsWizardStep
{
	StsVolumeDefinePanel panel;
	StsHeaderPanel header, info;

	/**
	 * StsEditVelocity creates or modifies the velocity model used for
	 * time-to-depth conversion.
	 *
	 * @param wizard StsVelocityWizard of which this is a step.
	 */
	public StsVolumeDefine(StsWizard wizard)
	{
		super(wizard);
		panel = new StsVolumeDefinePanel(this);
		header = new StsHeaderPanel();
		setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 400));
		header.setTitle("Velocity PostStack3d Definition");
		header.setSubtitle("Define Selected PostStack3d");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Velocity");
		header.setInfoText(wizardDialog,"(1) Specify type (average or instantaneous) of velocity in volume.\n"
                    + "(2) Review the minimum and maximum velocity detected in the volume.\n"
                    + "(3) Select the appropriate velocity units for volume data.\n"
                    + "(4) Specify whether one-way or two-way travel time.\n"
                    + "(5) Press Next> Button to proceed to model construction.");
	}
	public boolean start()
    {
        panel.initializePanel();
        return true;
    }

    public boolean end()
    {
        return true;
    }
}
