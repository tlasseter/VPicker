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
public class StsFromSonics extends StsWizardStep
{
	StsFromSonicsPanel panel;
	StsHeaderPanel header, info;

	/**
	 * StsimportProfiles imports ASCII profiles and defines volume for
	 * time-to-depth conversion.
	 *
	 * @param wizard StsVelocityWizard of which this is a step.
	 */
	public StsFromSonics(StsWizard wizard)
	{
		super(wizard);
		panel = new StsFromSonicsPanel(this);
		header = new StsHeaderPanel();
		setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 400));
		header.setTitle("Velocity Model Construction");
		header.setSubtitle("Use Time/Depth functions");
		header.setInfoText(wizardDialog,"Press Next> Button to construct model.\n");

		//header.setInfoText(wizardDialog,"(1) Select HANDVEL formatted velocity profile files\n" +
		//				   "(2) Define the volume to construct.\n" +
        //                   "(3) Press Next> Button to construct model.\n");
	}

    public boolean start()
	  {
		  return true;
	  }

	  public boolean end()
     {
		 return true;
	 }
}
