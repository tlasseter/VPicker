package com.Sts.Actions.Wizards.Velocity;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;

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
public class StsAssignVelocity extends StsWizardStep
{
	StsAssignVelocityPanel panel;
	StsHeaderPanel header, info;

	/**
	 * StsEditVelocity creates or modifies the velocity model used for
	 * time-to-depth conversion.
	 *
	 * @param wizard StsVelocityWizard of which this is a step.
	 */
	public StsAssignVelocity(StsWizard wizard)
	{
		super(wizard);
		panel = new StsAssignVelocityPanel(this);
		header = new StsHeaderPanel();
		setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 400));
		header.setTitle("Velocity Model Construction");
		header.setSubtitle("Define / Override Intervals");
		header.setInfoText(wizardDialog,"(1) Select a horizon from dropdown list.\n" +
						   "(2) Specify the interval velocity above the selected horizon.\n" +
                           "(3) Specify the interval velocity below the last horizon.\n" +
                           "(4) Specify the desired vertical sample rate.\n" +
                           " (5) Press Next> Button to construct model.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Velocity");

	}
	public boolean start()
	  {
          StsModelSurface[] surfaces = ((StsVelocityWizard)wizard).getSelectedSurfaces();
          panel.setSurfaces(surfaces);
		  return true;
	  }

	  public boolean end()
     {
		 return true;
	 }
}
