//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.VolumeInterpolation;

import com.Sts.Actions.Wizards.SeismicAttribute.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.SeismicAttributes.*;
import com.Sts.UI.Progress.*;

import java.awt.*;

public class StsProcessInterpolation extends StsWizardStep implements Runnable
{
	private StsProgressPanel panel;
	private StsHeaderPanel header;
	private StsSeismicVolume volume;
	private boolean isDataFloat;

	public StsProcessInterpolation(StsWizard wizard)
	{
		super(wizard);
		header = new StsHeaderPanel();
		panel = StsProgressPanel.constructor(5, 50);
		setPanels(panel, header);
		panel.setPreferredSize(new Dimension(400, 300));
		header.setTitle("Seismic Attribute Calculation");
		header.setSubtitle("Interpolate Attribute Volume Nulls");
		header.setInfoText(wizardDialog, "(1) Press the Finish Button to produce the desired cubes.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VolumeInterpolation");
	}

	public boolean start()
	{
		disablePrevious();
//		wizard.enableFinish();
		processVolume();
		return true;
	}

	public void run()
	{
		start();
	}

    public StsSeismicVolume getVolume() { return volume; }

    public void processVolume()
	{
		volume = ( (StsSeismicAttributeWizard) wizard).getSeismicVolume();
		isDataFloat = ( (StsSeismicAttributeWizard) wizard).isDataFloat();
//		StsVolumeConstructor.createAttributeVolumes(model, volume,
//		( (StsSeismicAttributeWizard) wizard).getAttributes(), isDataFloat, panel);

		Runnable runCreateVolume = (new Runnable()
		{
			public void run()
			{
				StsSeismicVolumeConstructor.createAttributeVolumes(model, volume,
					( (StsSeismicAttributeWizard) wizard).getAttributes(), isDataFloat, panel);
	            success = true;
                wizard.enableFinish();
                model.win3dDisplay();
			}
		});

		Thread createVolumeThread = new Thread(runCreateVolume);
		createVolumeThread.start();
	}

	public boolean end()
	{
        return success;
	}
}