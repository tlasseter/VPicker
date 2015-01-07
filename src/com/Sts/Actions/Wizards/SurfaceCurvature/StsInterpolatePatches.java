//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SurfaceCurvature;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.SeismicAttributes.*;
import com.Sts.UI.Progress.*;

import java.awt.*;

public class StsInterpolatePatches extends StsWizardStep implements Runnable
{
	private StsProgressPanel panel;
	private StsHeaderPanel header;
	private StsPatchVolume volume;

	public StsInterpolatePatches(StsWizard wizard)
	{
		super(wizard);
		header = new StsHeaderPanel();
		panel = StsProgressPanel.constructor(5, 50);
		setPanels(panel, header);
		panel.setPreferredSize(new Dimension(400, 300));
		header.setTitle("Seismic Attribute Calculation");
		header.setSubtitle("Calculate Attribute Volumes");
		header.setInfoText(wizardDialog, "(1) Press the Finish Button to produce the desired cubes.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VolumeCurvature");

	}

	public boolean start()
	{
		volume = ( (StsVolumeCurvatureWizard) wizard).getPatchVolume();
		if (volume == null)
		{
			panel.appendLine("Patch Volume is NULL.");
			return false;
		}
		disablePrevious();
		// wizard.disableFinish();
		processVolume();
		return true;
	}

	public void run()
	{
		start();
	}

    public StsPatchVolume getVolume() { return volume; }

    public void processVolume()
	{
		final int maxInterpolationSize =((StsVolumeCurvatureWizard)wizard).processVolume.panel.interpSize;
        final byte curveType =((StsVolumeCurvatureWizard)wizard).processVolume.panel.attributeType;

        Runnable runCreateVolume = (new Runnable()
		{
			public void run()
			{
				StsCurvatureVolumeConstructor.createInterpolatedVolume(model, volume, panel, curveType, maxInterpolationSize);
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
