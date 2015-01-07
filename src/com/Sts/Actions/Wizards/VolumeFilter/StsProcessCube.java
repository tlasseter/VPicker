//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.VolumeFilter;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.SeismicAttributes.*;
import com.Sts.UI.Progress.*;

import java.awt.*;

public class StsProcessCube extends StsWizardStep implements Runnable
{
    private StsVolumeFilterWizard wizard;

	private StsProgressPanel panel;
	private StsHeaderPanel header;
	private StsSeismicVolume volume;
	private boolean isDataFloat;

	public StsProcessCube(StsWizard wizard)
	{
		super(wizard);
        this.wizard = (StsVolumeFilterWizard)wizard;
		header = new StsHeaderPanel();
		panel = StsProgressPanel.constructor(5, 50);
		setPanels(panel, header);
		panel.setPreferredSize(new Dimension(400, 300));
		header.setTitle("Edge Preserving Volume Filter");
		header.setSubtitle("Defining Filtered Volume");
	    header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/VirtualVolumes.html");
	    header.setInfoText(wizardDialog,"(1) Select a volume to apply the smoothing filter.\n" +
	                           "(2) Select the type of filter to apply.\n" +
	                           "(3) Define the filter parameters.\n" +
	                           "(4) Press the Finish>> Button to create the volume.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VolumeFilter");

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
 		volume = ( (StsVolumeFilterWizard) wizard).getSeismicVolume();
		isDataFloat = ( (StsVolumeFilterWizard) wizard).isDataFloat();

		Runnable runCreateVolume = (new Runnable()
		{
			public void run()
			{
				int xSize = wizard.defineType.getXSize();
		    	int ySize = wizard.defineType.getYSize();
		    	int zSize = wizard.defineType.getZSize();
				byte analysis = wizard.defineType.panel.getAnalysisType();
				StsVolumeFilterConstructor.constructor(model, volume, true, panel, xSize, ySize, zSize, analysis);
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
