//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.PreStack3d;

import com.Sts.Actions.Wizards.PreStack2d.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.MVC.*;
import com.Sts.Types.PreStack.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;


public class StsPreStackProcess extends StsWizardStep implements Runnable
{
	public StsProgressPanel panel;
	public StsHeaderPanel header;
	public StsSeismicBoundingBox[] segyVolumes;
	StsSEGYFormatRec[] selectedAttributes;
    private StsPreStackSegyLine currentVolume = null;
    public int nSuccessfulVolumes = 0;

	public StsPreStackWizard wizard;

    public boolean canceled = false;
    public String processString;

	public StsPreStackProcess(StsPreStackWizard wizard)
	{
		super(wizard);
		this.wizard = wizard;
        processString = "pre-stack 3d";
    }

    public StsPreStackProcess(StsPreStack2dWizard wizard)
	{
		super(wizard);
		this.wizard = wizard;
        processString = "pre-stack 2d";
    }

    public boolean start()
	{
		run();
		return true;
	}

	public void constructPanel()
	{
		segyVolumes = wizard.getSegyVolumesToProcess();
		if (segyVolumes == null) return;

		panel = StsProgressPanel.constructorWithCancelButton(20, 50);
		header = new StsHeaderPanel();
		setPanels(panel, header);
		header.setTitle("Processing " + processString);
        header.setSubtitle("Processing Volumes");
        header.setInfoText(wizardDialog, "(1) Once complete, press the Finish Button to dismiss the screen\n" +
								 "     ***** Proceed to the pre-stack seismic load workflow step to load and view the data *****\n" +
								 "     ***** Cancel will dismis the screen and stop the processing after the current line is complete *****\n");
        wizard.dialog.setTitle("Process " + processString);
	}

	public void run()
	{
        boolean someExist = false;
		try
		{
			if (segyVolumes == null) return;
			if (panel.isCanceled())
			{
                setCancelStatus(null);
				return;
			}
            String existingFiles = new String( "Some of the volumes have already been processed: \n    Do you wish to continue?");

            for (int i = 0; i < segyVolumes.length; i++)
            {
                StsPreStackSegyLine volume = (StsPreStackSegyLine)segyVolumes[i];
                if(volume.overwriteFiles())
                {
                   someExist = true;
                   continue;
                }
            }
            if (someExist)
            {
                if (!StsYesNoDialog.questionValue(wizard.frame, existingFiles))
                {
                   wizard.gotoFirstStep();
                   return;
                }
            }
			// Process the volumes
			disableFinish();
			disablePrevious();
			enableCancel();

            nSuccessfulVolumes = 0;
            boolean ignoreMultiVolume = wizard.getIgnoreMultiVolume();
			for (int i = 0; i < segyVolumes.length; i++)
			{
				if (panel.isCanceled())
				{
					setCancelStatus(null);
					break;
				}
				currentVolume = (StsPreStackSegyLine)segyVolumes[i];

                panel.appendLine("Processing " + processString + " for " + currentVolume.getName());

				selectedAttributes = currentVolume.getSEGYFormat().getRequiredTraceRecords();
                for(int j=0; j<currentVolume.getNumVolumes(); j++)
                    success = currentVolume.readWritePreStackLines(panel, selectedAttributes, wizard.getOverrideGeometry(), ignoreMultiVolume, j);

                if (panel.isCanceled())
				{
					setCancelStatus(currentVolume);
					break;
				}

				Main.logUsage();
				if (success)
				{
                    panel.setLevel( StsProgressBar.INFO);
					panel.appendLine(processString + " " + currentVolume.getName() + " successfully processed.");
                    nSuccessfulVolumes++;
				}
				else
				{
                    panel.setLevel( StsProgressBar.ERROR);
					panel.appendLine(processString + " " + currentVolume.getName() + " processing failed.");
				}
				wizard.removeSegyVolumeToProcess(currentVolume);
			}
            if(nSuccessfulVolumes > 0)
			if (!canceled)
			{
                wizard.enableFinish();
                wizard.disableCancel();
			}
			if (nSuccessfulVolumes == segyVolumes.length)
            {
                panel.setLevel( StsProgressBar.INFO);
				panel.appendLine("All volumes completed OK.");
            }
            else if(nSuccessfulVolumes > 0)
            {
                panel.setLevel( StsProgressBar.WARNING);
                panel.appendLine("Some but not all volumes completed.");
            }
            panel.finished();
        }
		catch (Exception e)
		{
			StsException.outputException("StsPreStackProcess.run() failed.", e, StsException.WARNING);
			String message = new String("Failed to load volume " + currentVolume.getName() + ".\n" + "Error: " + e.getMessage());
            panel.progressBar.finished();
            panel.progressBar.setDescriptionAndLevel("Exception thrown.", StsProgressBar.ERROR);
            panel.appendLine(message);
			new StsMessage(wizard.frame, StsMessage.WARNING, message);
		}
	}

    public int getNSuccessfulVolumes() { return nSuccessfulVolumes;}

    public void setCancelStatus(StsPreStackSegyLine volume)
    {
        panel.progressBar.setDescriptionAndLevel("Cancelled", StsProgressBar.WARNING);
        if(volume == null)
            panel.appendLine("User has cancelled processing.");
        else
            panel.appendLine("User has cancelled processing for volume " + volume.getName());
        enablePrevious();
    }

	public boolean end()
	{
		return true;
	}
}
