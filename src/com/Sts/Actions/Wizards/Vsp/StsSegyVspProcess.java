package com.Sts.Actions.Wizards.Vsp;

import com.Sts.Actions.Wizards.PostStack.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsSegyVspProcess extends StsPostStackProcess implements Runnable
{

	public StsSegyVspProcess(StsSegyVspWizard wizard)
	{
		super(wizard);
    }

	public void run()
	{
        StsSegyVsp currentVolume = null;
        try
        {
            // Process the volumes
            disableFinish();
            disablePrevious();
            enableCancel();
            int nSuccessfulVolumes = 0;
            for (int i = 0; i < segyVolumes.length; i++)
			{
				if (panel.isCanceled())
				{
					setCancelStatus(null);
					break;
				}
				currentVolume = (StsSegyVsp)segyVolumes[i];

                panel.appendLine("Processing Segy VSP for " + currentVolume.getName());
				panel.appendLine("\tReading SEGY VSP: " + currentVolume.getSegyFilename());

				//selectedAttributes = ((StsPostStackWizard)wizard).getSelectedAttributes();
				selectedAttributes = currentVolume.getSegyFormat().getRequiredTraceRecords();
                //success = segyVolumes[i].readWriteVolume(panel, selectedAttributes);
                boolean success = currentVolume.readWriteVolume(panel, selectedAttributes);
             	if (panel.isCanceled())
				{
					setCancelStatus(currentVolume);
					break;
				}

				Main.logUsage();
				if (success)
				{
					panel.appendLine("PostStack3d " + currentVolume.getName() + " successfully processed.");
                    nSuccessfulVolumes++;
                }
				else
				{
					panel.appendLine("PostStack3d " + currentVolume.getName() + " processing failed.");
				}
				wizard.removeSegyVolumeToProcess(currentVolume);
			}
            if(nSuccessfulVolumes > 0)
            {
                wizard.enableFinish();
                wizard.disableCancel();
            }
			if (nSuccessfulVolumes == segyVolumes.length)
				panel.appendLine("All volumes completed OK.");
			else if(nSuccessfulVolumes > 0)
				panel.appendLine("Some but not all volumes completed.");
		}
		catch (Exception e)
		{
			StsException.outputException("StsSegyVspProcess.run() failed.", e, StsException.WARNING);
			String message = new String("Failed to load volume " + currentVolume.getName() + ".\n" + "Error: " + e.getMessage());
            panel.progressBar.finished();
            panel.progressBar.setDescriptionAndLevel("Exception thrown.", StsProgressBar.ERROR);
            panel.appendLine(message);
            new StsMessage(wizard.frame, StsMessage.WARNING, message);
		}
	}
}
