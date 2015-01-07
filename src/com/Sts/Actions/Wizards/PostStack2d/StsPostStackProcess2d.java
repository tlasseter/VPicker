package com.Sts.Actions.Wizards.PostStack2d;

import com.Sts.Actions.Wizards.PostStack.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsPostStackProcess2d extends StsPostStackProcess implements Runnable
{
    private StsSegyLine2d currentVolume = null;

    public StsPostStackProcess2d(StsPostStack2dWizard wizard)
    {
        super( wizard);
    }

    public void run()
    {
        try
        {
            // Process the volumes
            disableFinish();
            disablePrevious();
            enableCancel();
            nSuccessfulVolumes = 0;
            for (int i = 0; i < segyVolumes.length; i++)
            {
                if (panel.isCanceled())
                {
                    setCancelStatus(null);
                    break;
                }
                currentVolume = (StsSegyLine2d)segyVolumes[i];

                panel.appendLine("Processing Post-Stack2d " + currentVolume.getName());

                selectedAttributes = currentVolume.getSegyFormat().getRequiredTraceRecords();
                boolean success = currentVolume.readWriteLines(panel, selectedAttributes);
                 if (panel.isCanceled())
                {
                    setCancelStatus(currentVolume);
                    break;
                }

                Main.logUsage();
                if (success)
                {
                    panel.setLevel( StsProgressBar.INFO);
                    panel.appendLine("Post-Stack2d " + currentVolume.getName() + " successfully processed.");
                    nSuccessfulVolumes++;
                }
                else
                {
                    panel.setLevel( StsProgressBar.ERROR);
                    panel.appendLine("Post-Stack2d " + currentVolume.getName() + " processing failed.");
                }
                wizard.removeSegyVolumeToProcess(currentVolume);
            }
            if(nSuccessfulVolumes > 0)
            {
                wizard.enableFinish();
                wizard.disableCancel();
            }
            if (nSuccessfulVolumes == segyVolumes.length)
            {
                panel.setLevel( StsProgressBar.INFO);
                panel.appendLine("All lines completed OK.");
            }
            else if(nSuccessfulVolumes > 0)
            {
                panel.setLevel( StsProgressBar.WARNING);
                panel.appendLine("Some but not all lines completed.");
            }
            panel.finished();
        }
        catch (Exception e)
        {
            StsException.outputException("StsPostStack2dProcess.run() failed.", e, StsException.WARNING);
            String message = new String("Failed to load line " + currentVolume.getName() + ".\n" + "Error: " + e.getMessage());
            panel.progressBar.finished();
            panel.progressBar.setDescriptionAndLevel("Exception thrown.", StsProgressBar.ERROR);
            panel.appendLine(message);
            new StsMessage(wizard.frame, StsMessage.WARNING, message);
        }
    }
}

