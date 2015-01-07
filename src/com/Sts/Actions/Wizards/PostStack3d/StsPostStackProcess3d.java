package com.Sts.Actions.Wizards.PostStack3d;

import com.Sts.Actions.Wizards.PostStack.*;
import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsPostStackProcess3d extends StsPostStackProcess implements StsProgressRunnable
{
    private StsSegyVolume currentVolume = null;

    public StsPostStackProcess3d( StsSeismicWizard wizard)
    {
        super( wizard);
    }

    public void run()
    {
        try
        {
            if (segyVolumes == null) return;
            if (panel.isCanceled())
            {
                setCancelStatus(null);
                return;
            }
            // Process the volumes
            disableFinish();
            disablePrevious();
            enableCancel();
            nSuccessfulVolumes = 0;
            for (int i = 0; i < segyVolumes.length; i++)
            {
                if (panel.isCanceled())
                {
                    doCancel();
                    break;
                }
                currentVolume = (StsSegyVolume)segyVolumes[i];

                panel.appendLine("Processing Post-Stack3d for " + currentVolume.getName());
                panel.appendLine("\tReading SEGY volume: " + currentVolume.getSegyFilename());

                selectedAttributes = currentVolume.getSegyFormat().getRequiredTraceRecords();
                boolean success = currentVolume.readWriteVolume(model, panel, selectedAttributes);
                 if (panel.isCanceled())
                {
                    setCancelStatus(currentVolume);
                    break;
                }

                Main.logUsage();
                if (success)
                {
                    panel.setLevel( StsProgressBar.INFO);
                    panel.appendLine("Post-Stack3d " + currentVolume.getName() + " successfully processed.");
                    nSuccessfulVolumes++;
                }
                else
                {
                    panel.setLevel( StsProgressBar.ERROR);
                    panel.appendLine("Post-Stack3d " + currentVolume.getName() + " processing failed.");
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
            StsException.outputException("StsPostStackProcess.run() failed.", e, StsException.WARNING);
            String message = new String("Failed to load volume " + currentVolume.getName() + ".\n" + "Error: " + e.getMessage());
            panel.progressBar.finished();
            panel.progressBar.setDescriptionAndLevel("Exception thrown.", StsProgressBar.ERROR);
            panel.appendLine(message);
            new StsMessage(wizard.frame, StsMessage.WARNING, message);
        }
    }

    /** User may have hit cancel on progressBar in which case, running process was terminated and we just need to clean up.
     *  Otherwise user may have hit cancel on wizard which has to be sent to the panel which will terminate the process.
     *  If the first case, the call to panel.cancel() here is redundant but has no impact (it only sets the flag).
     */
    public void doCancel()
    {
        setCancelStatus(currentVolume);
        model.enableDisplay();
        panel.cancel();
    }
}
