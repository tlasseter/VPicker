
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.PostStack3dLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsVolumeLoad extends StsWizardStep implements Runnable
{
    private StsProgressPanel panel;
    private StsHeaderPanel header;
    StsSeismicVolume[] seismicVolumes;
    private StsVolumeWizard wizard = null;
    private boolean isDone = false;
	private boolean canceled = false;

    public StsVolumeLoad(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsVolumeWizard)wizard;
    }


    public void constructPanel()
    {
        panel = StsProgressPanel.constructorWithCancelButton();
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Post-Stack3d Selection");
        header.setSubtitle("Load Post-Stack3d File(s)");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#PostStack3dLoad");                
        header.setInfoText(wizardDialog,"(1) Once complete, press the Finish Button to dismiss the screen");
    }

    public boolean start()
    {
        run();
        return true;
    }

    public void run()

    {
        if (canceled)
       {
           success = false;
           return;
       }
       StsModel model = wizard.getModel();
       StsProject project = model.getProject();
       int nExistingVolumes = model.getObjectList(StsSeismicVolume.class).length;
       int n = -1;
       String name = "null";
       int nSuccessfulFiles = 0;

       try
       {
           StsFile[] selectedFiles = wizard.getSelectedFiles();
           String[] names = wizard.getFilenameEndings();
           disablePrevious();
           model.disableDisplay();

           panel.initialize(selectedFiles.length);
           for(n = 0; n < selectedFiles.length; n++)
           {
               panel.appendLine("Loading seismic volume: " + names[n]);
               StsSeismicVolume seismicVolume = StsSeismicVolume.constructor(selectedFiles[n], model, true);
               if(seismicVolume == null)
               {
                   panel.appendLine("Failed to add volume to project: " + selectedFiles[n].getPathname());
                   panel.setDescriptionAndLevel("Failed to load " + selectedFiles[n].getPathname(), StsProgressBar.ERROR);
                   StsToolkit.sleep(500);
                   panel.setDescriptionAndLevel("Continuing...", StsProgressBar.INFO);
                //   new StsMessage(model.win3d, StsMessage.WARNING, "Failed to add volume to project: " + selectedFiles[n].getPathname());
                   continue;
               }

                if((wizard).getSpectrum() != null)
                    seismicVolume.getColorscale().setSpectrum(wizard.getSpectrum());

                panel.setValue(n+1);
                panel.appendLine("Post-Stack3d " + n + " of " + selectedFiles.length + " loaded.");
                statusArea.setText("Post-Stack3d " + seismicVolume.getName() + " load successful.");
                nSuccessfulFiles++;
                checkSetZDomain(seismicVolume);
            }
            isDone = true;
            panel.finished();
            panel.appendLine("Loading Complete");
            if(nSuccessfulFiles != selectedFiles.length)
                panel.setDescriptionAndLevel("Some files not loaded", StsProgressBar.WARNING);
//            wizard.enableFinish();
        }
        catch(Exception e)
        {
            String message;
            if(n == -1) message = new String("Failed to load any seismic volumes.");
            else        message = new String("Failed to load volume " + name + ".\n" +
                                             "Error: " + e.getMessage());
            new StsMessage(wizard.frame, StsMessage.WARNING, message);
            panel.setDescriptionAndLevel(message , StsProgressBar.ERROR);
            success = false;
            return;
        }
        try
        {
            if(nSuccessfulFiles > 0) project.runCompleteLoading();
            success = true;
            return;
        }
        catch(Exception e)
        {
            StsException.outputException("StsVolumeLoad.start() failed.", e, StsException.WARNING);
            panel.setDescriptionAndLevel("Exception: " + e.getMessage(), StsProgressBar.ERROR);
            success = false;
            return;
        }
        finally
        {
             wizard.completeLoading(success);

        }
    }

	public void checkSetZDomain(StsSeismicVolume seismicVolume)
	{
		if(seismicVolume.getZDomain() != StsParameters.TD_NONE) return;
        new StsComboBoxDialog(wizard.frame, "Select z domain", "Select time or depth domain.", StsParameters.TD_STRINGS, true, seismicVolume, "zDomainString");
		seismicVolume.writeHeaderFile();
	}

    public boolean end()
    {
        return true;
    }
}
