
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2006
//Author:       TJLasseter
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.VspLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsVspLoad extends StsWizardStep implements Runnable
{
    public StsProgressPanel panel;
    private StsHeaderPanel header;
    StsVsp[] seismicVolumes;
    StsVspLoadWizard wizard;
    private boolean isDone = false;

    public StsVspLoad(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsVspLoadWizard)wizard;
    }

    public void constructPanel()
    {
        panel = StsProgressPanel.constructorWithCancelButton();
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("PostStack3d Selection");
        header.setSubtitle("Load PostStack3d(s)");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VspLoad");                
        header.setInfoText(wizardDialog,"(1) Once complete, press the Finish Button to dismiss the screen");
    }

    public boolean start()
    {
        run();
        return true;
    }

    public void run()

    {
        if (panel.isCanceled())
       {
           success = false;
           return;
       }
        StsModel model = wizard.getModel();
        StsProject project = model.getProject();
        int nExistingVolumes = model.getObjectList(StsVsp.class).length;
        int n = -1;
        String name = "null";

        try
        {
            StsAbstractFile[] selectedFiles = wizard.getSelectedFiles();
            String[] names = wizard.getFilenameEndings();

//            seismicVolumes = ((StsVolumeWizard)wizard).getSeismicVolumes();
//            if(seismicVolumes == null || seismicVolumes.length == 0) return false;

            model.disableDisplay();
            panel.initialize(selectedFiles.length);
            for(n = 0; n < selectedFiles.length; n++)
            {
                panel.appendLine("Loading vsp dataset: " + names[n]);

                StsVsp seismicVolume = StsVsp.constructor((StsFile)selectedFiles[n], model);
                 if(seismicVolume == null)
                 {
                     panel.appendLine("Failed to add VSP to project: " + selectedFiles[n].getPathname());
                     new StsMessage(model.win3d, StsMessage.WARNING, "Failed to add VSP to project: " + selectedFiles[n].getPathname());
                     continue;
                 }

                if((wizard).getSpectrum() != null)
                    seismicVolume.getColorscale().setSpectrum(((StsVspLoadWizard)wizard).getSpectrum());
                statusArea.setText("VSP " + seismicVolume.getName() + " load successful.");
                panel.setValue(n+1);
                panel.setDescription("VSP " + n + " of " + selectedFiles.length + " loaded.");
                checkSetZDomain(seismicVolume);
            }
            isDone = true;
            panel.appendLine("VSP data loaded successfully.\n");
            panel.finished();
            panel.setDescription("Loading Complete");
        }
        catch(Exception e)
        {
            String message;
            if(n == -1) message = new String("Failed to load any seismic volumes.");
            else        message = new String("Failed to load volume " + name + ".\n" +
                                             "Error: " + e.getMessage());
            panel.appendLine(message);
            new StsMessage(wizard.frame, StsMessage.WARNING, message);
            panel.setDescriptionAndLevel(message + "\n", StsProgressBar.WARNING);
            panel.finished();
            success = false;
            return;
        }
        try
        {
            if(nExistingVolumes > 0) project.runCompleteLoading();
            success = true;
            return;
        }
        catch(Exception e)
        {
            panel.setDescriptionAndLevel("Exception thrown." + e.getMessage(), StsProgressBar.ERROR);
            StsException.outputException("StsVspLoad.start() failed.", e, StsException.WARNING);
            success = false;
            return;
        }
        finally
        {
            wizard.completeLoading(success);
        }
    }

	public void checkSetZDomain(StsVsp seismicVolume)
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
