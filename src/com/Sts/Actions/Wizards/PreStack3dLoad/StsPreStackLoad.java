
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System



package com.Sts.Actions.Wizards.PreStack3dLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsPreStackLoad extends StsWizardStep implements Runnable
{
    private StsProgressPanel panel;
    private StsHeaderPanel header;
    private StsPreStackLoadWizard wizard;
    public StsPreStackLineSet3d preStackLineSet3d;
    
    public StsPreStackLoad(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsPreStackLoadWizard)wizard;
        constructPanel();
    }

    private void constructPanel()
    {
        panel = StsProgressPanel.constructorWithCancelButton();
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("PreStack-3d Loading");
        header.setSubtitle("Load 3D Line(s)");
        header.setInfoText(wizardDialog,"(1) Once loading is complete, press the Finish Button to dismiss the screen.");
    }

    public boolean start()
    {
        wizard.dialog.setTitle("Load 3D Line(s)");
        run();
        return true;
    }

    public void run()
	{
        if (panel.isCanceled())
        {
            panel.setDescriptionAndLevel("Seismic volume load canceled.", StsProgressBar.WARNING);
            return;
        }
        StsModel model = wizard.getModel();
        String name = wizard.getVolumeName();
        try
        {
            StsFile[] selectedFiles = wizard.getSelectedFiles();

            disablePrevious();

            model.disableDisplay();
            panel.appendLine("Loading Pre-Stack3d " + name);
            preStackLineSet3d = StsPreStackLineSet3d.constructor(name, selectedFiles, model, panel);

            if (preStackLineSet3d == null)
            {
                panel.setDescriptionAndLevel("Failed to load one or more seismic volumes.", StsProgressBar.ERROR);
                return;
            }

            preStackLineSet3d.setPredecessor(wizard.getPredecessor());

            if (wizard.getSpectrum() != null)
                preStackLineSet3d.getSeismicColorscale().setSpectrum(((StsPreStackLoadWizard)wizard).getSpectrum());

            panel.setDescription("Loading Complete.");
        }
        catch (Exception e)
        {
            panel.setDescriptionAndLevel("Exception thrown", StsProgressBar.ERROR);
            panel.appendLine("Exception thrown: " + e.getMessage());
            StsException.outputWarningException(this, "run", e);
            new StsMessage(wizard.frame, StsMessage.WARNING, "Exception thrown loading prestack lines.");
            return;
        }
        try
        {
       //     model.viewObjectChanged(this, preStackLineSet3d); //lineSet not initialized yet, if certain views are up causes NPE
//            model.project.runCompleteLoading();
            enableFinish();
            success = true;
            return;
        }
        catch (Exception e)
        {
            panel.appendLine("EsceptionThrown completing load: " + e.getMessage());
            StsException.outputException("StsPreStackLoad.start() failed.", e, StsException.WARNING);
            success = false;
            return;
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
