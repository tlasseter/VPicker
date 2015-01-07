
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System



package com.Sts.Actions.Wizards.PreStack2dLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsPreStackLoad2d extends StsWizardStep implements Runnable
{
    private StsProgressPanel panel;
    private StsHeaderPanel header;
    private StsPreStackLoad2dWizard wizard;
    protected StsPreStackLineSet2d preStackLineSet2d;

    public StsPreStackLoad2d(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsPreStackLoad2dWizard)wizard;
        constructPanel();
    }

    public void constructPanel()
    {
        panel = StsProgressPanel.constructorWithCancelButton();
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Pre-Stack 2D Loading");
        header.setSubtitle("Load 2D Line(s)");
        header.setInfoText(wizardDialog,"(1) Once loading is complete, press the Finish Button to dismiss the screen.");
    }

    public boolean start()
    {
        wizard.dialog.setTitle("Load 2D Line(s)");
        run();
        return true;
    }

    public void run()
	{
        if(panel.isCanceled())
        {
            panel.setDescriptionAndLevel("Seismic 2d load canceled.", StsProgressBar.WARNING);
            return;
        }
        StsModel model = wizard.getModel();
        String name = wizard.getLineName();
        try
        {
            StsFile[] selectedFiles = wizard.getSelectedFiles();

            disablePrevious();

            model.disableDisplay();

            preStackLineSet2d = StsPreStackLineSet2d.constructor(name, selectedFiles, model, panel, wizard.getBoundingBox());
            wizard.addVolume(preStackLineSet2d);
            //	if (((StsPreStackLoad2dWizard)wizard).getSpectrum() != null)
            //		line2d.getSeismicColorscale().setSpectrum(((StsPreStackLoadWizard)wizard).getSpectrum());

            preStackLineSet2d.setPredecessor(wizard.getPredecessor());

            panel.appendLine("Loading Complete");
            wizard.enableFinish();
            //			checkSetZDomain(line2d);
        }
        catch (Exception e)
        {
            panel.finished();
            panel.appendLine("Exception thrown: " + e.getMessage());
            panel.setDescriptionAndLevel("Exception thrown", StsProgressBar.ERROR);
            new StsMessage(wizard.frame, StsMessage.WARNING, "Failed to load any seismic volumes.");
            return;
        }
        try
        {
            model.getProject().runCompleteLoading();
            return;
        }
        catch (Exception e)
        {
            panel.appendLine("Exception thrown: " + e.getMessage());
            panel.setDescriptionAndLevel("Exception thrown", StsProgressBar.ERROR);
            StsException.outputException("StsPreStackLoad.start() failed.", e, StsException.WARNING);
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
