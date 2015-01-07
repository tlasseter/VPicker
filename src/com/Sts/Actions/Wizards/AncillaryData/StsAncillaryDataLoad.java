package com.Sts.Actions.Wizards.AncillaryData;

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;

public class StsAncillaryDataLoad extends StsWizardStep implements Runnable
{
    public StsProgressPanel panel;
    private StsHeaderPanel header;

    private boolean isDone = false;
    private boolean canceled = false;

    public StsAncillaryDataLoad(StsWizard wizard)
    {
        super(wizard);
    }

    public void constructPanel()
    {
        panel = StsProgressPanel.constructorWithCancelButton();
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Ancillary Data File Load");
        header.setSubtitle("Load File(s)");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#AncillaryData");
        header.setInfoText(wizardDialog,"(1) Once loading is complete, press the Finish Button to dismiss the screen");
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
        try
        {
            disablePrevious();
            disableCancel();
            panel.appendLine("Starting ancillary data loading...");
            success = ((StsAncillaryDataWizard)wizard).loadAncillaryData(panel);
            panel.setDescription("Completed loading ancillary data.");
        }
        catch(Exception e)
        {
            success = false;
            panel.setDescriptionAndLevel("Error loading ancillary data.", StsProgressBar.ERROR);
            StsException.outputException("StsAncillaryDataLoad.run() failed.", e, StsException.WARNING);
        }
    }

    public boolean end()
    {
        return true;
    }
}
