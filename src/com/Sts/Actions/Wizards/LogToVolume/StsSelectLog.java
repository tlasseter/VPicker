
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.LogToVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;

public class StsSelectLog extends StsWizardStep
{
    StsSelectLogPanel panel;
    StsHeaderPanel header;

    public StsSelectLog(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSelectLogPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Create Volume from Logs");
        header.setSubtitle("Select Log for Volume Computation");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#LogToVolume");
        header.setInfoText(wizardDialog,"(1) Select a log to be used to compute volume.\n" +
                "(2) If you wish to resample the log to a courser interval select the Resample checkbox.\n" +
                "(3) Specify the sample interval in number of samples.\n" +
                "(4) Select the resampling method.\n" +
                "    **** Instant - selects the sample at the interval location ****\n" +
                "    **** Average - averages the samples surrounding the interval sample ****\n" +
        		"(2) Press the Next>> Button");
    }

    public boolean start()
    {
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        return true;
    }

    public StsLogCurve getSelectedLog() { return panel.getSelectedLog(); }
}