
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Monitor;

import java.awt.Dimension;
import java.io.File;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.Actions.Wizards.WizardHeaders.StsHeaderPanel;
import com.Sts.DBTypes.StsSensor;
import com.Sts.DBTypes.StsMonitor;
import com.Sts.IO.StsFile;

import javax.swing.*;

public class StsMonitorSource extends StsWizardStep
{
    StsMonitorSourcePanel panel;
    StsHeaderPanel header;

    public StsMonitorSource(StsWizard wizard)
    {
        super(wizard);
        panel = new StsMonitorSourcePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Create a Real-time Monitor");
        header.setSubtitle("Create Manually or Define from File");
        header.setInfoText(wizard.dialog,
                "(1) Select Manual definition or define From File.\n" +
                "\nFile Header:\n" +
                "   Source Replace Reload Compute PollBy Name isStatic isRelative\n" +
                "Sample Data Row\n" +
                "   C:/S2S/Events true false true Time Stage01 false false\n\n" +
                "(2) Press the Next Button\n" +
                "(3) If Define from File is selected, a file selection dialog will be shown.\n" +
                "(4) Select the desired file (.txt) and press the Open button."
        );
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Monitor");
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
}