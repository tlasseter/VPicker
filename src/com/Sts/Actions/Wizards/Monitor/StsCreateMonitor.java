
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Monitor;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.UI.*;

public class StsCreateMonitor extends StsWizardStep
{
    StsWatchMonitorPanel panel;
    StsHeaderPanel header;

    public StsCreateMonitor(StsWizard wizard)
    {
        super(wizard);
        panel = new StsWatchMonitorPanel(model, null, false);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Create a Real-time Monitor");
        header.setSubtitle("Create the Monitor Object");
        header.setInfoText(wizardDialog, "(1)Press the <Back Button to create another monitor\n" +
                           "(2) By default all monitors are active. Use Enable and Disable to customize.\n" +
                           "       While operating the eyeglass button on the time control toolbar will provide similar capabilities.\n" +
                           "(3) Press the Finish Button to dismiss Wizard.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Monitor");                
    }

    public boolean start()
    {
        panel.initialize();
        wizard.rebuild();
        wizard.disableCancel();
        wizard.enableFinish();
        return true;
    }

    public boolean end()
    {
        return true;
    }

}