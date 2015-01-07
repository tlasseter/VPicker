
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

public class StsDefineMonitor extends StsWizardStep
{
    StsDefineMonitorPanel panel;
    StsHeaderPanel header;

    public StsDefineMonitor(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineMonitorPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Create a Real-time Monitor");
        header.setSubtitle("Define the Sensor Monitor");
        header.setInfoText(wizardDialog, "(1) Select the type of source.\n" +
                           "(2) Specify the source location.\n" +
                           "(3) Specify whether total events and cummlative amplitude attributes are computed.\n" +
                           "(4) Specify whether to reload entire file each time new data is detected.\n" +
                           "(5) Specify whether events are replaced when duplicate times are detected\n" +
                           "(6) Decide whether to poll on last modified time or on file size change.\n" +
                           "(7) Press the Next>> Button.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Monitor");
    }

    public boolean start()
    {
        wizard.enableNext();
        panel.initialize();
        return true;
    }

    public boolean end()
    {
    	if(!panel.verifySource(panel.getSourceLocation()))
    	{
        	new StsMessage(wizard.getModel().win3d, StsMessage.ERROR, "Source location not valid for this type, monitor not saved.");
        	wizard.gotoStep(this);
            return false;
        }   		
        return true;
    }

}

