
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Monitor;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

public class StsDefineSensor extends StsWizardStep
{
    StsDefineSensorPanel panel;
    StsHeaderPanel header;

    public StsDefineSensor(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineSensorPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Create a Real-time Monitor");
        header.setSubtitle("Define a New Sensor to Monitor");
        header.setInfoText(wizardDialog, "(1) Supply the sensor name.\n" +
                           "(2) Specify whether sensor is static or dynamic.\n" +
                           "(3) If static, specify the location.\n" +
                           "(4) If dynamic, specify whether the positions be absolute or relative.\n" +
                           "(5) If relative, specify an absolute origin.\n" +
                           "(6) Select whether both a time and date are provided.\n" +
                           "(7) If only time, specify a start date.\n" +
                           "(8) Press the Next>> Button.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Monitor");
    }

    public boolean start()
    {
        wizard.enableFinish();        
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        return true;
    }

}

