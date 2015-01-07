
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Monitor;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

public class StsSelectSensor extends StsWizardStep
{
    StsSelectSensorPanel panel;
    StsHeaderPanel header;

    public StsSelectSensor(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSelectSensorPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Create a Real-time Monitor");
        header.setSubtitle("Select Sensor to Monitor");
        header.setInfoText(wizard.dialog, "(1) Select one sensor from the list.\n" +
        		"                 OR\n" +
        		"(2) Press the New Sensor... Button to create a new sensor.\n" +
            "(3) Press the Next Button.");
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

