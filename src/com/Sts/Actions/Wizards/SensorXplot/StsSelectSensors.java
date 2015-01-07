
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SensorXplot;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

public class StsSelectSensors extends StsWizardStep
{
    StsSelectSensorsPanel panel;
    StsHeaderPanel header;

    public StsSelectSensors(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSelectSensorsPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Crossplot Sensors");
        header.setSubtitle("Select Sensor to Plot");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#SensorXplot");
        header.setInfoText(wizardDialog,"(1) Select a sensor for crossplot analysis.\n" +
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
}

