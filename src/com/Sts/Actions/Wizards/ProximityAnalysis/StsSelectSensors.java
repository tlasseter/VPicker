
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.ProximityAnalysis;

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
        header.setTitle("Select Sensors");
        header.setSubtitle("Selecting Sensors for Analysis");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#ProximityAnalysis");
        header.setInfoText(wizardDialog,"(1) Select a well to analyze against the sensor events.\n" +
        		"(2) Select one or more sensors for analysis.\n" +
        		"(3) Press the Next>> Button");
    }

    public boolean start()
    {
        panel.initialize();
        return true;
    }
    
    public Object[] getSelectedSensors()
    {
    	return panel.getSelectedSensors();
    }   

    public boolean end()
    {
        return true;
    }
}

