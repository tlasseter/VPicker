
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.ClusterAnalysis;

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
        header.setSubtitle("Selecting Sensors for Cluster Analysis");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#ClusterAnalysis");
        header.setInfoText(wizardDialog,"(1) Select one or more sensors for analysis.\n" +
        		"  **** Use the shift or ctrl keys to select multiple sensors ****\n" +
        		"(2) Press the Next>> Button");
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

