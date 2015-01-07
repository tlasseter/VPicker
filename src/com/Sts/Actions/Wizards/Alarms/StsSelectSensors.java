
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Alarms;

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
        header.setTitle("Define Alarms");
        header.setSubtitle("Selecting Sensors");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Alarms");
        header.setInfoText(wizardDialog,"(1) Select one or more sensors for alarm.\n" +
                                        "(2) Select the alarm type.\n" +
                                        "(3) Select the type of queue wanted when alarm is triggered.\n" +
                                        "(4) Select the sound to be played when alarm is triggered. \n" +
        		                        "(5) Press the Next>> Button");
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