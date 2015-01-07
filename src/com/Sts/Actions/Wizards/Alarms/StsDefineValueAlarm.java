
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Alarms;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.MicroseismicCorrelation.StsSelectSensorsPanel;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.StsSurface;
import com.Sts.DBTypes.StsWell;

public class StsDefineValueAlarm extends StsWizardStep
{
    StsDefineValueAlarmPanel panel;
    StsHeaderPanel header;

    public StsDefineValueAlarm(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineValueAlarmPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Define Alarms");
        header.setSubtitle("Define Value Alarm");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Alarms");
        header.setInfoText(wizardDialog,"(1) Specify the value for alarm.\n" +
                                        "(2) Specify a percentage threshold around value.\n" +
                                        "(3) Select whether inclusive or exclusive of defined value.\n" +
        		                        "(4) Press the Next>> Button");
    }

    public boolean start()
    {
        panel.initialize();
        wizard.enableFinish();
        wizard.disableNext();
        return true;
    }

    public boolean end()
    {
        return true;
    }
}