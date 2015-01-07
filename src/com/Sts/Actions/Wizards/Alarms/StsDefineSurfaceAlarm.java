
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

public class StsDefineSurfaceAlarm extends StsWizardStep
{
    StsDefineSurfaceAlarmPanel panel;
    StsHeaderPanel header;

    public StsDefineSurfaceAlarm(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineSurfaceAlarmPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Define Alarms");
        header.setSubtitle("Define Surface Alarm");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Alarms");
        header.setInfoText(wizardDialog,"(1) Select one surface for alarm.\n" +
                                        "(2) Specify whether depth is above or below reference surface.\n" +
                                        "(3) Specify the offset from the selected surface.\n" +
        		                        "(4) Press the Next>> Button");
    }

    public boolean start()
    {
        panel.initialize();
        wizard.enableFinish();
        wizard.disableNext();        
        return true;
    }

    public StsSurface getSelectedSurface() { return (StsSurface)panel.surfaceListBean.getSelectedObjects()[0]; }
    public boolean end()
    {
        return true;
    }
}