
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

public class StsDefineWellAlarm extends StsWizardStep
{
    StsDefineWellAlarmPanel panel;
    StsHeaderPanel header;

    public StsDefineWellAlarm(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineWellAlarmPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Define Alarms");
        header.setSubtitle("Define Well Alarm");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Alarms");
        header.setInfoText(wizardDialog,"(1) Select one well for alarm.\n" +
                                        "(2) Specify the radius from the selected well.\n" +
                                        "(3) Select whether inclusive or exclusive of defined zone.\n" +
        		                        "(4) Press the Next>> Button");
    }

    public boolean start()
    {
        panel.initialize();
        wizard.enableFinish();
        wizard.disableNext();
        return true;
    }

    public StsWell getSelectedWell() { return (StsWell)panel.wellListBean.getSelectedObjects()[0]; }
    public boolean end()
    {
        return true;
    }
}