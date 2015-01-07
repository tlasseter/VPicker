
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SensorLoad;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.Actions.Wizards.WizardHeaders.StsHeaderPanel;

public class StsSensorDefineRows extends StsWizardStep
{
    StsSensorDefineRowsPanel panel;
    StsHeaderPanel header;

    public StsSensorDefineRows(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSensorDefineRowsPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Sensor Definition");
        header.setSubtitle("Define the Number of Rows in the Header");
        header.setInfoText(wizardDialog,"(1) Select the desired files from the list.\n" +
                            "(2) Increase the number of header rows.\n" +
                           " **** Header rows will disappear as number is increased. ****\n" +
                           "(3) Once all header rows have disappeared, press the Next>> Button.");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Sensors");
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