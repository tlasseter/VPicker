
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

public class StsSensorDefineColumns extends StsWizardStep
{
    StsSensorDefineColumnsPanel panel;
    StsHeaderPanel header;

    public StsSensorDefineColumns(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSensorDefineColumnsPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Sensor Definition");
        header.setSubtitle("Define the Selected File Columns");
        header.setInfoText(wizardDialog,"(1) Select the desired files from the list.\n" +
                            "(2) Increase the number of columns to appropriate number.\n" +
                           " **** # columns is the number of values in file per event node. ****\n" +
                           "(3) Adjust required column assignments as required.\n" +
                            " **** If there is a header row, column numbers will be automatically assigned.\n" +
                            " **** All other columns in file will be loaded as attributes.\n" +
                           "(4) Once columns have been defined, press the Next>> Button.");
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