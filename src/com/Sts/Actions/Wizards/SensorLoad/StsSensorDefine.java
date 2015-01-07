
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SensorLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

import java.awt.*;

public class StsSensorDefine extends StsWizardStep
{
    StsSensorDefinePanel panel;
    StsHeaderPanel header;

    public StsSensorDefine(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSensorDefinePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 400));
        header.setTitle("Sensor Load");
        header.setSubtitle("Import Ascii Sensor File");
        header.setInfoText(wizardDialog, "(1) Press the Select button to find the sensor file to import.\n" +
                           " **** A column named Time is required, if dynamic positioned X, Y and Z are also required.\n" +
                           "(2) If the time column is not absolute, enter the starting time (dd-mm-yy hh:mm:ss.s).\n" +
                           "(3) Enter the column numbers of the attributes you wish to import. \n" +
                           " **** First column is 1. Attributes are delimited with commas ****\n" +
                           "(4) Select whether position is Dynamic (from file) or Static (user supplied).\n" +
                           "(5) If Static, specify the absolute coordiantes of the sensor. \n" +
                           "(6) Press the Next>> Button.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#SensorLoad");
    }

    public boolean start()
    {
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        panel.determineTimeType();
        return true;
    }

}

