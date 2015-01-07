
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Surfaces;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

public class StsDefineColumns extends StsWizardStep
{
    StsDefineColumnsPanel panel;
    StsHeaderPanel header;

    public StsDefineColumns(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineColumnsPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
//        panel.setPreferredSize(new Dimension(500, 400));
        header.setTitle("Surface Definition");
        header.setSubtitle("Define the Selected File Columns");
        header.setInfoText(wizardDialog,"(1) Increase the number of columns to appropriate number.\n" +
                           " **** # columns is the number of values in file per grid node. ****\n" +
                           "(2) Adjust column assignments as required.\n" +
                           "(3) Once columns have been defined, press the Next>> Button.");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Surfaces");
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

