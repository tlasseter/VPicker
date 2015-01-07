
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SuperGather;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

import java.awt.*;

public class StsStackParameter extends StsWizardStep
{
    StsStackParameterPanel panel;
    StsHeaderPanel header;

    public StsStackParameter(StsWizard wizard)
    {
        super(wizard);
        panel = new StsStackParameterPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 500));
        header.setTitle("Super Gather Definition");
        header.setSubtitle("Define Stacking Parameters");
        header.setInfoText(wizard.dialog,"(1) Specify the maximum Z data range to include in a stack.\n" +
                                         "(2) Specify the offset limit of CDP to be considered a neighbor of current gather.\n" +
                                         "      **** User will be allowed to stack entire line or neighbors as defined by this parameter. ****\n" +
                                         "(3) Press the Next>> Button to start picking velocity profiles.\n" +
                                         "      **** Velocity analysis toolbar will be added to the main window and used for velocity profile picking.\n");
    }

    public boolean start()
    {
        wizard.enableFinish();
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        panel.checkSaveSuperGatherProperties();
        return true;
    }
}

