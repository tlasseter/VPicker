
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.VelocityAnalysis;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

import java.awt.*;

public class StsIdentifyVelocityModels extends StsWizardStep
{
    StsIdentifyVelocityModelsPanel panel;
    StsHeaderPanel header;

    public StsIdentifyVelocityModels(StsWizard wizard)
    {
        super(wizard);
        panel = new StsIdentifyVelocityModelsPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 500));
        header.setTitle("Velocity Analysis Definition");
        header.setSubtitle("Identify Velocity Models");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VelocityAnalysis");
        header.setInfoText(wizard.dialog,"(1) Select a velocity model or imported seismic volume if dataset is migrated.\n" +
                                         "      **** Default is none, should select previous velocity model if residual picking. ****\n" +
                                         "(2) Specify the units of the selected velocity model or seismic volume.\n" +
                                         "(3) Select a velocity model to use/edit for NMO.\n" +
                                         "      **** Default is a new velocity model. ****\n" +
                                         "(4) Press the Next>> Button.");
    }

    public boolean start()
    {
       // wizard.enableFinish();  //force user to use "Next" button so that velocity model gets created
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        return true;
    }
}

