
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Collaboration;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

import java.awt.*;

public class StsCollaborationMode extends StsWizardStep
{
    StsCollaborationModePanel panel;
    StsHeaderPanel header;

    public StsCollaborationMode(StsWizard wizard)
    {
        super(wizard, new StsCollaborationModePanel(), null, new StsHeaderPanel());
        panel = (StsCollaborationModePanel) getContainer();
        panel.setPreferredSize(new Dimension(350,300));

        header = (StsHeaderPanel) getHdrContainer();
        header.setTitle("Collaboration Configuration");
        header.setSubtitle("Operating Mode");
        header.setLogo("AdvancedRealityLogo.gif");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Collaboration");
        header.setInfoText(wizardDialog,"(1) Select the desired mode of operation.\n" +
                           "(2) Press the Next>> Button to accept the mode setting\n");
    }

    public boolean start()
    {
        disableFinish();
        disablePrevious();
        enableNext();
        return true;
    }

    public boolean end()
    {
        wizard.dialog.setVisible(false);
        wizard.dialog.dispose();
        return true;
    }
}

