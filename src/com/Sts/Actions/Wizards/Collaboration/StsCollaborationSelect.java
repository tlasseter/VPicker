
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

public class StsCollaborationSelect extends StsWizardStep
{
    StsCollaborationSelectPanel panel;
    StsHeaderPanel header;

    public StsCollaborationSelect(StsWizard wizard)
    {
      super(wizard, new StsCollaborationSelectPanel(), null, new StsHeaderPanel());
      panel = (StsCollaborationSelectPanel) getContainer();
      panel.setPreferredSize(new Dimension(350,300));

      header = (StsHeaderPanel) getHdrContainer();
      header.setTitle("Collaboration Configuration");
      header.setSubtitle("Select Active Session");
      header.setLogo("AdvancedRealityLogo.gif");
      header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Collaboration");
      header.setInfoText(wizardDialog,"(1) Select the session that you wish to join.\n" +
                           "(2) Press the Next >> Button to classInitialize the session and join it.\n");

    }

    public boolean newSession()
    {
        return panel.newSessionChk.isSelected();
    }

    public boolean start()
    {
        enableNext();
        disableFinish();
        disablePrevious();
        return true;
    }

    public boolean end()
    {
//        ((StsCollaborationWizard)wizard).setCollaborationEntry(panel.currentEntry);
//        ((StsCollaborationWizard)wizard).gotoCollabStart();
//        wizard.dialog.setVisible(false);
//        wizard.dialog.dispose();
        return true;
    }
}

