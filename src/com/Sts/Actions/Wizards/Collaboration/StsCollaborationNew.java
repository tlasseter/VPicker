
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

public class StsCollaborationNew extends StsWizardStep
{
    StsCollaborationNewPanel panel;
    StsHeaderPanel header;

    public StsCollaborationNew(StsCollaborationWizard wizard)
    {
      super(wizard, new StsCollaborationNewPanel(wizard), null, new StsHeaderPanel());
      panel = (StsCollaborationNewPanel) getContainer();
      panel.setPreferredSize(new Dimension(350,300));

      header = (StsHeaderPanel) getHdrContainer();
      header.setTitle("Collaboration Configuration");
      header.setSubtitle("New Session");
      header.setLogo("AdvancedRealityLogo.gif");
      header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Collaboration");
      header.setInfoText(wizardDialog,"(1) Specify the collaboration session name.\n" +
                           "(2) Specify the collaboration session description so that others can identify it.\n" +
                           "(3) Input the session password which will need to be input by all participants.\n" +
                           "(4) Press the Next >> Button to proceed to next screen.\n");
    }

    public boolean start()
    {
//		StsCollaboration.getInstance(model, entry);
        return true;
    }

    public boolean end()
    {
        return true;
    }

}

