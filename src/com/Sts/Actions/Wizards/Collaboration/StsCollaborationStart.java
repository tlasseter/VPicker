//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Collaboration;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.Collaboration.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;

import java.awt.*;

public class StsCollaborationStart extends StsWizardStep implements Runnable
{
    private StsStatusPanel status;
    private StsHeaderPanel header;

    public StsCollaborationStart(StsCollaborationWizard wizard)
    {
        super(wizard, new StsStatusPanel(), null, new StsHeaderPanel());
        status = (StsStatusPanel) getContainer();
        status.setPreferredSize(new Dimension(350, 300));

        header = (StsHeaderPanel) getHdrContainer();
        header.setTitle("Collaboration Configuration");
        header.setSubtitle("Start Collaboration");
        header.setLogo("AdvancedRealityLogo.gif");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Collaboration");
        header.setInfoText(wizardDialog, "    **** No user input required. ****\n");
    }

    public boolean start()
    {
        status.setTitle("Starting Collaboration:");
		StsCollaborationEntry collaborationEntry = ((StsCollaborationWizard)wizard).collaborationEntry;
		StsCollaboration collaboration = StsCollaboration.construct(model, collaborationEntry);
        return true;
    }

    public void run()
    {

        startCollaboration();
        if (success)
        {
            status.setText("Collaboration started.");
            logMessage("Collaboration started.");
        }
        else
        {
            status.setText("Unable to start collaboration.");
            logMessage("Unable to start collaboration.");
        }

        status.sleep(1000);
        wizard.enableFinish();
        disablePrevious();
//        wizard.finish();
    }

    public boolean end()
    {
        wizard.dialog.setVisible(false);
        wizard.dialog.dispose();
        return success;
    }

    private void startCollaboration()
    {
		StsCollaborationEntry entry = ((StsCollaborationWizard)wizard).collaborationEntry;
		if(entry == null)
		{
			StsMessageFiles.errorMessage("No collaboration session is available.");
			return;
		}
		StsCollaboration collaboration = StsCollaboration.construct(model, entry);
        if(collaboration != null)
            success = true;
    }
}
