
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Collaboration;

import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;

public class StsCollaborationWizard extends StsWizard
{
	public String name = "test";
	public String ipAddress = "maui";
	public String password = "password";
	public int port = 9010;

	public StsCollaborationEntry collaborationEntry = null;
	private StsCollaborationNew newCollab = new StsCollaborationNew(this);
    private StsCollaborationSelect selectCollab = new StsCollaborationSelect(this);
    private StsCollaborationMode collabMode = new StsCollaborationMode(this);
    private StsCollaborationStart collabStart = new StsCollaborationStart(this);

    private StsWizardStep[] mySteps =
    {
        newCollab,
        selectCollab,
        collabMode,
        collabStart
    };

    public StsCollaborationWizard(StsActionManager actionManager)
    {
        super(actionManager);

        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
    	dialog.setTitle("Collaboration Configuration");
		collaborationEntry = new StsCollaborationEntry("test", "testDescription", null);
    	return super.start();
    }

    public void previous()
    {
        if(currentStep == collabMode)
        {
            if(selectCollab.newSession())
				gotoStep(newCollab);
            else
				gotoPreviousStep();
        }
        else
            gotoPreviousStep();
    }

    public void next()
    {
        if(currentStep == selectCollab)
        {
			collaborationEntry = selectCollab.panel.currentEntry;
            gotoStep(collabStart);
        }
        else if(currentStep == newCollab)
		{
//			collaborationEntry = newCollab.panel.currentEntry;
			gotoStep(collabStart);
		}
		else
            gotoNextStep();
    }

    public void gotoSelectCollab()
    {
        gotoStep(selectCollab);
    }

    public void gotoNewCollab()
    {
        gotoStep(newCollab);
    }

    public void gotoCollabMode()
    {
        gotoStep(collabMode);
    }

    public void gotoCollabStart()
    {
        gotoStep(collabStart);
    }

    public void setCollaborationEntry(StsCollaborationEntry entry) { collaborationEntry = entry; }
	public void setName(String name) {this.name = name; }
	public String getName() { return name; }

	public void setIpAddress(String ipAddress) { this.ipAddress = ipAddress; }
	public String getIpAddress() { return ipAddress; }

	public void setPassword(String password) { this.password = password; }
	public String getPassword() { return password; }

	public void setPort(int port) { this.port = port; }
	public int getPort() { return port; }

    public void finish()
    {
//SAJ    	next();
        super.finish();
    }

    public boolean end()
    {
        return super.end();
    }
}
