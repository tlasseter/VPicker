
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI.Toolbars;

import com.Sts.Actions.Wizards.Collaboration.*;
import com.Sts.Collaboration.*;
import com.Sts.DB.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.*;

import javax.swing.*;
import javax.swing.border.*;

public class StsCollaborationToolbar extends StsToolbar implements StsSerializable
{
    transient private StsModel model = null;
	transient private StsCollaboration collaboration = null;

    public static final String NAME = "Collaboration Toolbar";

    /** button filenames (also used as unique identifier button names) */
    public static final String CREATE_COLLAB = "startCollab";
    public static final String JOIN_COLLAB = "joinCollab";
    public static final String LEAVE_COLLAB = "leaveCollab";
    public static final String MODE_COLLAB = "modeCollab";

    public static final boolean defaultFloatable = true;

    public StsCollaborationToolbar()
    {
        super(NAME);
    }

    public StsCollaborationToolbar(StsWin3dBase win3d)
    {
        super(NAME);
        initialize(win3d);
    }

    public boolean getDefaultVisible() { return false; }

    public boolean initialize(StsWin3dBase win3d)
    {
        this.model = win3d.getModel();
        Border border = BorderFactory.createEtchedBorder();
        setBorder(border);

        add(new StsButton(CREATE_COLLAB,"Create a new collaboration session.", this, "createSession"));
        add(new StsButton(JOIN_COLLAB,"Join an existing collaboration.", this, "joinSession"));
        add(new StsButton(LEAVE_COLLAB,"Leave the current collaboration.", this, "leaveSession"));
        add(new StsButton(MODE_COLLAB,"Select collaboration mode.", this, "setMode"));

        addSeparator();
        addCloseIcon(win3d);

        setMinimumSize();
        revalidate();
        return true;
    }

	public void createSession()
	{
		StsCollaborationWizard wizard = new StsCollaborationWizard(model.mainWindowActionManager);
		wizard.start();
	}

    public void joinSession()
    {
		StsCollaborationWizard wizard = new StsCollaborationWizard(model.mainWindowActionManager);
		wizard.start();
        wizard.gotoSelectCollab();
    }
    public void leaveSession()
    {
        StsMessageFiles.infoMessage("Leaving the current collaboration session. - Not Enabled");
        return;
    }

    public void setMode()
    {
        StsMessageFiles.infoMessage("Change the mode of the current user in the current collaboration session");
        StsCollaborationWizard wizard = new StsCollaborationWizard(model.mainWindowActionManager);
        wizard.start();
        wizard.gotoCollabMode();
        return;
    }
}

