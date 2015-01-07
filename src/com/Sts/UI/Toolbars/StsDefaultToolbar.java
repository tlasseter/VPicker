
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI.Toolbars;

import com.Sts.Actions.Wizards.Color.*;
import com.Sts.DB.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Icons.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.border.*;

public class StsDefaultToolbar extends StsToolbar implements StsSerializable
{
    public static final String NAME = "Default Toolbar";
    transient StsModel model = null;

    public static final String TIME_DEPTH = "timeDepth";

    private static final String SAVE_HELP_TARGET = "ui.toolbars.common";
    private static final String SAVEAS_HELP_TARGET = "ui.toolbars.common";
    private static final String CURSOR3D_HELP_TARGET = "ui.toolbars.common";
    private static final String HELPBUTTON_HELP_TARGET = "ui.toolbars.common";

    public static final String UNDO = "undo";

    transient private StsComboBoxFieldBean timeDepthBean = null;
    transient private StsButton spectrumButton = null;
    transient private StsButton preferencesButton = null;
    transient Icon timeIcon = StsIcon.createIcon("time.gif");
    transient Icon depthIcon = StsIcon.createIcon("depth.gif");

    public static final boolean defaultFloatable = true;

    public StsDefaultToolbar()
    {
        super(NAME);
    }

    public StsDefaultToolbar(StsWin3dBase win3d)
    {
        super(NAME);
        initialize(win3d);
    }

    public boolean initialize(StsWin3dBase win3d)
    {
        StsActionManager actionManager = win3d.getActionManager();
        Border border = BorderFactory.createEtchedBorder();
        setBorder(border);
        this.model = win3d.getModel();

        StsButton button;

        preferencesButton = new StsButton("userPreferences", "Set User Preferences", model, "savePreferences");
        add(preferencesButton);

        timeDepthBean = new StsComboBoxFieldBean(StsProject.class, "zDomainString", null, StsParameters.TD_STRINGS);
        timeDepthBean.setName(TIME_DEPTH);
//        td = new StsButton("time", "View Time or Depth", model.project, "toggleZDomain");
        add(timeDepthBean);

        spectrumButton = new StsButton("spectrumDefine", "Define a New Spectrum", actionManager, StsSpectrumAction.class);
        add(spectrumButton);
//        HelpManager.setContextSensitiveHelp(b, SPECTRUM_HELP_TARGET, model.win3d);
//        b = StsButton.constructHelpButton(win3d);
//        add(b);

        add(new StsButton(UNDO, "Undo Previous Edit Action", actionManager, "undoAction"));
        addSeparator();
        addCloseIcon(win3d);

        setMinimumSize();

        initializeTimeDepth(model.getProject());
        return true;
    }

    public StsComboBoxFieldBean getTimeDepthBean() { return timeDepthBean; }

    public void initializeTimeDepth(StsProject project)
    {
        timeDepthBean.setBeanObject(project);
        timeDepthBean.setSelectedItem(project.getZDomainString());
        timeDepthBean.setEditable(project.getZDomainSupported() == StsProject.TD_TIME_DEPTH);
    }


    public boolean forViewOnly()
    {
        remove(spectrumButton);
        remove(preferencesButton);
        return true;
    }
/*
    public void changeTimeDepthIcon(byte type)
    {
        if(type == StsProject.TD_DEPTH)
            td.setIcon(depthIcon);
        else
            td.setIcon(timeIcon);
    }
 */
/*
    public void setButtonsEnabled()
    {
        boolean enabled = true;

        StsList buttonList = getButtonList();
        int nButtons = buttonList.getSize();
        for (int i=0; i<nButtons; i++)
        {
			StsButton b = (StsButton)buttonList.getElement(i);
			b.setEnabled(enabled);
        }

        super.setButtonsEnabled();
    }
*/
}

