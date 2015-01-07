
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI.Toolbars;

import com.Sts.DB.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.*;

import javax.swing.*;
import javax.swing.border.*;

public class StsTestingToolbar extends StsToolbar implements StsSerializable
{
    transient StsWin3dBase window;
    transient StsToggleButton uiTestButton, loadTestButton, getLocationButton;

    public static final String NAME = "Intra-Family Toolbar";
    public static final boolean defaultFloatable = true;

    public static final String UI_TEST = "uiTest";
    public static final String LOAD_TEST = "loadTest";
    public static final String LOCATION_QUERY = "getMouseLocation";

    transient private StsModel model;

    public StsTestingToolbar()
     {
         super(NAME);
     }

    public StsTestingToolbar(StsWin3dBase win3d)
    {
        super(NAME);
        initialize(win3d);
    }

    public boolean initialize(StsWin3dBase win3d)
    {
        this.window = win3d;

        Border border = BorderFactory.createEtchedBorder();
        setBorder(border);

        uiTestButton = new StsToggleButton(UI_TEST, "Test User Interfaces", this, UI_TEST);
        loadTestButton = new StsToggleButton(LOAD_TEST, "Test Data Loading", this, LOAD_TEST);
        getLocationButton = new StsToggleButton(LOCATION_QUERY, "Get Mouse Location", this, LOCATION_QUERY);

        add(uiTestButton);
        add(loadTestButton);
        add(getLocationButton);


        addSeparator();
        addCloseIcon(window);

        model = window.getModel();

        setMinimumSize();
        return true;
    }

    public void uiTest()
    {
        model.getWorkflowPlugIn().testWorkflow(model);
    }

    public void loadTest()
    {
        ;
    }

    public void getMouseLocation()
    {
        ;
    }
}
