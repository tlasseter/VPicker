
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

public class StsSeismicFlatteningToolbar extends StsToolbar implements StsSerializable
{
    transient private StsWin3dBase win3d;
    transient private StsModel model;
    transient private StsGLPanel3d glPanel3d;
    transient private ButtonGroup buttonGroup = new ButtonGroup();

    public static final String NAME = "Seismic Flattening Toolbar";
    public static final boolean defaultFloatable = true;

    /** button filenames (also used as unique identifier button names) */
    public static final String FLATTEN = "seismicFlatten";
    public static final String UNFLATTEN = "seismicUnflatten";

    public StsSeismicFlatteningToolbar()
     {
         super(NAME);
     }

	public StsSeismicFlatteningToolbar(StsWin3dBase win3d)
    {
        super(NAME);
       initialize(win3d);
    }

    public boolean initialize(StsWin3dBase win3d)
    {
        this.win3d = win3d;
        this.model = win3d.getModel();
        this.glPanel3d = win3d.getGlPanel3d();
 
        addButton(FLATTEN, "Flatten the 2D cursor", "flattenSeismic", true);
        addButton(UNFLATTEN, "Unflatten the 2D cursor", "unflattenSeismic", true);

        addSeparator();
        addCloseIcon(win3d);

        setMinimumSize();
        return true;
    }

    private AbstractButton addButton(String name, String tip, String methodName, boolean addToGroup)
    {
        StsToggleButton toggleButton = new StsToggleButton(name, tip, this, methodName);
        buttonGroup.add(toggleButton);
        add(toggleButton);
        return toggleButton;
    }

    /** Flatten the current seismic slice */
    public void flattenSeismic()
    {
        glPanel3d.setMouseMode(StsCursor.ZOOM);
        return;
    }

    /** Unflatten the current seismic slice */
    public void unflattenSeismic()
    {
        glPanel3d.setMouseMode(StsCursor.PAN);
        return;
    }


}

