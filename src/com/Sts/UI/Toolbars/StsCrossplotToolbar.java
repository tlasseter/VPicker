
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI.Toolbars;

import com.Sts.Actions.Crossplot.*;
import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import javax.swing.*;
import javax.swing.border.*;

public class StsCrossplotToolbar extends StsToolbar implements StsSerializable
{
    public static final String NAME = "Cross Plot Toolbar";

    /** button filenames (also used as unique identifier button names) */
//    public static final String NEW_CROSSPLOT = "newCrossplot";
    public static final String CREATE_POLYGON = "newPolygons";
    public static final String ATTRIBUTE_DENSITY = "attributeOrDensity";
    public static final String EDIT_COLORS = "editColors";
    public static final String SELECT_POINT_OR_POLYGON = "pointOrPolygon";

    public static final boolean defaultFloatable = true;

    transient StsModel model;
    transient StsToggleButton densityToggleButton;
    
    public StsCrossplotToolbar()
    {
        super(NAME);
    }

    public StsCrossplotToolbar(StsWin3dBase win3d)
	{
		super(NAME);
        initialize(win3d);
    }

    public boolean initialize(StsWin3dBase win3d)
    {
        StsActionManager actionManager = win3d.getActionManager();
        this.model = win3d.getModel();
        Border border = BorderFactory.createEtchedBorder();
        setBorder(border);

        add(new StsToggleButton(CREATE_POLYGON, "Create Polygons", actionManager, StsXPolygonAction.class));
        densityToggleButton = new StsToggleButton(ATTRIBUTE_DENSITY, "Color by density or average sample value.", model, "win3dDisplayAll");
        densityToggleButton.addIcons("attributeXP", "densityXP");
        add(densityToggleButton);

        StsToggleButton toggleButton = new StsToggleButton(SELECT_POINT_OR_POLYGON, "Select point (on 2D View only) or polygon creation (on 2D or XPlot views).");
        toggleButton.addIcons("point", "polygon");
        add(toggleButton);
        toggleButton.setEnabled(false);
        StsColorItemComboBoxFieldBean colorItemComboBoxBean = new StsColorItemComboBoxFieldBean(StsCrossplot.class, "typeItem", null);
        colorItemComboBoxBean.setName(EDIT_COLORS);
        add(colorItemComboBoxBean);
//        StsEditColorListPanel colorListPanel = new StsEditColorListPanel(EDIT_COLORS, StsXPolygonAction.class, "setPolygonColor", new Class[] { Color.class, String.class } );
//        add(colorListPanel);
//        colorListPanel.setEnabled(false);
        addSeparator();
        addCloseIcon(model.win3d);

        setMinimumSize();
        return true;
    }
    public void enableDensityToggle(boolean val)
    {
    	densityToggleButton.setEnabled(val);
    }

    public boolean validateState()
    {
        return model.getNObjects(StsCrossplot.class) > 0;
    }
}

