
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI.Toolbars;

import com.Sts.Actions.Boundary.*;
import com.Sts.DB.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;

import java.awt.*;

public class StsBuildBoundaryToolbar extends StsToolbar implements StsSerializable
{
    public static final String NAME = "Boundary Construction Toolbar";

    /** button filenames (also used as unique identifier button names) */
    public static final String AUTO_BOUNDARY = "autoBoundary";
    public static final String RECT_BOUNDARY = "rectangularBoundary";
    public static final String POLY_BOUNDARY = "polygonBoundary";

//    public static final String DELETE_BOUNDARY = "deleteBoundary";

    public StsBuildBoundaryToolbar()
    {
        super(NAME);
    }

    public StsBuildBoundaryToolbar(StsActionManager actionManager)
    {
        super(NAME);
        initialize(actionManager);
    }

    public boolean getDefaultFloatable() { return true; }

    public boolean initialize(StsActionManager actionManager)
    {
        Component[] components = getComponents();
        if(components != null && components.length > 0) return true; 
        add(new StsButton(AUTO_BOUNDARY, "Automatically build boundary around full grid.", actionManager, StsAutoBoundary.class));
        add(new StsButton(RECT_BOUNDARY, "Manually build row/col aligned rectangular boundary.", actionManager, StsRectangularBoundary.class));
        add(new StsButton(POLY_BOUNDARY, "Manually build arbitrary polygon boundary.", actionManager, StsPolygonBoundary.class));
 //       addSeparator();
//		add(new StsButton(DELETE_BOUNDARY, "Delete boundary", actionManager, StsDeleteBoundary.class));
        setMinimumSize();
        return true;
    }
}

