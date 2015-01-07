
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI.Toolbars;

import com.Sts.Actions.Build.*;
import com.Sts.Actions.Edit.*;
import com.Sts.DB.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.*;

public class StsBuildFrameToolbar extends StsToolbar implements StsSerializable
{
    public static final String NAME = "Fault Section Toolbar";

    /** button filenames (also used as unique identifier button names) */
    public static final String IMPORT_FAULT_CUTS = "importFaultCuts";
    public static final String EXPORT_FAULT_CUTS = "exportFaultCuts";
    public static final String PICK_FAULT_ON_CURSOR = "faultOnCursor";
    public static final String PICK_FAULT_ON_SECTION = "faultOnSection";
    public static final String PICK_FAULT_ON_SURFACES = "faultOnSurfaces";
    public static final String BUILD_FAULT_FROM_LINES = "faultSectionFrom2Lines";
    public static final String BUILD_FAULT_FROM_SURFACE_AND_LINES = "sectionEdgeOnSurface";
    public static final String INTERSECT_SECTION_EDGES = "intersectSectionEdges";
    public static final String PICK_EDGE_ON_SURFACE = "edgeOnSurface";
     public static final String PICK_EDGE_ON_CURSOR = "edgeOnCursor";
    public static final String BUILD_FAULTS_FROM_SURFACE_CUTS = "faultSection";

    public static final String EDIT_FAULT_ON_CURSOR = "editFaultOnCursor";
    public static final String EDIT_FAULT_ON_SECTION = "editFaultOnSection";
    public static final String EDIT_FAULT_EDGE = "editFaultEdge";
    //public static final String PICK_INTERIOR_FAULT_LINE = "pickInteriorFaultLine";
    public static final String ROTATE_FAULT_SECTION = "rotateFaultLine";
    public static final String PICK_FAULT_TRIM_EDGES = "trimFaultSection";

    public static final String DELETE_WELL_LINE = "deleteWellLine";
    public static final String DELETE_SECTION = "deleteSection";
    public static final String UNDO = "undo";

    public static final boolean defaultFloatable = true;

//    private static final Object[] FAULT_WELL = new Object[] { new Integer(StsWell.FAULT) };
//    private static final Object[] FAULT_SECTION = new Object[] { new Integer(StsSection.FAULT) };

    public StsBuildFrameToolbar()
    {
        super(NAME);
    }

    public StsBuildFrameToolbar(StsWin3dBase win3d)
    {
        super(NAME);
        initialize(win3d);
    }

    public boolean initialize(StsWin3dBase win3d)
    {
        StsActionManager actionManager = win3d.getActionManager();
        // icons
        add(new StsButton(PICK_FAULT_ON_CURSOR, "Pick fault section line on X or Y cursor plane", actionManager, StsLineOnCursor.class));
        add(new StsButton(PICK_FAULT_ON_SECTION, "Pick fault section line", actionManager, StsWellOnSection.class));
//        add(new StsButton(PICK_FAULT_ON_SURFACES, "Pick fault from points on surfaces", actionManager, StsFaultOnSurfaces.class));
//        add(new StsButton(BUILD_FAULT_FROM_SURFACE_AND_LINES, "Digitize fault cut on horizon between fault lines", actionManager, StsSectionEdgeOnSurface.class));
 		add(new StsButton(PICK_EDGE_ON_SURFACE, "Digitize fault cuts on horizon", actionManager, StsEdgeOnSurface.class));
        add(new StsButton(PICK_EDGE_ON_CURSOR, "Digitize fault cuts on slice cursor", actionManager, StsEdgeOnCursor.class));
        add(new StsButton(BUILD_FAULT_FROM_LINES, "Pick fault section from sticks", actionManager, StsAddSection.class));
        add(new StsButton(INTERSECT_SECTION_EDGES, "Automatically intersect/split sectionEdges", actionManager, StsIntersectSectionEdges.class));
        addSeparator();

		add(new StsButton(BUILD_FAULTS_FROM_SURFACE_CUTS, "Build fault sections from horizon cuts", actionManager, StsBuildSectionsFromEdges.class));
/*
        if (Main.NON_VERTICAL_FAULTS)
        {
            addSeparator();

            addButton(EDIT_FAULT_ON_CURSOR, "Edit fault section lines on X, Y or Z cursor plane", StsEditFaultLineOnCursor.class);
            addButton(EDIT_FAULT_ON_SECTION, "Edit fault section lines", StsEditFaultLineOnSection.class);
            addButton(EDIT_FAULT_EDGE, "Edit fault cut line", StsEditFaultEdge.class);
        }

        addSeparator();

        if (Main.NON_VERTICAL_FAULTS)
        {
            addButton(PICK_INTERIOR_FAULT_LINE, "Pick interior fault lines on current fault section", StsPickInteriorFaultSectionLine.class);
            addButton(ROTATE_FAULT_SECTION, "Rotate current fault section", StsRotateFaultLine.class);
        }
		    addButton(PICK_FAULT_TRIM_EDGES, "Trim current fault section", StsPickFaultSectionTrimEdge.class);
*/
        addSeparator();

//		add(new StsButton(DELETE_WELL_LINE, "Delete fault line", actionManager, StsDeleteLine.class));
		add(new StsButton(DELETE_SECTION, "Delete fault section", actionManager, StsDeleteSection.class));
/*
        addSeparator();

		addButton(UNDO, "Undo previous edit action", actionManager, "undoAction", null);
*/
        setMinimumSize();
        return true;
    }
}

