package com.Sts.UI.Toolbars;

import com.Sts.Actions.Edit.*;
import com.Sts.Actions.Sections.*;
import com.Sts.DB.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.*;

public class StsBuildVerticalSectionsToolbar extends StsToolbar implements StsSerializable
{
    public static final String NAME = "Fault Section Toolbar";

    /** button filenames (also used as unique identifier button names) */
    public static final String BUILD_VERTICAL_SECTIONS_ON_CURSORS = "verticalSectionOnCursors";
    public static final String BUILD_VERTICAL_SECTIONS_ON_SURFACE = "verticalSectionsOnSurface";

    public static final String DELETE_SECTION = "deleteSection";
    public static final String UNDO = "undo";

    public static final boolean defaultFloatable = true;

    public StsBuildVerticalSectionsToolbar()
    {
        super(NAME);
    }

    public StsBuildVerticalSectionsToolbar(StsWin3dBase win3d)
    {
        super(NAME);
        initialize(win3d);
    }

    public boolean initialize(StsWin3dBase win3d)
    {
        StsActionManager actionManager = win3d.getActionManager();
        // icons
//        add(new StsButton(BUILD_VERTICAL_SECTIONS_ON_CURSORS, "Build vertical section between vertical lines", actionManager, StsVerticalSectionsFromLines.class));
        add(new StsButton(BUILD_VERTICAL_SECTIONS_ON_SURFACE, "Build vertical sections on surface", actionManager, StsVerticalSectionsOnSurface.class));
        addSeparator();
        add(new StsButton(DELETE_SECTION, "Delete section", actionManager, StsDeleteSection.class));
        setMinimumSize();
        return true;
    }
}
