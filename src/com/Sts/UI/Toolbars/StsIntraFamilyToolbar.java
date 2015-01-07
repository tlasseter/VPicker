
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI.Toolbars;

import com.Sts.DB.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.*;

import javax.swing.*;
import javax.swing.border.*;

public class StsIntraFamilyToolbar extends StsToolbar implements StsSerializable
{
    transient StsWin3dBase window;
    transient StsToggleButton mouseTrackButton/*, cursorTrackButton*/, userTrackButton;

    public static final String NAME = "Intra-Family Toolbar";
    public static final String MOUSE_TRACK = "trackMouse";
//    public static final String CURSOR_TRACK = "trackCursor";
    public static final String USER_TRACK = "trackUser";

    public static final boolean defaultFloatable = true;

    public StsIntraFamilyToolbar()
     {
         super(NAME);
     }

    public StsIntraFamilyToolbar(StsWin3dBase win3d)
    {
        super(NAME);
        initialize(win3d);
    }

    public boolean initialize(StsWin3dBase win3d)
    {
        this.window = win3d;

        Border border = BorderFactory.createEtchedBorder();
        setBorder(border);

        mouseTrackButton = new StsToggleButton(MOUSE_TRACK, "Track the Mouse Position between windows", this, "trackMouse");
        mouseTrackButton.addIcons(MOUSE_TRACK + "Select", MOUSE_TRACK + "Deselect");
//        cursorTrackButton = new StsToggleButton(CURSOR_TRACK, "Track the Cursor Position in all windows", this, "trackCursor");
        userTrackButton = new StsToggleButton(USER_TRACK, "Lock 3d views or center on cursors for 2d.", this, "trackUser");
        userTrackButton.addIcons(USER_TRACK + "Select", USER_TRACK + "Deselect");
        add(mouseTrackButton);
//        add(cursorTrackButton);
        add(userTrackButton);

        userTrackButton.setSelected(true);
        addSeparator();
        addCloseIcon(window);

        setMinimumSize();
        return true;
    }
/*
    public void trackCursor()
    {
 //       boolean isCursorTracked = cursorTrackButton.isSelected();
        StsWindowFamily windowFamily =
        StsWindowFamily.ViewIterator = StsWindowFamily.getViewIterator(windowFamily);
        StsWin3dBase[] familyWindows = window.getFamilyWindows();
        int nWindows = familyWindows.length;
        for(int n = 0; n < nWindows; n++)
        {
            familyWindows[n].getSingleViewGlPanel3d().repaint();
        }
    }
*/
    public void trackUser()
    {
        boolean locked = userTrackButton.isSelected();

        StsWin3dBase[] familyWindows = window.getFamilyWindows();
        int nWindows = familyWindows.length;
        for(int n = 0; n < nWindows; n++)
        {
            familyWindows[n].isLocked = locked;
        }
    }

    public void trackMouse()
    {
        boolean isMouseTracking = mouseTrackButton.isSelected();

        StsWin3dBase[] familyWindows = window.getFamilyWindows();
        int nWindows = familyWindows.length;
        for(int n = 0; n < nWindows; n++)
        {
            familyWindows[n].setMouseTracking(isMouseTracking);
        }
    }
    public boolean forViewOnly() { return true; }
     
}
