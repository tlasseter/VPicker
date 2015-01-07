
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

public class StsMediaToolbar extends StsToolbar implements StsSerializable
{
    transient StsWin3dBase window;
    transient StsButton captureDesktopBtn, captureWindowBtn, captureGraphicBtn;
    transient StsToggleButton captureMovieBtn;

    public static final String NAME = "Media Toolbar";
    public static final boolean defaultFloatable = true;

    public StsMediaToolbar()
     {
         super(NAME);
     }

    public StsMediaToolbar(StsWin3dBase win3d)
    {
        super(NAME);
        initialize(win3d);
    }

    public boolean initialize(StsWin3dBase win3d)
    {
        this.window = win3d;

        Border border = BorderFactory.createEtchedBorder();
        setBorder(border);

        captureDesktopBtn = new StsButton("Desktop", "Capture an image of the Desktop", window, "outputImage", new Integer(StsWin3dBase.DESKTOP));
        captureWindowBtn = new StsButton("Window", "Capture an image of the Window", window, "outputImage", new Integer(StsWin3dBase.WINDOW));
        captureGraphicBtn = new StsButton("Graphic", "Capture an image of the Graphics", window, "outputImage", new Integer(StsWin3dBase.GRAPHIC));
        captureMovieBtn = new StsToggleButton("movie", "Capture movie frames", window, "outputMovie");
        captureMovieBtn.addIcons("movieSelect", "movieDeselect");

        add(captureDesktopBtn);
        add(captureWindowBtn);
        add(captureGraphicBtn);
        add(captureMovieBtn);

        addSeparator();
        addCloseIcon(window);

        setMinimumSize();
        return true;
    }

    public void resetMovieBtn()
    {
    	captureMovieBtn.setSelected(window.captureMovie());
    }

    public boolean forViewOnly() { return true; }
}
