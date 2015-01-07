
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI.Toolbars;

import com.Sts.Actions.*;
import com.Sts.DB.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.*;

import javax.swing.*;
import java.awt.*;

public class StsInformationToolbar extends StsToolbar implements StsSerializable
{
    public static final String NAME = "Information Toolbar";

    transient JLabel dirLabel = new JLabel();
    public static final String MALL_MAP = "mallMap";

    public static final boolean defaultFloatable = true;

    public StsInformationToolbar()
    {
        super(NAME);
    }

    public StsInformationToolbar(StsWin3dBase win3d)
    {
        super(NAME);
        initialize(win3d);
    }

    public boolean initialize(StsWin3dBase win3d)
    {
        StsMallMap mm = new StsMallMap();

        // Add Directional Text
        dirLabel.setHorizontalAlignment(SwingConstants.CENTER);
        dirLabel.setHorizontalTextPosition(SwingConstants.CENTER);
        dirLabel.setSize(new Dimension(30,20));
        dirLabel.setFont(new java.awt.Font("Dialog", 1, 12));
        dirLabel.setBorder(BorderFactory.createEtchedBorder());
        dirLabel.setText(" N ");
        this.add(dirLabel);

        addSeparator();
        add(new StsButton(MALL_MAP,"2D Map of the Project Area", mm, "init", win3d.getModel()));
        addSeparator();
        addCloseIcon(win3d);

        setMinimumSize();
        return true;
    }

    public void setText(String txt) {
        dirLabel.setText(txt);
        return;
    }
}

