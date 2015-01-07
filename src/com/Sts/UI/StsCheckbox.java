

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import javax.swing.*;
import java.awt.*;

public class StsCheckbox extends JCheckBox
{
	Font defaultFont = null;
	public StsCheckbox()
	{
    	super();
	}

    public StsCheckbox(String label, String tooltip)
    {
        super();
        setText(label);
        setToolTipText(tooltip);
	}

    public void setGrayed(boolean b)
    {
        if( defaultFont == null )
            defaultFont = getFont();

    	if( b )
	    	setFont(new Font("Dialog", Font.ITALIC, defaultFont.getSize()));
        else
	    	setFont(defaultFont);
    }
    public void setState(boolean selected, boolean grayed)
    {
    	setSelected(selected);
        setGrayed(grayed);
    }
}


