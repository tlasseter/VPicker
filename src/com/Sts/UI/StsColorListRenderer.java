
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

// class to use with color list item for building color icon lists

package com.Sts.UI;

import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;


public class StsColorListRenderer extends JLabel implements ListCellRenderer
{
    public Component getListCellRendererComponent(JList list, Object value, int index,  boolean isSelected, boolean hasFocus)
    {
        /*
        if(value == null) 
        {
            value = list.getModel().getElementAt(0);
        }
        */
        try
        {
            StsColorListItem item = (StsColorListItem)value;
            setOpaque(true);
            setIcon(item.getIcon());
            setText(item.getName());
            if (isSelected)
            {
                setBackground(Color.black);
                setForeground(Color.white);
            }
            else // not selected
            {
                setBackground(Color.white);
                setForeground(Color.black);
            }
            return this;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "getListCellRendererComponent", e);
            return null;
        }
    }
}