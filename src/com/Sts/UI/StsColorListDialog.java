
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

// dialog box to set a text field

package com.Sts.UI;

import com.Sts.DBTypes.*;

import javax.swing.*;

public class StsColorListDialog extends StsListDialog
{
    static public final StsColor PROTOTYPE_COLOR = StsColor.BLACK;

    private int iconHeight;
    private int iconWidth;
    private int borderWidth;

	public StsColorListDialog(JFrame parent, String title, String text, StsColor[] colors, String[] names)
	{
        this(parent, title, text, colors, 16, 16, 1, names, StsListDialog.DEFAULT_ROW_HEIGHT, true, false);
    }

	public StsColorListDialog(JFrame parent, String title, String text, StsColor[] colors,
            int iconWidth, int iconHeight, int borderWidth, String[] names, float rowHeight, boolean finishDialog, boolean modal)
	{
        super(parent, title, text, null, rowHeight, false, modal, false,
                    new StsColorListItem(PROTOTYPE_COLOR, StsListDialog.PROTOTYPE_STRING,
                            iconWidth, iconHeight, borderWidth), new StsColorListRenderer());
        this.iconWidth = iconWidth;
        this.iconHeight = iconHeight;
        this.borderWidth = borderWidth;
        if (finishDialog)
        {
            setItems(colors, names);
	    	pack();
            adjustSize();
        }
	}

    public void setItems(StsColor[] colors, String[] names)
    {
        if (colors==null || names==null) return;
        int nItems = Math.min(colors.length, names.length);
        StsColorListItem[] items = new StsColorListItem[nItems];
        for (int i=0; i<nItems; i++)
        {
            items[i] = new StsColorListItem(colors[i], names[i], iconWidth,
                    iconHeight, borderWidth);
        }
        list.setListData(items);
    }

    public static void main(String[] args)
    {
        String[] names = StsColor.colorNames8;
        StsColor[] colors = StsColor.colors8;
        StsColorListDialog d = new StsColorListDialog(null, "StsColorListDialog test",
                "list of colors:", colors, names);
        d.setModal(true);
        d.setVisible(true);
        System.exit(0);
    }

}


