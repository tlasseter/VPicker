
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

// dialog box to set a text field

package com.Sts.UI;

import com.Sts.DBTypes.*;

public class StsColorListPanel extends StsCellListPanel
{
    private int iconHeight;
    private int iconWidth;
    private int borderWidth;

	public StsColorListPanel(String[] names, StsColor[] colors)
	{
        this(names, colors, 16, 16, 1, StsListDialog.DEFAULT_ROW_HEIGHT);
    }

	public StsColorListPanel(String[] names, StsColor[] colors, int iconWidth, int iconHeight, int borderWidth, float rowHeight)
	{
        super(null, rowHeight, new StsColorListItem(StsColor.RED, "prototype string", iconWidth, iconHeight, borderWidth), new StsColorListRenderer());
        this.iconWidth = iconWidth;
        this.iconHeight = iconHeight;
        this.borderWidth = borderWidth;
        setItems(colors, names);
	}

    public void setItems(StsColor[] colors, String[] names)
    {
        if (colors==null || names==null) return;
        int nItems = Math.min(colors.length, names.length);
        StsColorListItem[] items = new StsColorListItem[nItems];
        for (int i=0; i<nItems; i++)
            items[i] = new StsColorListItem(colors[i], names[i], iconWidth, iconHeight, borderWidth);
        list.setListData(items);
    }

    public static void main(String[] args)
    {
        String[] names = StsColor.colorNames8;
        StsColor[] colors = StsColor.colors8;
        StsColorListPanel listPanel = new StsColorListPanel(names, colors);
        StsCellListPanel.createDialog(listPanel);
    }
}
