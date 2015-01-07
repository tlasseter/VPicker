
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import com.Sts.DBTypes.*;

import javax.swing.*;

public class StsColorListSelector extends StsListSelector
{
    private int iconHeight;
    private int iconWidth;
    private int borderWidth;
    private StsColorListItem[] items = null;

	public StsColorListSelector(JFrame parent, String title, String text,
            StsColor[] colors, String[] names)
	{
        this(parent, title, text, colors, 16, 16, 1, names,
                StsListDialog.DEFAULT_ROW_HEIGHT, true, true, null, null);
    }
	public StsColorListSelector(JFrame parent, String title, String text,
            StsColor[] colors, int iconWidth, int iconHeight, int borderWidth,
            String[] names, float rowHeight, boolean finishDialog)
    {
        this(parent, title, text, colors, iconWidth, iconHeight, borderWidth,
                names, rowHeight, finishDialog, false, null, null);
    }
	public StsColorListSelector(JFrame parent, String title, String text,
            StsColor[] colors, int iconWidth, int iconHeight, int borderWidth,
            String[] names, float rowHeight, boolean finishDialog,
            boolean useCancelBtn, String button1Text, String button2Text)
    {
	    super(parent, title, text, names, rowHeight, finishDialog, useCancelBtn,
                button1Text, button2Text,
                new StsColorListItem(StsColorListDialog.PROTOTYPE_COLOR,
                    StsDoubleListDialog.PROTOTYPE_STRING, iconWidth, iconHeight,
                    borderWidth),
                new StsColorListRenderer());
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
        items = new StsColorListItem[nItems];
        for (int i=0; i<nItems; i++)
        {
            items[i] = new StsColorListItem(colors[i], names[i], iconWidth, iconHeight, borderWidth);
        }
        list.setListData(items);
    }

    public String[] getSelectedItems()
    {
        int[] selectedIndices = getSelectedIndices();
        if (selectedIndices==null) return null;

        String[] selectedItems = new String[selectedIndices.length];
        for (int i=0; i<selectedIndices.length; i++)
        {
            selectedItems[i] = new String(items[selectedIndices[i]].getName());
        }
        return selectedItems;
    }
}


