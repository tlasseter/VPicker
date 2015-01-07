
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import javax.swing.*;

public class StsListSelector extends StsDoubleListSelector
{
    private int[] selectedIndices = null;

	public StsListSelector(JFrame parent, String title, String text,
            String[] items)
	{
        this(parent, title, text, items, StsListDialog.DEFAULT_ROW_HEIGHT, true,
                false, null, null);
    }
	public StsListSelector(JFrame parent, String title, String text,
            String[] items, float rowHeight, boolean finishDialog,
            boolean useCancelBtn, String button1Text, String button2Text)
	{
        this(parent, title, text, items, rowHeight, finishDialog, useCancelBtn,
            button1Text, button2Text, StsDoubleListDialog.PROTOTYPE_STRING, null);
    }
	public StsListSelector(JFrame parent, String title, String text,
            String[] items, float rowHeight, boolean finishDialog,
            boolean useCancelBtn, String button1Text, String button2Text,
            Object itemPrototype, ListCellRenderer cellRenderer)
	{
        super(parent, title, text, null, items, null, null, false, rowHeight,
                finishDialog, useCancelBtn, button1Text, button2Text, itemPrototype,
                cellRenderer);
    }

}

