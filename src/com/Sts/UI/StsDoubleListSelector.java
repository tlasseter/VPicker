
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import com.Sts.DBTypes.*;

import javax.swing.*;

public class StsDoubleListSelector extends StsDoubleListDialog
{
    private int[] selectedIndices = null;

	public StsDoubleListSelector(JFrame parent, String title, String text,
            String label1, String[] items1, String label2, String[] items2)
	{
        this(parent, title, text, label1, items1, label2, items2, true,
                StsListDialog.DEFAULT_ROW_HEIGHT, true,
                false, null, null, StsDoubleListDialog.PROTOTYPE_STRING, null);
    }
	public StsDoubleListSelector(JFrame parent, String title, String text,
            String label1, String[] items1, String label2, String[] items2,
            boolean offset, float rowHeight, boolean finishDialog,
            boolean useCancelBtn, String button1Text,
            String button2Text, Object itemPrototype, ListCellRenderer cellRenderer)
	{
        super(parent, title, text, label1, items1, label2, items2,
            offset, rowHeight, finishDialog, true, useCancelBtn,
            button1Text, button2Text, itemPrototype, cellRenderer);
    }

    /** list selection mode */
    public void setSingleSelectionMode(boolean singleMode)
    {
        ListSelectionModel model = list.getSelectionModel();
        if (singleMode)
        {
            model.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        }
        else
        {
            model.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        }
    }

    /** pre-selected indices/items */
    public void setSelectedIndices(int[] indices)
    {
        if (indices==null) return;
        list.setSelectedIndices(indices);
    }
    public void setSelectedIndex(int index)
    {
        list.setSelectedIndex(index);
    }
    public void setSelectedItems(String[] items)
    {
        if (items==null) return;
        list.setSelectedValue(items[0], true);
        for (int i=1; i<items.length; i++) list.setSelectedValue(items[i], false);
    }
    public void setSelectedItem(String item)
    {
        if (item==null) return;
        list.setSelectedValue(item, true);
    }

    /** return values */
    public int[] getSelectedIndices()
    {
        boolean debug = false; //StsTrace.getTrace();
        if (selectedIndices==null)
        {
            int[] selected = list.getSelectedIndices();
            if (selected==null) selectedIndices = null;
            else if (selected.length==0) selectedIndices = null;
            else
            {
                selectedIndices = new int[selected.length];
                System.arraycopy(selected, 0, selectedIndices, 0, selected.length);
            }
        }
        return selectedIndices;
    }
    public String[] getSelectedItems()
    {
        String[] items = getItems();
        if (items==null) return null;
        getSelectedIndices();
        String[] selectedItems = null;
        if (selectedIndices!=null)
        {
            selectedItems = new String[selectedIndices.length];
            for (int i=0; i<selectedIndices.length; i++)
            {
                selectedItems[i] = new String(items[selectedIndices[i]]);
            }
        }
        return selectedItems;
    }
    public String getSelectedItem()
    {
        String[] selectedItems = getSelectedItems();
        return (selectedItems==null) ? null : selectedItems[0];
    }

    public static void main(String[] args)
    {
        StsColor[] colors1 = { StsColor.RED, StsColor.GREEN, StsColor.BLUE, StsColor.YELLOW,
                StsColor.MAGENTA, StsColor.BLACK };
        String[] names1 = { "Red", "Green", "Blue", "Yellow", "Magenta", "Black" };
        StsColorListItem[] items1 = new StsColorListItem[colors1.length];
        for (int i=0; i<items1.length; i++)
        {
            items1[i] = new StsColorListItem(colors1[i], names1[i]);
        }

        StsColor[] colors2 = { StsColor.ORANGE, StsColor.WHITE, StsColor.CYAN, StsColor.PINK,
                StsColor.DARK_GRAY };
        String[] names2 = { "square", "circle", "triangle", "diamond", "hexagon" };
        StsColorListItem[] items2 = new StsColorListItem[colors2.length];
        for (int i=0; i<items2.length; i++)
        {
            items2[i] = new StsColorListItem(colors2[i], names2[i]);
        }
        StsDoubleListSelector d = new StsDoubleListSelector(new JFrame(),
                "StsDoubleListSelector test", "Double List Example", "Colors: ",
                names1, "Shapes: ", names2, true, 2, false, true, "1", "2",
                new StsColorListItem(StsColorListDialog.PROTOTYPE_COLOR,
                    StsDoubleListDialog.PROTOTYPE_STRING),
                new StsColorListRenderer());
        d.setSingleSelectionMode(false);
        d.getList().setListData(items1);
        d.getList2().setListData(items2);
        d.pack();
        d.adjustSize();
        d.setVisible(true);
        String[] selected = d.getSelectedItems();
        System.out.println("Items:");
        if (selected == null) System.out.println("\tNothing selected!");
        else for (int i=0; i<selected.length; i++) System.out.println("\t"+selected[i]);
        System.exit(0);
    }

}