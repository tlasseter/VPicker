
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System


package com.Sts.RFUI;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.event.*;

public class StsOrderedListDialog extends StsColorListSelector
{
    private StsWin3d win3d;
    protected StsModel model;
    private int lastSelectedIndex = 0;
    private String[] orderedNames = null;
    private Class objectClass;
    private String[] names = null;
    private StsColor[] colors = null;

    public StsOrderedListDialog(StsWin3d win3d, String title, String label,
            Class objectClass, String[] orderedNames)
    {
        this(win3d, title, label, objectClass, orderedNames, false, "Move Up",
                "Move Down");
    }
    public StsOrderedListDialog(StsWin3d win3d, String title, String label,
            Class objectClass, String[] orderedNames, boolean useCancelBtn,
            String button1Text, String button2Text)
    {
        super((JFrame)win3d, title, label, null, 16, 16, 1, null,
                StsListDialog.DEFAULT_ROW_HEIGHT, false, useCancelBtn,
                button1Text, button2Text);
        if (win3d==null) return;
        this.win3d = win3d;
        model = win3d.getModel();
        if (model==null) return;
        this.orderedNames = orderedNames;
        this.objectClass = objectClass;

        //okayButton.setText("Build");
        setSingleSelectionMode(true);
        if (!setListItems()) return;
        pack();
        adjustSize();
    }

    private boolean setListItems()
    {
        if (orderedNames==null) return false;
        if (orderedNames.length<1) return false;
        StsClass instList = model.getCreateStsClass(objectClass);
        if (instList==null) return false;
        names = new String[orderedNames.length];
        colors = new StsColor[orderedNames.length];
        for (int i=0; i<orderedNames.length; i++)
        {
            StsObject obj = instList.getElement(orderedNames[i]);
            if (obj==null) return false;
            names[i] = obj.getName();
            colors[i] = obj.getStsColor();
        }
        setItems(colors, names);
        return true;
    }

    protected void button1_actionPerformed(ActionEvent e)
    {
        // perform up action
        moveActionUp(true);
    }
    protected void button2_actionPerformed(ActionEvent e)
    {
        // perform down action
        moveActionUp(false);
    }

    private void moveActionUp(boolean up)
    {
        int[] index = list.getSelectedIndices();
        if (index==null) return;
        int i = index[0];
        if (!ableToMove(orderedNames[i]))
        {
            StsMessageFiles.logMessage("Sorry, unable to move:  " + orderedNames[i]);
            return;
        }
        if (i<0 || i>orderedNames.length-1) return;
        if (up)
        {
            if (i==0) return;
            StsMath.swap(orderedNames, i-1, i);
            StsMath.swap(names, i-1, i);
            StsMath.swap(colors, i-1, i);
            lastSelectedIndex = i-1;
        }
        else  // down
        {
            if (i==orderedNames.length-1) return;
            StsMath.swap(orderedNames, i, i+1);
            StsMath.swap(names, i, i+1);
            StsMath.swap(colors, i, i+1);
            lastSelectedIndex = i+1;
        }
        setItems(colors, names);
        list.setSelectedIndex(lastSelectedIndex);
    }

    protected String[] getOrderedNames() { return orderedNames; }

    // override in derived class
    protected boolean ableToMove(String name) { return true; }

    static public String parseName(String name)
    {
        String parsedName = null;
        if (name == null) return null;
        int index = name.indexOf("(");
        if (index < 1) parsedName = name;
        else parsedName = name.substring(0, index-1).trim();
        return parsedName;
    }

    static public String[] parseNames(String[] names)
    {
        if (names == null) return null;
        String[] parsedNames = new String[names.length];
        for (int i=0; i<names.length; i++)
        {
            parsedNames[i] = parseName(names[i]);
        }
        return parsedNames;
    }
}

