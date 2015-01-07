
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

public class StsRenameDialog extends StsColorListSelector
{
    private StsWin3d win3d;
    private StsModel model;
    private StsClass instList;
    private StsObject[] objects;
    private String[] names;
    private StsColor[] colors;
    private StsSetTextDialog textDialog = null;

    /** constructor for a sorted model instance list of objects */
    public StsRenameDialog(StsModel model, Class StsObjectClass,
            String title, String label)
        throws StsException
    {
        this(model, StsObjectClass, title, label, true);
    }

    /** constructor for a model instance list of objects with optional sort */
    public StsRenameDialog(StsModel model, Class StsObjectClass,
            String title, String label, boolean sort)
        throws StsException
    {
        super((JFrame)model.win3d, title, label, null, 16, 16, 1, null,
                StsListDialog.DEFAULT_ROW_HEIGHT, false, false,
                "Rename", null);
        if (StsObjectClass==null)
        {
            throw new StsException(StsException.FATAL,
                    "StsRenameDialog.StsRenameDialog:  Object class is null!");
        }
        if (!StsSelectStsObjects.StsObjectIsSelectable(StsObjectClass))
        {
            throw new StsException(StsException.FATAL,
                    "StsRenameDialog.StsRenameDialog:  Object isn't selectable.");
        }
        this.model = model;
        win3d = model.win3d;
        setLocationRelativeTo(win3d);
        instList = model.getCreateStsClass(StsObjectClass);
        if (instList == null || instList.getSize() <= 0)
        {
            StsMessageFiles.logMessage("No items to rename found!");
            return;
        }

        StsClass objectList = model.getCreateStsClass(StsObjectClass);
        int nObjects = (objectList==null) ? 0 : objectList.getSize();
        if (nObjects > 0)
        {
            objects = new StsObject[nObjects];
            names = new String[nObjects];
            colors = new StsColor[nObjects];
            for (int i=0; i<nObjects; i++)
            {
                objects[i] = objectList.getElement(i);
                names[i] = objects[i].getName();
                colors[i] = objects[i].getStsColor();
            }
            if (sort) StsMath.qsort(names, colors, objects);
            setItems(colors, names);
        }

        setModal(false);
        setSingleSelectionMode(true);
        pack();
        adjustSize();
    }

    /** constructor for a model instance list of objects with optional sort */
    public StsRenameDialog(StsModel model, StsObject[] objects,
            String title, String label, boolean sort)
        throws StsException
    {
        super((JFrame)model.win3d, title, label, null, 16, 16, 1, null,
                StsListDialog.DEFAULT_ROW_HEIGHT, false, false,
                "Set Name", null);
        if (objects == null || objects[0] == null)
        {
            throw new StsException(StsException.FATAL,
                    "StsRenameDialog.StsRenameDialog:  Object array is null!");
        }
        Class StsObjectClass = objects[0].getClass();
        if (StsObjectClass==null)
        {
            throw new StsException(StsException.FATAL,
                    "StsRenameDialog.StsRenameDialog:  Object class is null!");
        }
        if (!StsSelectStsObjects.StsObjectIsSelectable(StsObjectClass))
        {
            throw new StsException(StsException.FATAL,
                    "StsRenameDialog.StsRenameDialog:  Object isn't selectable.");
        }
        this.objects = objects;
        this.model = model;
        win3d = model.win3d;
        setLocationRelativeTo(win3d);
        instList = model.getCreateStsClass(StsObjectClass);
        if (instList == null || instList.getSize() <= 0)
        {
            StsMessageFiles.logMessage("No items to rename found!");
            return;
        }

        names = new String[objects.length];
        colors = new StsColor[objects.length];
        for (int i=0; i<objects.length; i++)
        {
            names[i] = objects[i].getName();
            colors[i] = objects[i].getStsColor();
        }
        if (sort) StsMath.qsort(names, colors, objects);
        setItems(colors, names);

        setModal(false);
        setSingleSelectionMode(true);
        pack();
        adjustSize();
    }

    protected void button1_actionPerformed(ActionEvent e)
    {
        int[] index = list.getSelectedIndices();
        if (index==null) return;
        int i = index[0];

        if (textDialog == null)
        {
            textDialog = new StsSetTextDialog(new JFrame(), "Set new name",
                    objects[i].getName());
        }
        else textDialog.setText(objects[i].getName());
        textDialog.setVisible(true);
        String newName = textDialog.getText();
        if (newName == null) return;
        if (instList.getElement(newName) != null)
        {
            StsMessageFiles.logMessage("The item can't be renamed to: " + newName +
                    ".  That name is already used.");
            return;
        }

        StsMainObject mainObject = (StsMainObject)objects[i];
        mainObject.setName(newName);
        names[i] = newName;
        setItems(colors, names);
        model.refreshObjectPanel();
    }

    protected void okayButton_actionPerformed(ActionEvent parm1)
    {
        if (textDialog != null) textDialog.setVisible(false);
        super.okayButton_actionPerformed( parm1);
    }
}

