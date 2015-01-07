
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.RFUI;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.event.*;

public class StsColorSelectDialog extends StsColorListSelector
{
    private StsModel model;
    private StsSpectrum spectrum;
    private StsObject[] objects;
    private String[] names;
    private StsColor[] colors;

    /** constructor for a sorted model instance list of objects */
    public StsColorSelectDialog(StsModel model, Class StsObjectClass, String title, String label)
        throws StsException
    {
        this(model, StsObjectClass, title, label, true);
    }

    /** constructor for a model instance list of objects with optional sort */
    public StsColorSelectDialog(StsModel model, Class StsObjectClass, String title, String label, boolean sort)
        throws StsException
    {
        super((JFrame)model.win3d, title, label, null, 16, 16, 1, null,
                StsListDialog.DEFAULT_ROW_HEIGHT, false, false,
                "Set", "Choose");
        if (StsObjectClass==null)
        {
            throw new StsException(StsException.WARNING,
                    "StsColorSelectDialog.StsColorSelectDialog:  Object class is null!");
        }
        if (!StsSelectStsObjects.StsObjectIsSelectable(StsObjectClass))
        {
            throw new StsException(StsException.WARNING,
                    "StsColorSelectDialog.StsColorSelectDialog:  Object isn't selectable.");
        }
        StsClass instList = model.getCreateStsClass(StsObjectClass);
        if (instList == null || instList.getSize() <= 0)
        {
            StsMessageFiles.logMessage("No items to change colors found!");
            return;
        }
        this.model = model;
        setLocationRelativeTo(model.win3d);
        spectrum = model.getSpectrum("Basic");

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

    protected void button1_actionPerformed(ActionEvent e)
    {
        int[] index = list.getSelectedIndices();
        if (index==null) return;
        int i = index[0];
        StsColor currentColor = spectrum.getCurrentColor();
        objects[i].setStsColor(currentColor);
        colors[i] = currentColor;
        setItems(colors, names);
        model.refreshObjectPanel();
    }
    protected void button2_actionPerformed(ActionEvent e)
    {
        spectrum.display();
    }

    protected void okayButton_actionPerformed(ActionEvent parm1)
    {
        spectrum.closeDisplay();
        super.okayButton_actionPerformed( parm1);
    }

    static public void main(String[] args)
    {
        try
        {
            StsModel model = new StsModel();
            int nPolygons = 5;
            for(int n = 0; n < nPolygons; n++)
            {
                StsType type = new StsType("red", StsColor.RED);
                com.Sts.DBTypes.StsXPolygon polygon = new com.Sts.DBTypes.StsXPolygon(type);
//                model.add(polygon);
            }
            StsColorSelectDialog d = new StsColorSelectDialog(model, StsXPolygon.class, "Title", "Label", true);
            d.pack();
            d.setVisible(true);
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }

}

