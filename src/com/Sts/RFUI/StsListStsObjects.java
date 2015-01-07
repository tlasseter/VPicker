
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

public class StsListStsObjects
{
    // constants
    static public final String EMPTY_LIST = "< none to list >";
    static public final int LIST_ORDER = 0;
    static public final int ALPHA_ASCENDING_ORDER = 1;
    static public final int ALPHA_DESCENDING_ORDER = 2;
    static public final int USER_ASCENDING_ORDER = 3;
    static public final int USER_DESCENDING_ORDER = 4;

    static public final StsColor[] NULL_COLORS = new StsColor[1];
    static public final String[] NULL_NAMES = { EMPTY_LIST };

    private StsModel model;
    private Class StsObjectClass;
    private String title;
    private String text;
    private int ordering;
    protected StsListDialog dialog = null;

    /** constructor for a model instance list of sorted objects */
    public StsListStsObjects(StsModel model, Class StsObjectClass,
            String title, String text)
        throws StsException
    {
        this(model, StsObjectClass, title, text, null, ALPHA_ASCENDING_ORDER, false);
    }

    /** constructor for a model instance list of objects with optional sort */
    public StsListStsObjects(StsModel model, Class StsObjectClass,
            String title, String text, String[] orderedNames, int ordering,
            boolean displayOrderingValue)
        throws StsException
    {
        if (model==null)
        {
            throw new StsException(StsException.FATAL,
                    "StsListStsObjects.StsListStsObjects:  Model is null!");
        }
        if (StsObjectClass==null)
        {
            throw new StsException(StsException.FATAL,
                    "StsListStsObjects.StsListStsObjects:  Object class is null!");
        }
        if (!StsSelectStsObjects.StsObjectIsSelectable(StsObjectClass))
        {
            throw new StsException(StsException.FATAL,
                    "StsListStsObjects.StsListStsObjects:  Objects can't be listed.");
        }
        this.model = model;
        this.StsObjectClass = StsObjectClass;
        this.title = title;
        this.text = text;
        this.ordering = ordering;

        StsClass instList = model.getCreateStsClass(StsObjectClass);
        int nObjects = (instList==null) ? 0 : instList.getSize();

        String[] names;
        StsColor[] colors;
        if (nObjects>0)
        {
            // get names and colors
            names = new String[nObjects];
            colors = new StsColor[nObjects];
            if (ordering==USER_DESCENDING_ORDER || ordering==USER_ASCENDING_ORDER)
            {
                String[] newOrderedNames = StsOrderedList.getOrderedNames(instList,
                        orderedNames, ordering==USER_ASCENDING_ORDER);
                for (int i=0; i<newOrderedNames.length; i++)
                {
                    StsObject obj = instList.getElement(newOrderedNames[i]);
                    names[i] = obj.getName();
                    if (displayOrderingValue)
                    {
                        names[i] += "  (" + (int)obj.getOrderingValue() + ")";
                    }
                    colors[i] = obj.getStsColor();
                }
            }
            else
            {
                for (int i=0; i<nObjects; i++)
                {
                    StsObject obj = instList.getElement(i);
                    names[i] = obj.getName();
                    if (displayOrderingValue)
                    {
                        names[i] += "  (" + (int)obj.getOrderingValue() + ")";
                    }
                    colors[i] = obj.getStsColor();
                }
                if (ordering==ALPHA_ASCENDING_ORDER) StsMath.qsort(names, colors);
            }
        }
        else  // empty list
        {
            colors = NULL_COLORS;
            names = NULL_NAMES;
        }

        // set the list dialog and show it
        setDialog(colors, names);
    }

    /** constructor for a sorted array of objects */
    public StsListStsObjects(StsModel model, StsObject[] objects, String title,
            String text)
        throws StsException
    {
        this(model, objects, title, text, ALPHA_ASCENDING_ORDER, false);
    }

    /** constructor for an array of objects with optional sort */
    public StsListStsObjects(StsModel model, StsObject[] objects, String title,
            String text, int ordering, boolean displayOrderingValue)
        throws StsException
    {
        if (model==null)
        {
            throw new StsException(StsException.FATAL,
                    "StsListStsObjects.StsListStsObjects:  Model is null!");
        }
        this.model = model;
        this.title = title;
        this.text = text;
        this.ordering = ordering;

        String[] names = null;
        StsColor[] colors = null;
        if (objects!=null)
        {
            // get names and colors
            names = new String[objects.length];
            colors = new StsColor[objects.length];
            if (ordering==USER_DESCENDING_ORDER || ordering==USER_ASCENDING_ORDER)
            {
                String[] orderedNames = StsOrderedList.getOrderedNames(objects,
                        null, ordering==USER_ASCENDING_ORDER);
                for (int i=0; i<orderedNames.length; i++)
                {
                    names[i] = orderedNames[i];
                    for (int j=0; j<objects.length; j++)
                    {
                        if (names[i].equals(objects[j].getName()))
                        {
                            if (displayOrderingValue)
                            {
                                names[i] += "  (" + (int)objects[j].getOrderingValue() + ")";
                            }
                            colors[i] = objects[j].getStsColor();
                            break;
                        }
                    }
                }
            }
            else
            {
                for (int i=0; i<objects.length; i++)
                {
                    names[i] = objects[i].getName();
                    if (displayOrderingValue)
                    {
                        names[i] += "  (" + (int)objects[i].getOrderingValue() + ")";
                    }
                    colors[i] = objects[i].getStsColor();
                }
                if (ordering==ALPHA_ASCENDING_ORDER) StsMath.qsort(names, colors);
            }
        }
        else  // empty list
        {
            colors = NULL_COLORS;
            names = NULL_NAMES;
        }

        // set the list dialog and show it
        setDialog(colors, names);
    }

    protected void setDialog(StsColor[] colors, String[] names)
    {
        // build and show dialog
        if (colors[0]==null) // no color
        {
            dialog = new StsListDialog(new JFrame(), title, text, names);
        }
        else
        {
            dialog = new StsColorListDialog(new JFrame(), title, text, colors, names);
        }
        dialog.setLocationRelativeTo(model.win3d);
        dialog.setVisible(true);
    }

}
