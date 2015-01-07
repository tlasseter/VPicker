
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

// select from an object instance list of an StsObject that implements StsSelectable

package com.Sts.RFUI;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsSelectStsObjects
{
    private StsModel model;
	private StsWin3d win3d;
    private Class objectClass;
    private String title;
    private String text;
    private boolean sort;
    protected StsObject[] objects = null;
    protected boolean singleSelect;
    protected StsListSelector selector = null;
    private boolean visibleOnly = false;
    private int nVisibleObjects = 0;

	public StsSelectStsObjects()
	{
	}

    static public StsSelectStsObjects constructor(StsModel model, Class objectClass,
            String title, String text, boolean singleSelect)
	{
		try
		{
			return new StsSelectStsObjects(model, objectClass, title, text, singleSelect);
		}
		catch(Exception e)
		{
		    StsException.outputException("StsSelectStsObjects constructor() failed.",
				e, StsException.WARNING);
			return null;
		}
	}

    /** constructor for a sorted model instance list of objects */
    private StsSelectStsObjects(StsModel model, Class objectClass,
            String title, String text, boolean singleSelect)
        throws StsException
    {
        this(model, objectClass, title, text, singleSelect, true);
    }

    /** constructor for a model instance list of objects with optional sort */
    private StsSelectStsObjects(StsModel model, Class objectClass,
            String title, String text, boolean singleSelect, boolean sort)
        throws StsException
    {
        if (model==null)
        {
            throw new StsException(StsException.WARNING,
                    "StsSelectStsObjects.StsSelectStsObjects:  Model is null!");
        }
        if (objectClass==null)
        {
            throw new StsException(StsException.WARNING,
                    "StsSelectStsObjects.StsSelectStsObjects:  Object class is null!");
        }
        if (!StsSelectStsObjects.StsObjectIsSelectable(objectClass))
        {
            throw new StsException(StsException.WARNING,
                    "StsSelectStsObjects.StsSelectStsObjects:  Object isn't selectable.");
        }
        this.model = model;
		this.win3d = model.win3d;
        this.objectClass = objectClass;
        this.title = title;
        this.text = text;
        this.singleSelect = singleSelect;
        this.sort = sort;
        objects = model.getObjectList(objectClass);
    }

    static public StsSelectStsObjects constructor(StsModel model, Class objectClass,
            String title, String text, boolean singleSelect, boolean sort)
	{
		try
		{
			return new StsSelectStsObjects(model, objectClass, title, text, singleSelect, sort);
		}
		catch(Exception e)
		{
		    StsException.outputException("StsSelectStsObjects constructor() failed.",
				e, StsException.WARNING);
			return null;
		}
	}

    /** constructor for an array of objects */
    private StsSelectStsObjects(StsModel model, StsObject[] objects, String title,
            String text, boolean singleSelect)
        throws StsException
    {
        if (model==null)
        {
            throw new StsException(StsException.WARNING,
                    "StsSelectStsObjects.StsSelectStsObjects:  Model is null!");
        }
        this.model = model;
        if (objects==null || objects.length == 0)
        {
            throw new StsException(StsException.WARNING,
                    "StsSelectStsObjects.StsSelectStsObjects:  Objects are null!");
        }
        this.objects = objects;
        if (objects!=null)
        {
            objectClass = objects[0].getClass();
            if (!StsObjectIsSelectable(objectClass))
            {
                throw new StsException(StsException.WARNING,
                    "StsSelectStsObjects.StsSelectStsObjects:  Objects can't be selected.");
            }
        }
        this.title = title;
        this.text = text;
        this.singleSelect = singleSelect;
    }

    static public StsSelectStsObjects constructor(StsModel model, StsObject[] objects, String title,
            String text, boolean singleSelect)
	{
		try
		{
			return new StsSelectStsObjects(model, objects, title, text, singleSelect);
		}
		catch(Exception e)
		{
		    StsException.outputException("StsSelectStsObjects constructor() failed.",
				e, StsException.WARNING);
			return null;
		}
	}

    public void setVisibleOnly(boolean visibleOnly) { this.visibleOnly = visibleOnly; }

    private StsObject[] getVisibleObjects()
    {
        // count visible objects
        int nObjects = (objects==null) ? 0 : objects.length;
        if (nObjects == 0) return null;

        StsObject[] visibleObjects = new StsObject[nObjects];

        nVisibleObjects = 0;
        for (int i=0; i<nObjects; i++)
        {
            StsMainObject object = (StsMainObject)objects[i];
            if (object.getIsVisible()) visibleObjects[nVisibleObjects++] = (StsObject)object;
        }

        return (StsObject[])StsMath.trimArray(visibleObjects, nVisibleObjects);
    }

    public int getNVisibleObjects() { return nVisibleObjects; }

    public String[] selectNames()
    {
        StsObject[] objs = visibleOnly ? getVisibleObjects() : objects;
        if (objs == null) return null;
        if (objs.length == 1)
        {
            String[] selectedNames = { objs[0].getName() };
            return selectedNames;
        }
        setSelector(objs);

        return selector.getSelectedItems();
    }

    protected void setSelector(StsObject[] objects)
    {
        String[] names;
        StsColor[] colors;

        if (objects!=null)
        {
            names = new String[objects.length];
            colors = new StsColor[objects.length];
            for (int i=0; i<objects.length; i++)
            {
                names[i] = objects[i].getName();
                colors[i] = objects[i].getStsColor();
            }
            if (sort) StsMath.qsort(names, colors, objects);
        }
        else  // empty list
        {
            colors = StsListStsObjects.NULL_COLORS;
            names = StsListStsObjects.NULL_NAMES;
        }

        setSelector(colors, names);
    }

    protected void setSelector(StsColor[] colors, String[] names)
    {
        if (colors[0]==null)  // use regular selector
        {
            if (selector==null)
            {
                selector = new StsListSelector(win3d, title, text,
                        names);
                selector.setSingleSelectionMode(singleSelect);
            }
            else
            {
                selector.setItems(names);
            }
        }
        else // use color selector
        {
            if (selector==null)
            {
                selector = new StsColorListSelector(win3d, title, text, colors, names);
                selector.setSingleSelectionMode(singleSelect);
            }
            else
            {
                ((StsColorListSelector)selector).setItems(colors, names);
            }
        }
        selector.setLocationRelativeTo(model.win3d);
        selector.setVisible(true);
    }
    public StsListSelector getSelector() { return selector; }

    public String selectName()
    {
        if (objects==null) return null;
        String[] selectedNames = selectNames();
        return (selectedNames==null) ? null : selectedNames[0];
    }

    public int[] selectIndices()
    {
        StsObject[] objs = visibleOnly ? getVisibleObjects() : objects;
        if (objs == null) return null;
        if (objs.length==1)
        {
            int index = 0;
            if (visibleOnly)
            {
                for (int i=0; i<objects.length; i++)
                {
                    if (objs[0]==objects[i])
                    {
                        index = i;
                        break;
                    }
                }
            }

            int[] selectedIndices = { index };
            return selectedIndices;
        }
        setSelector(objs);

        return selector.getSelectedIndices();
    }

    public StsObject[] selectObjects()
    {
        int[] selectedIndices = selectIndices();
        if (selectedIndices==null) return null;
        if (selectedIndices[selectedIndices.length-1]>=objects.length) return null;
        StsObject[] selectedObjects = new StsObject[selectedIndices.length];
        for (int i=0; i<selectedIndices.length; i++)
        {
            selectedObjects[i] = objects[selectedIndices[i]];
        }
        return selectedObjects;
    }

    public StsObject selectObject()
    {
        StsObject[] objects = selectObjects();
        return (objects==null) ? null : objects[0];
    }

    static boolean StsObjectIsSelectable(Class objectClass)
    {
        Class[] interfaces = objectClass.getInterfaces();
        if (interfaces==null) return false;
        for (int i=0; i<interfaces.length; i++)
        {
            if (interfaces[i]==StsSelectable.class) return true;
        }
        return false;
    }
}
