
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.RFUI;

import com.Sts.DBTypes.*;
import com.Sts.MVC.View3d.*;

import java.awt.event.*;

public class StsSurfaceLinkDialog extends StsOrderedListDialog//StsSelectStsObjects
{
    public StsSurfaceLinkDialog(StsWin3d win3d, String[] surfaceNames)
    {
        super(win3d, "Link surfaces with markers", "Surfaces (with marker links):",
                StsModelSurface.class, surfaceNames, true, null, null);
    }

    public void setItems(StsColor[] colors, String[] names)
    {
        int nItems = (names==null) ? 0 : names.length;
        if (nItems > 0)
        {
            names = parseNames(names);
            for (int i=0; i<nItems; i++)
            {
                StsModelSurface s = null;
                StsSurfaceClass surfaceClass = (StsSurfaceClass)model.getCreateStsClass(StsModelSurface.class);
                s = (StsModelSurface)surfaceClass.getObjectWithName(names[i]);
                if (s==null) continue;
                StsMarker marker = s.getMarker();
                if (marker == null) names[i] += " ( )";
                else names[i] += " (" + marker.getName() + ")";
            }
        }
        super.setItems(colors, names);
    }

    public StsModelSurface getSelectedSurface()
    {
        String name = parseName(getSelectedItem());
        if (name == null) return null;

        StsModelSurface s = null;
        StsSurfaceClass surfaceClass = (StsSurfaceClass)model.getCreateStsClass(StsModelSurface.class);
        return (StsModelSurface)surfaceClass.getObjectWithName(name);
    }

    protected void okayButton_actionPerformed(ActionEvent e)
    {
        getSelectedItems();
        super.okayButton_actionPerformed(e);
    }

    /*
    public StsSurfaceLinkDialog(StsModel model) throws StsException
    {
        super(model, model.getSurfaces(StsModelSurface.IMPORTED),
                "Link surfaces with markers",
                "Select a surface (marker links shown):", true);
    }

    // overrides StsSelectStsObjects method
    protected void setSelector(StsObject[] objects)
    {
        String[] names;
        Color[] colors;
        if (objects!=null)
        {
            names = new String[objects.length];
            colors = new Color[objects.length];
            for (int i=0; i<objects.length; i++)
            {
                StsSurface surface = (StsSurface)objects[i];
                names[i] = surface.getName();
                StsMarker marker = surface.getMarker();
                if (marker == null) names[i] += "   ( )";
                else names[i] += "   (" + marker.getName() + ")";
                colors[i] = surface.getColor();
            }
        }
        else  // empty list
        {
            colors = StsListStsObjects.NULL_COLORS;
            names = StsListStsObjects.NULL_NAMES;
        }

        setSelector(colors, names);
    }
    */
}
