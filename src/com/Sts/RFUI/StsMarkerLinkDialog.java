
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

public class StsMarkerLinkDialog extends StsOrderedListDialog //StsSelectStsObjects
{
    public StsMarkerLinkDialog(StsWin3d win3d, String title)
    {
        super(win3d, title, "Markers (with surface links):",
                StsMarker.class, ((StsWellMarkerClass)(win3d.getModel().getCreateStsClass(StsWellMarker.class))).getOrderedMarkerNames(),
                true, null, null);
    }

    public void setItems(StsColor[] colors, String[] names)
    {
        int nItems = (names==null) ? 0 : names.length;
        if (nItems > 0)
        {
            names = parseNames(names);
            for (int i=0; i<nItems; i++)
            {
                StsMarker marker = (StsMarker)model.getObjectWithName(StsMarker.class, names[i]);
                if (marker == null) continue;
                StsModelSurface s = marker.getModelSurface();
                if (s == null) names[i] += " ( )";
                else names[i] += " (" + s.getName() + ")";
            }
        }
        super.setItems(colors, names);
    }


    public StsMarker getSelectedMarker()
    {
        String name = parseName(getSelectedItem());
        if (name == null) return null;
        return (StsMarker)model.getObjectWithName(StsMarker.class, name);
    }

    protected void okayButton_actionPerformed(ActionEvent e)
    {
        getSelectedItems();
        super.okayButton_actionPerformed(e);
    }

    /*
    public StsMarkerLinkDialog(StsModel model, String title) throws StsException
    {
        super(model, StsMarker.class, title,
                "Select a marker (surface links shown):", true, false);
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
                StsMarker marker = (StsMarker)objects[i];
                names[i] = marker.getName();
                StsSurface surface = marker.getParentSurface();
                if (surface == null) names[i] += "   ( )";
                else names[i] += "   (" + surface.getName() + ")";
                colors[i] = marker.getColor();
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
