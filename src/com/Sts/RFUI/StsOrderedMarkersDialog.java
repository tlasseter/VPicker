
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

public class StsOrderedMarkersDialog extends StsOrderedListDialog
{
    public StsOrderedMarkersDialog(StsWin3d win3d)
    {
        super(win3d, "Reorder markers by stratigraphy",
                "Markers (with average depth & associated well zones):",
                StsMarker.class,
                ((StsWellMarkerClass)win3d.getModel().getCreateStsClass(StsWellMarker.class)).getOrderedMarkerNames(true));
    }

    protected void okayButton_actionPerformed(ActionEvent e)
    {
        ((StsWellMarkerClass)model.getCreateStsClass(StsWellMarker.class)).setOrderedMarkerNames(getOrderedNames());
        super.okayButton_actionPerformed(e);
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
                names[i] += "  (" + ((int)marker.getOrderingValue()) + ")";
                StsWellZoneSet wzs = marker.getWellZoneSetAbove();
                if (wzs != null) names[i] += "  :  base of " + wzs.getName();
                wzs = marker.getWellZoneSetBelow();
                if (wzs != null) names[i] += "  :  top of " + wzs.getName();
            }
        }
        super.setItems(colors, names);
    }

    // override base class
    protected boolean ableToMove(String name)
    {
        if (name == null) return false;
        try
        {
            StsMarker marker = (StsMarker)model.getObjectWithName(StsMarker.class, name);
            if (marker==null) return false;
            if (marker.getWellZoneSetAbove() != null || marker.getWellZoneSetBelow() != null) return false;
        }
        catch (Exception e) { return false; }
        return true;
    }
}

