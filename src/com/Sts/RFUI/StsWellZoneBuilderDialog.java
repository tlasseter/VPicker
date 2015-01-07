
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.RFUI;

import com.Sts.DBTypes.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.*;

import java.awt.event.*;

public class StsWellZoneBuilderDialog extends StsOrderedListDialog
{
    public StsWellZoneBuilderDialog(StsWin3d win3d)
    {
        super(win3d, "Build well zones",
                "Markers (with associated surface & well zones):",
                StsMarker.class, ((StsWellMarkerClass)(win3d.getModel().getCreateStsClass(StsWellMarker.class))).getOrderedMarkerNames(true),
                true, null, null);
        setSingleSelectionMode(false);
    }

    protected void okayButton_actionPerformed(ActionEvent e)
    {
        String[] selectedItems = getSelectedItems();
        if (selectedItems==null || selectedItems.length<2)
        {
            new StsMessage(null, StsMessage.WARNING,
                    "Please select at least two markers or cancel.");
            return;
        }
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
                if (marker==null) continue;
                StsModelSurface s = marker.getModelSurface();
                if (s == null) names[i] += "  ( )";
                else names[i] += "  (" + s.getName() + ")";
                StsWellZoneSet wzs = marker.getWellZoneSetAbove();
                if (wzs != null) names[i] += "  :  base of " + wzs.getName();
                wzs = marker.getWellZoneSetBelow();
                if (wzs != null) names[i] += "  :  top of " + wzs.getName();
            }
        }
        super.setItems(colors, names);
    }

}
