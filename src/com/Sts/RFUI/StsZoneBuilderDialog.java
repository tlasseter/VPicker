
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

public class StsZoneBuilderDialog extends StsOrderedListDialog
{
    private StsWin3d win3d;

    public StsZoneBuilderDialog(StsWin3d win3d, String[] surfaceNames)
    {
        super(win3d, "Define zones from  ", "Surfaces (with associated marker & zones):",
                StsModelSurface.class, surfaceNames, true, null, null);
        setSingleSelectionMode(false);
    }

    protected void okayButton_actionPerformed(ActionEvent e)
    {
        String[] selectedItems = getSelectedItems();
        if (selectedItems==null || selectedItems.length<2)
        {
            new StsMessage(null, StsMessage.WARNING,
                    "Please select at least two surfaces or cancel.");
            return;
        }
        super.okayButton_actionPerformed(e);
    }

    public void setItems(StsColor[] colors, String[] names)
    {
        int nItems = (names==null) ? 0 : names.length;
        if (nItems > 0)
        {
            StsSurfaceClass surfaceClass = (StsSurfaceClass)model.getCreateStsClass(StsModelSurface.class);
            names = parseNames(names);
            for (int i=0; i<nItems; i++)
            {
                StsModelSurface s = null;
                s = (StsModelSurface)surfaceClass.getObjectWithName(names[i]);
                if (s==null) continue;
                StsMarker m = s.getMarker();
                if (m == null) names[i] += "  ( )";
                else names[i] += "  (" + m.getName() + ")";
                StsZone z = s.getZoneAbove();
                if (z != null) names[i] += "  :  base of " + z.getName();
                z = s.getZoneBelow();
                if (z != null) names[i] += "  :  top of " + z.getName();
            }
        }
        super.setItems(colors, names);
    }
}
