
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
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.event.*;

public class StsSurfaceAssignmentDialog extends StsDoubleListSelector
{
    static public final String NO_SURFACE = " ";

    private StsModel model;
    private StsObject[] surfaces;
    private StsObject[] zones;
    private boolean[] fixedIndices;

    private StsColorListItem[] surfaceItems = null;
    private StsColorListItem[] zoneItems = null;

    private int lastSelectedIndex = 0;
    private String[] orderedNames = null;

    /** constructor */
    public StsSurfaceAssignmentDialog(StsWin3d win3d, StsObject[] allSurfaces,
            StsObject[] zones)
    {
        super((JFrame)win3d, "Reorder Surfaces & Assign to Zones",
            "Move surfaces Up & Down between zones:",
            "Surfaces:", null, "Zones:", null, true, 2.5f, false, true, "Move Up", "Move Down",
            new StsColorListItem(StsColorListDialog.PROTOTYPE_COLOR,
                StsDoubleListDialog.PROTOTYPE_STRING), new StsColorListRenderer());
        model = win3d.getModel();
        this.zones = zones;
        buildLists(allSurfaces);
        pack();
        adjustSize();
    }

    private void buildLists(StsObject[] allSurfaces)
    {
        try
        {
            if (zones == null || allSurfaces == null) return;

            // get only imported surfaces
            StsList surfaceList = new StsList(1, 1);
            for (int i=0; i<allSurfaces.length; i++)
            {
                if (allSurfaces[i].getName().startsWith(StsZone.ZONE_SURFACE_PREFIX)) continue;
                surfaceList.add(allSurfaces[i]);
            }
            surfaces = new StsModelSurface[surfaceList.getSize()];
            for (int i=0; i<surfaces.length; i++)
            {
                surfaces[i] = (StsModelSurface)surfaceList.getElement(i);
            }
            surfaceList = null;

            zoneItems = new StsColorListItem[zones.length];
            surfaceItems = new StsColorListItem[surfaces.length+zones.length+1];
            fixedIndices = new boolean[surfaceItems.length];
            for (int i=0; i<fixedIndices.length; i++) fixedIndices[i] = false;
            for (int i=0; i<surfaceItems.length; i++)
            {
                surfaceItems[i] = new StsColorListItem(StsColor.WHITE, NO_SURFACE);
            }
            StsList fixedSurfaces = new StsList(1, 1);
            StsSurface top = null;
            StsSurface base = null;
            for (int i=0; i<zoneItems.length; i++)
            {
                zoneItems[i] = new StsColorListItem(zones[i].getStsColor(), zones[i].getName(), 16, 16, 3);
                StsZone zone = (StsZone)zones[i];
                top = zone.getTopSurface();
                StsModelSurface topModel = zone.getTopModelSurface();
                if (top != null && top != base)
                {
                    if (topModel != null) surfaceItems[i] = new StsColorListItem(topModel.getStsColor(), topModel.getName(), 16, 16, 3);
                    else surfaceItems[i] = new StsColorListItem(top.getStsColor(), top.getName(), 16, 16, 3);
                    fixedIndices[i] = true;
                    if (!fixedSurfaces.contains(top)) fixedSurfaces.add(top);
                }
                base = zone.getBotSurface();
                StsModelSurface baseModel = zone.getBaseModelSurface();
                if (base != null)
                {
                    if (baseModel != null) surfaceItems[i+1] = new StsColorListItem(baseModel.getStsColor(), baseModel.getName(), 16, 16, 3);
                    else surfaceItems[i+1] = new StsColorListItem(base.getStsColor(), base.getName(), 16, 16, 3);
                    fixedIndices[i+1] = true;
                    if (!fixedSurfaces.contains(base)) fixedSurfaces.add(base);
                }

            }
            int j = zones.length+1;
            for (int i=0; i<surfaces.length; i++)
            {
                if (fixedSurfaces.contains(surfaces[i])) continue;
                surfaceItems[j] = new StsColorListItem(surfaces[i].getStsColor(), surfaces[i].getName());
                j++;
            }
            list.setListData(surfaceItems);
            list2.setListData(zoneItems);
        }
        catch(Exception e) { StsException.outputException(e, StsException.WARNING); }
    }

    protected void button1_actionPerformed(ActionEvent e)
    {
        // perform up action
        moveAction(true);
    }
    protected void button2_actionPerformed(ActionEvent e)
    {
        // perform down action
        moveAction(false);
    }

    private void moveAction(boolean up)
    {
        int[] index = list.getSelectedIndices();
        if (index==null) return;
        int i = index[0];
        if (i<0 || i>surfaceItems.length-1) return;
        if (fixedIndices[i]) return;
        if (surfaceItems[i].getName().equals(NO_SURFACE)) return;
        if (up)
        {
            int j = i - 1;
            for (; j>=0; j--)
            {
                if (!fixedIndices[j])
                {
                    StsMath.swap(surfaceItems, j, i);
                    lastSelectedIndex = j;
                    break;
                }
            }
            if (j < 0) return;
        }
        else  // down
        {
            int j = i + 1;
            for (; j<surfaceItems.length; j++)
            {
                if (!fixedIndices[j])
                {
                    StsMath.swap(surfaceItems, i, j);
                    lastSelectedIndex = j;
                    break;
                }
            }
            if (j >= surfaceItems.length) return;
        }
        list.setListData(surfaceItems);
        list.setSelectedIndex(lastSelectedIndex);
    }

    public String[] getSurfaceNames()
    {
        if (surfaceItems == null || zones == null) return null;
        int nSurfaces = zones.length + 1;
        String[] names = new String[nSurfaces];
        for (int i=0; i<nSurfaces; i++)
        {
            names[i] = surfaceItems[i].getName();
        }
        return names;
    }

    protected void cancelButton_actionPerformed(ActionEvent e)
    {
        surfaceItems = null;
        super.cancelButton_actionPerformed(e);
    }
}

