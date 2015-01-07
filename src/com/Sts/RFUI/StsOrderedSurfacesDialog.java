
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.RFUI;

import com.Sts.DBTypes.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Utilities.*;

import java.awt.event.*;

public class StsOrderedSurfacesDialog extends StsOrderedListDialog
{
    public StsOrderedSurfacesDialog(StsWin3d win3d, String[] surfaceNames)
    {
        super(win3d, "Reorder surfaces by stratigraphy", "Surfaces (with center depth & associated zones):",
                StsModelSurface.class, surfaceNames);
    }

    protected void okayButton_actionPerformed(ActionEvent e)
    {
        StsSurfaceClass surfaceClass = (StsSurfaceClass)model.getCreateStsClass(StsModelSurface.class);
        surfaceClass.reorderListByNames(getOrderedNames());
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
                StsModelSurface s = null;
                StsSurfaceClass surfaceClass = (StsSurfaceClass)model.getCreateStsClass(StsModelSurface.class);
                s = (StsModelSurface)surfaceClass.getObjectWithName(names[i]);
                if (s==null) continue;
                float centerZ = s.getOrderingValue();
                if (centerZ != StsParameters.nullValue)
                {
                    names[i] += "  (" + ((int)centerZ) + ")";
                }
                StsZone z = s.getZoneAbove();
                if (z != null) names[i] += "  :  base of " + z.getName();
                z = s.getZoneBelow();
                if (z != null) names[i] += "  :  top of " + z.getName();
                //StsSurface relatedSurface = s.getRelatedSurface();
                //if (relatedSurface == null ||
                //    !relatedSurface.isType(StsModelSurface.MODEL)) continue;
                //StsZone z = relatedSurface.getZoneAbove();
                //if (z != null) names[i] += "  :  base of " + z.getName();
                //z = relatedSurface.getZoneBelow();
                //if (z != null) names[i] += "  :  top of " + z.getName();
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
            StsSurfaceClass surfaceClass = (StsSurfaceClass)model.getCreateStsClass(StsModelSurface.class);
            StsModelSurface surface = (StsModelSurface)surfaceClass.getObjectWithName(name);
            if (surface==null) return false;
            if (surface.getZoneAbove() != null ||
                surface.getZoneBelow() != null) return false;
            /*
            StsSurface relatedSurface = surface.getRelatedSurface();
            if (relatedSurface==null ||
                !relatedSurface.isType(StsModelSurface.MODEL)) return true;
            if (relatedSurface.getZoneAbove() != null ||
                relatedSurface.getZoneBelow() != null) return false;
            */
        }
        catch (Exception e) { return false; }
        return true;
    }
}

