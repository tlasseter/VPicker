
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;
import com.Sts.MVC.View3d.*;

public class StsCursorSectionRowCol extends StsRowCol implements StsSerializable
{
    protected StsCursorSection section;

    public StsCursorSectionRowCol()
    {
    }

    public StsCursorSectionRowCol(float rowF, float colF, StsCursorSection section)
    {
        super(rowF, colF);
        this.section = section;
    }

    public StsCursorSectionRowCol(StsPoint point, StsCursorSection section)
    {
		this.section = section;
        recompute(point);
    }

    public void recompute(StsPoint point)
    {
        try
        {
            float[] xyz = point.getXYZorT();
            rowF = section.getRowCoor(xyz);
            colF = section.getColCoor(xyz);
            computeRowOrCol();
        }
        catch(Exception e)
        {
            StsException.outputException("StsSectionRowCol.recompute() failed.",
                e, StsException.WARNING);
        }
    }

    public String toString()
    {
        return new String("StsSectionRowCol for cursorSection dir: " + section.dir +
                          " at: rowF " + rowF + " colF: " + colF +
            " on: " + StsParameters.rowCol(rowOrCol));
    }
}