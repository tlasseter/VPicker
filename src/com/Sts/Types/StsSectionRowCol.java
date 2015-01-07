
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

public class StsSectionRowCol extends StsRowCol implements StsSerializable
{
    protected StsSection section;

    public StsSectionRowCol()
    {
    }

    public StsSectionRowCol(float rowF, float colF, StsSection section)
    {
        super(rowF, colF);
        this.section = section;
    }

    public StsSectionRowCol(StsSectionRowCol rowCol0, StsSectionRowCol rowCol1, double f)
    {
		super(rowCol0, rowCol1, f);
        this.section = rowCol0.getSection();
    }

    public StsSectionRowCol(StsPoint point, StsSection section)
    {
		this.section = section;
        recompute(point);
    }

    public StsSectionRowCol(StsGridSectionPoint gridPoint, StsSection section)
    {
        StsLine line = null;

        if(gridPoint == null || section == null) return;
        StsSurfaceVertex vertex = gridPoint.getVertex();
        if(vertex != null) line = vertex.getSectionLine();
        initialize(gridPoint, section, line);
    }

    private void initialize(StsGridSectionPoint gridPoint, StsSection section, StsLine line)
    {
        this.section = section;
        StsPoint point = gridPoint.getPoint();
        this.colF = section.getColF(point, line);
        float z = point.getZorT();
        z = StsMath.checkRoundOffInteger(z);
        point.setZorT(z);
        this.rowF = section.getRowF(z);
		computeRowOrCol();
    }

    public StsSectionRowCol(StsGridSectionPoint gridPoint, StsSection section, StsLine line)
    {
        initialize(gridPoint, section, line);
    }


    public StsSection getSection() { return section; }

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
        return new String("StsSectionRowCol for section: " + section.getIndex() +
                          " at: rowF " + rowF + " colF: " + colF + " on: " + StsParameters.rowCol(rowOrCol));
    }
}
