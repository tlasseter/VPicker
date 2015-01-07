
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.Utilities.*;

public class StsGridRowCol extends StsRowCol implements StsSerializable
{
//    static public StsXYSurfaceGridable horizonGrid;

    public StsGridRowCol()
    {
    }

    public StsGridRowCol(StsPoint point)
    {
//        if(horizonGrid == null) horizonGrid = currentModel.getGridDefinition();
        recompute(point);
    }

    public StsGridRowCol(float rowF, float colF)
    {
        this.rowF = rowF;
        this.colF = colF;
        computeRowOrCol();
    }

    public StsGridRowCol(int rowOrCol, float rowF, float colF)
    {
        this.rowOrCol = rowOrCol;
        this.rowF = rowF;
        this.colF = colF;
    }

    public StsGridRowCol(StsXYSurfaceGridable grid, float[] point)
    {
       rowF = grid.getRowCoor(point);
       colF = grid.getColCoor(point);
       computeRowOrCol();
    }

    public StsGridRowCol(StsGridRowCol rowCol0, StsGridRowCol rowCol1, double f)
    {
		super(rowCol0, rowCol1, f);
    }

    public float compare(int rowOrCol, StsRowCol otherRowCol)
    {
        float rowColF = getCrossingRowColF(rowOrCol);
        float otherRowColF = otherRowCol.getCrossingRowColF(rowOrCol);
        return rowColF - otherRowColF;
    }

    static public StsGridRowCol getCommonGridRowCol(StsRowCol rowCol, StsRowCol otherRowCol)
    {
        if(rowCol == null || otherRowCol == null) return null;
        int sameRowOrCol = rowCol.getSameRowOrCol(otherRowCol);
        if(sameRowOrCol == rowCol.rowOrCol)
            return new StsGridRowCol(sameRowOrCol, rowCol.rowF, rowCol.colF);
        else if(sameRowOrCol == otherRowCol.rowOrCol)
            return new StsGridRowCol(sameRowOrCol, otherRowCol.rowF, otherRowCol.colF);
        else if(sameRowOrCol == StsGridRowCol.NONE)
        {
            // StsException.systemError(StsGridRowCol.class, "getCommonGridRowCol", " gridRowCols ambiguous for " + rowCol.toString() + " and " + otherRowCol.toString());
            float rowF = (rowCol.rowF + otherRowCol.rowF)/2;
            float colF = (rowCol.colF + otherRowCol.colF)/2;
            return new StsGridRowCol(sameRowOrCol, rowF, colF);
        }
        else if(sameRowOrCol == StsGridRowCol.ROW)
        {
            // StsException.systemError(StsGridRowCol.class, "getCommonGridRowCol", " gridRowCols ambiguous for " + rowCol.toString() + " and " + otherRowCol.toString());
            float colF = (rowCol.colF + otherRowCol.colF)/2;
            return new StsGridRowCol(sameRowOrCol, rowCol.rowF, colF);
        }
        else if(sameRowOrCol == StsGridRowCol.COL)
        {
            float rowF = (rowCol.rowF + otherRowCol.rowF)/2;
            return new StsGridRowCol(sameRowOrCol, rowF, rowCol.colF);
        }
        else
        {
            StsException.systemError(StsGridRowCol.class, "getCommonGridRowCol");
            float rowF = (rowCol.rowF + otherRowCol.rowF)/2;
            float colF = (rowCol.colF + otherRowCol.colF)/2;
            return new StsGridRowCol(sameRowOrCol, rowF, colF);
        }
    }
/*
    public int[] getLowerLeftRowCol(StsRowCol nextRowCol)
    {
        float nextRowF, nextColF;

        int row = (int)rowF;
        int col = (int)colF;

        if(rowOrCol == NONE) return new int[] { row, col };

        if(nextRowCol == null) return new int[] { row, col };

        if(rowOrCol == ROW)
        {
            nextRowF = nextRowCol.getRowF();
            if(nextRowF < rowF) row--;
        }
        else if(rowOrCol == COL)
        {
            nextColF = nextRowCol.getColF();
            if(nextColF < colF) col--;
        }
        else // ROWCOL
        {
            nextRowF = nextRowCol.getRowF();
            nextColF = nextRowCol.getColF();
            if(nextRowF < rowF) row--;
            else if(nextRowF == rowF && nextColF > colF) row--;
            if(nextColF < colF) col--;
            else if(nextColF == colF && nextRowF < rowF) col--;
        }
        return new int[] { row, col };
    }
*/
    public void recompute(StsPoint point)
    {
        try
        {
            StsXYGridable horizonGrid = currentModel.getGridDefinition();
            if(horizonGrid == null) return;
            rowF = horizonGrid.getRowCoor(point.v);
            colF = horizonGrid.getColCoor(point.v);
            computeRowOrCol();
        }
        catch(Exception e)
        {
            StsException.outputException("StsGridRowCol.recompute() failed.",
                e, StsException.WARNING);
        }
    }

    public String toString()
    {
        return new String("gridRowCol at: rowF " + rowF + " colF: " + colF + " on: " + StsParameters.rowCol(rowOrCol));
    }
}
