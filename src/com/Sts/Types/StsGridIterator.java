
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.Utilities.*;

public class StsGridIterator
{
    StsEdgeLoopRadialLinkGrid linkedGrid;
    public int rowMin, rowMax;
    public int row, nextRow;
    public int colMin, colMax, nextColMin, nextColMax;

    static public final int nullInteger = StsParameters.nullInteger;

    public StsGridIterator(StsEdgeLoopRadialLinkGrid linkedGrid)
    {
        if(linkedGrid == null)
        {
            StsException.outputException(new StsException(StsException.WARNING,
                "StsGridIterator constructTraceAnalyzer failed: linkedGrid argument is null."));
            return;
        }

        this.linkedGrid = linkedGrid;

        rowMin = linkedGrid.getRowMin();
        rowMax = linkedGrid.getRowMax();

        nextRow = rowMin;
        nextColMin = linkedGrid.getColMin(nextRow);
        nextColMax = linkedGrid.getColMax(nextRow);
    }

    public StsGridIterator(StsXYSurfaceGridable grid)
    {
        if(grid == null)
        {
            StsException.outputException(new StsException(StsException.WARNING,
                "StsGridIterator constructTraceAnalyzer failed: grid argument is null."));
            return;
        }

        rowMin = 0;
        rowMax = grid.getNRows()-1;

        colMin = 0;
        colMax = grid.getNCols()-1;

        nextRow = 0;
        nextColMin = 0;
        nextColMax = colMax;
    }

    public int getNRows()
    {
        return rowMax-rowMin+1;
    }

    public int getNextRow()
    {
        row = nextRow;
        nextRow++;
        if(nextRow > rowMax) return nullInteger;

        if(linkedGrid != null)
        {
            colMin = nextColMin;
            if(colMin == nullInteger) return nullInteger;
            colMax = nextColMax;
            if(colMax == nullInteger) return nullInteger;

            nextColMin = linkedGrid.getColMin(nextRow);
            nextColMax = linkedGrid.getColMax(nextRow);
        }

        return row;
    }

    public boolean hasRow(int row)
    {
        if(linkedGrid == null)
            return row >= rowMin && row <= rowMax;
        else
        {
            if(!linkedGrid.hasRow(row)) return false;

            colMin = linkedGrid.getColMin(row);
            colMax = linkedGrid.getColMax(row);
            return true;
        }
    }
}

