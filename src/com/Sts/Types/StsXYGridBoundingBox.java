//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

// add row and col limits to bounding box
public class StsXYGridBoundingBox extends StsSerialize implements Cloneable
{
    public int rowMin = largeInt; // row at yMin
    public int rowMax = -largeInt; // row at yMax
    public int colMin = largeInt; // col at xMin
    public int colMax = -largeInt; // col at xMax
	public int nRows = 0;
	public int nCols = 0;

	static public final int largeInt = StsParameters.largeInt;

	public StsXYGridBoundingBox()
    {
    }

    public StsXYGridBoundingBox(int rowMin, int rowMax, int colMin, int colMax)
    {
        this.rowMin = rowMin;
        this.rowMax = rowMax;
        this.colMin = colMin;
        this.colMax = colMax;
    }

    public StsGridBoundingBox getGridBoundingBox()
    {
        return new StsGridBoundingBox(this);
    }
    public StsXYGridBoundingBox cloneBox()
    {
        try
        {
            return (StsXYGridBoundingBox)clone();
        }
        catch (Exception e)
        {
            StsException.outputException("StsGridBoundingBox.getClone(grid) failed.", e, StsException.WARNING);
            return null;
        }
    }

    public StsXYGridBoundingBox unionBox(StsXYGridBoundingBox otherBox)
    {
        StsXYGridBoundingBox combinedBox = cloneBox();
        combinedBox.addBoundingBox(otherBox);
        return combinedBox;
    }

    public void initialize(StsXYGridBoundingBox boundingBox)
    {
        rowMin = boundingBox.rowMin;
        rowMax = boundingBox.rowMax;
        colMin = boundingBox.colMin;
        colMax = boundingBox.colMax;
    }

    public void initialize()
    {
        rowMin = largeInt;
        rowMax = -largeInt;
        colMin = largeInt;
        colMax = -largeInt;
    }

    public int getNRows()
    {
        if(rowMin == largeInt) return 0;
        return rowMax - rowMin + 1;
    }

    public int getNCols()
    {
        if(colMin == largeInt) return 0;
        return colMax - colMin + 1;
    }

    public int getRowMin() { return rowMin; }
    public int getRowMax() { return rowMax; }
    public int getColMin() { return colMin; }
    public int getColMax() { return colMax; }

    public boolean isInsideRowCol(float rowF, float colF)
    {
        return rowF >= rowMin && rowF <= rowMax && colF >= colMin && colF <= colMax;
    }

    public boolean isInsideRowCol(int row, int col)
    {
        return row >= rowMin && row <= rowMax && col >= colMin && col <= colMax;
    }

    public boolean isInsideGridRowCol(int row, int col)
    {
        return row >= 0 && row < nRows && col >= 0 && col < nCols;
    }

    public boolean addPoint(int row, int col)
    {
        if(rowMin == largeInt)
        {
            rowMin = row;
            rowMax = row;
            colMin = col;
            colMax = col;
            return true;
        }
        boolean changed = false;
        if(row < rowMin)
        {
            changed = true;
            rowMin = row;
            if(rowMax == -largeInt) rowMax = rowMin;
        }
        else if(row > rowMax)
        {
            changed = true;
            rowMax = row;
        }
        if(col < colMin)
        {
            changed = true;
            colMin = col;
            if(colMax == -largeInt) colMax = colMin;
        }
        else if(col > colMax)
        {
            changed = true;
            colMax = col;
        }
        return changed;
    }

    public boolean addPoint(float rowF, float colF)
    {
        if(rowMin == largeInt)
        {
            rowMin = StsMath.floor(rowF);
            rowMax = StsMath.ceiling(rowF);
            colMin = StsMath.floor(colF);
            colMax = StsMath.ceiling(colF);
            return true;
        }
        boolean changed = false;
         if(StsMath.floor(rowF) < rowMin)
         {
             changed = true;
             rowMin = StsMath.ceiling(rowF);
             if(rowMax == -largeInt) rowMax = rowMin;
         }
         else if(StsMath.ceiling(rowF) > rowMax)
         {
             changed = true;
             rowMax = StsMath.ceiling(rowF);
         }
         if(StsMath.floor(colF) < colMin)
         {
             changed = true;
             colMin = StsMath.floor(colF);
             if(colMax == -largeInt) colMax = colMin;
         }
         else if(StsMath.ceiling(colF) < colMax)
         {
             changed = true;
             colMax = StsMath.ceiling(colF);
         }
         return changed;
     }

    public boolean addBoundingBox(StsXYGridBoundingBox box)
    {
        boolean changed = false;
        if(box.rowMin < rowMin)
        {
            changed = true;
            rowMin = box.rowMin;
        }
        if(box.rowMax > rowMax)
        {
            changed = true;
            rowMax = box.rowMax;
        }
        if(box.colMin < colMin)
        {
            changed = true;
            colMin = box.colMin;
        }
        if(box.colMax > colMax)
        {
            changed = true;
            colMax = box.colMax;
        }
        return changed;       
    }

    public void intersectBoundingBox(StsXYGridBoundingBox box)
    {
        rowMin = Math.max(rowMin, box.rowMin);
        rowMax = Math.min(rowMax, box.rowMax);
        colMin = Math.max(colMin, box.colMin);
        colMax = Math.min(colMax, box.colMax);
    }

	public int getAbsRow(int relativeRow)
	{
		return relativeRow + rowMin;
	}

	public int getAbsCol(int relativeCol)
	{
		return relativeCol + colMin;
	}

    public String toString()
    {
        return new String("rowMin: " + rowMin + " rowMax: " + rowMax + " colMin: " + colMin + " colMax: " + colMax);
    }

    public static void main(String[] args)
    {
        StsXYGridBoundingBox box = new StsXYGridBoundingBox(10, 20, 30, 40);
        System.out.println("box: " + box.toString());
        StsXYGridBoundingBox cloneBox = box.cloneBox();
        System.out.println("cloned box: " + cloneBox.toString());
    }
}