
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

abstract public class StsRowCol extends StsSerialize implements StsSerializable
{
    public int rowOrCol = NONE;      /** @param rowOrCol point is on row or col or both */
    public float rowF = nullValue;   /** @param rowF row coordinate value    */
    public float colF = nullValue;   /** @param colF col coordinate value    */

    // State flags for connect
    static public final int MINUS = StsParameters.MINUS;
    static public final int PLUS = StsParameters.PLUS;

    // State flags for rowOrCol
    static public final int ROW = StsParameters.ROW;
    static public final int COL = StsParameters.COL;
    static public final int ROWCOL = StsParameters.ROWCOL;
    static public final int NONE = StsParameters.NONE;

	static final float nullValue = StsParameters.nullValue;

    public StsRowCol()
    {
    }

	public StsRowCol(float rowF, float colF)
	{
		this.rowF = rowF;
		this.colF = colF;
		computeRowOrCol();
	}

    public StsRowCol(StsRowCol rowCol0, StsRowCol rowCol1, double f)
    {
        float rowF0 = rowCol0.getRowF();
        float rowF1= rowCol1.getRowF();
        rowF = rowF0 + (float)f*(rowF1 - rowF0);
        float colF0 = rowCol0.getColF();
        float colF1= rowCol1.getColF();
        colF = colF0 + (float)f*(colF1 - colF0);
        computeRowOrCol();
    }

    public int getRowOrCol() { return rowOrCol; }
    public float getRowF() { return rowF; }
    public float getColF() { return colF; }

    public boolean isRowOrCol() { return rowOrCol != NONE; }
    public boolean isRowAndCol() { return rowOrCol == ROWCOL; }
    public boolean isRow() { return rowOrCol == ROW || rowOrCol == ROWCOL; }
    public boolean isCol() { return rowOrCol == COL || rowOrCol == ROWCOL; }

    public int getRow() { return StsMath.roundOffInteger(rowF); }
    public int getCol() { return StsMath.roundOffInteger(colF); }

    public boolean isRowOrCol(int matchRowOrCol)
    {
        if(matchRowOrCol == rowOrCol) return true;

        if(matchRowOrCol == ROW)
            return rowOrCol == ROWCOL;
        else if(matchRowOrCol == COL)
            return rowOrCol == ROWCOL;
        else
            return false;
    }

    static public boolean isRow(int rowOrCol) { return rowOrCol == ROW || rowOrCol == ROWCOL; }
    static public boolean isCol(int rowOrCol) { return rowOrCol == COL|| rowOrCol == ROWCOL; }

    public boolean isThisRow(int row)
    {
        return isRow() && rowF == row;
    }

    public boolean isThisCol(int col)
    {
        return isCol() && colF == col;
    }

    public boolean isSameRowAndCol(StsRowCol otherRowCol)
    {
        if(otherRowCol == null) return false;
        if(rowF != otherRowCol.getRowF()) return false;
        return colF == otherRowCol.getColF();
    }

    public boolean isSameRowAndCol(int row, int col)
    {
        if(rowF != (float)row) return false;
        return colF == (float)col;
    }

    public int getSameRowOrCol(StsRowCol otherRowCol)
    {
        if(otherRowCol == null) return NONE;
        boolean isRow = rowF == otherRowCol.getRowF() && isRow();
        boolean isCol = colF == otherRowCol.getColF() && isCol();
        if(isRow && isCol)
            return ROWCOL;
        else if(isRow)
            return ROW;
        else if(isCol)
            return COL;
        else
            return NONE;
    }

    public boolean isSameRowOrCol(StsRowCol otherRowCol)
    {
        return getSameRowOrCol(otherRowCol) != NONE;
    }

    public float getRowColIndexF(int rowOrCol)
    {
        if(rowOrCol == ROW) return rowF;
        else if(rowOrCol == COL) return colF;
        else return nullValue;
    }

     public int getRowCol(int rowOrCol)
     {
         if (rowOrCol == ROW) return (int)rowF;
         else if (rowOrCol == COL) return (int)colF;
         else return -1;
     }

    public float getCrossingRowColF()
    {
        return getCrossingRowColF(rowOrCol);
    }

    public float getCrossingRowColF(int rowOrCol)
    {
        if(rowOrCol == ROW) return colF;
        else if(rowOrCol == COL) return rowF;
        else return nullValue;
    }

    public void computeRowOrCol()
    {
        boolean isRow = StsMath.isExactIntegral(rowF);
        boolean isCol = StsMath.isExactIntegral(colF);

//        if(isRow) rowF = StsMath.roundOffInteger(rowF);
//        if(isCol) colF = StsMath.roundOffInteger(colF);

        if(isRow)
        {
            if(isCol) rowOrCol = ROWCOL;
            else      rowOrCol = ROW;
        }
        else if(isCol)
            rowOrCol = COL;
        else
            rowOrCol = NONE;
    }

    public void setRowColF(float rowF, float colF)
    {
        this.rowF = rowF;
        this.colF = colF;
        computeRowOrCol();
    }

    public void setRowColF(int rowOrCol, float rowF, float colF)
    {
        this.rowOrCol = rowOrCol;
        this.rowF = rowF;
        this.colF = colF;
    }

    public void setRowF(float rowF)
    {
        this.rowF = rowF;
        computeRowOrCol();
    }

    public void setColF(float colF)
    {
        this.colF = colF;
        computeRowOrCol();
    }

    public void setRowOrCol(int setRowOrCol)
    {
        if(rowOrCol == NONE)
            rowOrCol = setRowOrCol;
        else if(rowOrCol != setRowOrCol)
            rowOrCol = ROWCOL;

        if(setRowOrCol == ROW || setRowOrCol == ROWCOL)
            rowF = StsMath.roundOffInteger(rowF);
        else if(setRowOrCol == COL || setRowOrCol == ROWCOL)
            colF = StsMath.roundOffInteger(colF);
    }

	public void setRowOrColIndex(int setRowOrCol, int index)
	{
		setRowOrColIndex(setRowOrCol, index, true);
	}

	public void setRowOrColIndex(int setRowOrCol, int index, boolean checkRedundant)
	{
		// mainDebug check
//		if(checkRedundant) debugCheckIsRedundant(setRowOrCol, index);
        setRowOrColIndexF(setRowOrCol, (float)index);
    }

    public void setRowOrColIndexF(int setRowOrCol, float indexF)
    {
//		if(rowOrCol == setRowOrCol) return;
        if(rowOrCol == NONE) rowOrCol = setRowOrCol;
//        if(rowOrCol == NONE || rowOrCol == ROWCOL) rowOrCol = setRowOrCol;
        else if(rowOrCol != setRowOrCol) rowOrCol = ROWCOL;

        if(setRowOrCol == ROW) rowF = indexF;
        else if(setRowOrCol == COL) colF = indexF;
    }

    private void debugCheckIsRedundant(int setRowOrCol, int index)
    {

		boolean isRedundant;

		isRedundant = isRowOrCol(setRowOrCol);

		if(setRowOrCol == ROW || setRowOrCol == ROWCOL)
		{
			if(rowF != (float)index) isRedundant = false;
		}
		if(setRowOrCol == COL || setRowOrCol == ROWCOL)
		{
			if(colF != (float)index) isRedundant = false;
		}
		if(isRedundant)
		{
			String type;
			if(this instanceof StsSectionRowCol) type = "Section";
			else type = "Grid";

			System.out.println("Redundant: " + type +
				" rowOrCol: " + StsParameters.rowCol(rowOrCol) + " row: " + rowF + " col: " + colF +
				" setRowOrCol: " + StsParameters.rowCol(rowOrCol) + " value: " + index);
		}
    }

    public boolean matchesRowCol(int[] rowColIndexes)
    {
		return rowF == (float)rowColIndexes[0] && colF == (float)rowColIndexes[1];
    }

    public float getGridDistance(int rowOrCol, StsRowCol nextRowCol)
    {
        float rowColF0 = getRowColIndexF(rowOrCol);
        float rowColF1 = nextRowCol.getRowColIndexF(rowOrCol);
        return Math.abs(rowColF1 - rowColF0);
    }

	public float[] getRowAndColF()
	{
		return new float[] { rowF, colF };
	}

	public int[] getRowAndCol()
	{
		return new int[] { (int)rowF, (int)colF };
	}
    
    public float getGridDistance(StsRowCol nextRowCol)
    {
        float dRow = nextRowCol.getRowF()- rowF;
        float dCol = nextRowCol.getColF() - colF;
        return (float)Math.sqrt(dRow*dRow + dCol*dCol);
    }

    public int[] getLowerLeftRowCol(StsRowCol nextRowCol)
    {
        int row = (int)rowF;
        int col = (int)colF;

        if(rowOrCol == NONE) return new int[] { row, col };
        if(nextRowCol == null) return new int[] { row, col };

        float nextRowF = nextRowCol.getRowF();
        float nextColF = nextRowCol.getColF();

        if(isRow())
        {
            if(nextRowF < rowF) row--;
            else if(nextRowF == rowF && nextColF > colF) row--;
        }
        if(isCol())
        {
            if(nextColF < colF) col--;
            else if(nextColF == colF && nextRowF < rowF) col--;
        }
        return new int[] { row, col };
    }

    static public int[] getLowerLeftRowCol(StsRowCol[] rowCols)
    {
        int nPoints = rowCols.length;
        int row = (int)rowCols[0].rowF;
        int col = (int)rowCols[0].colF;
        for(int n = 1; n < nPoints; n++)
        {
            row = Math.min(row, (int)rowCols[n].rowF);
            col = Math.min(col, (int)rowCols[n].colF);
        }
        return new int[] { row, col };
    }

    // Compute distance from point to line; if intersection off end of line, return false;
	// if distance/lineLength > fracDistance, return false; otherwise return true.
	public boolean between(StsRowCol gridRowCol0, StsRowCol gridRowCol1, float fracDistance)
	{
		float rowF0, rowF1;
		float colF0, colF1;
		float[] v; // vector from gridPoints 0 to 1;
		float[] vp;  // vector from gridPoint 0 to p;

		rowF0 = gridRowCol0.getRowF();
		colF0 = gridRowCol0.getColF();
		rowF1 = gridRowCol1.getRowF();
		colF1 = gridRowCol1.getColF();
		v = new float[] { rowF1-rowF0, colF1-colF0 };
		vp =  new float[] { rowF-rowF0, colF-colF0 };
		return StsMath.isOnLine(vp, v, fracDistance);
	}

    abstract public void recompute(StsPoint point);
	abstract public String toString();
}
