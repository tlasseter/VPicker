
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Utilities.*;

public class StsGridPoint extends StsSerialize implements StsSerializable
{
    public StsPoint point; // elements are X, Y, time or approx-depth, depth
    public float[] nxyz;
    public float rowF, colF;
    public int row, col;
    public byte nullType = SURF_PNT;
    public StsXYGridable grid;
    public int rowOrCol;

   /* null type constants */
    static public final byte SURF_PNT = StsParameters.SURF_PNT;
    static public final byte SURF_BOUNDARY = StsParameters.SURF_BOUNDARY;
    static public final byte SURF_GAP = StsParameters.SURF_GAP;
    static public final byte SURF_GAP_SET = StsParameters.SURF_GAP_SET;
    static public final byte SURF_GAP_FILLED = StsParameters.SURF_GAP_FILLED;
//    static public final byte NULL_GAP_OR_BOUNDARY = StsParameters.NULL_GAP_OR_BOUNDARY;

    protected static final int NONE = StsParameters.NONE;
    protected static final int ROW = StsParameters.ROW;
    protected static final int COL = StsParameters.COL;
    protected static final int ROWCOL = StsParameters.ROWCOL;

    static public final String[] nullNames = new String[]
        {
        "none", "not", "gap set", "gap", "boundary", "gap filled"};

    public StsGridPoint()
    {
    }

    public StsGridPoint(StsXYSurfaceGridable grid)
    {
        this.grid = grid;
        point = new StsPoint(4);
    }

    public StsGridPoint(StsPoint point, float rowF, float colF)
    {
        this.point = point;
        setRowF(rowF);
        setColF(colF);
        setRowOrCol();
    }

    public StsGridPoint(double[] xyz, StsXYGridable grid)
    {
        point = new StsPoint(xyz);
        this.grid = grid;
        setGridRow();
        setGridCol();
        setRowOrCol();
    }

    public StsGridPoint(StsPoint point, StsXYGridable grid)
    {
        this.point = point;
        this.grid = grid;
        setGridRow();
        setGridCol();
        setRowOrCol();
    }

    public StsGridPoint(float rowF, float colF, StsXYGridable grid)
    {
        this.grid = grid;
        point = new StsPoint(4);
        setRowF(rowF);
        setColF(colF);
        setRowOrCol();
    }

    public StsGridPoint(int row, int col, StsXYGridable grid)
    {
        this.grid = grid;
        point = new StsPoint(4);
        setRow(row);
        setCol(col);
        rowOrCol = ROWCOL;
    }
/*
    public StsGridPoint(StsPoint point, float rowF, float colF, StsXYGridable grid)
    {
        this.grid = grid;
        pxyz = new float[4];
        pxyz[0] = point.getX();
        pxyz[1] = point.getY();
        pxyz[2] = point.getT();
        pxyz[3] = point.getZ();
        this.rowF = rowF;
        this.colF = colF;
        row = Math.round(rowF);
        col = Math.round(colF);
        setRowOrCol();
    }
*/
    public StsGridPoint(StsGridCrossingPoint gridCrossingPoint)
    {
        point = gridCrossingPoint.point;
        setRowF(gridCrossingPoint.iF);
        setColF(gridCrossingPoint.jF);
        setRowOrCol();
    }

    public void setRowCol(int row, int col)
    {
        setRow(row);
        setCol(col);
        rowOrCol = ROWCOL;
    }

    private void setGridRow()
    {
        if (point == null)return;
        if (grid == null)return;
        setRowF(grid.getRowCoor(getPoint().v));
//        row = StsMath.floor(rowF);
//        row = StsMath.roundOffInteger(rowF);
    }

    public void setY(float y)
    {
        if (point == null)return;
        getPoint().setY(y);
        setGridRow();
        setRowOrCol();
    }

    public float getY()
    {
        return point.getY();
    }

    // If we have changed the grid basis, then reset grid for this point
    public void setY(float y, StsXYGridable grid)
    {
        this.grid = grid;
        if(point == null) return;
        point.setY(y);
        setGridRow();
        setRowOrCol();
    }

    private void setGridCol()
    {
        if (point == null)return;
        if (grid == null)return;
        setColF(grid.getColCoor(point.v));
//        col = StsMath.floor(colF);
//        col = StsMath.roundOffInteger(colF);
    }

    public void setX(float x)
    {
        if (point == null)return;
        point.setX(x);
        setGridCol();
        setRowOrCol();
    }

    // If we have changed the grid basis, then reset grid for this point
    public void setX(float x, StsXYGridable grid)
    {
        this.grid = grid;
        if (point == null)return;
        point.setX(x);
        setGridCol();
        setRowOrCol();
    }

    public float getX()
    {
        return point.getX();
    }

    public void setRow(int row)
    {
        this.row = row;
        rowF = (float) row;
        if (grid != null) point.setY(grid.getYCoor(rowF, colF));
        setRowOrCol();
    }

    public void setCol(int col)
    {
        this.col = col;
        colF = (float)col;
        if(grid != null) point.setX(grid.getXCoor(rowF, colF));
        setRowOrCol();
    }

    public void setRowF(float rowF)
    {
        this.row = (int)rowF;
		if(StsMath.sameAs((float)row, rowF))
			this.rowF = row;
		else
			this.rowF = rowF;

        if(grid != null) point.setY(grid.getYCoor(rowF, colF));
    }

    public void setColF(float colF)
    {
        this.col = (int)colF;
		if(StsMath.sameAs((float)col, colF))
			this.colF = col;
		else
			this.colF = colF;
        if(grid != null) point.setX(grid.getXCoor(rowF, colF));
    }

    public void setRowColF(float rowF, float colF)
    {
        setRowF(rowF);
        setColF(colF);
        setRowOrCol();
    }

    public void setZ(float z) { point.setZ(z); }
    public float getZ() { return point.getZ(); }
    public void setT(float t) { point.setT(t); }
    public float getT() { return point.getT(); }
    public void setZorT(boolean isDepth, float zt)
    {
        point.setZorT(zt, isDepth);
    }

    public float[] getXYZorT()
    {
        return getXYZorT(isDepth);
    }

    public float[] getXYZorT(boolean isDepth)
    {
        return point.getXYZorT(isDepth);
    }

    public StsPoint getXYZorTPoint()
    {
        return getXYZorTPoint(isDepth);
    }

    public StsPoint getXYZorTPoint(boolean isDepth)
    {
        if(point.length() < 4) return point;
        if(isDepth)
            return new StsPoint( getX(), getY(), getZ() );
        else
            return new StsPoint( getX(), getY(), getT() );
    }

    public float getZorT() { return getZorT(isDepth); }

    public float getZorT(boolean isDepth)
    {
        return getPoint().getZorT(isDepth);
    }

    public void setNormal(float[] normal) { nxyz = normal; }

    public void setNullType(byte nullType) { this.nullType = nullType; }

    public boolean notNull() { return nullType == SURF_PNT; }
/*
    public float[] setGetLeftNormal(float[] vA, float[] vB)
    {
        float[] normal = new float[3];

        normal[0] = -vA[1]*vB[2] + vA[2]*vB[1];
        normal[1] = -vA[2]*vB[0] + vA[0]*vB[2];
        normal[2] = -vA[0]*vB[1] + vA[1]*vB[0];

        return normal;
    }
*/
    public void interpolate(StsGridPoint start, StsGridPoint end, float f)
    {
        setPoint(StsPoint.staticInterpolatePoints(start.point, end.point, f));
        setRowF(start.rowF + f*(end.rowF - start.rowF));
        setColF(start.colF + f*(end.colF - start.colF));
        if(start.grid == end.grid) this.grid = start.grid;
        setRowOrCol();
    }

    public void adjustToNearestRowCol()
    {
        colF = (float)col;
        getPoint().setX(grid.getXCoor(rowF, colF));
        rowF = (float)row;
        getPoint().setY(grid.getYCoor(rowF, colF));
        rowOrCol = ROWCOL;
    }

	public final StsPoint getPoint()
	{
		return point;
	}

    void setRowOrCol()
    {
        boolean isRow = StsMath.isExactIntegral(rowF);
        boolean isCol = StsMath.isExactIntegral(colF);

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

    public int[] getLowerLeftRowCol(StsGridPoint prevGridPoint)
     {
		 int rowLL = Math.min(row, prevGridPoint.row);
		 int colLL = Math.min(col, prevGridPoint.col);
		 return new int[] { rowLL, colLL };
	}

    public void display(StsGLPanel3d glPanel3d)
    {
		StsGLDraw.drawPoint(getPoint().getXYZorT(), StsColor.WHITE, glPanel3d,
                            3, StsGraphicParameters.vertexShift);
    }

    public void setPoint(StsPoint point)
    {
        this.point = point;
    }

    static public int getSameRow(StsGridPoint g0, StsGridPoint g1)
    {
        if(g0.row == g1.row)
            return g0.row;
        else
            return -1;
    }

    static public int getSameCol(StsGridPoint g0, StsGridPoint g1)
    {
        if(g0.col == g1.col)
            return g0.col;
        else
            return -1;
    }

	static public boolean isSameZ(StsGridPoint g0, StsGridPoint g1)
    {
		return g0.getZ() == g1.getZ();
    }
    public String toString()
    {
        if(point != null)
            return point.toString() + " rowF: " + rowF + " colF: " + colF;
        else
            return "undefined";
    }
}
