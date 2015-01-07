package com.Sts.MVC.View3d;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.1
 */

import com.Sts.DBTypes.*;
import com.Sts.Types.*;

public class StsCrossplotPoint
{
    public StsCursorPoint cursorPoint;
    public StsColor stsColor;
    public int crossplotRow; // row in XP display
    public int crossplotCol; // col in XP display

    public StsCrossplotPoint(StsCursorPoint cursorPoint, int dataRow, int dataCol, StsColor color, byte colValue, byte rowValue)
    {
        this.cursorPoint = cursorPoint;
        stsColor = color;
        this.crossplotCol = dataCol;
        this.crossplotRow = dataRow;
    }

    public StsCrossplotPoint(StsColor color, int row, int col, byte colValue, byte rowValue)
    {
        stsColor = color;
        this.crossplotCol = col;
        this.crossplotRow = row;
    }

    public int getCursorRow() { return Math.round(cursorPoint.rowNum); }
    public int getCursorCol() { return Math.round(cursorPoint.colNum); }

    public float getCursorRowF() { return cursorPoint.rowNum; }
    public float getCursorColF() { return cursorPoint.colNum; }

    public float[] getVolumeXYZ() { return cursorPoint.point.v; }
}