package com.Sts.Types;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.View3d.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsCursorPoint extends StsSerialize implements StsSerializable
{
    public int dirNo = -1;
    public StsPoint point;
    public float rowNum;
    public float colNum;

    public StsCursorPoint()
    {
    }

    public StsCursorPoint(int dirNo, StsPoint point)
    {
        this.dirNo = dirNo;
        this.point = point;
    }

    public StsCursorPoint(int dirNo, float dirCoordinate, boolean axesFlipped, StsPoint point)
    {
        this.dirNo = dirNo;
        this.point = point;
        if(dirNo == StsCursor3d.XDIR)
        {
            point.v[2] = point.v[1];
            point.v[1] = point.v[0];
            point.v[0] = dirCoordinate;
        }
        else if(dirNo == StsCursor3d.YDIR)
        {
            point.v[2] = point.v[1];
            point.v[1] = dirCoordinate;
        }
        else if(dirNo == StsCursor3d.ZDIR)
        {
            if(axesFlipped)
            {
                float temp = point.v[0];
                point.v[0] = point.v[1];
                point.v[1] = temp;
            }
            point.v[2] = dirCoordinate;
        }
    }

    public StsCursorPoint(int dirNo, float[] rowColF, StsPoint point)
    {
        this.dirNo = dirNo;
        this.rowNum = rowColF[0];
        this.colNum = rowColF[1];
        this.point = point;
    }

    public int[] getRowCol()
    {
        int row = Math.round(rowNum);
        int col = Math.round(colNum);
        return new int[] { row, col };
    }
}
