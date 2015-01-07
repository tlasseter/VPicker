package com.Sts.SeismicAttributes;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.Utilities.Seismic.*;
import com.Sts.Utilities.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jun 12, 2009
 * Time: 9:24:52 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsPatchPoint extends StsSerialize implements StsSerializable, Comparable<StsPatchPoint>, Cloneable
{
    public float value;
    public float z = StsParameters.nullValue;
    public byte pointType;
    public int patchID = -1;
    public int row;
    public int col;
    public int slice;
    public float rowCorrel;
    public float colCorrel;

    public StsPatchPoint()
    {
    }

    public StsPatchPoint(int row, int col, int slice, float z, float value, byte pointType)
    {
        this.row = row;
        this.col = col;
        this.slice = slice;
        this.z = z;
        this.value = value;
        this.pointType = pointType;
    }

    public int compareTo(StsPatchPoint otherPoint)
    {
        if(slice > otherPoint.slice) return 1;
        else if (slice < otherPoint.slice) return -1;
        else return 0;
    }

    public StsPatchPoint clone()
    {
        try
        {
            return(StsPatchPoint)super.clone();
        }
        catch(Exception e)
        {
            StsException.systemError(this, "clone");
            return null;
        }
    }

    public StsPatchPoint cloneClearCorrels()
    {
        StsPatchPoint point = clone();
        point.rowCorrel = 0.0f;
        point.colCorrel = 0.0f;
        return point;
    }

    public Integer hashCode(int nVolumeCols)
    {
        return new Integer(row*nVolumeCols + col);       
    }

    public int getSlice() { return slice; }

    public int getIndex(int nVolumeCols)
    {
        return col + row*nVolumeCols;
    }

    public String toString()
    {
        return "id: " + patchID + " row: " + row + " col: " + col + " slice: " + slice + " value: " + value +
                " type: " + StsTraceUtilities.typeStrings[pointType];
    }
}
