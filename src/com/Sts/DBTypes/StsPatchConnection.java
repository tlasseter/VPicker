package com.Sts.DBTypes;

import com.Sts.SeismicAttributes.*;
import com.Sts.Utilities.*;

/**
 * Created by IntelliJ IDEA.
* User: Tom Lasseter
* Date: Sep 8, 2009
* Time: 2:55:17 PM
* To change this template use File | Settings | File Templates.
*/
public class StsPatchConnection implements Comparable<StsPatchConnection>
{
    StsPatchPoint newPoint;
    StsPatchPoint otherPoint;
    int nSheet = -1;
    double correl;

    public StsPatchConnection(StsPatchPoint newPoint, StsPatchPoint otherPoint, double correl, double distance)
    {
        initialize(newPoint, otherPoint, correl, distance);
    }

    void initialize(StsPatchPoint newPoint, StsPatchPoint otherPoint, double correl, double distance)
    {
        this.newPoint = newPoint;
        this.otherPoint = otherPoint;
        this.correl = correl;
    }

    void addCorrelation()
    {
         if(isRow())
            otherPoint.rowCorrel = (float)correl;
        else
            otherPoint.colCorrel = (float)correl;

    }

    public void clearCorrelation()
    {
        if(isRow())
            otherPoint.rowCorrel = 0.0f;
        else
            otherPoint.colCorrel = 0.0f;
    }

    public String toString()
    {
        float cor = 0.0f;
        if(isRow())
            return "newPoint: " + newPoint.toString() + " otherPoint: " + otherPoint.toString() + " sheet: " + nSheet + " rowCorrel: " + correl;
        else
            return "newPoint: " + newPoint.toString() + " otherPoint: " + otherPoint.toString() + " sheet: " + nSheet + " colCorrel: " + correl;
    }

    boolean isRow() { return otherPoint.row == newPoint.row; }

    public int compareTo(StsPatchConnection otherConnection)
    {
        StsPatchPoint otherPoint1 = otherConnection.newPoint;
        StsPatchPoint thisPoint1 = this.newPoint;
        // sort first by row & col of newPoint
        if(thisPoint1.row < otherPoint1.row) return -1;
        if(thisPoint1.row > otherPoint1.row) return 1;
        if(thisPoint1.col < otherPoint1.col) return -1;
        if(thisPoint1.col > otherPoint1.col) return 1;
        // sort by z of newPoint
        double thisZ1 = thisPoint1.slice;
        double otherZ1 = otherPoint1.slice;
        if(thisZ1 < otherZ1) return -1;
        if(thisZ1 > otherZ1) return 1;

        // if here, the first points of both connections are identical
        // sort by col connection to previous row and then by row connection to previous column
        boolean thisIsRow = isRow();
        boolean otherIsRow = otherConnection.isRow();
        if(!thisIsRow && otherIsRow) return -1;
        else if(thisIsRow && !otherIsRow) return 1;
        // sort by z of otherPoint
        thisZ1 = this.otherPoint.slice;
        otherZ1 = otherConnection.otherPoint.slice;
        if(thisZ1 < otherZ1) return -1;
        if(thisZ1 > otherZ1) return 1;
        StsException.systemError(this, "compareTo", "Two connections are equal. " + toString() + " and " + otherConnection.toString());
        return 0;
    }

    public void clearPointPatchIDs()
    {
        otherPoint.patchID = -1;
        newPoint.patchID = -1;
    }

    public boolean sameRowCol(StsPatchConnection otherConnection)
    {
        if(otherConnection == null) return false;
        if(newPoint.row != otherConnection.newPoint.row) return false;
        if(newPoint.col != otherConnection.newPoint.col) return false;
        if(otherPoint.row != otherConnection.otherPoint.row) return false;
        if(otherPoint.col != otherConnection.otherPoint.col) return false;
        return true;
    }
}
