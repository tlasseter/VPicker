package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Utilities.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Mar 8, 2008
 * Time: 3:29:08 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsLineSegment implements StsInstance3dDisplayable
{
    private StsPoint firstPoint;
    private StsPoint lastPoint;
    private StsColor stsColor = StsColor.DARK_YELLOW;

    public StsLineSegment()
    {
    }

    public StsLineSegment(StsColor color)
    {
        stsColor = color;
    }

    public void clearPoints()
    {
        firstPoint = null;
        lastPoint = null;
    }

    public void display(StsGLPanel3d glPanel3d)
    {
        if(firstPoint == null || lastPoint == null) return;
        StsGLDraw.drawLineSegment(glPanel3d.getGL(), stsColor, firstPoint, lastPoint, 4);
    }

    public StsPoint getFirstPoint()
    {
        return firstPoint;
    }

    public void setFirstPoint(double[] xyz)
    {
        this.firstPoint = new StsPoint(xyz);
    }

    public void setFirstPoint(float[] xyz)
    {
        this.firstPoint = new StsPoint(xyz);
    }

    public void setFirstPoint(StsPoint firstPoint)
    {
        this.firstPoint = firstPoint;
        lastPoint = null;
    }

    public StsPoint getLastPoint()
    {
        return lastPoint;
    }

    public void setLastPoint(StsPoint lastPoint)
    {
        this.lastPoint = lastPoint;
    }

    public void setLastPoint(double[] xyz)
    {
        this.lastPoint = new StsPoint(xyz);
    }

    public void setLastPoint(float[] xyz)
    {
        this.lastPoint = new StsPoint(xyz);
    }

    public float getDistance()
    {
        if(firstPoint == null || lastPoint == null) return 0.0f;
        return firstPoint.distance(lastPoint);
    }
    
    public double getAngle()
    {
        if(firstPoint == null || lastPoint == null) return 0.0f;
        return lastPoint.fromAngle(firstPoint);
    }

    public float getDip()
    {
        if(firstPoint == null || lastPoint == null) return 0.0f;
        return firstPoint.fromDip(lastPoint);
    }

    public float getHorizontalDistance()
    {
        if(firstPoint == null || lastPoint == null) return 0.0f;
        return firstPoint.horizontalDistance(lastPoint);
    }
    
    public float getVerticalDistance()
    {
        if(firstPoint == null || lastPoint == null) return 0.0f;
        return firstPoint.verticalDistance(lastPoint);
    }
    
    public StsColor getStsColor()
    {
        return stsColor;
    }

    public void setStsColor(StsColor color)
    {
        stsColor = color;
    }
}
