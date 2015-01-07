package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Aug 31, 2007
 * Time: 10:08:18 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsCursorZSurface extends StsRotatedGridBoundingBox implements  StsXYSurfaceGridable
{
    float zt;
    byte zDomain;

    static float[] verticalNormal = new float[] { 0.0f, 0.0f, 1.0f };

    public StsCursorZSurface()
    {
    }

    public StsCursorZSurface(boolean persistent)
    {
        super(persistent);
    }

    public StsCursorZSurface(float z)
    {
 //       super(false);
        this.zt = z;
        zDomain = currentModel.getProject().getZDomain();
    }

    // StsXYSurfaceGridable methods

    public float getZMin() { return StsParameters.nullValue; }
    public float getZMax() { return StsParameters.nullValue; }
    public float getZ() { return zt; }
    public float getZInc() { return 0.0f; }
    public float[][] getPointsZ() { return null; }

	public float interpolateBilinearZ(StsGridPoint gridPoint, boolean computeIfNull, boolean setPoint) 
    {
        float z =  interpolateBilinearZ(gridPoint.getX(), gridPoint.getY());
        if(setPoint) gridPoint.setZ(z);
        return z;
    }

    public float interpolateBilinearZ(StsPoint point, boolean computeIfNull, boolean setPoint)
    {
        return interpolateBilinearZ(point.getX(), point.getY());
    }

    public StsPoint getXYZTPoint(StsPoint cursorPoint)
    {
        float x = cursorPoint.getX();
        float y = cursorPoint.getY();
        if(zDomain == StsProject.TD_DEPTH)
        {
            float t;

            StsSeismicVelocityModel velocityModel = currentModel.getProject().velocityModel;
            if(velocityModel == null)
                t = zt;
            else
            {
                try
                {
                    t = (float)velocityModel.getT(x, y, zt, zt);
                }
                catch(Exception e)
                {
                    t = zt;
                }
            }
            return new StsPoint(x, y, zt, 0.0, t);
        }
        else // zDomain == StsProject.TD_TIME
        {
            float t = zt;
            float zz;

            StsSeismicVelocityModel velocityModel = currentModel.getProject().velocityModel;
            if(velocityModel == null)
                zz = t;
            else
            {
                try
                {
                    zz = (float)velocityModel.getZ(x, y, t);
                }
                catch(Exception e)
                {
                    zz = t;
                }
            }
            return new StsPoint(x, y, zz, 0.0, t);
        }
    }

    public float interpolateBilinearZ(float x, float y)
    {
        if(zDomain == StsProject.TD_DEPTH)
        {
            if(isDepth) return zt;
            StsSeismicVelocityModel velocityModel = currentModel.getProject().velocityModel;
            if(velocityModel == null) return zt;
            try
            {
                return (float)velocityModel.getT(x, y, zt, zt);
            }
            catch(Exception e)
            {
                return zt;
            }
        }
        else // TD_TIME
        {
            if(!isDepth) return zt;
            StsSeismicVelocityModel velocityModel = currentModel.getProject().velocityModel;
            if(velocityModel == null) return zt;
            try
            {
                return (float)velocityModel.getZ(x, y, zt);
            }
            catch(Exception e)
            {
                return zt;
            }
        }
    }

    public float getComputePointZ(int row, int col) { return zt; }

    public boolean toggleSurfacePickingOn() { return true; }
    public void toggleSurfacePickingOff() {}
    public String getName() { return Float.toString(zt); }
    public StsGridPoint getSurfacePosition(StsMouse mouse, boolean display, StsGLPanel3d glPanel3d) { return null; }
    public void setIsVisible(boolean isVisible) { }
    public boolean getIsVisible() { return true; }

    // StsSurfaceGridable methods

    public float getRowCoor(float[] xyz) { return 0; }
    public float getColCoor(float[] xyz) { return 0; }

    public StsPoint getPoint(int row, int col) { return null; }
    public StsPoint getPoint(float rowF, float colF) { return null; }
    public float[] getXYZorT(int row, int col) { return null; }
    public float[] getXYZorT(float rowF, float colF) { return null; }
    public float[] getNormal(int row, int col) { return verticalNormal; }
    public float[] getNormal(float rowF, float colF) { return verticalNormal; }
    public void checkConstructGridNormals() { }


}
