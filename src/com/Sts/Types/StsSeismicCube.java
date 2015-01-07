
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Stuart A. Jackson
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;

public class StsSeismicCube extends StsSerialize implements Cloneable, StsSerializable
{
    /** bit type combinations used in calling distanceSquaredType method */
    static public final int INLINE =   1;
    static public final int XLINE =   2;
    static public final int TD =   3;

    protected int nXLines;                            // number of crosslines
    protected int nLines;                             // number of lines
    protected int nSamples;                           // number of samples per trace
    protected double xOrigin;                         // x at first line and crossline
    protected double yOrigin;                         // y at first line and crossline
    protected float lineSpacing;                      // spacing between lines
    protected float xLineSpacing;                     // spacing between crosslines
    protected float sampleSpacing;                    // spacing between samples (msecs or secs)
    protected float zMin;                             // minimum time (or depth if depth cube)
    protected float zMax;                             // maximum time (or depth if depth cube)
    protected float lineAngle;                        // angle from +x axis to a line in the cube
    protected boolean isXLineCCW;                     // true if crosslineAngle is 90 degrees CCW from line angle
    protected float lineNumMin;                       // line number for first line
    protected float lineNumMax;                       // line number for last line
    protected float xLineNumMin;                      // crossline number for first crossline
    protected float xLineNumMax;                      // crossline number for last crossline
    protected float lineNumInc;                       // interval between lineNumbers;
    protected float xLineNumInc;                      // interval between xlineNumbers;
    protected float dataMin = Float.MAX_VALUE;        // minimum data value
    protected float dataMax = -Float.MAX_VALUE;       // maximum data value

    public StsSeismicCube()
    {

    }

    public int getNumberXlines() { return nXLines; }
    public int getNumberInlines() { return nLines; }
    public int getNumberSamples() { return nSamples; }
    public int[] getVolumeStats()
    {
        int[] vol = new int[3];
        vol[INLINE] = nLines;
        vol[XLINE] = nXLines;
        vol[TD] = nSamples;
        return vol;
    }
    public double getXOrigin() { return xOrigin; }
    public double getYOrigin() { return yOrigin; }
    public double[] getOrigin()
    {
        double[] origin = new double[2];
        origin[0] = xOrigin;
        origin[1] = yOrigin;
        return origin;
    }
    public float getInlineSpacing() { return lineSpacing; }
    public float getXlineSpacing() { return xLineSpacing; }
    public float getSampleSpacing() { return sampleSpacing; }
    public float getZMin() { return zMin; }
    public float getZMax() { return zMax; }
    public float[] getZRange()
   {
       float[] z = new float[2];
       z[0] = zMin;
       z[1] = zMax;
       return z;
   }
    public float getInlineAngle() { return lineAngle; }
    public float getInlineMin() { return lineNumMin; }
    public float getInlineMax() { return lineNumMax; }
    public float getInlineIncrement() { return lineNumInc; }
    public float[] getInlineStats()
   {
       float[] i = new float[3];
       i[0] = lineNumMin;
       i[1] = lineNumMax;
       i[2] = lineNumInc;
       return i;
   }
    public float getXlineMin() { return xLineNumMin; }
    public float getXlineMax() { return xLineNumMax; }
    public float getXlineIncrement() { return xLineNumInc; }
    public float[] getXlineStats()
    {
        float[] x = new float[3];
        x[0] = xLineNumMin;
        x[1] = xLineNumMax;
        x[2] = xLineNumInc;
        return x;
    }
    public float getDataMin() { return dataMin; }
    public float getDataMax() { return dataMax; }
    public float[] getDataRange()
    {
        float[] d = new float[2];
        d[0] = dataMin;
        d[1] = dataMax;
        return d;
    }

    // Set Values
    public void setNumberXlines(int nx) { nXLines = nx; }
    public void setNumberInlines(int ni) { nLines = ni; }
    public void setNumberSamples(int ns) { nSamples = ns; }
    public boolean setVolumeStats(int[] vstats)
    {
        if(vstats.length != 3) return false;
        nLines = vstats[INLINE];
        nXLines = vstats[XLINE];
        nSamples = vstats[TD];
        return true;
    }
    public void setInlineAngle(float iangle) { lineAngle = iangle; }
    public void setInlineMin(float imin) { lineNumMin = imin; }
    public void setInlineMax(float imax) { lineNumMax = imax; }
    public void setInlineIncrement(float iinc) { lineNumInc = iinc; }
    public boolean setInlineStats(float[] istats)
    {
        if(istats.length != 3) return false;
        lineNumMin = istats[0];
        lineNumMax = istats[1];
        lineNumInc = istats[2];
        return true;
    }
    public void setXlineMin(float xmin) { xLineNumMin = xmin; }
    public void setXlineMax(float xmax) { xLineNumMax = xmax; }
    public void setXlineIncrement(float xinc) { xLineNumInc = xinc; }
    public boolean setXlineStats(float[] xstats)
    {
        if(xstats.length != 3) return false;
        lineNumMin = xstats[0];
        lineNumMax = xstats[1];
        lineNumInc = xstats[2];
        return true;
    }
    public void setDataMin(float dmin) { dataMin = dmin; }
    public void setDataMax(float dmax) { dataMax = dmax; }
    public boolean setDataRange(float[] dstats)
    {
        if(dstats.length != 2) return false;
        dataMin = dstats[0];
        dataMax = dstats[1];
        return true;
    }

}




