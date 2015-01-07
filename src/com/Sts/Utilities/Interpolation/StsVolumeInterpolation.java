//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Utilities.Interpolation;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Utilities.*;

import java.nio.*;
import java.util.*;

public class StsVolumeInterpolation extends StsInterpolation
{
     /** volume to be interpolated */
    StsSeismicVolume volume;
    /** horizontal slice of volume to be interpolated */
    int nSlice;
    /** byteBuffer for volume plane being interpolated */
    ByteBuffer planeBuffer;
   /** byteBuffer for volume slice thru point or points currently being interpolated. */
    byte[] sliceBytes;
    int nRows, nCols, nSlices;
    float determinant;

    static final byte nullByte = StsParameters.nullByte;
    static final boolean debug = false;

    static final int XDIR = StsCursor3d.XDIR;
    static final int YDIR = StsCursor3d.YDIR;
    static final int ZDIR = StsCursor3d.ZDIR;

    public StsVolumeInterpolation()
    {
    }

    public StsVolumeInterpolation(StsSeismicVolume volume)
    {
        initializeGrid(volume);
    }

    static public StsVolumeInterpolation getInstance(StsSeismicVolume volume)
    {
        if (interpolator == null) interpolator = new StsVolumeInterpolation(volume);
        interpolator.useGradient = false;
        return (StsVolumeInterpolation)interpolator;
    }

    public void initializeGrid(StsSeismicVolume volume)
	{
		if(this.volume == volume) return;
        this.volume = volume;
        nRows = volume.nRows;
        nCols = volume.nCols;
        nSlices = volume.nSlices;
        rowMin = 0;
		rowMax = nRows-1;
		colMin = 0;
		colMax = nCols-1;
    }

    /** included for abstract class compatability. */
    boolean isPointGrid(int row, int col) { return false; }

    public ByteBuffer interpolatePlane(int dir, int nPlane, ByteBuffer plane)
    {
        int capacity = plane.capacity();
        byte[] filledBytes = new byte[capacity];
        plane.get(filledBytes);

        switch(dir)
        {
            case XDIR:
                interpolateColPlane(nPlane, nRows, nSlices, filledBytes);
                break;
            case YDIR:
                interpolateRowPlane(nPlane, nCols, nSlices, filledBytes);
                break;
            case ZDIR:
                interpolateSlicePlane(nPlane, nRows, nCols, filledBytes);

        }
        ByteBuffer byteBuffer = ByteBuffer.allocateDirect(capacity);
        byteBuffer.put(filledBytes);
        byteBuffer.rewind();
        return byteBuffer;
    }

    private void interpolateRowPlane(int row, int nCols, int nSlices, byte[] filledBytes)
    {
        int i;
        for(int n = 0; n < nSlices; n++)
        {
            setSlice(n);
            i = n;
            for(int col = 0; col < nCols; col++, i += nSlices)
            {
                byte byteValue = filledBytes[i];
                if(byteValue != nullByte) continue;
                int value = (int)interpolate(row, col, false, false);
                byteValue = StsMath.unsignedIntToUnsignedByte(value);
                filledBytes[i] = byteValue;
            }
        }
    }

    private void interpolateColPlane(int col, int nRows, int nSlices, byte[] filledBytes)
    {
        int i;
        for(int n = 0; n < nSlices; n++)
        {
            setSlice(n);
            i = n;
            for(int row = 0; row < nRows; row++, i += nSlices)
            {
                byte byteValue = filledBytes[i];
                if(byteValue != nullByte) continue;
                int value = (int)interpolate(row, col, false, false);
                byteValue = StsMath.unsignedIntToUnsignedByte(value);
                filledBytes[i] = byteValue;
            }
        }
    }

    private void interpolateSlicePlane(int nPlane, int nRows, int nCols, byte[] filledBytes)
    {
        setSlice(nPlane, filledBytes);
        int i = 0;
        for(int row = 0; row < nRows; row++)
        {
            for(int col = 0; col < nCols; col++, i++)
            {
                byte byteValue = filledBytes[i];
                if(byteValue != nullByte) continue;
                int value = (int)interpolate(row, col, false, false);
                byteValue = StsMath.unsignedIntToUnsignedByte(value);
                filledBytes[i] = byteValue;
            }
        }
    }

    public void setSlice(int slice, byte[] filledBytes)
    {
        this.nSlice = slice;
        int nBytes = filledBytes.length;
        this.sliceBytes = new byte[nBytes];
        System.arraycopy(filledBytes, 0, sliceBytes, 0, nBytes);
    }

    public void setSlice(int slice)
    {
        this.nSlice = slice;
        sliceBytes = volume.readBytePlaneData(StsCursor3d.ZDIR, nSlice);
    }

    public float interpolate(int iCenter, int jCenter, boolean useGradient, boolean debugInterpolate)
    {
        interpolator.initialize(iCenter, jCenter, useGradient, debugInterpolate);
        return interpolator.interpolatePoint();
    }

    public float interpolatePoint(int i, int j)
    {
        iCenter = i;
        jCenter = j;
        float z = nullValue;
        return interpolatePoint();
    }

    private int getByteBufferIntValue(int i, int j)
    {
        byte byteValue =  getByteBufferValue(i, j);
        return StsMath.signedByteToUnsignedInt(byteValue);
    }

    private byte getByteBufferValue(int i, int j)
    {
        int n = i*nCols + j;
        return  sliceBytes[n];
    }

    /** return true if this point at ij in the source surface can be used for interpolation */
    protected boolean addPointOK(int i, int j)
    {
        return getByteBufferValue(i, j) != nullByte;
    }


    void addPoint(int i, int j)
    {
        int dx, dy, distSq;
        float wt;

        byte byteValue = getByteBufferValue(i, j);
        if(byteValue == nullByte) return;
        int z = StsMath.signedByteToUnsignedInt(byteValue);
        dx = j - jCenter;
        dy = i - iCenter;
        distSq = dx*dx + dy*dy;
        wt = 1.0f/distSq ;

        weightSum += wt;
        weightMinimum = wt < weightMinFactor*weightSum;

        // when first point found, use it to determine minNPoints required for termination
        if(nPoints++ == 0) setMinNPoints(distSq);

        Point point = new Point(i, j, z, dx, dy, wt);
        sortedPoints.add(point);
        spiralHasGoodPoints = true;
    }

    float getZ(ArrayList sortedPoints)
    {
        Point point;
        int row, col;
        float dZdX, dZdY;
        float z0, w, z, wz, wt, dx, dy;

        try
        {
            int nTotalPoints = sortedPoints.size();
            if (nTotalPoints == 0) return nullValue;

            nPointsUsed = Math.min(nTotalPoints, currentMinNPoints);
            point = (Point) sortedPoints.get(nPointsUsed - 1);
            float weightMin = point.wt;
            for (int n = nPointsUsed + 1; n < nTotalPoints; n++)
            {
                point = (Point) sortedPoints.get(n);
                if (point.wt < weightMin) break;
                nPointsUsed++;
            }

            // remove the points we are not going to use
            for (int n = nTotalPoints - 1; n >= nPointsUsed; n--)
                sortedPoints.remove(n);

            // useGradient is ignored

            w = 0.0f;
            wz = 0.0f;

            for (int n = 0; n < nPointsUsed; n++)
            {
                point = (Point) sortedPoints.get(n);
                row = point.row;
                col = point.col;
                dx = -point.dx;
                dy = -point.dy;
                z = point.z;

                wt = point.wt;
                w += wt;
                wz += wt * z;

                if (debugInterpolate)
                {
                    StsMessageFiles.infoMessage(" point " + n + " row: " + row + " col: " + col + " z: " + z + " wt: " + wt +
                        " z: " + z + " dx: " + dx + " dy: " + dy);
                }
            }
            z = wz / w;
            if (debugInterpolate) StsMessageFiles.infoMessage("    z: " + z);
            return z;
        }
        catch (Exception e)
        {
            StsException.outputException("StsBlockGridInterpolationWeightedPlane.getZ() failed.",
                    e, StsException.WARNING);
            return nullValue;
        }
    }

    public String getGridName() { return volume.getName(); }
/*
    protected boolean isConverged()
    {
        if(spiralHasGoodPoints) nGoodSpirals++;
        if(nGoodSpirals < minNGoodSpirals) return false;
        return weightMinimum || nPoints >= maxNPoints;
    }
*/
}