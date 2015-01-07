package com.Sts.Types;

import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.util.*;

/** Primarily contains attribute data and methods for all seismic lines.
 * All prestack and postack lines, 2d and 3d, defined in both processing and in loading
 * are subclassed from this abstract class through intermediate abstract subclasses.
 * Poststack 2d seismic lines are also subclassed: StsSegyVsp.
 */
abstract public class StsSeismicLine extends StsSeismicBoundingBox
{
    /** min values array for attributes */
    public double[] attributeMinValues;
    /** max values array for attributes */
    public double[] attributeMaxValues;
    /** index of this line in the set of lines */
     public int lineIndex = -1;
   /** zero value for byte data */
    public int signedIntDataZero = 0;

    public String handVelName = null;
    public transient StsFile[] handVelFiles = null;

    public float[] cdpX; /** array of y values for a 2d line in rotated local coordinate system */
    public float[] cdpY; /** array of cdp numbers for a 2d line */
    public int[] cdp;
    transient protected float[][] normals;

    static final float nullValue = StsParameters.nullValue;
    static public final byte nullByte = StsParameters.nullByte;
    static final double roundOff = StsParameters.roundOff;

    static final boolean debug = false;

    static public final String floatFormat = "floats";
    static public final String byteFormat = "bytes"; /** array of x values for a 2d line in rotated local coordinate system */

    public StsSeismicLine()
    {
    }

    public StsSeismicLine(boolean persistent)
    {
        super(persistent);
    }

	public StsSeismicLine(StsFile file, String stsDirectory)
	{
		super(file, stsDirectory);
	}

    public String createHeaderFilename()
    {
        return createFilename(headerFormat);
    }

    public String createFloatDataFilename()
    {
        return createFilename(floatFormat);
    }

    public String createByteDataFilename()
    {
        return createFilename(byteFormat);
    }

    public String getAttributeDataFilename()
    {
        return createFilename(attributeFormat);
    }

    public File getHeaderFile()
    {
        return new File(stsDirectory, createHeaderFilename());
    }

    public String toString()
    {
        return stemname;
    }

    public double[] getAttributeMinValues()
    {
        return attributeMinValues;
    }

    public double[] getAttributeMaxValues()
    {
        return attributeMaxValues;
    }

    public double getAttributeMinValue(int index)
    {
        return attributeMinValues[index];
    }

    public double getAttributeMaxValue(int index)
    {
        return attributeMaxValues[index];
    }

    public double[] getAttributeArray(String attributeName, int nLineTraces) throws FileNotFoundException
    {
		StsMappedDoubleBuffer doubleBuffer = getAttributeArrayBuffer(attributeName, nLineTraces);
        if(doubleBuffer == null)
        {
            StsMessageFiles.infoMessage("Unable to load attribute named: " + attributeName);
            return null;
        }
		double[] attributes = new double[nLineTraces];
		doubleBuffer.get(attributes);
		doubleBuffer.close();
		return attributes;
    }

	public StsMappedDoubleBuffer getAttributeArrayBuffer(String attributeName, int nLineTraces) throws FileNotFoundException
	{
		int attributeIndex = this.getAttributeIndex(attributeName);
		if (attributeIndex == -1)
		{
//			new StsMessage(currentModel.win3d, StsMessage.WARNING, "Failed to find " + attributeName + " attribute in attribute names list.");
			return null;
		}
		return getAttributeArrayBuffer(attributeIndex, nLineTraces);
	}
    
    public StsMappedDoubleBuffer getAttributeArrayBuffer(int attributeIndex, int nLineTraces) throws FileNotFoundException
	{
		StsMappedDoubleBuffer doubleBuffer;
		try
		{
			doubleBuffer = StsMappedDoubleBuffer.constructor(stsDirectory, attributeFilename, "r");
		}
		catch(FileNotFoundException fnfe)
		{
			new StsMessage(currentModel.win3d, StsMessage.WARNING, "Failed to find attribute file " + attributeFilename);
			return null;
		}
		if (attributeIndex == -1)
		{
//			new StsMessage(currentModel.win3d, StsMessage.WARNING, "Failed to find " + attributeName + " attribute in attribute names list.");
			return null;
		}
		doubleBuffer.map(attributeIndex*nLineTraces, nLineTraces);
		return doubleBuffer;
	}

    public StsMappedDoubleBuffer getAllAttributesArrayBuffer(int nLineTraces)
	{
		StsMappedDoubleBuffer doubleBuffer;
		try
		{
			doubleBuffer = StsMappedDoubleBuffer.constructor(stsDirectory, attributeFilename, "r");
		}
		catch(FileNotFoundException fnfe)
		{
			new StsMessage(currentModel.win3d, StsMessage.WARNING, "Failed to find attribute file " + attributeFilename);
			return null;
		}
		doubleBuffer.map(0, nAttributes*nLineTraces);
		return doubleBuffer;
	}

    public String getAssocLineName() { return assocLineName;}
    public void setAssocLineName( String assocLineName) { this.assocLineName = assocLineName;}
    public String getHandVelName() { return handVelName;}
    public void setHandVelName( String handVelName) { this.handVelName = handVelName;}
    public StsFile[] getHandVelFiles() { return handVelFiles;}
    public void setHandVelFiles( ArrayList<StsFile> files) 
    { 
        this.handVelFiles = new StsFile[files.size()];
        for (int i=0; i < files.size(); i++) handVelFiles[i] = files.get(i);
    }

/*
        public void actionPerformed(ActionEvent e)
        {
            if (!(e.getSource() instanceof StsColorscale))
                return;
            StsColorscale colorscale = (StsColorscale) e.getSource();
            if(seismicLineSet.colorscale == colorscale)
            {
                seismicColorList.setColorListChanged(true);
                currentModel.viewObjectRepaint(seismicColorList);
                textureChanged = true;
            }
        }
    */
    public boolean setLineXYsFromCDPs()
    {
        try
        {
            double[] lineCdpX = getAttributeArray(StsSEGYFormat.CDP_X, nCols);
            double[] lineCdpY = getAttributeArray(StsSEGYFormat.CDP_Y, nCols);
            double[] lineCdp = getAttributeArray(StsSEGYFormat.CDP, nCols);
            cdpX = new float[nCols];
            cdpY = new float[nCols];
            cdp = new int[nCols];
            for(int n = 0; n < nCols; n++)
            {
                float[] xy = currentModel.getProject().getRotatedRelativeXYFromUnrotatedAbsoluteXY(lineCdpX[n], lineCdpY[n]);
                cdpX[n] = xy[0];
                cdpY[n] = xy[1];
                cdp[n] = (int)lineCdp[n];
            }
            computeNormals();
            return true;
        }
        catch(FileNotFoundException fnfe)
        {
            new StsMessage(currentModel.win3d, StsMessage.WARNING, "Failed to find attribute file " + attributeFilename);
            return false;
        }
    }

    public boolean setLineXYsFromShtRecs()
    {
        try
        {
            double[] lineShtX = getAttributeArray(StsSEGYFormat.SHT_X, nCols);
            double[] lineShtY = getAttributeArray(StsSEGYFormat.SHT_Y, nCols);
            double[] lineRcvX = getAttributeArray(StsSEGYFormat.REC_X, nCols);
            double[] lineRcvY = getAttributeArray(StsSEGYFormat.REC_Y, nCols);
            double[] lineCdp = getAttributeArray(StsSEGYFormat.CDP, nCols);
            cdpX = new float[nCols];
            cdpY = new float[nCols];
            cdp = new int[nCols];
            int nLastTrace = -1;
            for(int n = 0; n < nCols; n++)
            {
                double x = lineRcvX[n];
                double y = lineRcvY[n];
                float[] xy = currentModel.getProject().getRotatedRelativeXYFromUnrotatedAbsoluteXY(x, y);
                cdpX[n] = xy[0];
                cdpY[n] = xy[1];
                cdp[n] = (int)lineCdp[n];
            }
            computeNormals();
            return true;
        }
        catch(FileNotFoundException fnfe)
        {
            new StsMessage(currentModel.win3d, StsMessage.WARNING, "Failed to find attribute file " + attributeFilename);
            return false;
        }
    }

    protected void computeNormals()
    {
        if(cdpX == null) return;
        normals = new float[nCols][];
        if(nCols < 2)
        {
            normals[0] = new float[]{1.0f, 0.0f, 0.0f};
        }
        else if(nCols == 2)
        {
            normals[0] = new float[]{cdpY[0] - cdpY[1], cdpX[1] - cdpX[0], 0.0f};
            normals[1] = normals[0];
        }
        else
        {
            float x0, y0, x1, y1, x2, y2;
            x1 = cdpX[0];
            y1 = cdpY[0];
            x2 = cdpX[1];
            y2 = cdpY[1];
            for(int n = 1; n < nCols - 1; n++)
            {
                x0 = x1;
                y0 = y1;
                x1 = x2;
                y1 = y2;
                x2 = cdpX[n + 1];
                y2 = cdpY[n + 1];
                normals[n] = new float[]{y0 - y1, x1 - x0, 0.0f};
            }
            normals[0] = normals[1];
            normals[nCols - 1] = normals[nCols - 2];
        }
    }}
