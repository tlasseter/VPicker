
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Actions.Import.*;
import com.Sts.DB.*;
import com.Sts.IO.*;
import com.Sts.MVC.Main;
import com.Sts.Types.StsPoint;
import com.Sts.Utilities.*;

import java.io.*;
import java.net.*;
import java.util.ArrayList;

public class StsLogVector extends StsSerialize implements Cloneable, StsSerializable, Serializable
{
    // constants - type

    static public final String[] types = new String[] { "X", "Y", "DEPTH", "MDEPTH", "SUBSEA", "TIME", "LOG", "DRIFT", "AZIMUTH", "TWT", "OWT"};
    static public final String[] depth_types = new String[] {"DEPTH", "DEPT", "TVD" };
    static public final byte X = 0;
    static public final byte Y = 1;
    static public final byte DEPTH = 2;
    static public final byte MDEPTH = 3;
    static public final byte SUBSEA = 4;
    static public final byte TIME = 5;
    static public final byte LOG = 6;
    static public final byte DRIFT = 7;
    static public final byte AZIMUTH = 8;
    static public final byte TWT = 9;
    static public final byte OWT = 10;
    static public final byte UNDEFINED = -1;

    static public final String WELL_DEV_PREFIX = "well-dev";
    static public final String WELL_LOG_PREFIX = "well-logs";
	static public final String WELL_TD_PREFIX = "well-td";
    static public final String SENSOR_PREFIX = "sensor";
    static public final String WELL_TIME_LOG_PREFIX = "well-timeLogs";

    static public final String FORMAT_BIN = "bin";
    static public final String FORMAT_TXT = "txt";

    static public final float largeFloat = Float.MAX_VALUE;

    // instance fields

    protected String asciiFilename;
    protected String binaryFilename;
    protected int asciiFileColumn;
    protected String name;
    protected byte type = LOG;
    protected byte units = StsParameters.DIST_NONE;
    protected int version = 0;

    protected float minCutoff;
    protected float maxCutoff;

	public float nullValue = StsParameters.nullValue;
	protected boolean isNull = true;

	protected double origin = 0.0;
	protected boolean offsetFromOrigin = false;

	transient double originAdjust = 0.0;
    transient protected StsFloatVector values = null;
    transient String stsDirectory = ".";
    /** Histogram of the data distribution */
    public float[] dataHist = null;

    /** Total samples in each of 255 steps */
    transient private int dataCnt[] = new int[255];
    transient private int ttlHistogramSamples = 0;

	public StsLogVector()
   {
   }

   public StsLogVector(String name, StsFloatVector values)
    {
        this.name  = name;
        setValues(values);
//        if (values != null) setCutoffs(values.getMinValue(), values.getMaxValue());
    }

    public StsLogVector(String binaryFilename)
    {
        this.binaryFilename = binaryFilename;
        StsKeywordIO.parseBinaryFilename(binaryFilename);
        name = StsKeywordIO.subname;
        setType();
    }
/*
    public StsLogVector(String asciiFilename, String name, int asciiFileColumn)
    {
        this.asciiFilename = asciiFilename;
        setName(name);
        this.asciiFileColumn = asciiFileColumn;
        setType();
    }
*/
    public StsLogVector(String asciiFilename, String binaryFilename, String name)
    {
		this(asciiFilename, binaryFilename, name, 0, -1);
    }

    public StsLogVector(String asciiFilename, String binaryFilename, String name, int version, int asciiFileColumn)
    {
        this.asciiFilename = asciiFilename;
        binaryFilename = StsStringUtils.cleanString(binaryFilename);
        this.binaryFilename = binaryFilename;
        this.version = version;
        this.name = name;
        this.asciiFileColumn = asciiFileColumn;
        setType();
    }

    public StsLogVector(byte type, StsFloatVector values)
    {
        this.type = type;
        setNameFromType(type);
        this.values = values;
    }

    public StsLogVector(byte type, ArrayList<StsPoint> points, int index)
    {
        this.type = type;
        setNameFromType(type);
        int nValues = points.size();
        float[] values = new float[nValues];
        for(int n = 0; n < nValues; n++)
            values[n] = points.get(n).v[index];
        this.values = new StsFloatVector(values);
    }

    public StsLogVector(String name, ArrayList<StsPoint> points, int index)
    {
        this.name = name;
        int nValues = points.size();
        float[] values = new float[nValues];
        for(int n = 0; n < nValues; n++)
            values[n] = points.get(n).v[index];
        this.values = new StsFloatVector(values);
    }

    public StsLogVector(byte type, float[] values)
    {
        StsFloatVector floatVector = new StsFloatVector(values);
        this.type = type;
        setNameFromType(type);
        this.values = floatVector;
    }

    public StsLogVector(float[] values)
    {
        StsFloatVector floatVector = new StsFloatVector(values);
        setNameFromType(type);
        this.values = floatVector;
    }

    /** Set Log Name */
    public void setNameAndType(String name, String binaryFilename)
    {
        this.name = name;
        this.binaryFilename = binaryFilename;
        setType();
    }

    public void setName(String name) { this.name = name; }
    public void setBinaryFileName(String name) { this.binaryFilename = name; }

	/** values are checked against this null; if null or out of range, value is skipped */
	public void setNullValue(float nullValue) { this.nullValue = nullValue; }

    /** Set the native units of the ASCII log vector */
	public void setUnits(byte units) { this.units = units; }

	/** isNull indicates vector is all nulls */
	public boolean isNull() { return isNull; }

    private void setType()
    {
        String clippedType = null;
        int nTypes = types.length;

        for(int n = 0; n < nTypes; n++)
        {
            if(name.length() < types[n].length())
                clippedType = types[n].substring(0,name.length());
            else
                clippedType = types[n];

            if (name.equalsIgnoreCase(clippedType))
                type = (byte) n;
        }
    }

    private void setNameFromType(byte type) { name = getStringFromType(type); }

    static public String getStringFromType(byte type)
    {
        int nTypes = types.length;
        if(type < 0 || type >= nTypes) return null;
        return types[type];
    }

    // Accessors
    public String getName(){ return name; }
    public byte getType() { return type; }
    public StsFloatVector getValues() { return values; }
	public int getAsciiFileColumn() { return asciiFileColumn; }
//    public float getNullValue() { return nullValue; }
    public float getMinCutoff(){ return minCutoff; }
    public float getMaxCutoff(){ return maxCutoff; }
    public float getMinValue(){ return values==null ? StsParameters.nullValue : values.getMinValue(); }
    public float getMaxValue(){ return values==null ? StsParameters.nullValue : values.getMaxValue(); }
    public void setMinValue(float minValue){ values.setMinValue(minValue); }
    public void setMaxValue(float maxValue){ values.setMaxValue(maxValue); }
    public double getOrigin() { return origin; }
    public void setVersion(int version) { this.version = version; }
    public int getVersion() { return version; }

	public void resetMinMaxValue(float val)
	{
    	if(val > getMaxValue())
    		setMaxValue(val);
    	if(val < getMinValue())
    		setMinValue(val);
	}

	public void resetMinMaxValue(double[] vals)
	{
        for(int i=0; i<vals.length; i++)
        {
    	    if(vals[i] > getMaxValue())
    		    setMaxValue((float)vals[i]);
    	    if(vals[i] < getMinValue())
    		    setMinValue((float)vals[i]);
        }
	}

	public void setOrigin(double origin)
	{
        if(origin == 0.0) return;
		offsetFromOrigin = true;
		this.origin = origin;
	}

    public boolean isType(byte type) { return this.type == type; }

    public String getTypeString() { return types[type]; }

    static public StsLogVector getVectorOfType(StsLogVector[] logVectors, byte type)
    {
        if(logVectors == null) return null;
        int nLogVectors = logVectors.length;
        for(int n = 0; n < nLogVectors; n++)
        {
        	if(logVectors[n] != null && logVectors[n].isType(type))
        		return logVectors[n];

        }
        return null;
    }

    static public byte getTypeFromString(String typeString)
    {
        int nTypes = types.length;

        for(int n = 0; n < nTypes; n++)
            if(typeString.equals(types[n])) return (byte)n;

        return UNDEFINED;
    }

    static public float[] getValuesArrayOfType(StsLogVector[] logVectors, byte type)
    {
        StsLogVector logVector = getVectorOfType(logVectors, type);
        if(logVector == null) return null;
        return logVector.getFloats();
    }

//    public void setLogCurveType(StsLogCurveType curveType) { this.logCurveType = curveType; }

    public float[] getValuesArray()
    {
        if(values == null) return null;
        return values.getValues();
    }

	public float[] getFloats()
	{
        if(!checkLoadVector()) return null;
		return getValuesArray();
	}

    public void checkAdjustOrigin(double modelOrigin)
    {
        float adjustOrigin = (float)(origin - modelOrigin);
        if(!StsMath.sameAs(adjustOrigin, 0.0f)) offsetValues(adjustOrigin);
        origin = modelOrigin;
        values.adjustRange(adjustOrigin);
    }

    public void offsetValues(float offset)
    {
        if(values == null) return;
        float[] array = values.getValues();
        int nValues = array.length;
        for(int n = 0; n < nValues; n++)
            array[n] += offset;
    }

    public void setValues(float[] vals)
    {
    	StsFloatVector vector = new StsFloatVector(vals);
        this.values = vector;
        if(values == null) return;
        //values.setMinMax();
        for(int i=0; i<values.getValues().length; i++)
        {
            float value = values.getElement(i);
            if (value == nullValue)
                values.setElement(i, StsParameters.nullValue);
        }
        values.setMinMax();
    }

    public void setValues(StsFloatVector vals)
    {
        this.values = vals;
        if(values == null) return;
        values.setMinMax();
        for(int i=0; i<values.getValues().length; i++)
        {
            float value = values.getElement(i);
            if (value == nullValue)
                values.setElement(i, StsParameters.nullValue);
        }
    }

    /** validate and set min/max cutoffs */
    public void setCutoffs(float minCutoff, float maxCutoff)
    {
	    if (minCutoff > maxCutoff)
    	{
            StsException.systemError("StsLogVector.setCutoffs failed." +
                		 " min cutoff must be less than max cutoff");
            return;
 	    }
    	this.minCutoff = minCutoff;
    	this.maxCutoff = maxCutoff;
    }

    static public boolean deviationVectorsOK(StsLogVector[] vectors)
    {
        StsLogVector xVector, yVector, zVector, mVector;
        xVector = getVectorOfType(vectors, X);
        yVector = getVectorOfType(vectors, Y);
        zVector = getVectorOfType(vectors, DEPTH);
        mVector = getVectorOfType(vectors, MDEPTH);
        return xVector != null && yVector != null && (zVector != null || mVector != null);
    }

    /** null value */
/*
    public void setNullValue(float nullValue)
    {
        this.nullValue = nullValue;
        if (values != null) values.setNullValue(nullValue);
    }
*/
    public void clearArray() { values = null; }

    public void resetVector()
    {
        values.resetVector();
    }

	public boolean hasBinaryFile(String binaryFileDir)
	{
		if(binaryFilename == null) return false;
		String fullFilename = binaryFileDir + File.separator + StsStringUtils.cleanString(binaryFilename);
		File file = new File(fullFilename);
		return file.exists();
	}

	public boolean deleteWriteBinaryFile()
	{
        String binaryDataDirectory = currentModel.getProject().getBinaryFullDirString();
        return checkWriteBinaryFile(binaryDataDirectory, true, true);
	}

    public boolean checkWriteBinaryFile(String binaryDataDirectory)
	{
		return checkWriteBinaryFile(binaryDataDirectory, true);
	}

    public boolean checkWriteBinaryFile(String binaryDataDirectory, boolean append)
	{
        return checkWriteBinaryFile(binaryDataDirectory, append, false);
    }

    public boolean checkWriteBinaryFile(String binaryDataDirectory, boolean append, boolean delete)
	{
		StsBinaryFile binaryFile = null;
		String fullFilename;

		try
		{
//			if(values == null || values.getValues() == null) return false;
			if(binaryFilename == null) return false;

            StsFile file = StsFile.constructor(binaryDataDirectory, binaryFilename);
            if(delete && file.exists()) file.delete();
//			fullFilename = binaryFileDir + File.separator + binaryFilename;
//            URL url = new URL("file:" + fullFilename);

			binaryFile = new StsBinaryFile(file);
//			binaryFile = new StsBinaryFile(url);
			if(!binaryFile.openWrite(append)) return false;
            binaryFile.setByteValues(new byte[] {units} );
            binaryFile.setDoubleValues( new double[] { origin } );
			binaryFile.setFloatVector(values);
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsLogVector.writeBinaryFile() failed.",
					e, StsException.WARNING);
			return false;
		}
		finally
		{
			try { if(binaryFile != null) binaryFile.close(); }
			catch(Exception e) { }
		}
	}

    public void convertToProjectUnits(float scalar)
    {
       // Convert values to project units
       if(values != null)
       {
           origin = origin * scalar;
           if(scalar != 1.0)
               for(int i=0; i<values.getSize(); i++)
                   values.setElement(i, values.getElement(i) * scalar);
       }
    }

    public int getSize()
    {
        checkLoadVector();
        if(values == null) return 0;
        return values.getSize();
    }

    public boolean readBinaryFile(StsAbstractFile file, boolean loadValues)
	{
		StsBinaryFile binaryFile = null;

		try
		{
			binaryFile = new StsBinaryFile(file);
			if(!binaryFile.openRead()) return false;
            units = binaryFile.getByteValues()[0];
            double[] origins = binaryFile.getDoubleValues();
            origin = origins[0];
			values = new StsFloatVector();
			binaryFile.getFloatVector(values, loadValues);
			values.setGrowIncrement(1);
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsLogVector.readBinaryFile() failed.",
					e, StsException.WARNING);
			return false;
		}
		finally
		{
			try { if(binaryFile != null) binaryFile.close(); }
			catch(Exception e) { }
		}
	}
/*
	public boolean readBinaryFile(URL url, boolean loadValues)
	{
		StsBinaryFile binaryFile = null;
		String fullFilename;

		try
		{
			binaryFile = new StsBinaryFile(url);
			if(!binaryFile.openReadAndCheck()) return false;
            double[] origins = binaryFile.getDoubleValues();
            origin = origins[0];
			values = new StsFloatVector();
			binaryFile.getFloatVector(values, loadValues);
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsLogVector.readBinaryFile() failed.",
					e, StsException.WARNING);
			return false;
		}
		finally
		{
			try { if(binaryFile != null) binaryFile.close(); }
			catch(Exception e) { }
		}
	}
*/
   public boolean checkVector()
   {
       StsFloatVector floatVector = getValues();
	   if(floatVector != null && floatVector.getValues() != null)
           return true;
       else
           return false;
   }

   public boolean checkLoadVector()
   {
	   StsFloatVector floatVector = getValues();
	   if(floatVector != null && floatVector.getValues().length != 0)
       {
           return true;
       }
//	   String binaryDataDir = stsDirectory + currentModel.getProject().getBinaryDirString();
	   String binaryDataDir = currentModel.getProject().getBinaryFullDirString();
//       System.out.println("Vector Directory=" + binaryDataDir);
       if(Main.isJarDB)
           return readBinaryFile(Main.jarDBFilename, true);
       else
          return readBinaryFile(binaryDataDir, true);
   }

   public boolean checkLoadVector(String directory)
   {
       stsDirectory = directory;
       return checkLoadVector();
   }

	public boolean readBinaryFile(String binaryFileDir, boolean loadValues)
	{
        URL url;
		StsBinaryFile binaryFile = null;
		String fullFilename;
        StsAbstractFile file = null;

		try
		{
            if(Main.isJarDB)
            {
                String[] filenames = Main.jar.getFilenames();
                for(int i=0; i<filenames.length; i++)
                {
                    if(filenames[i].endsWith(binaryFilename))
                    {
                        file = Main.jar.getFile(i);
                        break;
                    }
                }
                if(file == null) return false;
            }
            else
            {
                file = StsFile.constructor(binaryFileDir, binaryFilename);
                if(file == null) return false;
                if(!file.exists()) return false;
            }
            return readBinaryFile(file, loadValues);
        }
		catch(Exception e)
		{
			StsException.outputException("StsLogVector.readBinaryFile() failed.",
					e, StsException.WARNING);
			return false;
		}
		finally
		{
			try { if(binaryFile != null) binaryFile.close(); }
			catch(Exception e) { }
		}
	}
/*
	public boolean readBinaryFile(DataInputStream dis, boolean loadValues)
	{
		StsBinaryFile binaryFile = null;
		String fullFilename;

		try
		{
			binaryFile = new StsBinaryFile(dis);
            double[] origins = binaryFile.getDoubleValues();
            origin = origins[0];
			values = new StsFloatVector();
			binaryFile.getFloatVector(values, loadValues);
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsLogVector.readBinaryFile() failed.",
					e, StsException.WARNING);
			return false;
		}
	}
*/
    private InputStream getInputStream(URL url)
    {
        try
        {
            URLConnection conn = url.openConnection();
            return conn.getInputStream();
        }
        catch(Exception e)
        {
            StsException.outputException("StsSelectCurveTableModel.getInputStream() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

	public void appendRT(double[] vals)
	{
        for(int i=0; i<vals.length; i++)
        {
		    if ((float)vals[i] == nullValue || vals[i] <= -Float.MAX_VALUE || vals[i] >= Float.MAX_VALUE)
			    values.append(StsParameters.nullValue);
		    else
		    {
			    if(offsetFromOrigin)
				    values.append((float)(vals[i] - origin));
                else
                    values.append((float)vals[i]);
		    }
        }
	}

	public void appendRT(double value)
	{
		if ((float)value == nullValue || value <= -Float.MAX_VALUE || value >= Float.MAX_VALUE)
			values.append(StsParameters.nullValue);
		else
		{
			if(offsetFromOrigin)
				values.append((float)(value - origin));
            else
                values.append((float)value);
		}
	}
	/** Check if value is null because 1) it is out of reasonable range, or 2) it is the defined nullValue.
	 *  On the first nonNull, reset the origin (if origin is being used: offsetFromOrigin = true).
	 *  All values are offsets from origin and first value is 0, i.e., absolute coordinate is origin value.
	 */
	public void append(double value)
	{
		if ((float)value == nullValue || value <= -Float.MAX_VALUE || value >= Float.MAX_VALUE)
			values.append(StsParameters.nullValue);
		else
		{
			// on first non-null, reset origin to fileHeaderOrigin + thisValue
			// adjust subsequent values by the first value
			if(isNull)
            {
                if(offsetFromOrigin)
                {
                    originAdjust = value;
                    origin += value;
                }
                isNull = false;
            }

			if(offsetFromOrigin)
				values.append((float)(value - originAdjust));
            else
                values.append((float)value);
		}
	}

	public void writeAppend(double value, String binaryDir)
	{
		values.setGrowIncrement(1);
		appendRT(value);
        values.setMinMax();
		checkWriteBinaryFile(binaryDir, false);
	}

	public void writeReplaceAt(double value, int index, String binaryDir)
	{
        values.setElement(index, (float)value);
        values.setMinMax();
		checkWriteBinaryFile(binaryDir, false);
	}

	public void writeAppend(double[] vals, String binaryDir)
	{
		values.setGrowIncrement(1);
		appendRT(vals);
        values.setMinMax();
		checkWriteBinaryFile(binaryDir, false);
	}
	public void append(double[] vals)
	{
		values.setGrowIncrement(1);
		appendRT(vals);
        values.setMinMax();
    }
	public void checkMonotonic()
	{
		if(values == null) return;
		values.checkMonotonic();
	}

	public void setMinMaxAndNulls(float logCurveNull)
	{
		if(values == null) return;
		nullValue = logCurveNull;
		isNull = values.setMinMaxAndNulls(logCurveNull);
	}

    /** Assume values are monotonic: increasing if values[1] >= values[0]; decreasing otherwise */
    public float getIndexF(float value)
    {
        float prevValue, nextValue;

        if(!checkLoadVector()) return StsParameters.nullValue;
        float[] values = getValuesArray();
        int nValues = values.length;
        if(nValues < 2) return StsParameters.nullValue;

        if(values[1] >= values[0])
        {
            nextValue = values[0];
            for (int n = 1; n < nValues; n++)
            {
                prevValue = nextValue;
                nextValue = values[n];
                if (value <= nextValue)
                    return n - 1 + (value - prevValue) / (nextValue - prevValue);
            }
        }
        else
        {
            prevValue = values[nValues-1];
            for(int n = nValues-2; n == 0; n--)
            {
                nextValue = prevValue;
                prevValue = values[n];
                if (value >= nextValue)
                  return n - 1 + (value - prevValue) / (nextValue - prevValue);

            }
        }
        return StsParameters.nullValue;
    }

    public float getValue(float indexF)
    {
        if(!checkLoadVector()) return StsParameters.nullValue;
        float[] values = getValuesArray();
        int nValues = values.length;
        if(nValues < 2) return StsParameters.nullValue;
        int index = (int)indexF;
        index = StsMath.minMax(index, 0, nValues-2);
        float f = indexF - index;
        float prevValue = values[index];
        float nextValue = values[index+1];
        return prevValue*(1.0f - f) + nextValue*f;
    }

    public byte getUnits() { return units; }

    public float[] getHistogram()
    {
        if(dataHist != null)
            return dataHist;
        dataHist = new float[255];
        clearHistogram();
        for(int i=0; i<getValuesArray().length; i++)
        {
            accumulateHistogram(getValuesArray()[i]);
        }
        calculateHistogram();
        return dataHist;
    }

    public void accumulateHistogram(int bindex)
    {
        if (bindex > 254)
        {
            bindex = 254;
        }
        if (bindex < 0)
        {
            bindex = 0;
        }
        dataCnt[bindex] = dataCnt[bindex] + 1;
        ttlHistogramSamples++;
    }

    private void accumulateHistogram(float value)
    {
        float scaledFloat = 254.0f * (float)(((double)value - getMinValue()) / (getMaxValue() - getMinValue()));
        int scaledInt = StsMath.minMax(Math.round(scaledFloat), 0, 254);
        accumulateHistogram(scaledInt);
    }

    public void calculateHistogram()
    {
        for (int i = 0; i < 255; i++)
        {
            dataHist[i] = (float) ( (float)dataCnt[i] / (float)ttlHistogramSamples) * 100.0f;
        }
    }

    public void clearHistogram()
    {
        for (int i = 0; i < 255; i++)
        {
            dataCnt[i] = 0;
            dataHist[i] = 0.0f;
        }
        ttlHistogramSamples = 0;
    }

    static public float interpolateValue(float inputValue, StsLogVector inputVector, StsLogVector outputVector)
    {
        if(!inputVector.checkLoadVector()) return StsParameters.nullValue;
        if(!outputVector.checkLoadVector()) return StsParameters.nullValue;
        float[] inputValues = inputVector.getValuesArray();
        float[] outputValues = outputVector.getValuesArray();
        return StsMath.interpolateValue(inputValue, inputValues, outputValues);
    }
}

































