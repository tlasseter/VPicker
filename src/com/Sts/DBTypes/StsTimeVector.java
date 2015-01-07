
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Actions.Import.*;
import com.Sts.IO.*;
import com.Sts.MVC.Main;
import com.Sts.Utilities.*;

import java.io.*;
import java.net.*;


public class StsTimeVector extends StsSerialize implements Cloneable, Serializable
{
    // constants - type

    static public final String[] types = new String[] { "X", "Y", "DEPTH", "MDEPTH", "SUBSEA", "TIME", "LOG", "DRIFT", "AZIMUTH", "TWT", "OWT"};
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

	static public final long largeLong = StsParameters.largeLong;

    // instance fields

    protected String asciiFilename;
    protected String binaryFilename;
    protected int asciiFileTimeColumn = -1;
    protected int asciiFileDateColumn = -1;    
    protected String name;
    protected byte type = LOG;
    protected byte units = StsParameters.DIST_NONE;
    protected int version = 0;
    protected int timeShift = 0;

    protected long minCutoff;
    protected long maxCutoff;

	protected long nullValue = StsParameters.nullLongValue;
	protected boolean isNull = true;

	protected double origin = 0.0;
	protected boolean offsetFromOrigin = false;

	transient double originAdjust = 0.0;
    transient protected StsLongVector values = null;
    transient String stsDirectory = ".";

	public StsTimeVector()
   {
   }

   public StsTimeVector(String name, StsLongVector values)
    {
        this.name  = name;
        setValues(values);
//        if (values != null) setCutoffs(values.getMinValue(), values.getMaxValue());
    }

    public StsTimeVector(String binaryFilename)
    {
        this.binaryFilename = binaryFilename;
        StsKeywordIO.parseBinaryFilename(binaryFilename);
        name = StsKeywordIO.subname;
        setType();
    }

    public StsTimeVector(String asciiFilename, String binaryFilename, String name)
    {
		this(asciiFilename, binaryFilename, name, 0, -1);
    }
    public StsTimeVector(String asciiFilename, String binaryFilename, String name, int version, int asciiFileTimeColumn)
    {
        this(asciiFilename, binaryFilename, name, version, asciiFileTimeColumn, -1);
    }
    
    public StsTimeVector(String asciiFilename, String binaryFilename, String name, int version, int asciiFileTimeColumn,
    		int asciiFileDateColumn)
    {
        this.asciiFilename = asciiFilename;
        binaryFilename = StsStringUtils.cleanString(binaryFilename);
        this.binaryFilename = binaryFilename;
        this.version = version;
        this.name = name;
        this.asciiFileTimeColumn = asciiFileTimeColumn;
        this.asciiFileDateColumn = asciiFileDateColumn;        
        setType();
    }

    public StsTimeVector(byte type, StsLongVector values)
    {
        this.type = type;
        setNameFromType(type);
        this.values = values;
    }

    public StsTimeVector(byte type, long[] values)
    {
        StsLongVector floatVector = new StsLongVector(values);
        this.type = type;
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
	/** values are checked against this null; if null or out of range, value is skipped */
	public void setNullValue(long nullValue) { this.nullValue = nullValue; }

    /** Set the native units of the ASCII log vector */
	public void setUnits(byte units) { this.units = units; }

    /** Set the native units of the ASCII log vector */
	public void setTimeShift(int shift)
    {
        if(shift == timeShift) return;
        timeShift = shift;
        clearArray();
        checkLoadVector();
    }
    public int getTimeShift() { return timeShift; }
    
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

            if (name.equals(clippedType))
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
    public StsLongVector getValues() { return values; }
	public int getAsciiFileTimeColumn() { return asciiFileTimeColumn; }
	public int getAsciiFileDateColumn() { return asciiFileDateColumn; }	
//    public float getNullValue() { return nullValue; }
    public long getMinCutoff(){ return minCutoff; }
    public long getMaxCutoff(){ return maxCutoff; }
    public long getMinValue(){ return values==null ? nullValue : values.getMinValue(); }
    public long getMaxValue(){ return values==null ? nullValue : values.getMaxValue(); }
    public void setMinValue(long minValue){ values.setMinValue(minValue); }
    public void setMaxValue(long maxValue){ values.setMaxValue(maxValue); }
    public double getOrigin() { return origin; }
    public void setVersion(int version) { this.version = version; }
    public int getVersion() { return version; }

	public void setOrigin(double origin)
	{
		offsetFromOrigin = true;
		this.origin = origin;
	}

    public boolean isType(byte type) { return this.type == type; }

    public String getTypeString() { return types[type]; }

    static public StsTimeVector getVectorOfType(StsTimeVector[] logVectors, byte type)
    {
        if(logVectors == null) return null;
        int nLogVectors = logVectors.length;
        for(int n = 0; n < nLogVectors; n++)
            if(logVectors[n] != null && logVectors[n].isType(type)) return logVectors[n];
        return null;
    }

    public void resetVector()
    {
        values.resetVector();
    }

    static public byte getTypeFromString(String typeString)
    {
        int nTypes = types.length;

        for(int n = 0; n < nTypes; n++)
            if(typeString.equals(types[n])) return (byte)n;

        return UNDEFINED;
    }

    static public long[] getValuesArrayOfType(StsTimeVector[] logVectors, byte type)
    {
        StsTimeVector logVector = getVectorOfType(logVectors, type);
        if(logVector == null) return null;
        return logVector.getLongs();
    }

//    public void setLogCurveType(StsLogCurveType curveType) { this.logCurveType = curveType; }

    public long[] getValuesArray()
    {
        if(values == null) return null;
        return values.getValues();
    }

	public long[] getLongs()
	{
        if(!checkLoadVector()) return null;
		return getValuesArray();
	}
	
	public float[] getRelativeFloats()
	{
		if(values == null) return null;
		float[] floatValues = values.getValuesAsRelativeFloats(values.getMinValue());
		if(floatValues != null) return floatValues;
		if(!checkLoadVector()) return null;
		return values.getValuesAsRelativeFloats(values.getMinValue());
	}
	
    public void checkAdjustOrigin(double modelOrigin)
    {
        long adjustOrigin = (long)(origin - modelOrigin);
        if(!StsMath.sameAs(adjustOrigin, 0.0f)) offsetValues(adjustOrigin);
        origin = modelOrigin;
        values.adjustRange(adjustOrigin);
    }

    public void offsetValues(long offset)
    {
        if(values == null) return;
        long[] array = values.getValues();
        int nValues = array.length;
        for(int n = 0; n < nValues; n++)
            array[n] += offset;
    }

    public void setValues(StsLongVector values) { this.values = values; }

    /** validate and set min/max cutoffs */
    public void setCutoffs(long minCutoff, long maxCutoff)
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

    /** null value */
/*
    public void setNullValue(float nullValue)
    {
        this.nullValue = nullValue;
        if (values != null) values.setNullValue(nullValue);
    }
*/
    public void clearArray() { values = null; }

	public boolean hasBinaryFile(String binaryFileDir)
	{
		if(binaryFilename == null) return false;
		String fullFilename = binaryFileDir + File.separator + binaryFilename;
		File file = new File(fullFilename);
		return file.exists();
	}
	
	public boolean checkWriteBinaryFile(String binaryFileDir)
	{
		return checkWriteBinaryFile(binaryFileDir, true);
	}
	
	public boolean checkWriteBinaryFile(String binaryFileDir, boolean append)
	{
		StsBinaryFile binaryFile = null;
		String fullFilename;

		try
		{
//			if(values == null || values.getValues() == null) return false;
			if(binaryFilename == null) return false;

            StsFile file = StsFile.constructor(binaryFileDir, binaryFilename);

			binaryFile = new StsBinaryFile(file);
			if(!binaryFile.openWrite(append)) return false;
            binaryFile.setByteValues(new byte[] {units} );
            binaryFile.setDoubleValues( new double[] { origin } );
			binaryFile.setLongVector(values);
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
                   values.setElement(i, values.getElement(i) * (long)scalar);
       }
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
            origin = origins[0] + timeShift;
			values = new StsLongVector();
			binaryFile.getLongVector(values, loadValues);

            if(timeShift != 0)
            {
                values.setMaxValue(values.getMaxValue() + timeShift);
                values.setMinValue(values.getMinValue() + timeShift);
                long[] lvalues = values.getValues();
                for(int i=0; i<lvalues.length; i++)
                    lvalues[i] = lvalues[i] + timeShift;
                values.setValues(lvalues);
            }
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

   public boolean checkVector()
   {
       StsLongVector longVector = getValues();
	   if(longVector != null && longVector.getValues() != null)
           return true;
       else
           return false;
   }

   public boolean checkLoadVector()
   {
	   StsLongVector longVector = getValues();
	   if(longVector != null && longVector.getValues() != null)
       {
           return true;
       }
	   String binaryDataDir = currentModel.getProject().getBinaryFullDirString();
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
                file = Main.jar.getFileEndingWith(binaryFilename);
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

	public void appendRT(long[] times)
	{
        for(int i=0; i<times.length; i++)
        {
		    if (times[i] == nullValue)
			    values.append(StsParameters.nullLongValue);
		    else
		    {
			    if(offsetFromOrigin)
				    values.append((long)(times[i] - origin));    // Origin already has time shift applied
                else
                    values.append(times[i] + timeShift);
		    }
        }
	}

	public void appendRT(long value)
	{
		if ((long)value == nullValue)
			values.append(StsParameters.nullLongValue);
		else
		{
			if(offsetFromOrigin)
				values.append((long)(value - origin));   // Origin already has time shift applied
            else
                values.append((long)value - timeShift);
		}
	}

	/** Check if value is null because 1) it is out of reasonable range, or 2) it is the defined nullValue.
	 *  On the first nonNull, reset the origin (if origin is being used: offsetFromOrigin = true).
	 *  All values are offsets from origin and first value is 0, i.e., absolute coordinate is origin value.
	 */ 
	public void append(long value)
	{
		if (value <= -largeLong || value >= largeLong)
        {
            values.append(StsParameters.nullLongValue);
            return;
        }
		float floatValue = (float) value;
		if (floatValue == nullValue)
			values.append(StsParameters.nullLongValue);
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
				values.append((long)(value - originAdjust));
            else
                values.append((long)value);
		}
	}

	public void writeAppend(long[] vals, String binaryDir)
	{
		values.setGrowIncrement(1);
		appendRT(vals);
        values.setMinMax();
        // Not the best approach but since RT is being re-written it will do for now.
        // Must back out the timeshift before writting the values and then re-set it
        int timeShiftTemp = timeShift;
        setTimeShift(0);
		checkWriteBinaryFile(binaryDir, false);
        setTimeShift(timeShiftTemp);
	}

	public void writeAppend(long value, String binaryDir)
	{
		values.setGrowIncrement(1);		
		appendRT(value);
        values.setMinMax();

        // Not the best approach but since RT is being re-written it will do for now.
        // Must back out the timeshift before writting the values and then re-set it
        int timeShiftTemp = timeShift;
        setTimeShift(0);
		checkWriteBinaryFile(binaryDir, false);
        setTimeShift(timeShiftTemp);         
	}
	
	public void resetMinMaxValue(long time)
	{
    	if(time > getMaxValue())
    		setMaxValue(time);
    	if(time < getMinValue())
    		setMinValue(time);		
	}

	public void resetMinMaxValue(long[] times)
	{
        for(int i=0; i<times.length; i++)
        {
    	    if(times[i] > getMaxValue())
    		    setMaxValue(times[i]);
    	    if(times[i] < getMinValue())
    		    setMinValue(times[i]);
        }
	}

	public int checkMonotonic()
    {
		return values.checkMonotonic();
	}

	public void setMinMaxAndNulls(long logCurveNull)
	{
		if(values == null) return;
		nullValue = logCurveNull;
		isNull = values.setMinMaxAndNulls(logCurveNull);
	}

    /** Assume values are monotonic: increasing if values[1] >= values[0]; decreasing otherwise */
    public int getIndexNearest(long value)
    {
        long prevValue, nextValue;

        if(!checkLoadVector()) return 0;
        long[] values = getValuesArray();
        int nValues = values.length;
        if(nValues < 2) return 0;

        if(values[1] >= values[0])
        {
            nextValue = values[0];
            for (int n = 1; n < nValues; n++)
            {
                prevValue = nextValue;
                nextValue = values[n];
                if (value <= nextValue)
                    return n - 1;
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
                  return n - 1;

            }
        }
        return nValues;
    }

    /** Assume values are monotonic: increasing if values[1] >= values[0]; decreasing otherwise */
    public float getIndexF(long value)
    {
        float prevValue, nextValue;

        if(!checkLoadVector()) return nullValue;
        long[] values = getValuesArray();
        int nValues = values.length;
        if(nValues < 2) return nullValue;

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
        return nullValue;
    }

    public long getValue(float indexF)
    {
        if(!checkLoadVector()) return nullValue;
        long[] values = getValuesArray();
        int nValues = values.length;
        if(nValues < 2) return nullValue;
        int index = (int)indexF;
        index = StsMath.minMax(index, 0, nValues-2);
        float f = indexF - index;
        long prevValue = values[index];
        long nextValue = values[index+1];
        return prevValue*(1 - (long)f) + nextValue*(long)f;
    }

    public int getIndex(long time)
    {
        int idx = -1;
        if(!checkLoadVector()) return idx;
        long[] values = getValuesArray();
        int nValues = values.length;
        for(int i=0; i<nValues; i++)
        {
            if(values[i] == time)
                return i;
        }
        return -1;
    }

    public byte getUnits() { return units; }
}

































