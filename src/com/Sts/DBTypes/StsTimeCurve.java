
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;


public class StsTimeCurve extends StsMainObject implements StsSelectable
{
	protected String wellname;
	protected String logname;
	//protected String name;
    protected int version = 0;
    protected float curveMin, curveMax;
    protected StsTimeCurveType timeCurveType = null;  // Not used yet but required for backwards compatibility with dbs from prior releases

    /** directory where the data is stored */
	public String stsDirectory;

    protected StsTimeVector timeVector = null;
    protected StsLogVector mDepthVector = null;
    protected StsLogVector valueVector = null;

    transient boolean debug = false;
    static public final double doubleNullValue = StsParameters.doubleNullValue;

	public StsTimeCurve()
	{
	}

    public StsTimeCurve(boolean persistent)
    {
        super(persistent);
    }

    public StsTimeCurve(StsTimeVector timeVector, StsLogVector valueVector, int version)
    {
        this(timeVector, null, valueVector, version);
    }
	public StsTimeCurve(StsTimeVector timeVector, StsLogVector mDepthVector, StsLogVector valueVector, int version)
	{
		super(false);
		if (timeVector==null && valueVector==null) return;

		this.timeVector = timeVector;
		this.valueVector = valueVector;
        this.mDepthVector = mDepthVector;
		this.version = version;

        this.stsDirectory = currentModel.getProject().getDataFullDirString();


		setName(valueVector.getName());
		//if(timeCurveType == null) timeCurveType = getTimeCurveTypeFromName(getName());
		addToModel();
		//timeCurveType.addLogCurve(this);
	}

	static public StsTimeCurve constructTimeCurve(StsTimeVector timeVector, StsLogVector vector, int version)
	{
		try
		{
			if (vector == null) return null;
			if (timeVector == null) return null;

			vector.setMinMaxAndNulls(StsParameters.largeFloat);
			if(vector.isNull)
				return null;
			StsTimeCurve logCurve = new StsTimeCurve(timeVector, vector, version);

			return logCurve;
		}
		catch (Exception e)
		{
			StsException.outputException("StsTimeCurve:constructTimeCurve() failed.", e, StsException.WARNING);
			return null;
		}
	}

	static public StsTimeCurve[] constructTimeCurves(StsTimeVector timeVector, StsLogVector[] timeVectors, int version)
	{
		try
		{
			if (timeVectors == null)return null;
			StsTimeCurve[] timeCurves = new StsTimeCurve[0];

			if (timeVector == null)
                return null;

			int nTimeCurveVectors = timeVectors.length;
			for (int n = 0; n < nTimeCurveVectors; n++)
			{
				timeVectors[n].setMinMaxAndNulls(StsParameters.largeFloat);
				if(timeVectors[n].isNull)
					continue;
				StsTimeCurve logCurve = new StsTimeCurve(timeVector, timeVectors[n], version);
				timeCurves = (StsTimeCurve[])StsMath.arrayAddElement(timeCurves, logCurve);
			}
			return timeCurves;
		}
		catch (Exception e)
		{
			StsException.outputException("StsTimeCurve:constructTimeCurves() failed.", e, StsException.WARNING);
			return null;
		}
	}

    public boolean initialize(StsModel model)
    {
        return initialize();
    }

    public boolean initialize()
    {
        return true;
    }

	public void setWell(StsWell well)
	{
		wellname = well.getName();
		this.dbFieldChanged("wellname", wellname);
	}

    public String getWellname() { return wellname; }

    //public void setTimeCurveType(StsTimeCurveType curveType) { this.timeCurveType = curveType; }
    //public StsTimeCurveType getTimeCurveType() { return timeCurveType; }

    //public StsColor getStsColor() { return timeCurveType.getStsColor(); }

    public StsTimeVector getTimeVector()
    {
        checkLoadVectors();
        return timeVector;
    }
    public StsLogVector getValueVector()
    {
        checkLoadVectors();
        return valueVector;
    }
    public StsLogVector getMDepthVector() { return mDepthVector; }

    public StsLongVector getTimeLongVector()
    {
        if (timeVector == null) return null;
		timeVector.checkLoadVector(stsDirectory);
        return timeVector.getValues();
    }

    public StsFloatVector getValuesFloatVector()
    {
        if (valueVector == null) return null;
		valueVector.checkLoadVector(stsDirectory);
        return valueVector.getValues();
    }

    public float[] getVectorHistogram() { return valueVector.getHistogram(); }
    public StsFloatVector getMDepthFloatVector()
    {
        if (mDepthVector == null) return null;
        mDepthVector.checkLoadVector(stsDirectory);
        return mDepthVector.getValues();
    }

    public float[] getMDepthVectorFloats()
    {
        return getVectorFloats(mDepthVector);
	}

	public long[] getTimeVectorLongs()
	{
		return getVectorLongs(timeVector);
	}

	public long[] getVectorLongs(StsTimeVector timeVector)
	{
		if(timeVector == null) return null;
		return timeVector.getLongs();
	}

    public float[] getVectorFloats(StsLogVector logVector)
    {
        if(logVector == null) return null;
        return logVector.getFloats();
	}

	public float[] getValuesVectorFloats()
	{
		StsFloatVector floatVector = getValuesFloatVector();
		if(floatVector == null) return null;
		return floatVector.getValues();
	}

	public boolean checkLoadVectors()
	{
		if(timeVector == null || !timeVector.checkLoadVector(stsDirectory)) return false;
		if(valueVector == null || !valueVector.checkLoadVector(stsDirectory)) return false;
		return true;
	}

    public void convertTimeUnits(float scalar)
    {
        if(timeVector != null)
            timeVector.convertToProjectUnits(scalar);
        return;
    }

    public boolean deletePoint(int idx)
    {
        // Not implemented yet since it will require a re-write of all the binary files.
        return true;
    }
    public void setValueVector(StsLogVector valueVector){ this.valueVector = valueVector; }

    //public void setName(String name) { this.name = name; }
    //public String getName() { return name; }
    public String getLogname() { return logname; }
    public String getCurvename() { return getName(); }
    public String toString() { return getName(); }
    public int getNumValues()
    {
        float[] valVector = this.getValuesVectorFloats();
        if(valVector == null) return 0;
        else return valVector.length;
    }

    /** clear floatArrays to reduce memory requirements */
    public void clearFloatArrays()
    {
        if(valueVector != null) valueVector.clearArray();
        if(timeVector != null) timeVector.clearArray();
        if(mDepthVector != null) mDepthVector.clearArray();
    }

    public void setCurveMin(float curveMin)
	{
		fieldChanged("curveMin", curveMin);
	}

    public void setCurveMax(float curveMax)
	{
		fieldChanged("curveMax", curveMax);
	}

    public float getCurveMin()
    {
        if(valueVector == null)
        {
        	return curveMin;
        	//return StsParameters.largeFloat;
        }
        StsFloatVector floatVector = valueVector.getValues();
        if(floatVector == null)
        {
        	return curveMin;
        	//return StsParameters.largeFloat;
        }
        return floatVector.getMinValue();
    }

    public float getCurveMax()
    {
        if(valueVector == null)
        {
        	return curveMax;
        	//return StsParameters.largeFloat;
        }
        StsFloatVector floatVector = valueVector.getValues();
        if(floatVector == null)
        {
        	return curveMax;
        	//return -StsParameters.largeFloat;
        }
        return floatVector.getMaxValue();
    }

    public void convertDepthUnits(float scalar)
    {
        if(mDepthVector != null)
            mDepthVector.convertToProjectUnits(scalar);
        return;
    }

    public StsTimeCurveType getTimeCurveTypeFromName(String name)
    {
		StsTimeCurveType timeCurveType;
		timeCurveType = (StsTimeCurveType)getStsClassObjectWithName(StsTimeCurveType.class, name);
		if(timeCurveType == null)
            timeCurveType = new StsTimeCurveType(this, null);
		return timeCurveType;
    }

    public boolean retrieveValuesFromDataStore()
    {
        return true;
    }

    public float getMinValue()
    {
        if(valueVector == null) return StsParameters.largeFloat;
        return valueVector.getMinValue();
    }

    public float getMaxValue()
    {
        if(valueVector == null) return -StsParameters.largeFloat;
        return valueVector.getMaxValue();
    }

    public float getValueAt(int index, long time)
    {
        int idx = getTimeIndexGE(time);
        return valueVector.getValue(idx);
    }

    public float getValueAt(long time)
    {
        int idx = getTimeIndexGE(time);
        return valueVector.getValue(idx);
    }

    public int getTimeIndexGE(long time)
    {
        long timeF = (long)time;

        long[] timeArray = getTimeVector().getValues().getValues();
        return StsMath.arrayIndexAbove(timeF, timeArray);
    }

    public String getStsDirectory() { return stsDirectory; }

    public int getTimeIndexLE(long time)
    {
        long timeF = (long)time;
        long[] timeArray = getTimeVector().getValues().getValues();
        return StsMath.arrayIndexBelow(timeF, timeArray);
    }

    public boolean addValueToCurve(Long time, Double value)
    {
    	return addValueToCurve(value.doubleValue());
    }

    public boolean addValueToCurve(double value)
    {
    	String binaryDataDir = currentModel.getProject().getBinaryFullDirString();
    	getValueVector().writeAppend(value, binaryDataDir);
    	return true;
    }
    public boolean replaceValueInCurve(long time, double value)
    {
    	String binaryDataDir = currentModel.getProject().getBinaryFullDirString();
        int index = getTimeVector().getIndex(time);
    	getValueVector().writeReplaceAt(value, index, binaryDataDir);
    	return true;
    }
    public boolean addTimeToCurve(long time)
    {
    	String binaryDataDir = currentModel.getProject().getBinaryFullDirString();
    	getTimeVector().writeAppend(time, binaryDataDir);
    	return true;
    }
    public boolean addValuesToCurve(double[] values)
    {
    	String binaryDataDir = currentModel.getProject().getBinaryFullDirString();
    	getValueVector().writeAppend(values, binaryDataDir);
    	return true;
    }
    public boolean addTimesToCurve(long[] times)
    {
    	String binaryDataDir = currentModel.getProject().getBinaryFullDirString();
    	getTimeVector().writeAppend(times, binaryDataDir);
    	return true;
    }
    public void resetVectors()
    {
        getValueVector().resetVector();
        getTimeVector().resetVector();
    }
}















