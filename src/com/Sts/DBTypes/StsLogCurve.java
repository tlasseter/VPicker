
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.io.*;
import java.util.ArrayList;

public class StsLogCurve extends StsMainTimeObject implements StsSelectable
{
//    protected StsWell well;
    /** Name of the well this logCurve belongs to */
    protected String wellname;
    /** A specific name for this logCurve, different from its logCurveType name (curvename); not currently used */
    protected String logname;
    /** this is the logCurveType name; it could be the same as or an alias of an existing logCurveType name */
    protected String curvename;
    /** If edited, this log curve will created a new instance with the same name, different values, and an incremented version number */
    protected int version = 0;
    /** The type or family this log curve belongs to (GR, NPHI, etc). */
    protected StsLogCurveType logCurveType = null;
//    protected float curveMin, curveMax;
    /** directory where the data is stored */
	public String stsDirectory;

    transient StsWellClass wellClass;

    protected StsLogVector mdepthVector = null;
    protected StsLogVector depthVector = null;
    protected StsLogVector timeVector = null;
    protected StsLogVector valueVector = null;

    static public final double doubleNullValue = StsParameters.doubleNullValue;

	public StsLogCurve()
	{

	}

    public StsLogCurve(boolean persistent)
    {
        super(persistent);

        stsDirectory = currentModel.getProject().getDataFullDirString();
        wellClass = (StsWellClass)currentModel.getStsClass(StsWell.class);
    }

    public StsLogCurve(String name, boolean persistent)
    {
        super(persistent);
        setName(name);
    }

    public StsLogCurve(StsLogVector mdepthVector, StsLogVector depthVector, StsLogVector valueVector, int version)
	{

		super(false);

		if (depthVector==null && valueVector==null) return;

		this.mdepthVector = mdepthVector;
		this.depthVector = depthVector;
		this.valueVector = valueVector;
		this.version = version;

        stsDirectory = currentModel.getProject().getDataFullDirString();

		curvename = valueVector.getName();
		logCurveType = getLogCurveTypeFromName(curvename);
		addToModel();
		logCurveType.addLogCurve(this);
        wellClass = (StsWellClass)currentModel.getStsClass(StsWell.class);
    }

	public StsLogCurve( String name )
		{

			super(false);

			//if (depthVector==null && valueVector==null) return;

			//this.mdepthVector = mdepthVector;
			//this.depthVector = depthVector;
			//this.valueVector = valueVector;
			//this.version = version;

			stsDirectory = currentModel.getProject().getDataFullDirString();

			curvename = name;
			logCurveType = getLogCurveTypeFromName(curvename);

			addToModel();
			//if (logCurveType != null)
			   logCurveType.addLogCurve(this);
			wellClass = (StsWellClass)currentModel.getStsClass(StsWell.class);
			StsLogCurve[] logCurves = new StsLogCurve[0];
			logCurves = (StsLogCurve[])StsMath.arrayAddElement(logCurves, this);
		}

    static public StsLogCurve nullLogCurveConstructor(String name)
    {
        return new StsLogCurve(name, false);
    }

	static public StsLogCurve[] constructLogCurves(StsWell well, StsLogVector[] logVectors, float curveNullValue, int version)
	{
		StsLogVector xVector = null, yVector = null, depthVector = null, mdepthVector = null;
		try
		{
			if (logVectors == null)return null;
			int nLogVectors = logVectors.length;
			StsLogCurve[] logCurves = new StsLogCurve[0];

			StsLogVector[] logCurveVectors = new StsLogVector[0];
			for (int n = 0; n < nLogVectors; n++)
			{
                if(logVectors[n].values == null) continue;
                if (logVectors[n].isType(StsLogVector.DEPTH))
					depthVector = logVectors[n];
				else if (logVectors[n].isType(StsLogVector.MDEPTH))
					mdepthVector = logVectors[n];
				else if (logVectors[n].isType(StsLogVector.X) || logVectors[n].isType(StsLogVector.Y))
					continue;
				else
					logCurveVectors = (StsLogVector[]) StsMath.arrayAddElement(logCurveVectors, logVectors[n]);
			}

			if (depthVector == null && mdepthVector == null)
            {
                return null;
            }

			if (mdepthVector == null && well != null)
			{
				mdepthVector = well.getMDepthsVectorFromDepths(depthVector.getFloats());
				if (mdepthVector != null) logVectors = (StsLogVector[]) StsMath.arrayAddElement(logVectors, mdepthVector);
			}

			//		if(!setVectors(logVectors)) return false;
			//		StsLogCurve depthLogCurve = new StsLogCurve(well, mdepthVector, depthVector, null, 0);
			//		well.setDepthLogCurve(depthLogCurve);

			if (depthVector != null)
                depthVector.checkMonotonic();

			int nLogCurveVectors = logCurveVectors.length;
			for (int n = 0; n < nLogCurveVectors; n++)
			{
				logCurveVectors[n].setMinMaxAndNulls(curveNullValue);
				StsLogCurve logCurve = new StsLogCurve(mdepthVector, depthVector, logCurveVectors[n], version);
				logCurves = (StsLogCurve[])StsMath.arrayAddElement(logCurves, logCurve);
			}
			return logCurves;
		}
		catch (Exception e)
		{
			StsException.outputException("StsWellKeywordIO.constructLogCurves() failed.",
										 e, StsException.WARNING);
			return null;
		}
	}

    public boolean initialize(StsModel model)
    {
        return initialize();
    }

    public boolean initialize()
    {
        StsObject object = currentModel.getObjectWithName(StsWell.class, wellname);
        if(object != null)
            wellClass = (StsWellClass)currentModel.getStsClass(StsWell.class);
        else
            wellClass = (StsWellClass)currentModel.getStsClass(StsLiveWell.class);

        return true;
    }

	public void setWell(StsWell well)
	{
//		this.well = well;
		wellname = well.getName();
		this.dbFieldChanged("wellname", wellname);
	}

    public String getWellname() { return wellname; }

//    public String getName() { return parentSet==null ? null : parentSet.getName(); }
//    public StsColor getStsColor() { return (parentSet==null) ? null : parentSet.getStsColor(); }
//    public Color getColor() { return (parentSet==null) ? null : parentSet.getColor(); }
/*
	public boolean checkAddMDepth(StsWell well)
	 {
		 if (mdepthVector != null) return true;
		 mdepthVector = well.getMDepthsFromDepthsUsingWellDev(depthVector.getFloats());
		 return mdepthVector != null;
	 }
*/

    public void setCurvename(String curvename)
    {

        this.curvename = curvename;
        logCurveType = getLogCurveTypeFromName(curvename);

    }

//    public void setLogCurveType(StsLogCurveType curveType) { this.logCurveType = curveType; }
    public StsLogCurveType getLogCurveType() { return logCurveType; }

    public StsColor getStsColor() { return logCurveType.getStsColor(); }

    public StsLogVector getMDepthVector() { return mdepthVector; }
    public StsLogVector getDepthVector() { return depthVector; }
    public StsLogVector getTimeVector() { return timeVector; }
    public StsLogVector getValueVector() { return valueVector; }

    public StsFloatVector getMDepthFloatVector()
    {
        if (mdepthVector == null) return null;
		mdepthVector.checkLoadVector(stsDirectory);
        return mdepthVector.getValues();
    }

    public StsFloatVector getDepthFloatVector()
    {
        if (depthVector == null) return null;
        depthVector.checkLoadVector(stsDirectory);
        return depthVector.getValues();
    }

    public StsFloatVector getTimeFloatVector()
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

	public float[] getMDepthVectorFloats()
	{
		StsFloatVector floatVector = getMDepthFloatVector();
		if(floatVector == null) return null;
		return floatVector.getValues();
	}

	public float[] getDepthVectorFloats()
	{
		StsFloatVector floatVector = getDepthFloatVector();
		if(floatVector == null) return null;
		return floatVector.getValues();
	}

	public float[] getTimeVectorFloats()
	{
		return getVectorFloats(timeVector);

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
		if(depthVector == null || !depthVector.checkLoadVector(stsDirectory)) return false;
		if(valueVector == null || !valueVector.checkLoadVector(stsDirectory)) return false;
        if(mdepthVector != null)
            mdepthVector.checkLoadVector(stsDirectory);
		return true;
	}

    public void convertDepthUnits(float scalar)
    {
        if(depthVector != null)
            depthVector.convertToProjectUnits(scalar);
        if(mdepthVector != null)
            mdepthVector.convertToProjectUnits(scalar);
        return;
    }

    public void convertTimeUnits(float scalar)
    {
        if(timeVector != null)
            timeVector.convertToProjectUnits(scalar);
        return;
    }

    public void setValueVector(StsLogVector valueVector){ this.valueVector = valueVector; }

    public String getName() { return curvename; }
    public String getLogname() { return logname; }
    public String getCurvename() { return curvename; }
    public String toString() { return curvename; }

    /** clear floatArrays to reduce memory requirements */
    public void clearFloatArrays()
    {
        if(valueVector != null) valueVector.clearArray();
        if(mdepthVector != null) mdepthVector.clearArray();
        if(depthVector != null) depthVector.clearArray();
        if(timeVector != null) timeVector.clearArray();
    }

    /** set/get logCurveType */
    public StsLogCurveType getLogCurveTypeFromName(String name)
    {
		StsLogCurveType logCurveType;
		logCurveType = (StsLogCurveType)getStsClassObjectWithName(StsLogCurveType.class, name);
		if(logCurveType == null) logCurveType = new StsLogCurveType(this);
		return logCurveType;
    }

    public String getLogCurveTypeName()
    {
        return logCurveType.getName();
    }

    public boolean logCurveTypeNameMatches(String name)
    {
        return logCurveType.getName().equals(name);
    }

    /** get interpolated value for one array given another */
	public float getInterpolatedValue(float depth)
	{
		float[] depths = getDepthVectorFloats();
		if(depths == null) return StsParameters.nullValue;
		float[] values = getValuesVectorFloats();
		if(values == null) return StsParameters.nullValue;
		return getInterpolatedValue(depths, values, depth);
	}

	static public float getInterpolatedValue(float[] depths, float[] values, float depth)
	{
		int indexBelow = StsMath.arrayIndexBelow(depth, depths);
		int nValues = depths.length;
		indexBelow = StsMath.minMax(indexBelow, 0, nValues-2);
		float f = (depth - depths[indexBelow])/(depths[indexBelow+1] - depths[indexBelow]);
		return values[indexBelow] + f*(values[indexBelow+1] - values[indexBelow]);
	}

    public boolean matchesName(String name)
    {
        if(name.equals(curvename)) return true;
        if(!logCurveType.hasAlias()) return false;
        return name.equals(logCurveType.aliasToType.toString());
    }

    public boolean retrieveValuesFromDataStore()
    {
        return true;
    /*
        try
        {
            // get the HDF values
            StsProject p = currentModel.getProject();
            String wellName = well.getName();
            boolean zDomain = HDFWellFile.DEPTH_TAG.equals(p.getZDomain());
            File file = new File(p.getDataFullDirString(), wellName);
            SDSFile hdfFile = new SDSFile(file, SDSFile.ACCESS_READ);
            HDFWellFile wellFile = new HDFWellFile(hdfFile, zDomain);
            HDFWellLog wellLog = zDomain ? wellFile.getZLog(logName)
                    : wellFile.getTLog(logName);
            if (wellLog==null)  // try the S2S created HDF well file
            {
                hdfFile.end();
                file = new File(p.getBinaryFullDirString(), wellName);
                hdfFile = new SDSFile(file, SDSFile.ACCESS_READ);
                wellFile = new HDFWellFile(hdfFile, zDomain);
                wellLog = zDomain ? wellFile.getZLog(logName)
                        : wellFile.getTLog(logName);
            }
            if (wellLog==null) { hdfFile.end(); return false; }
            float[] depthsArray = wellLog.getDepthArray();
            float[] valuesArray = wellLog.getValueArray();
            hdfFile.end();

            // build the log curve

            StsLogVector zVector = new StsLogVector(currentModel, wellName, depthName,
                        new StsFloatVector(depthsArray, depthsArray.length, 1));
            depthDomain = zDomain ? StsLogVector.MEASURED_Z : StsLogVector.TIME;
            zVector.setDepthDomain(depthDomain);
            if (zDomain) mdepthVector = zVector;
            else timeVector = zVector;
            values = new StsLogVector(currentModel, wellName, logName,
                        new StsFloatVector(valuesArray, valuesArray.length, 1));
        }
        catch (Exception e) { return false; }
        return true;
    */
    }


/*
    static public boolean constructLogCurves(StsWell well, StsLogVector[] logVectors, float curveNullValue, String prefix)
    {
        StsLogVector xVector = null, yVector = null, depthVector = null, mdepthVector = null;

        try
        {
            if(logVectors == null) return false;
            int nLogVectors = logVectors.length;

            StsLogVector[] logCurveVectors  = new StsLogVector[0];
            for(int n = 0; n < nLogVectors; n++)
            {
                if(logVectors[n].isType(StsLogVector.MDEPTH))
                    mdepthVector = logVectors[n];
                else if(logVectors[n].isType(StsLogVector.X))
				    xVector = logVectors[n];
			    else if(logVectors[n].isType(StsLogVector.Y))
					yVector = logVectors[n];
				else if(logVectors[n].isType(StsLogVector.DEPTH))
					depthVector = logVectors[n];

              else
                    logCurveVectors = (StsLogVector[])StsMath.arrayAddElement(logCurveVectors, logVectors[n]);
            }

            if(depthVector == null && mdepthVector == null) return false;

			if(mdepthVector == null)
			{
				mdepthVector = well.getMDepthsFromDepths(xVector, yVector, depthVector, prefix);
				if(mdepthVector == null) return false;
				logVectors = (StsLogVector[])StsMath.arrayAddElement(logVectors, mdepthVector);
			}

			StsLogCurve depthLogCurve = new StsLogCurve(well, mdepthVector, depthVector, null, 0);
			well.setDepthLogCurve(depthLogCurve);

			if(depthVector != null) depthVector.checkMonotonic();

            int nLogCurveVectors = logCurveVectors.length;
            for(int n = 0; n < nLogCurveVectors; n++)
            {
				logCurveVectors[n].setMinMaxAndNulls(curveNullValue);
                StsLogCurve logCurve = new StsLogCurve(well, mdepthVector, depthVector, logCurveVectors[n], 0);
                well.addLogCurve(logCurve);
            }
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsWellKeywordIO.constructLogCurves() failed.",
                e, StsException.WARNING);
            return false;
        }
    }
*/

    /** return average of values over a depth range */
    public float getAverageOverZRange(float zTop, float zBase)
    {
        boolean debug = false;
        StsFloatVector values = getValuesFloatVector();
        int[] indexRange = getDepthFloatVector().getIndicesInValueRange(zTop, zBase);
        if (indexRange==null)
        {
            if (debug) System.out.println("No index range returned for "
                    + "curve: " + getName() + " with zTop = " + zTop
                    + " & zBase = " + zBase);
            return StsParameters.nullValue;
        }
        float total = 0.0f;
        int nValues = indexRange[1] - indexRange[0] + 1;
        if (debug)
        {
            System.out.println("\nStsLogCurves.getAverageOverZRange:");
            System.out.println("curve:  " + getName());
            System.out.println("Z range:  " + zTop + " - " + zBase);
            System.out.print("Values: ");
        }
        for (int i=indexRange[0]; i<=indexRange[1]; i++)
        {
            float value = values.getElement(i);
            if (debug) System.out.print(value + " ");
            if (value == StsParameters.nullValue) nValues--;
            else total += value;
        }
        if (debug) System.out.println(" ");
        float average = total / (float)nValues;
        if (debug) System.out.println("Average: " + average);
        return average;
    }

    /** return most common categorical facies over a depth range */
    public float getCategoricalValueOverZRange(StsCategoricalFacies categoricalFacies,
            float zTop, float zBase)
    {
        try
        {
            StsFloatVector valueVector = getValuesFloatVector();
            float[] allValues = valueVector.getValues();
            int[] indexRange = getDepthFloatVector().getIndicesInValueRange(zTop, zBase);

            int nValues = indexRange[1] - indexRange[0] + 1;
            float[] values = new float[nValues];
            System.arraycopy(allValues, indexRange[0], values, 0, nValues);
            return (float)categoricalFacies.getMostCommonFaciesCategory(values);
        }
        catch (NullPointerException e) { return StsParameters.nullValue; }
    }

    /** print out log curve values & their average over a depth range
        Note: this method doesn't do any tvd-measured depth conversion. */
    public void displayRange(float zMin, float zMax)
    {
        try
        {
            StringBuffer buffer = new StringBuffer(" \n");

            String depthName = depthVector.getName();
            String curveName = getName();
            int nValues = getValuesFloatVector().getSize();
            buffer.append("Curve: " + curveName +
                    "\nNumber of values: " + nValues + "\n");
            buffer.append(" \n");
            buffer.append("index\t" + depthName + "\tValue\n");
            StsFloatVector depthValues;
            depthValues = getDepthFloatVector();
            StsFloatVector curveValues = getValuesFloatVector();
            for (int j=0; j<nValues; j++)
            {
                float depth = depthValues.getElement(j);
                if (depth>=zMin && depth<=zMax)
                {
                    buffer.append(j + "\t" + depthValues.getElement(j) +
                        "\t" + curveValues.getElement(j) + "\n");
                    if (depth>zMax) break;
                }
            }
            buffer.append(" \n");

            // print average value
            float average = getAverageOverZRange(zMin, zMax);
            buffer.append(" \nAverage over range: " + average);

            // print out depth and curve names
            PrintWriter out = new PrintWriter(System.out, true); // needed for correct formatting
            out.println(buffer.toString());

            // display dialog box
            StsTextAreaDialog dialog = new StsTextAreaDialog(null, "Log Curve Range Listing for "
                    + curveName, buffer.toString(), 30, 40);
            dialog.setVisible(true);
        }
        catch(Exception e)
        {
            StsException.outputException("StsLogCurve.displayRange() failed.",
                e, StsException.WARNING);
        }
    }
/*
    public double getMinDepth()
    {
        if(depths == null) return doubleNullValue;
		checkLoadVector(depths);
        return (double)depths.getMinValue();
    }

    public double getMaxDepth()
    {
        if(depths == null) return doubleNullValue;
		checkLoadVector(depths);
        return (double)depths.getMaxValue();
    }
*/
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

    public int getDepthIndexGE(double depth)
    {
        float depthF = (float)depth;

        float[] depthArray = getDepthVector().getValues().getValues();
        return StsMath.arrayIndexAbove(depthF, depthArray);
    }

    public String getStsDirectory() { return stsDirectory; }

    public int getDepthIndexLE(double depth)
    {
        float depthF = (float)depth;
        float[] depthArray = getDepthVector().getValues().getValues();
        return StsMath.arrayIndexBelow(depthF, depthArray);
    }

    public double getValueFromPanelXFraction(double fraction)
    {
		if (logCurveType == null) return StsParameters.nullValue;
		if (logCurveType.name.equals("TIME/DEPTH")) return StsParameters.nullValue;

		if (logCurveType.name.equals("Interval Velocity"))
			return StsParameters.nullValue;
        double[] scale = logCurveType.getScale();
        double gridMin = scale[0];
        double gridMax = scale[1];
        if (logCurveType.isLinear())
            return gridMin + fraction * (gridMax - gridMin);
        else
        {
            double logGridMin = StsMath.log10(gridMin);
            double logGridMax = StsMath.log10(gridMax);
            return Math.pow(10, logGridMin + fraction * (logGridMax - logGridMin));
        }
    }

    /** draws two logs next to well, one on left, one on right.  From the viewer's perspective,
     *  logs are drawn from min to the left side to max at the right side.
     *
     * @param glPanel3d
     * @param well
     * @param origin zero origin of log axis: -1 for left and 0 for right.
     */
    public void display3d(StsGLPanel3d glPanel3d, StsWell well, float origin)
    {
        display3d(glPanel3d, well, origin, well.getMaxMDepth());
    }

    public void display3d(StsGLPanel3d glPanel3d, StsWell well, float origin, float mdepth)
    {
        displayLog3d(glPanel3d, well, origin, mdepth);
    }

    protected void displayLog3d(StsGLPanel3d glPanel3d, StsWell well, float origin)
    {
        displayLog3d(glPanel3d, well, origin, well.getMaxMDepth());
    }
    protected void displayLog3d(StsGLPanel3d glPanel3d, StsWell well, float origin, float mdLimit)
    {
        if(mdepthVector == null) return;
        float mdepths[] = mdepthVector.getFloats();
        if(mdepths == null) return;
        StsLogVector valueVector = getValueVector();
        if(valueVector == null) return;
        float values[] = valueVector.getFloats();
        if(values == null) return;

        float curveMin = valueVector.getMinValue();
        float curveMax = valueVector.getMaxValue();
        if(curveMin == curveMax) return;

        int logCurveWidth = wellClass.getLogCurveDisplayWidth();
        int logLineWidth = wellClass.getLogCurveLineWidth();

        int nValues = values.length;
        float displayCurveMin = logCurveType.getDisplayCurveMin();
        float displayCurveMax = logCurveType.getDisplayCurveMax();

        int scaleType = logCurveType.getScaleType();
        if(scaleType == StsLogCurveType.LOG)
        {
            displayCurveMin = (float)Math.log10(displayCurveMin);
            displayCurveMax = (float)Math.log10(displayCurveMax);
            curveMin = (float)Math.log10(curveMin);
            curveMax = (float)Math.log10(curveMax);
        }
        float scale = 1.0f/(displayCurveMax - displayCurveMin);
        float offset = -displayCurveMin*scale;
        boolean clip = curveMin < displayCurveMin || curveMax > displayCurveMax;
        GL gl = glPanel3d.getGL();
        logCurveType.getStsColor().setGLColor(gl);
        gl.glLineWidth((float)logLineWidth);
        gl.glDisable(GL.GL_LIGHTING);
        gl.glBegin(GL.GL_LINE_STRIP);
        boolean draw = true;
        for(int n = 0; n < nValues; n++)
        {
            float value = values[n];
            if(value == StsParameters.nullValue)
            {
                if(draw)
                {
                    gl.glEnd();
                    draw = false;
                }
                continue;
            }
            else if(!draw)
            {
                gl.glBegin(GL.GL_LINE_STRIP);
                draw = true;
            }
            float mdepth = mdepths[n];
            if(mdepth > mdLimit)
                break;
            StsPoint point = well.getPointAtMDepth(mdepth, true);
            double[] screenPoint = glPanel3d.getScreenCoordinates(point);
            StsPoint slope = well.getSlopeAtMDepth(mdepth);
            StsPoint slopePoint = StsPoint.multByConstantAddPointStatic(slope, 1000.0f, point);
            double[] screenSlopePoint = glPanel3d.getScreenCoordinates(slopePoint);

            // screen normal is negative reciprocal of screen slope
            double dsx = -(screenSlopePoint[1] - screenPoint[1]);
		    double dsy =  (screenSlopePoint[0] - screenPoint[0]);

            double s = Math.sqrt(dsx * dsx + dsy * dsy );

            if(s == 0.0)
            {
                dsx = 1.0;
                dsy = 0.0;
            }
            else
            {
                dsx /= s;
                dsy /= s;
            }
        /*
            else if(dsx < 0.0)
            {
                dsx = -dsx/s;
                dsy = -dsy/s;
            }
            else if(dsx > 0.0)
            {
                dsx = dsx/s;
                dsy = dsy/s;
            }
            else
            {
                dsx = 0.0;
                dsy = 1.0;
            }
        */
            if(scaleType == StsLogCurveType.LOG)
                value = (float)Math.log10(value);
            float scaledValue = value*scale + offset;
            if(clip) scaledValue = StsMath.minMax(scaledValue, 0.0f, 1.0f);
            screenPoint[0] += logCurveWidth*dsx*(origin + scaledValue);
            screenPoint[1] += logCurveWidth*dsy*(origin + scaledValue);
            double[] logPoint = glPanel3d.getWorldCoordinates(screenPoint);
            gl.glVertex3dv(logPoint, 0);
        }
        if(draw) gl.glEnd();
        gl.glEnable(GL.GL_LIGHTING);
    }

    public boolean hasColorscale() { return false; }
    public void resetVectors()
    {
        getValueVector().resetVector();
        getTimeVector().resetVector();
        getDepthVector().resetVector();
        getMDepthVector().resetVector();
    }

    public void applyPoints(ArrayList<StsPoint> points, StsWell well)
    {
        mdepthVector = new StsLogVector(StsLogVector.MDEPTH, points, 1);
        depthVector = well.getDepthsFromMDepths(mdepthVector);
        valueVector = new StsLogVector(name, points, 0);
    }
    public boolean addValuesToCurve(double[] mdepths, double[] values)
    {
    	String binaryDataDir = currentModel.getProject().getBinaryFullDirString();
    	getValueVector().writeAppend(values, binaryDataDir);
        getMDepthVector().writeAppend(mdepths, binaryDataDir);
    	return true;
    }

    public boolean addValuesToMDepthCurveWrite(double[] values)
    {
    	String binaryDataDir = currentModel.getProject().getBinaryFullDirString();
    	getMDepthVector().writeAppend(values, binaryDataDir);
    	return true;
    }
    public boolean addValuesToMDepthCurve(double[] values)
    {
    	getMDepthVector().append(values);
    	return true;
    }
    public boolean addValuesToCurveWrite(double[] values)
    {
    	String binaryDataDir = currentModel.getProject().getBinaryFullDirString();
    	getValueVector().writeAppend(values, binaryDataDir);
    	return true;
    }
    public boolean addValuesToCurve(double[] values)
    {
    	getValueVector().append(values);
    	return true;
    }
    public boolean replaceValueInCurve(double mdepth, double value, int index)
    {
    	String binaryDataDir = currentModel.getProject().getBinaryFullDirString();
    	getValueVector().writeReplaceAt(value, index, binaryDataDir);
        getMDepthVector().writeReplaceAt(mdepth, index, binaryDataDir);
    	return true;
    }
    public boolean replaceValueInCurve(double value, int index)
    {
    	String binaryDataDir = currentModel.getProject().getBinaryFullDirString();
    	getValueVector().writeReplaceAt(value, index, binaryDataDir);
    	return true;
    }
    public boolean replaceValueInMDepthCurve(double value, int index)
    {
    	String binaryDataDir = currentModel.getProject().getBinaryFullDirString();
    	getMDepthVector().writeReplaceAt(value, index, binaryDataDir);
    	return true;
    }
    /** draw a log curve */
/*
    public void display(StsLogTrack track) throws StsException
    {
        StsTrace.methodIn(this, "display");

        if (getValues()==null || getDepthValues()==null) return;

        StsWellWinContainer wellWinCntr = track.getWellDisplay().getWellWinContainer();
        GL gl = wellWinCntr.getGL();
        GLComponent glc = wellWinCntr.getGLComponent();
        GLU glu = wellWinCntr.getGLU();

        // set bounds
        Rectangle bounds = track.getBounds();
        gl.glPushMatrix();
        System.out.println("Viewport:  origin = " + bounds.x + ", " + bounds.y +
                            ", size = " + bounds.getSize().width + ", " +
                           bounds.getSize().height);
		gl.glViewport(glc, bounds.x, bounds.y, bounds.getSize().width,
                       bounds.getSize().height);
        gl.glScissor(bounds.x, bounds.y, bounds.getSize().width,
                     bounds.getSize().height);
        int mm[] = new int[1];
	    gl.glGetIntegerv(GL.GL_MATRIX_MODE, mm);
	    gl.glMatrixMode(GL.GL_PROJECTION);
	    gl.glLoadIdentity();
        System.out.println("Ortho:  xMin = " + values.getMinValue() +
                       ", xMax = " + values.getMaxValue() +
                       ", yMin = " + depths.getMaxValue() +
                       ", yMax = " + depths.getMinValue());
	    glu.gluOrtho2D((double)values.getMinValue(), (double)values.getMaxValue(),
                   (double)depths.getMaxValue(), (double)depths.getMinValue());
	    gl.glMatrixMode(mm[0]);

        // draw curve
        int nValues = values.getValues().getSize();
        System.out.println("\nlog = " + values.getName());
        System.out.println("number of values = " + nValues);
        if (nValues>0)
        {
            float[] zSamples = depths.getValues().getValues();
            float[] samples = values.getValues().getValues();
            StsColor color = this.getStsColor();
            gl.glColor3f(color.red, color.green, color.blue);
            float nullValue = values.getNullValue();
            gl.glBegin(GL.GL_LINE_STRIP);
                for (int i=0; i<nValues; i++)
                {
                    if (i%50==0) System.out.println(i + ": x = " + samples[i] + ", y = "
                                        + zSamples[i] + " ");
                    if (samples[i] != nullValue) gl.glVertex2f(samples[i], zSamples[i]);
                }
            gl.glEnd();
            zSamples=null;  // ok to garbage collect
            samples=null;   // ok to garbage collect
        }
        gl.glPopMatrix();
        Rectangle winBounds = wellWinCntr.getBounds();
		gl.glViewport(glc, winBounds.x, winBounds.y, winBounds.getSize().width,
                       winBounds.getSize().height);
        gl.glScissor(winBounds.x, winBounds.y, winBounds.getSize().width,
                     winBounds.getSize().height);

        StsTrace.methodOut(this, "display");
    }
*/
}












