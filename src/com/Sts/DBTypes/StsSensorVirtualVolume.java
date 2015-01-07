package com.Sts.DBTypes;

/**
 * <p>Title: Workflow development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2001</p>
 * <p>Company: 4D Systems LLC</p>
 * @author unascribed
 * @version 1.0
 */

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.Utilities.*;
import java.lang.reflect.*;

import java.text.*;
import java.util.Arrays;

public class StsSensorVirtualVolume extends StsVirtualVolume implements StsTreeObjectI
{
    private Object[] sensors = null;
    private StsTimeCurve curve = null;
    private byte accumMethod = ONOFF;
    private byte shapeType = SPHERE;
    private float xyOffset = 100.0f;
    private float zOffset = 100.0f;
    private float azimuth = 0.0f;
    private float majorMinorRatio = 1;
    private int areaScale = 1;
    private boolean scaleByAttribute = true;
    private boolean drawEllipses = false;
    private float sensorZMin = 0.0f;
    private float sensorZMax = 1000.0f;
    private float sensorZInc = 10.0f;
    private int totalEvents = 0;
    private float[][] sensorXyRange = new float[2][2];
    private float[] ampMax = new float[ACCUM_METHOD.length];
	private float[] ampMin = new float[ACCUM_METHOD.length];

    transient public static final byte ONOFF = 0;
    transient public static final byte CNT = 1;
    transient public static final byte SUM = 2;
    transient public static final byte AVG = 3;
    transient public static final byte MIN = 4;
    transient public static final byte MAX = 5;
    transient public static final byte DENSITY = 6;
    transient public static final byte AZIMUTH = 7;
    transient public static String[] ACCUM_METHOD = {"On/Off", "Count", "Sum", "Average", "Minimum", "Maximum", "Density (sqft)"};
	transient public static byte SPHERE = 0;
    transient public static byte CYLINDER = 1;
    transient public static String[] shapes = {"Sphere", "Cylinder"};
    transient int planeDir = ZDIR;
	transient float planeCoordinate = 0.0f;
	transient boolean rangeSet = false;
    transient StsProgressPanel ppanel = null;

    static StsObjectPanel virtualVolumeObjectPanel = null;
    static StsEditableColorscaleFieldBean vcolorscaleBean = new StsEditableColorscaleFieldBean(StsSensorVirtualVolume.class,"colorscale");

    static public final StsFieldBean[] virtualDisplayFields =
    {
        new StsBooleanFieldBean(StsVirtualVolume.class, "isVisible", "Enable"),
        new StsBooleanFieldBean(StsVirtualVolume.class, "readoutEnabled", "Mouse Readout"),
        vcolorscaleBean
    };
    
    static public final StsFieldBean[] virtualPropertyFields = new StsFieldBean[]
    {
        new StsStringFieldBean(StsSensorVirtualVolume.class, "name", true, "Name"),
        new StsComboBoxFieldBean(StsSensorVirtualVolume.class, "accumMethodString", "Accumulation Method:", ACCUM_METHOD),
        new StsComboBoxFieldBean(StsSensorVirtualVolume.class, "shape", "Shape:", shapes),
        new StsFloatFieldBean(StsSensorVirtualVolume.class, "xyOffset", false, "Horzontal:", false),
        new StsFloatFieldBean(StsSensorVirtualVolume.class, "zOffset", false, "Vertical:", false),
        new StsIntFieldBean(StsSeismicVolume.class, "nRows", false, "Number of Lines"),
        new StsIntFieldBean(StsSeismicVolume.class, "nCols", false, "Number of Crosslines"),
        new StsDoubleFieldBean(StsSeismicVolume.class, "xOrigin", false, "X Origin"),
        new StsDoubleFieldBean(StsSeismicVolume.class, "yOrigin", false, "Y Origin"),
        new StsFloatFieldBean(StsSeismicVolume.class, "xInc", false, "X Inc"),
        new StsFloatFieldBean(StsSeismicVolume.class, "yInc", false, "Y Inc"),
        new StsFloatFieldBean(StsSeismicVolume.class, "zInc", false, "Z Inc"),
        new StsFloatFieldBean(StsSeismicVolume.class, "xMin", false, "X Loc Min"),
        new StsFloatFieldBean(StsSeismicVolume.class, "yMin", false, "Y Loc Min"),
        new StsFloatFieldBean(StsSeismicVolume.class, "zMin", false, "Min Z or T"),
        new StsFloatFieldBean(StsSeismicVolume.class, "zMax", false, "Max Z or T"),
        new StsFloatFieldBean(StsSeismicVolume.class, "angle", false, "Angle to Line Direction"),
        new StsFloatFieldBean(StsSeismicVolume.class, "dataMin", false, "Data Min"),
        new StsFloatFieldBean(StsSeismicVolume.class, "dataMax", false, "Data Max"),
        new StsFloatFieldBean(StsSeismicVolume.class, "dataAvg", false, "Data Avg"),
    };
    
    public StsSensorVirtualVolume()
    {
        //System.out.println("SensorVirtualVolume constructor called.");
    }
    
    public StsSensorVirtualVolume(boolean persistent)
    {
        super(persistent);    	
    }
    
    public static StsSensorVirtualVolume constructor(String name, Object[] sensors, StsTimeCurve selectedCurve, byte accumType, float zStep, float offset)
    {
    	StsSeismicVolume[] volumes = (StsSeismicVolume[])currentModel.getCastObjectList(StsSeismicVolume.class);
        if(volumes.length == 0)
        {
        	new StsMessage(currentModel.win3d, StsMessage.ERROR, "Must have one seismic volume to get volume geometry.");
        	return null;
        }
        else
        {
        	return new StsSensorVirtualVolume(name, sensors, null, volumes[0], selectedCurve, accumType, zStep, offset, offset, SPHERE);
        }
    }

    public static StsSensorVirtualVolume constructor(String name, Object[] sensors, StsProgressPanel panel, StsSeismicBoundingBox volume, StsTimeCurve selectedCurve, byte accumType, float zStep, float xyoffset, float zoffset, byte shape)
    {
        return new StsSensorVirtualVolume(name, sensors, panel, volume, selectedCurve, accumType, zStep, xyoffset, zoffset, shape);
    }

    public StsSensorVirtualVolume(String name, Object[] sensors, StsProgressPanel panel, StsSeismicBoundingBox volume, StsTimeCurve selectedCurve, byte accumType, float zStep, float xyoffset, float zoffset, byte shape)
    {
        super(false);

        ppanel = panel;
        StsToolkit.copyDifferentStsObjectAllFields(volume, this);
		clearNonRelevantMembers();
		
    	setName(name);
    	this.sensors = sensors;
    	this.curve = selectedCurve;
    	this.accumMethod = accumType;
        this.xyOffset = xyoffset;
        this.zOffset = zoffset;
        this.shapeType = shape;

		getSensorRanges();
        ppanel.appendLine("Computed Z range " + sensorZMin + " to " + sensorZMax);

    	// compute z range from sensor data if volume only supports time.
    	if(volume.getZDomain() == StsProject.TD_TIME)
    	{
            // TODO: this volume needs to be congruent with existing project orthogonal grid

            this.zMin = sensorZMin;
            this.zMax = sensorZMax;
            this.zInc = zStep;
            initializedZ = true;
    	}
    	setZDomain(StsProject.TD_DEPTH);
        computeHalo();
        computeGridRanges();

        setName(name);
        colorscale = null;
        initializeColorscale();
        ppanel.appendLine("Computing data range....this may take a minute.");
        initialize(currentModel);
        isVisible = true;
        ppanel = null;
    }
    
    static public StsSensorVirtualVolumeClass getSensorVirtualVolumeClass()
    {
        return (StsSensorVirtualVolumeClass)currentModel.getCreateStsClass(StsSensorVirtualVolumeClass.class);
    }
    
    public void setDrawStimulated(boolean draw)
    {
    	drawEllipses = draw;
    }

    public float getSensorZMin() { return sensorZMin; }
    public float getSensorZMax() { return sensorZMax; }
    public void recomputeSensorRanges()
    {
        int count = 0;
		for(int i=0; i<sensors.length; i++)
		{
			StsDynamicSensor sensor = (StsDynamicSensor)sensors[i];
            count = count + sensor.getNumValues();
        }
        if(count != totalEvents)
        {
            StsCursor3d cursor = currentModel.getCursor3d();
            if(cursor != null)
                cursor.clearTextureDisplays();
            seismicColorList.setColorListChanged(true);
            totalEvents = count;
            getSensorRanges();
            computeHalo();
        }
    }

	public void getSensorRanges()
	{
		sensorZMin = StsParameters.largeFloat;
		sensorZMax = -StsParameters.largeFloat;
        sensorXyRange[0][0] = StsParameters.largeFloat;
        sensorXyRange[0][1] = -StsParameters.largeFloat;
        sensorXyRange[1][0] = StsParameters.largeFloat;
        sensorXyRange[1][1] = -StsParameters.largeFloat;
		for(int i=0; i<sensors.length; i++)
		{
			StsDynamicSensor sensor = (StsDynamicSensor)sensors[i];
            for(int j=0; j<sensor.getNumValues(); j++)
            {
                float value = sensor.getZorT(j);
                float[] xy = sensor.getXYZ(j);
                if(sensorZMax < value)
                    sensorZMax = value;
                if(sensorZMin > value)
                    sensorZMin = value;
                if(sensorXyRange[0][0] > xy[0])
                    sensorXyRange[0][0] = xy[0];
                if(sensorXyRange[0][1] < xy[0])
                    sensorXyRange[0][1] = xy[0];
                if(sensorXyRange[1][0] > xy[1])
                    sensorXyRange[1][0] = xy[1];
                if(sensorXyRange[1][1] < xy[1])
                    sensorXyRange[1][1] = xy[1];
            }
		}

    }

    public void computeHalo()
    {
        // Need to add the radius of the search area
        sensorZMin = sensorZMin - zOffset;
        if(sensorZMin < zMin) sensorZMin = zMin;
        sensorZMax = sensorZMax + zOffset;
        if(sensorZMax > zMax) sensorZMax = zMax;

        sensorXyRange[0][0] = sensorXyRange[0][0] - xyOffset;
        if(sensorXyRange[0][0] < xMin) sensorXyRange[0][0] = xMin;
        sensorXyRange[0][1] = sensorXyRange[0][1] + xyOffset;
        if(sensorXyRange[0][1] > xMax) sensorXyRange[0][1] = xMax;

        sensorXyRange[1][0] = sensorXyRange[1][0] - xyOffset;
        if(sensorXyRange[1][0] < yMin) sensorXyRange[1][0] = yMin;

        sensorXyRange[1][1] = sensorXyRange[1][1] + xyOffset;
        if(sensorXyRange[1][1] > yMax) sensorXyRange[1][1] = yMax;

	}
    public void computeGridRanges()
    {
        xMin = StsMath.intervalRoundDown(xMin, xInc);
        xMax = StsMath.intervalRoundUp(xMax, xInc);
        nCols = Math.round((xMax - xMin)/xInc) + 1;
        setColNumMax(getColNumMin() + nCols - 1);
        setColNumInc(1.0f);

        yMin = StsMath.intervalRoundDown(yMin, yInc);
        yMax = StsMath.intervalRoundUp(yMax, yInc);
        nRows = Math.round((yMax - yMin)/yInc) + 1;
        setRowNumMax(getRowNumMin() + nRows - 1);
        setRowNumInc(1.0f);

        zMin = StsMath.intervalRoundDown(zMin, zInc);
        zMax = StsMath.intervalRoundUp(zMax, zInc);
        nSlices = Math.round((zMax - zMin)/zInc) + 1;
    }

    public void computeGridRangesOld()
    {
        nCols = (int)((xMax - xMin)/xInc) + 1;
        setColNumMax(getColNumMin() + nCols);
        setColNumInc(1.0f);
        nRows = (int)((yMax - yMin)/yInc) + 1;
        setRowNumMax(getRowNumMin() + nRows);
        setRowNumInc(1.0f);
        nSlices = (int)((zMax - zMin)/zInc) + 1;
    }

    public boolean initialize(StsModel model)
    {
    	boolean result = super.initialize(model);
    	//getSensorRanges();
    	initializeMinMax();
    	return result;
    }

 	public void initializeColorscale()
	{
		try
		{
			if (colorscale == null)
			{
				StsSpectrumClass spectrumClass = currentModel.getSpectrumClass();
			    colorscale = new StsColorscale("Seismic", spectrumClass.getSpectrum(spectrumClass.SPECTRUM_BWYR), dataMin, dataMax);
				colorscale.setEditRange(dataMin, dataMax);
			}
			seismicColorList = new StsColorList(colorscale);
			colorscale.addActionListener(this);
		}
		catch (Exception e)
		{
			StsException.outputException("StsSensorVirtualVolume.initializeColorscale() failed.", e, StsException.WARNING);
		}
	}

    public void initializeMinMax()
    {
        for(int i=0; i<ACCUM_METHOD.length; i++)
        {
        	ampMax[i] = -StsParameters.largeFloat;
        	ampMin[i] = StsParameters.largeFloat;
        }    	
        ampMax[ONOFF] = 2.0f;
        ampMin[ONOFF] = 0.0f;	
        ampMax[CNT] = 1.0f;
        ampMin[CNT] = 0.0f;  
        
    	float z = sensorZMin + (sensorZMax - sensorZMin)/2.0f;
        readBytePlaneData(ZDIR, z);

    	resetDataRange(ampMin[getAccumMethod()], ampMax[getAccumMethod()], false);
    }
    
    public static byte getAccumMethodFromString(String aType)
    {
    	for(int i=0; i<ACCUM_METHOD.length; i++)
    	{
    		if(aType.equalsIgnoreCase(ACCUM_METHOD[i]))
    			return (byte)i;
    	}
    	return -1;
    }
    
    public StsTimeCurve[] getProperties()
    {
    	if(sensors == null)
    		return null;
    	else
    		return ((StsSensor)sensors[0]).getPropertyCurvesExcludingXYZ();
    }
    
    public byte getAccumMethod() { return accumMethod; }
    public void setAccumMethod(byte method) { accumMethod = method; }
    public void setAccumMethodString(String sMethod) 
    {
    	setAccumMethod(getAccumMethodFromString(sMethod));
    	resetDataRange(ampMin[getAccumMethod()], ampMax[getAccumMethod()], false);
        dbFieldChanged("accumMethod", accumMethod);
		currentModel.viewObjectChanged(this, this);
		currentModel.viewObjectRepaint(this, this);
    }
    public float getAzimuth() { return azimuth; }
    public void setAzimuth(float az) 
    {	
    	azimuth = az; 
    	dbFieldChanged("azimuth", azimuth);
    }
    public float getMajorMinorRatio() { return majorMinorRatio; }
    public void setMajorMinorRatio(float ratio) 
    { 
    	majorMinorRatio = ratio; 
    	dbFieldChanged("majorMinorRatio", majorMinorRatio);
    }
    public boolean getScaleByAttribute() { return scaleByAttribute; }
    public void setScaleByAttribute(boolean scaleBy) 
    { 
    	scaleByAttribute = scaleBy; 
    	dbFieldChanged("scaleByAttribute", scaleByAttribute);
    }       
    public void setAreaScale(int scale) 
    { 
    	areaScale = scale; 
    	dbFieldChanged("areaScale", areaScale);
    }
    public int getAreaScale() { return areaScale; }    
    public void resetDataRange()
    {    	   	
    	resetDataRange(ampMin[getAccumMethod()], ampMax[getAccumMethod()], false);
		currentModel.viewObjectRepaint(this, this);
    }
    
    public String getAccumMethodString()
    {
    	return ACCUM_METHOD[accumMethod];
    }

    public byte[] readBytePlaneData(int dir, float dirCoordinate)
    {
        try
        {
        	int nPlanePoints = 0;
        	planeDir = dir;
        	planeCoordinate = dirCoordinate;
        	switch(dir)
        	{
        	case ZDIR:
                if((dirCoordinate > zMax) || (dirCoordinate < zMin))
                    return null;
        		nPlanePoints = nRows * nCols;
        		break;
        	case XDIR:
                if((dirCoordinate > xMax) || (dirCoordinate < xMin))
                    return null;
        		nPlanePoints = nRows * nSlices;
        		break;
        	case YDIR:
                if((dirCoordinate > yMax) || (dirCoordinate < yMin))
                    return null;
        		nPlanePoints = nCols * nSlices;
        		break;
        	default:
        		return null;
        	}
        	byte[][] planes = new byte[1][nPlanePoints];
            return processPlaneData(planes);
        }
        catch (Exception e)
        {
            StsException.outputException("StsSeismicVolume.readPlaneData() failed.",
                                         e, StsException.WARNING);
            return null;
        }
    }

    public byte[] processPlaneData(byte[][] planes)
    {
        try
        {
        	float value = 0.0f;
            float min = StsParameters.largeFloat;
            float max = -StsParameters.largeFloat;
            float[] xyz = new float[3];
        	int nPlanePoints = planes[0].length;
        	byte[] planePoints = new byte[nPlanePoints];
            float[] values = new float[nPlanePoints];
            float volume = (float)(4/3 * Math.PI * xyOffset * xyOffset * xyOffset);

            switch(planeDir)
            {
                case XDIR:
                    xyz[0] = planeCoordinate;
                    xyz[1] = yMin;
                    if((planeCoordinate > sensorXyRange[0][1]) || (planeCoordinate < sensorXyRange[0][0]))
                    {
                        Arrays.fill(planePoints, StsParameters.nullByte);
                        return planePoints;
                    }
                    int cnt = 0;
                    for(int row = 0; row < nRows; row++, xyz[1] += yInc)
                    {
                        xyz[2] = zMin;
                        for(int slice = 0; slice < nSlices; slice++, xyz[2] += zInc)
                        {
                            if((xyz[1] < sensorXyRange[1][1]) && (xyz[1] > sensorXyRange[1][0]) && (xyz[2] < sensorZMax) && (xyz[2] > sensorZMin))
                            {
                                values[cnt] = computeValue(xyz, volume);
                                if(values[cnt] != StsParameters.nullValue)
                                {
                                    if(min > values[cnt]) min = values[cnt];
                                    if(max < values[cnt]) max = values[cnt];
                                }
                                cnt++;
                            }
                            else
                            {
                                values[cnt++] = StsParameters.nullValue;
                                continue;
                            }
                        }
                    }
                    break;
                case YDIR:
                    xyz[1] = planeCoordinate;
                    xyz[0] = xMin;
                    if((planeCoordinate > sensorXyRange[1][1]) || (planeCoordinate < sensorXyRange[1][0]))
                    {
                        Arrays.fill(planePoints, StsParameters.nullByte);
                        return planePoints;
                    }
                    cnt = 0;
                    for(int col = 0; col < nCols; col++, xyz[0] += xInc)
                    {
                        xyz[2] = zMin;
                        for(int slice = 0; slice < nSlices; slice++, xyz[2] += zInc)
                        {
                            if((xyz[0] < sensorXyRange[0][1]) && (xyz[0] > sensorXyRange[0][0]) && (xyz[2] < sensorZMax) && (xyz[2] > sensorZMin))
                            {
                                values[cnt] = computeValue(xyz, volume);
                                if(values[cnt] != StsParameters.nullValue)
                                {
                                    if(min > values[cnt]) min = values[cnt];
                                    if(max < values[cnt]) max = values[cnt];
                                }
                                cnt++;
                            }
                            else
                            {
                                values[cnt++] = StsParameters.nullValue;
                                continue;
                            }
                        }
                    }
                    break;
                 case ZDIR:
                    xyz[2] = planeCoordinate;
                    xyz[1] = yMin;
                    if((planeCoordinate > sensorZMax) || (planeCoordinate < sensorZMin))
                    {
                        return null;
                    }
                    cnt = 0;
                    for(int row = 0; row < nRows; row++, xyz[1] += yInc)
                    {
                        if((planeCoordinate > sensorZMax) || (planeCoordinate < sensorZMin))
                        {
                            Arrays.fill(planePoints, StsParameters.nullByte);
                            return planePoints;
                        }
                        xyz[0] = xMin;
                        if(ppanel != null) ppanel.appendLine("   Processing row #" + row + " for z= " + planeCoordinate);
                        for(int col = 0; col < nCols; col++, xyz[0] += xInc)
                        {
                            if((xyz[0] < sensorXyRange[0][1]) && (xyz[0] > sensorXyRange[0][0]) && (xyz[1] < sensorXyRange[1][1]) && (xyz[1] > sensorXyRange[1][0]))
                            {
                                values[cnt] = computeValue(xyz, volume);
                                if(values[cnt] != StsParameters.nullValue)
                                {
                                    if(min > values[cnt]) min = values[cnt];
                                    if(max < values[cnt]) max = values[cnt];
                                }
                                cnt++;
                            }
                            else
                            {
                                values[cnt++] = StsParameters.nullValue;
                                continue;
                            }
                        }
                    }
                    break;
            }
            // Adjust the min and max if required.
            adjustMinMax(min,max);

            float rangeMin = ampMin[accumMethod];
            float rangeMax = ampMax[accumMethod];
            float scale = StsMath.floatToUnsignedByteScale(rangeMin, rangeMax);
            float scaleOffset = StsMath.floatToUnsignedByteScaleOffset(scale, rangeMin);
            int num = 0;
            for(int k=0; k<nPlanePoints; k++)
    		{
    			if(values[k] != StsParameters.nullValue)
                {
    				planePoints[k] = StsMath.floatToUnsignedByte254WithScale(values[k], scale, scaleOffset);
                    num++;
                }
                else
    				planePoints[k] = (byte)255;
    		}
            return planePoints;
        }
        catch (Exception e)
        {
            StsException.outputException("StsSeismicVolume.processBytePlaneData() failed.",
                                         e, StsException.WARNING);
            return null;
        }
    }

    private void adjustMinMax(float min, float max)
    {
         if((min < ampMin[accumMethod]) && (min != StsParameters.largeFloat))
            ampMin[accumMethod] = min;
         if((max > ampMax[accumMethod]) && (max != -StsParameters.largeFloat))
        	ampMax[accumMethod] = max;

        if(ampMin[accumMethod] == ampMax[accumMethod])
        {
            switch(accumMethod)
            {
        		case CNT:
                    ampMax[accumMethod]++;
                    break;
        	    case SUM:
                    ampMax[accumMethod] = ampMax[accumMethod] + ampMin[accumMethod];
                    break;
        		case AVG:
                    ampMax[accumMethod] = ampMax[accumMethod] + ampMin[accumMethod];
                    break;
                case MIN:
                    ampMax[accumMethod] = ampMax[accumMethod] + ampMin[accumMethod]*.10f;
                    break;
        		case MAX:
                    ampMax[accumMethod] = ampMax[accumMethod] + ampMin[accumMethod]*.10f;
                    break;
        		case DENSITY:
                    ampMax[accumMethod] = ampMax[accumMethod] + ampMin[accumMethod]*.10f;
        			break;
        		default:
        			break;
        	}
        }
        resetDataRange();
    }

    public float computeValue(float[] xyzort, float  cellVolume)
    {
        float[] att = null;
        float value = 0.0f;
        // Use only XY to compute volume since Z is highly prone to errors

            // Find all the events within the user specified offset of the cell
            float[] xyz = xyzort;
            if(getZDomain() == StsProject.TD_TIME)
                xyz = ((StsDynamicSensor)sensors[0]).getZfromT(xyzort);

            float[] newAtt = null;
            for(int j=0; j<sensors.length; j++)
        	{
                StsDynamicSensor sensor = (StsDynamicSensor)sensors[j];
                if(shapeType == SPHERE)
                    newAtt = sensor.getAttributeInRadius(curve, xyz , xyOffset);
                else
                    newAtt = sensor.getAttributeInCylinder(curve, xyz , xyOffset, zOffset);

                att = (float[])StsMath.arrayAddArray(att, newAtt);
            }
            if(att == null) return StsParameters.nullValue;
            if(att.length < 1) return StsParameters.nullValue;
             // Compute the value to place in the cell
            switch(accumMethod)
            {
        	    case ONOFF:
                    if(att.length > 0) return 1.0f;
                    else return StsParameters.nullValue;
        		case CNT:
                    return att.length;
        	    case SUM:
                    value = 0.0f;
                    for(int i=0; i<att.length; i++)
                        value += att[i];
                    return value;
        		case AVG:
                    value = 0.0f;
                    for(int i=0; i<att.length; i++)
                        value += att[i];
                    return value/att.length;
                case MIN:
                    value = StsParameters.largeFloat;
                    for(int i=0; i<att.length; i++)
                        if(value > att[i]) value = att[i];
                    if(value == StsParameters.largeFloat)
                        value = StsParameters.nullValue;
                    return value;
        		case MAX:
                    value = -StsParameters.largeFloat;
                    for(int i=0; i<att.length; i++)
                        if(value < att[i]) value = att[i];
                    if(value == -StsParameters.largeFloat)
                        value = StsParameters.nullValue;
                    return value;
        		case DENSITY:
        			// Compute the density of points around the surrounding points
                    value = att.length/cellVolume;
                    return value;
        		case AZIMUTH:
        		    // Compute the azimuth through the surrounding points
                    System.out.println("Azimuth - Not implemented yet.");
                    return StsParameters.nullValue;
        		default:
        			break;
        	}
        return StsParameters.nullValue;
    }

	public void display(StsGLPanel glPanel)
	{
		super.display(glPanel);
		if(drawEllipses) drawEllipses((StsGLPanel3d)glPanel);
	}
	
    public void drawEllipses(StsGLPanel3d glPanel3d)
    {
    	float zCoor = currentModel.win3d.getCursor3d().getCurrentDirCoordinate();
    	byte[] plane = readBytePlaneData(ZDIR, zCoor);
    	float[] xyz = new float[3];
    	float[] xy = null;
    	int row, col;
    	float xRadius = 1.0f * getAreaScale();
    	float yRadius = getMajorMinorRatio() * getAreaScale();
    	float scale = getAreaScale();
    	float range = ampMax[accumMethod] - ampMin[accumMethod];
    	for(int i=0; i<plane.length; i++)
    	{
    		if(plane[i] == -2) 
    			continue; 
    		
    		row = i/nCols;
    		col = i - (row*nCols);
    		xy = getXYCoors(row, col);
    		if(scaleByAttribute)
    			scale = getScaledValue(plane[i])/range * getAreaScale();
    		else
    			scale = getAreaScale();
        	xRadius = 1.0f * scale;
        	yRadius = getMajorMinorRatio() * scale;    		
    		xyz[0] = xy[0];
    		xyz[1] = xy[1];
    		xyz[2] = zCoor;
    		StsGLDraw.drawEllipse(xyz, StsColor.BLUE, glPanel3d, (int)xRadius, (int)yRadius, getAzimuth(), 2.0f);
    		StsGLDraw.drawFilledEllipse(xyz, StsColor.BLUE, glPanel3d, (int)xRadius, (int)yRadius, getAzimuth(), 2.0f);
    	}
    }
    
    private int computePlaneIndex(StsPoint sensorPt)
    {
    	int dimOne = 0, dimTwo = 0;
    	// Determine the index of the sensor point on the plane
    	switch(planeDir)
    	{
    		case ZDIR:
    			// Verify the point is on the plane
    			if(getNearestSliceCoor(sensorPt.getZorT()) != getNearestSliceCoor(planeCoordinate))
    			{
    				return -1;
    			}
    			
    			// Determine the position of the point on the plane
    			dimOne = getNearestRowCoor(sensorPt.getY());
    			dimTwo = getNearestColCoor(sensorPt.getX());
    			return (dimOne * nCols) + dimTwo;
    		case XDIR:
    			if(getNearestColCoor(sensorPt.getX()) != getNearestColCoor(planeCoordinate))
    				return -1;
    			
    			dimOne = getNearestRowCoor(sensorPt.getY());
    			dimTwo = getNearestSliceCoor(sensorPt.getZorT());
    			return (dimOne * nSlices) + dimTwo;	
    		case YDIR:
    			if(getNearestRowCoor(sensorPt.getY()) != getNearestRowCoor(planeCoordinate))
    				return -1;
    			
    			dimOne = getNearestColCoor(sensorPt.getX());
    			dimTwo = getNearestSliceCoor(sensorPt.getZorT());
    			return (dimOne * nSlices) + dimTwo;   					
    		default:
    			return -1;
    	}
    }
    
    public String computeVolume(StsProgressPanel panel)
    {
    	DecimalFormat volumeFormat = new DecimalFormat("#,###."); 
    	String unitString = "cuft";
    	if(currentModel.getProject().getDepthUnits() == StsParameters.DIST_METER)
    		unitString = "cumtr";
    	byte[] plane = null;
    	float slice = 0.0f;
    	float incVol = xInc * yInc * zInc;
    	float volume = 0.0f;
    	panel.setMaximum(nSlices);    	
    	panel.setCount(0);
    	for(int i=0; i<nSlices; i++)
    	{
    		slice = i*zInc;
    	    plane = readBytePlaneData(ZDIR, slice);
    	    if(plane == null)
    	    	continue;
    	    for(int j=0; j<plane.length; j++)
    	    {
    	    	if(plane[j] != -2)
    	    	{
    	    		volume += incVol;   	    		
    	    	}
    	    }
    	    panel.incrementCount();
    	    panel.setDescription(" Volume= " + volumeFormat.format(volume) + " " + unitString);
    	}
    	panel.finished();
    	return volumeFormat.format(volume) + " " + unitString;
    }

    public void setShape(String shape)
    {
        for(byte i = 0; i < 2; i++)
            if(shape == shapes[i])
                shapeType = i;
    }

    public String getShape() { return shapes[shapeType]; }
    
    public StsFieldBean[] getDisplayFields()
    {
        return virtualDisplayFields;
    }

    public StsFieldBean[] getPropertyFields()
    {
        return virtualPropertyFields;
    }

    public StsObjectPanel getObjectPanel()
    {
        if (virtualVolumeObjectPanel == null)
        {
            virtualVolumeObjectPanel = StsObjectPanel.constructor(this, true);
        }
        return virtualVolumeObjectPanel;
    }

}