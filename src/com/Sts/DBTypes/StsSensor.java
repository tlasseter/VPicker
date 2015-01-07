//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.SensorLoad.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.Sounds.*;
import com.Sts.UI.*;
import com.Sts.Utilities.Shaders.*;
import com.Sts.Utilities.*;
import info.monitorenter.gui.chart.*;

import javax.media.opengl.glu.*;
import javax.media.opengl.*;

import javax.swing.*;
import java.awt.event.*;
import java.awt.*;
import java.io.*;
import java.util.*;
import java.text.SimpleDateFormat;

public class StsSensor extends StsMainTimeObject implements StsSelectable, StsTreeObjectI, StsTimeEnabled, StsMonitorable // , StsViewable
{
    transient StsProgressPanel progressPanel = StsProgressPanel.constructor();
    transient JDialog progressDialog = null;
    /** Processing realtimne data
    transient boolean processing = false;
    /** List of time dependent curves or attributes */
    protected StsObjectRefList timeCurves = null;
    /** Strings for time and depth domains */
    static final String[] TD_BOTH_STRINGS = StsParameters.TD_BOTH_STRINGS;
    /** X shift for all X coordinate values */
    protected float xShift = 0.0f;
    /** Y shift for all Y coordiante values */
    protected float yShift = 0.0f;
    /** Z shift for all Z values */
    protected float zShift = 0.0f;
    /** Color vector name used to persist the color by setting */
    protected String colorVectorName = "None";
    /** Property vector name used to persist the size by setting */
    protected String propertyVectorName = "None";
    /** X, Y origin used to convert relative coordinates to actual */
    protected double xOrigin, yOrigin;
    /** Well associated with this sensor data */
    protected StsWell well = null;
    /** Alarms */
    protected StsAlarm[] alarms = null;
    /** Hand Edited Events */
    protected byte[] handEdits = null;
    /** Coordinate type - dynamic or static */
    protected byte coordinateType = NONE;
    /** Are the supplied coordiantes relative */
    protected boolean isRelative = false;
    /** Currently supported Z domains (time, depth) */
    protected byte zDomainSupported = StsParameters.TD_TIME_DEPTH;
    /** The Z domain of the original sensor data */
    protected byte zDomainOriginal = StsParameters.TD_NONE;
    /** The horizontal units of the original sensor data */
    protected byte nativeHorizontalUnits = StsParameters.DIST_FEET;
    /** The vertical units of the original sensor data */
    protected byte nativeVerticalUnits = StsParameters.DIST_FEET;
    /** Display type - All events prior to current time, a range of events or a single event closest to the current time */
    protected byte displayType = SHOW_ALL;
    /** The range of events to show when displayType is set to INCREMENT */
    protected int displayDuration = 3600000; // One Hour
    /** Display symbol */
    protected byte pointType = CUBE;
    /** The data minimum of the currently selected sensor attribute */
    protected double dataMin = StsParameters.largeDouble;
    /** The data maximum of the currently selected sensor attribute */
    protected double dataMax = -StsParameters.largeDouble;
    /** The scaled minimum of the currently selected property attribute */
    protected float scaleMin = StsParameters.largeFloat;
    /** The scaled maximum of the currently selected property attribute */
    protected float scaleMax = -StsParameters.largeFloat;
    /** The scale linear or log10 */
    protected boolean scaleLinear = true;
    /** The scaled minimum of the currently selected color attribute */
    protected float colorMin = StsParameters.largeFloat;
    /** The scaled maximum of the currently selected color attribute */
    protected float colorMax = -StsParameters.largeFloat;
    /** The color scale linear or log10 */
    protected boolean colorLinear = true;
    /** Clip the current data range */
    protected boolean clipRange = true;
    /** The maximum time duration (seconds) of real-time time series plots of this sensor */
    protected int timeSeriesDuration = 18000;  // default is 30 minutes
    /** A scale factor used to adjust the relative size of sensor events */
    protected int numberBins = 100;
    /** The default spectrum used to color sensor events */
    protected String spectrumName = "Basic";
    /** The selected colorscale */
    protected StsColorscale colorscale;
    /** Does the provided sensor data have a date or just times */
    protected boolean hasDate = true;
    /** A string containing the date the data was collected when it is not provided in the file. */
    protected String dateString = null;
    /** Last modified date when sensor is realtime accumulated */
    protected long lastModified = 0;
    /** Y Axis type (log or linear) for time series plots. */
    protected byte yAxisType = LINEAR;
    /** Area computation method for time series plots */
    protected byte areaComputationMethod = Chart2D.NONE;
    /** Compute additional attributes on load */
    protected boolean computeCurves = false;
    /** Enable verbose Realtime Logging*/
    protected boolean isRealtimeLogging = false;
    /** Enable best fit plane through events */
    boolean bestFitPlane = false;

    /** Currently selected attribute vector to size events */
    protected StsTimeCurve propertyVector = null;
    /** Currently selected attribute vector to color events */
    protected StsTimeCurve colorVector = null;
    /** Current filter     */
    protected StsObject primaryFilter = null;
    protected StsObject secondaryFilter = null;

    //public StsTimeSeriesDisplayProperties timeSeriesDisplayProperties = null;
    static public final byte TWO_D = 0;
    static public final byte THREE_D = 1;
    static public final byte BOTH_D = 2;
    /** Coordinate types */
    static public final byte NONE = 0;
    static public final byte STATIC = 1;
    static public final byte DYNAMIC = 2;
    /** Display Types */
    static public final byte SHOW_SINGLE = 0;
    static public final byte SHOW_INCREMENT = 1;
    static public final byte SHOW_ALL = 2;
    static public final String[] displayTypeStrings = new String[] {"Single", "Increment", "All"};
    static public final String[] AMPLITUDE_KEYWORDS = new String[] {"Amplitude", "Amp", "Energy", "Ampl", "Magnitude", "Mag", "Magnevt"};

    /** Filename or exported file */
    transient String exportViewFilename = null;
    /** Export sensor events in view (true) or all events (false). */
    transient boolean exportView = false;
    /** Export progress bar dialog */
    transient StsProgressBarDialog progressBarDialog = null;
    /** Has the sensor been initialed. Used on database load */
    transient boolean initialized = false;
    /** Graphics */
    transient GLU glu;
    /** Array containing user selected size values */
    transient float[] sizeFloats = null;
    /** Array containing user selected color values */
    transient float[] colorFloats = null;
    /** Array containing time values */
    transient long[] timeLongs = null;
    /** Arrays containing beachball values */
    transient float[] dipFloats = null;
    transient float[] strikeFloats = null;
    transient float[] rakeFloats = null;

    /** Used to turn on debug messages */
    transient boolean debug = false;
    /** Accent time - accents event closest to this time */
    transient long durationAccentTime = 0l;
    transient long singleAccentTime = 0l;
    transient int accented = 0;
    transient boolean reloadRequired = false;

    /** Current event index */
    transient int sidx = 0;
    /** Is the current point and realtime loaded point */
    transient boolean realtimePoint = false;
    /** List of indices of the highlighted points */
    transient int[] highlightedPoints = null;
    /** Index of the current graphically picked point */
    transient int pickedIdx = -1;
    /** Array of values for a new attribute */
    transient float[] attribute = null;
    /** Name of the new attribute to add to sensor */
    transient String attributeName = null;
    transient boolean debugRealtime = false;
    transient int lastNumLines = 0;

    /** Cluster Analysis used to keep track of events in view to allow export of view. If this flag is on, clustering
    //  will impact display and export. */
    transient boolean clustering = false;
    /** Array containing values used to display and possibly export events in view. Values are assigned
    //  by the cluster analysis, proximity analysis and crossplot analysis workflow steps and may
    //  represent the number of events in the cluster, whether it is in proximity or by the polygon
    //  it is assigned to, respectively. */
    transient int[] clusters = null;

    /** Event symbols   */
    public static final byte SQUARE = 0;
    public static final byte CUBE = 1;
    public static final byte SPHERE = 2;
    public static final byte CYLINDER = 3;
    public static final byte DISK = 4;
    public static final byte TRIANGLE = 5;
    public static final byte STAR = 6;
    public static final byte DIAMOND = 7;
    public static final byte CIRCLE = 8;
    public static final byte NBSQUARE = 9;
    public static final byte BEACHBALL = 10;
    static public final String[] SYMBOL_TYPE_STRINGS = new String[] { "Square", "Cube", "Sphere", "Cylinder", "Disk", "Triangle", "Star", "Diamond", "Circle", "NoBorderSquare", "Beachballs"};
    static public final byte[] SYMBOL_TYPES = new byte[] { SQUARE, CUBE, SPHERE, CYLINDER, DISK, TRIANGLE, STAR, DIAMOND, CIRCLE, NBSQUARE, BEACHBALL };
    /** Time series Y Axis types */
    static public final byte LINEAR = 0;
    static public final byte LOG10 = 1;
    static public final byte LOGE = 2;
    static public final String[] AXIS_STRINGS = new String[] { "Linear", "Log10", "LogE"};
    static public final byte[] AXIS_TYPES = new byte[] { LINEAR, LOG10, LOGE };

    /** Required beachball fields */
    static public final String DIP = "DIP";
    static public final String STRIKE = "STRIKE";
    static public final String RAKE = "RAKE";

    /** Reference to the object panel */
    static protected StsObjectPanel objectPanel = null;

    /** List of all possible filters */
    transient StsObject[] filterList = new StsObject[0];
    /** List of all possible size attributes */
    transient StsTimeCurve[] propertyList = null;
    /** List of all possible color attributes including constant colors */
    transient StsTimeCurve[] colorByList;
    /** List of indices of points to export - used when exporting current view */
    transient int[] exportPoints = null;

    /** User interface beans to alow user to set object properties */
    static protected StsComboBoxFieldBean filterListBean;
    static protected StsComboBoxFieldBean sfilterListBean;
    static protected StsComboBoxFieldBean propertyListBean;
    static protected StsFloatFieldBean scaleMinBean;
    static protected StsFloatFieldBean scaleMaxBean;
    static protected StsBooleanFieldBean scaleLinearBean;
    static protected StsComboBoxFieldBean colorByBean;
    static protected StsFloatFieldBean colorMinBean;
    static protected StsFloatFieldBean colorMaxBean;
    static protected StsBooleanFieldBean colorLinearBean;
    static protected StsComboBoxFieldBean symbolListBean;
    static protected StsComboBoxFieldBean axisTypeListBean;
    static protected StsComboBoxFieldBean displayTypeBean;
    static protected StsEditableColorscaleFieldBean colorscaleBean;
    static protected StsDoubleFieldBean xBean;
    static protected StsDoubleFieldBean yBean;
    static protected StsDoubleFieldBean zBean;

    static public StsFieldBean[] displayFields = null;
    static public StsFieldBean[] propertyFields = null;

    /** default constructor */
    public StsSensor()
    {
    }

    public StsSensor(StsWell well, String name)
    {
        setName(name);
        numberBins = getSensorClass().getDefaultSize();
        pointType = getSensorClass().getDefaultPointType();
        setIsVisible(true);
        setWell(well);
        nativeHorizontalUnits = currentModel.getProject().getXyUnits();
        nativeVerticalUnits = currentModel.getProject().getDepthUnits();
        initialized = true;
    }

    static public StsSensor nullSensorConstructor(String name)
    {
        return new StsSensor(null, name);
    }

    /**
     * Set the well associated with this sensor (optional)
     * @param well object
     */
    public void setWell(StsWell well)
    {
        this.well = well;
    }
    
    /**
     * Get the Time shift
     * @return time shift
     */
    public int getTimeShift()
    {
        if(timeCurves != null)
        {
            StsTimeCurve timeCurve = (StsTimeCurve)timeCurves.getElement(0);
            return timeCurve.getTimeVector().getTimeShift() / 1000; // Convert user specified seconds to milli-seconds
        }
        return 0;
    }

    /**
     * Get the X shift
     * @return x shift
     */
    public float getXShift()
    {
        return xShift;
    }

    /**
     * Get the Y shift
     * @return y shift
     */
    public float getYShift()
    {
        return yShift;
    }
    /**
     * Get the Z shift
     * @return z shift
     */
    public float getZShift()
    {
        return zShift;
    }
    /**
     * Set a user supplied time shift for all time values
     * @param time - time shift
     */
    public void setTimeShift(int time)
    {
        if(timeCurves != null)
        {
            for(int i=0; i<getNTimeCurves(); i++)
            {
                StsTimeCurve timeCurve = (StsTimeCurve)timeCurves.getElement(i);
                timeCurve.getTimeVector().setTimeShift(time * 1000); // Convert user specified seconds to millseconds
                currentModel.viewObjectRepaint(this, this);                
            }
        }
        currentModel.viewObjectRepaint(this, this);
    }
    /**
     * Set a user supplied X shift for all X coordinates
     * @param x - X shift
     */
    public void setXShift(float x)
    {
        xShift = x;
        dbFieldChanged("xShift", xShift);
        currentModel.viewObjectRepaint(this, this);
    }
    /**
     * Set the born and death dates based on the time vector.
     */
    public void setBornDeathWithVector()
    {
        setBornDate(getTimeMin());
        setDeathDate(getTimeMax());
    }

    /**
     * Set compute curves flag. If set, cumulative amplitude and total events attributes will be added and computed automatically.
     * @param compute
     */
    public void setComputeCurves(boolean compute)
    {
        computeCurves = compute;
    }
    /**
     * Set the minimum value to scale the color range
     * used to scale colors the same across multiple sensors
     * @param min
     */
    public void setColorMin(float min)
    {
        colorMin = min;
        if(colorscale == null)
            return;
        colorscale.setRange((float)colorMin, (float)colorMax);
        dbFieldChanged("colorMin", colorMin);
        currentModel.viewObjectRepaint(this, this);
    }

    /**
     * Set the maximum value to scale the color range
     * used to scale colors the same across multiple sensors
     * @param max
     */
    public void setColorMax(float max)
    {
        colorMax = max;
        if(colorscale == null)
            return;
        colorscale.setRange((float)colorMin, (float)colorMax);
        dbFieldChanged("colorMax", colorMax);
        currentModel.viewObjectRepaint(this, this);
    }
    /**
      * Set the color of the events from supplied name
      * @param colorname - name of supported color
      */
     public void setColorFromString(String colorname)
     {
         setColorBy((StsTimeCurve)colorByBean.fromString(colorname));
         currentModel.refreshObjectPanel();
     }

    /**
     * Set the scaling type (log/linear) for color
     * @param val
     */
    public void setColorLinear(boolean val)
    {
        colorLinear = val;
        dbFieldChanged("colorLinear", colorLinear);
        currentModel.viewObjectRepaint(this, this);
    }
    public boolean getColorLinear() { return colorLinear; }

    /**
     * Set the minimum value to scale the event size
     * used to scale event size the same across multiple sensors
     * @param min
     */
    public void setScaleMin(float min)
    {
        scaleMin = min;
        dbFieldChanged("scaleMin", scaleMin);
        currentModel.viewObjectRepaint(this, this);
    }
    /**
     * Set the maximum value to scale the event size
     * used to scale event size the same across multiple sensors
     * @param max
     */
    public void setScaleMax(float max)
    {
        scaleMax = max;
        dbFieldChanged("scaleMax", scaleMax);
        currentModel.viewObjectRepaint(this, this);
    }
    /**
     * Set the scaling type (log/linear) for size
     * @param val
     */
    public void setScaleLinear(boolean val)
    {
        scaleLinear = val;
        dbFieldChanged("scaleLinear", scaleLinear);
        currentModel.viewObjectRepaint(this, this);
    }
    public boolean getScaleLinear() { return scaleLinear; }
    /**
     * Set a user supplied Y shift for all Y coordinates
     * @param y - Y shift
     */
    public void setYShift(float y)
    {
        yShift = y;
        dbFieldChanged("yShift", yShift);
        currentModel.viewObjectRepaint(this, this);
    }

    /**
     * Set a user supplied Z shift for all Z coordinates
     * @param z - Z shift
     */
    public void setZShift(float z)
    {
        zShift = z;
        dbFieldChanged("zShift", zShift);
        currentModel.viewObjectRepaint(this, this);
    }

    /**
     * A string containing the name of the sensor object
     * @return name of object
     */
    public String toString() { return getName().toString(); }

    /**
     * Get a string with the minimum time from the time vector
     * this method overides StsMainTimeObject method
     * @return minimum time in time vector
     */
    public String getBornDate()
    {
        String timeStg = StsDateFieldBean.convertToString(System.currentTimeMillis());
        if(timeCurves.getSize() > 0)
            timeStg = StsDateFieldBean.convertToString(((StsTimeCurve)timeCurves.getElement(0)).getTimeVector().getMinValue());
        return timeStg;
    }

    public long getBornDateLong()
    {
        long ldate = System.currentTimeMillis();
        if(timeCurves.getSize() > 0)
            ldate = ((StsTimeCurve)timeCurves.getElement(0)).getTimeVector().getMinValue();
        return ldate;
    }

    public long getDeathDateLong()
    {
        StsTimeCurve curve = (StsTimeCurve)timeCurves.getElement(0);
        if(curve == null)
            return super.getDeathDateLong();
        return curve.getTimeVector().getMaxValue();
    }

    public void setDisplayTypeString(String stype)
    {
        for(int i=0; i<displayTypeStrings.length; i++)
        {
            if(stype == displayTypeStrings[i])
            {
                if(displayType == i) return;
                displayType = (byte)i;
                /*
                if(displayType == SHOW_INCREMENT)
                    durationBean.setEditable(true);
                else
                    durationBean.setEditable(false);
                */
                dbFieldChanged("displayType", displayType);
                currentModel.viewObjectRepaint(this, this);
                return;
            }
        }
    }
    public void setAxisTypeString(String stype)
    {
        for(int i=0; i<AXIS_STRINGS.length; i++)
        {
            if(stype.equals(AXIS_STRINGS[i]))
                yAxisType = (byte)i;
        }
        dbFieldChanged("yAxisType", yAxisType);
        currentModel.viewObjectChanged(this, this);
    }
    //public void setAreaComputeTypeString(String stype)
    //{
    //    for(int i=0; i<Chart2D.timeEpochStrings.length; i++)
    //    {
    //        if(stype.equals(Chart2D.timeEpochStrings[i]))
    //            areaComputationMethod = (byte)i;
    //    }
    //    dbFieldChanged("areaComputationMethod", areaComputationMethod);
    //    currentModel.viewObjectChanged(this, this);
    //}
    public byte getAxisType() { return yAxisType; }
    //public byte getAreaComputeType() { return areaComputationMethod; }
    public String getDisplayTypeString()
    {
        return displayTypeStrings[displayType];
    }

    public String getAxisTypeString()
    {
        return AXIS_STRINGS[yAxisType];
    }
    //public String getAreaComputeTypeString()
    //{
    //    return Chart2D.timeEpochStrings[areaComputationMethod];
    //}
    public int getTimeSeriesDuration() { return timeSeriesDuration; }
    public void setTimeSeriesDuration(int numPts)
    {

    	timeSeriesDuration = numPts;
        dbFieldChanged("timeSeriesDuration", timeSeriesDuration);
        currentModel.viewObjectRepaint(this, this);
    }

    /**
     * Does nothing...required for born date bean
     * @param born
     */
    public void setBornDate(String born)
    {
        System.out.println("Setting");
    }

    /**
     * Get a string with the maximum time from the time vector
     * this method overides StsMainTimeObject.getDeathDate()
     * @return
     */
    public String getDeathDate()
    {
        StsTimeCurve curve = getTimeCurve(StsLogVector.types[StsLogVector.TIME]);
        return StsDateFieldBean.convertToString(curve.getTimeVector().getMaxValue());
    }
    /**
     * Does nothing...required for death date bean
     * @param death
     */
    public void setDeathDate(String death) { }

    /**
     * Get the current symbol type setting
     * @return symbol type as string
     */
    public String getSymbolString()
    {
        return SYMBOL_TYPE_STRINGS[pointType];
    }
    /**
     * Get the realtime logging value
     * @return is logging of realtime data on or off
     */
    public boolean getIsRealtimeLogging() { return isRealtimeLogging; }
    /**
     * Set the realtime logging flag
     * @param val
     */
    public void setIsRealtimeLogging(boolean val)
    {
        isRealtimeLogging = val;
    }
    /**
     * Set the current symbol type
     * @param symbolString
     */
    public void setSymbolString(String symbolString)
    {
        for(int i=0; i<SYMBOL_TYPE_STRINGS.length; i++)
        {
            if (symbolString.equals(SYMBOL_TYPE_STRINGS[i]))
            {
                // Verify Dip, Strike and Rake Vectors Exist
                if(i == BEACHBALL)
                {
                   if(!canDisplayBeachballs())
                   {
                       symbolListBean.setSelectedItem(SYMBOL_TYPE_STRINGS[pointType]);
                       return;
                   }
                }
                pointType = SYMBOL_TYPES[i];
                break;
            }
        }
        dbFieldChanged("pointType", pointType);
        currentModel.viewObjectRepaint(this, this);
        return;
    }

    /**
     * Verify Beachballs can be run
     */
    public boolean canDisplayBeachballs()
    {
        if(getTimeCurve(DIP, false) == null)
            return false;
        if(getTimeCurve(STRIKE, false) == null)
            return false;
        if(getTimeCurve(RAKE, false) == null)
            return false;
        return true;
    }
    /**
     * Get the currently supported Z domain (time, depth or both)
     * Once a velocity model is constructed all Z domains are supported.
     * @return
     */
    public byte getZDomainSupported()
    {
        return zDomainSupported;
    }

    /**
     * Does this sensor contain X, Y and Z vectors and therefore is a dynamic sensor
     * @return true - dynamic, false - static
     */
    public boolean canBeDynamic()
    {
        boolean hasX = false;
        boolean hasY = false;
        boolean hasZ = false;

        if(timeCurves == null) return true;
        for(int i=0; i<timeCurves.getSize(); i++)
        {
            StsTimeCurve tc = (StsTimeCurve)timeCurves.getElement(i);
            if(tc.getValueVector().getType() == StsLogVector.X)
                hasX = true;
            if(tc.getValueVector().getType() == StsLogVector.Y)
                hasY = true;
            else
            {
                for(int j=0; j<StsSensorKeywordIO.Z_KEYWORDS.length; j++)
                {
                    if(tc.getName().equalsIgnoreCase(StsSensorKeywordIO.Z_KEYWORDS[j]))
                       hasZ = true;
                }
            }
        }
        if(!hasX || !hasY || !hasZ)
        {
            return false;
        }
        return true;
    }

    /**
     * Get the index of the first attribute that is not X, Y or Z
     * @return
     */
    public int firstNonPositionalCurve()
    {
        for(int i=0; i<timeCurves.getSize(); i++)
        {
            boolean isZ = false;
            StsTimeCurve tc = (StsTimeCurve)timeCurves.getElement(i);
            if((tc.getValueVector().getType() != StsLogVector.X) &&
               (tc.getValueVector().getType() != StsLogVector.Y))
            {
                for(int j=0; j<StsSensorKeywordIO.Z_KEYWORDS.length; j++)
                {
                    if(tc.getName().equalsIgnoreCase(StsSensorKeywordIO.Z_KEYWORDS[j]))
                    {
                        isZ = true;
                    }
                }
                if(!isZ)
                    return i;  // Include NONE
            }
        }
        return 0;
    }

    /**
     * Get the vertical units that the data was supplied as a string
     * @return
     */
    public String getNativeVerticalUnitsString()
    {
        return StsParameters.DIST_STRINGS[nativeVerticalUnits];
    }
    /**
     * Get the vertical units that the data was supplied
     * @return
     */
    public byte getNativeVerticalUnits()
    {
        return nativeVerticalUnits;
    }
    /**
     * Get the horizontal units that the data was supplied as a string
     * @return
     */
    public String getNativeHorizontalUnitsString()
    {
        return StsParameters.DIST_STRINGS[nativeHorizontalUnits];
    }
    /**
     * Get the horizontal units that the data was supplied
     * @return
     */
    public byte getNativeHorizontalUnits()
    {
        return nativeHorizontalUnits;
    }
    /**
     * Get an array of the user available field beans
     * @return
     */
    public StsFieldBean[] getDisplayFields()
    {
       try
       {
           if (displayFields == null)
           {
               propertyListBean = new StsComboBoxFieldBean(StsSensor.class, "property", "Size By:", "propertyList");
               filterListBean = new StsComboBoxFieldBean(StsSensor.class, "primaryFilter", "Primary Filter:", "filterList");
               sfilterListBean = new StsComboBoxFieldBean(StsSensor.class, "secondaryFilter", "Secondary Filter:", "filterList");
       		   scaleMinBean = new StsFloatFieldBean(StsSensor.class, "scaleMin", true, "Scale Min");
    		   scaleMaxBean = new StsFloatFieldBean(StsSensor.class, "scaleMax", true, "Scale Max");
               scaleLinearBean = new StsBooleanFieldBean(StsSensor.class, "scaleLinear", true, "Linear Size Scale");
               colorByBean = new StsComboBoxFieldBean(StsSensor.class, "colorBy", "Color By:", "colorByList");
       		   colorMinBean = new StsFloatFieldBean(StsSensor.class, "colorMin", true, "Color Min");
    		   colorMaxBean = new StsFloatFieldBean(StsSensor.class, "colorMax", true, "Color Max");
               colorLinearBean = new StsBooleanFieldBean(StsSensor.class, "colorLinear", true, "Linear Color Scale");
               symbolListBean = new StsComboBoxFieldBean(StsSensor.class, "symbolString", "Symbol:", SYMBOL_TYPE_STRINGS);
               displayTypeBean = new StsComboBoxFieldBean(StsSensor.class, "displayTypeString", "Show:", displayTypeStrings);
               axisTypeListBean = new StsComboBoxFieldBean(StsSensor.class, "axisTypeString", "Series - Y Axis Type:", AXIS_STRINGS);
               colorscaleBean = new StsEditableColorscaleFieldBean(StsSensor.class, "colorscale");
               StsFloatFieldBean xShiftBean = new StsFloatFieldBean(StsSensor.class, "xShift", true, "X Shift:", true);
               StsFloatFieldBean yShiftBean = new StsFloatFieldBean(StsSensor.class, "yShift", true, "Y Shift:", true);
               StsFloatFieldBean zShiftBean = new StsFloatFieldBean(StsSensor.class, "zShift", true, "Depth Shift:", true);
               StsIntFieldBean timeShiftBean = new StsIntFieldBean(StsSensor.class, "timeShift", -2592000, 2592000, "Time Shift (seconds):", false);   // Max is 30 days

               //xBean = new StsDoubleFieldBean(StsSensor.class, "xLoc", true, "X:", false);
               //yBean = new StsDoubleFieldBean(StsSensor.class, "yLoc", true, "Y:", false);
               //zBean = new StsDoubleFieldBean(StsSensor.class, "zLoc", true, "Depth:", false);
               StsDateFieldBean bornField = new StsDateFieldBean(StsSensor.class, "bornDate", false, "Born Date:");
               bornField.setFormat(currentModel.getProject().getTimeDateFormat());
               StsIntFieldBean timeSeriesDurationBean = new StsIntFieldBean(StsSensor.class, "timeSeriesDuration", true, "Series - Realtime Duration (secs):", true);
               StsIntFieldBean numBinsBean = new StsIntFieldBean(StsSensor.class, "numberBins", 1, 5000, "Scale:",true);
               displayFields = new StsFieldBean[]
               {
                    new StsBooleanFieldBean(StsSensor.class, "isVisible", "Enable"),
                    new StsBooleanFieldBean(StsSensor.class, "isLogging", "Enable Real-time Log"),
                    new StsBooleanFieldBean(StsSensor.class, "fitPlane", "Show Best Fit Plane"),
                    filterListBean,
                    sfilterListBean,
                    propertyListBean,
                    scaleMinBean,
                    scaleMaxBean,
                    scaleLinearBean,
                    symbolListBean,
                    colorByBean,
                    colorMinBean,
                    colorMaxBean,
                    numBinsBean,
                    displayTypeBean,
                    axisTypeListBean,
                    //new StsComboBoxFieldBean(StsSensor.class, "areaComputeTypeString", "Series - Compute Area By:", Chart2D.timeEpochStrings),
                    timeSeriesDurationBean,
                    bornField,
                    //xBean, yBean, zBean,
                    xShiftBean, yShiftBean, zShiftBean, timeShiftBean,
                    colorscaleBean
               };
           }
            return displayFields;
        }
        catch (Exception e)
        {
            StsException.outputException("StsSensor.getDisplayFields() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public StsFieldBean[] getPropertyFields()
    {
        if(propertyFields == null)
        {
            propertyFields = new StsFieldBean[]
            {
                    new StsStringFieldBean(StsSensor.class, "nativeHorizontalUnitsString", false, "Native Horizontal Units:"),
                    new StsStringFieldBean(StsSensor.class, "nativeVerticalUnitsString", false, "Native Vertical Units:"),
                    new StsIntFieldBean(StsSensor.class, "numValues", false, "Number of Events:"),
            };
        }
        return propertyFields;
    }

    public Object[] getChildren()
    {
        return new Object[0];
    }

    public void setIsLogging(boolean val)
    {
        isRealtimeLogging = val;
    }
    public boolean getIsLogging() { return isRealtimeLogging; }
    /**
     * Are there any dependencies on this object
     * @return
     */
    public boolean anyDependencies()
    {
        return false;
    }
    /**
     * Get a reference to the object panel
     * @return
     */
    public StsObjectPanel getObjectPanel()
    {
        if (objectPanel == null)
        {
            objectPanel = StsObjectPanel.constructor(this, true);
        }
        return objectPanel;
    }
    /**
     * Called when the object is selected on the object panel.
     * @return
     */
    public void treeObjectSelected()
    {
        getSensorClass().selected(this);
    }

    /**
     * Get the sensor class containing this object
     * @return
     */
    public StsSensorClass getSensorClass()
    {
        return (StsSensorClass)currentModel.getCreateStsClass(this);
    }

    /**
     * Draw an event in 2D view
     * @param pt2d - Position of the point
     * @param color - Color of the point
     * @param gl - Graphics handle
     * @param size - Size of the point
     */
    public void draw2DPoint(float[] pt2d, StsColor color, float size, float scaleFactor, StsGLPanel3d glPanel3d, GL gl)
    {
        if(color.alpha == 0.0) return;
        if(pointType == DIAMOND)
        {
        	if(size > 100)
        		StsGLDraw.drawDiamond2d(pt2d, color, gl, StsBitmap.BIG);
        	else
        		StsGLDraw.drawDiamond2d(pt2d, color, gl, StsBitmap.SMALL);
        }
        else if(pointType == STAR)
        {
        	if(size > 100)
        		StsGLDraw.drawStar2d(pt2d, color, gl, StsBitmap.BIG);
        	else
        		StsGLDraw.drawStar2d(pt2d, color, gl, StsBitmap.SMALL);
        }
        else if(pointType == TRIANGLE)
        {
        	if(size > 100)
        		StsGLDraw.drawTriangle2d(pt2d, color, gl, StsBitmap.BIG);
        	else
        		StsGLDraw.drawTriangle2d(pt2d, color, gl, StsBitmap.SMALL);
        }
        else if(pointType == CIRCLE)
        {
            gl.glDisable(GL.GL_LIGHTING);
            StsGLDraw.drawDisk2d(glPanel3d, pt2d, color, size * scaleFactor);
        }
        else if(pointType == NBSQUARE)  // Needed because when there are too many points the border takes all the pixels (solid black)
        {
        	StsGLDraw.drawPoint2d(pt2d, color, gl, (int)size/10);
        }
        else
        {
            gl.glDisable(GL.GL_LIGHTING);
            StsGLDraw.drawDisk2d(glPanel3d, pt2d, color, size * scaleFactor);
        }
    }

    /**
     * Get the selected property (size) attribute as an array of floats
     * @return attribute array for sizing events
     */
    public float[] getPropertyAsFloats()
    {
        if(getProperty() != null)
        {
        	if(getProperty().getName().equals("Relative Time"))
        		return ((StsTimeCurve)timeCurves.getElement(0)).getTimeVector().getRelativeFloats();
        	else
        		return getProperty().getValuesVectorFloats();
        }
        return null;
    }
    /**
     * Get the property (size) attribute as an array of floats for the supplied time curve
     * @param curve - a time curve from the available attribute curves
     * @return attribute array for sizing events
     */
    public float[] getPropertyAsFloats(StsTimeCurve curve)
    {
        if(curve.getName().equals("Relative Time"))
        	return ((StsTimeCurve)timeCurves.getElement(0)).getTimeVector().getRelativeFloats();
        else
        	return curve.getValuesVectorFloats();
    }
    /**
     * Get the selected color attribute as an array of floats
     * @return attribute array for coloring events
     */
    public float[] getColorByAsFloats()
    {
        if(getColorBy() != null)
        {
        	if(getColorBy().getName().equals("Relative Time"))
        		return ((StsTimeCurve)timeCurves.getElement(0)).getTimeVector().getRelativeFloats();
            else
        	{
        		float[] cFloats = getColorBy().getValuesVectorFloats();
        		return cFloats;
        	}
        }
        return null;
    }

    public float[] getColorFloats() { return colorFloats; }

    //TODO need to implement this if we want axes drawn on beachballs....very messy though
    public void drawBeachballAxes(StsGLPanel3d glPanel3d)
    {

    }

    /**
     * Determine which event should be accented either most current point or realtime point.
     * @param timeLongs - Array of time values for the events
     * @param index - Index of the current point being displayed.
     * @return - Is it to be accented or not?
     */
    public boolean determineDurationAccent(long[] timeLongs, int index)
    {
    	if ((index + 1) >= timeLongs.length)
    	{
    		if(!realtimePoint)
    			return false;
    	    else
    	    {
//    	    	realtimePoint = false;
    		    return true;
    	    }
        }
        // Is it time enabled, if not plot all points
        if(getSensorClass().getEnableTime())
        {
            if(getSensorClass().getEnableAccent() || getSensorClass().getEnableSound() || getSensorClass().getEnableGoTo())
            {
            	if (timeLongs[index + 1] > currentModel.getProject().getProjectTime())
            	{
            		if(timeLongs[index] != durationAccentTime)
            		{
                        accented++;
            			if(accented%getSensorClass().getAccentDuration() == 0)
                            durationAccentTime = timeLongs[index];
                        return true;
            		}
            		else
                    {
                        accented = 0;
            			return false;
                    }
            	}
            }
        }
        return false;
    }

    /**
     * Determine which event should be accented either most current point or realtime point.
     * @param timeLongs - Array of time values for the events
     * @param index - Index of the current point being displayed.
     * @return - Is it to be accented or not?
     */
    public boolean determineSingleAccent(long[] timeLongs, int index)
    {
    	if ((index + 1) >= timeLongs.length)
    	{
    		if(!realtimePoint)
    			return false;
    	    else
    	    {
    	    	realtimePoint = false;
    		    return true;
    	    }
        }
        // Is it time enabled, if not plot all points
        if(getSensorClass().getEnableTime())
        {
            if(getSensorClass().getEnableAccent() || getSensorClass().getEnableSound() || getSensorClass().getEnableGoTo())
            {
            	if (timeLongs[index + 1] > currentModel.getProject().getProjectTime())
            	{
            		if(timeLongs[index] != singleAccentTime)
            		{
                        singleAccentTime = timeLongs[index];
                        return true;
            		}
            		else
            			return false;
            	}
            }
        }
        return false;
    }

    public void accentEvents(int index, StsSensorClass sensorClass, StsGLPanel3d glPanel3d, GL gl, float[] xyz, float size, float origSize, StsColor aColor, byte domain)
    {
        if(determineDurationAccent(timeLongs, index) && !clustering)
        {
            if(sensorClass.getEnableAccent())
            {
                if(domain == THREE_D)                
                    draw3DPoint(glPanel3d, xyz, size, aColor, false);
                else
                    draw2DPoint(xyz, aColor, size, 1.0f, glPanel3d, gl);
            }
        }
        if(determineSingleAccent(timeLongs, index) && !clustering)
        {
            if(sensorClass.getEnableGoTo())
                goToRelativePoint(xyz, sensorClass.getGoToOffset(), BOTH_D);
            if(sensorClass.getEnableSound())
               	StsSound.play(sensorClass.getDefaultSound());
        }
        if(handEdits != null)
        {
            if(handEdits[index] == 0)
            {
                if(domain == THREE_D)
                    draw3DPoint(glPanel3d, xyz, origSize * 1.05f, aColor, false);
                else
                    draw2DPoint(xyz, aColor, origSize * 1.05f, 1.0f, glPanel3d, glPanel3d.getGL());
            }
        }
    }

    /**
     * Determine if this sensor is displayable. May not be displayable because it does not support
     * the current domain or it is turned off.
     * @return
     */
    public boolean verifyDisplayable()
    {
        if(getNTimeCurves() == 0) return false;     // May be null if setup for real-time
        if(!currentModel.getProject().canDisplayZDomain(zDomainSupported))
            return false;

        if(zDomainOriginal != currentModel.getProject().getZDomain() &&
           currentModel.getProject().velocityModel == null)
            return false;

        if(isVisible == false)
            return false;

        return true;
    }

    /**
     * check and load the positioning, time and required attribute arrays.
     * @return - success
     */
    public boolean checkLoadVectors()
    {
        // If processing realtime data
        //if(processing)
        //    return false;

        if(pointType == BEACHBALL)
        {
            dipFloats = getTimeCurve(DIP, false).getValuesVectorFloats();
            strikeFloats =  getTimeCurve(STRIKE, false).getValuesVectorFloats();
            rakeFloats = getTimeCurve(RAKE, false).getValuesVectorFloats();
        }
        if(timeLongs == null)
        {
        	StsMessageFiles.infoMessage("Failed to load value vectors for sensor: " + this.getName());
        	return false;
        }
        return true;
    }

    public float[] getXYZValue(float x, float y, float z)
    {
        float zVal = z;
        try
        {
          if(zDomainOriginal != currentModel.getProject().getZDomain())
          {
              if(zDomainOriginal == StsProject.TD_DEPTH)
                  zVal = (float) currentModel.getProject().velocityModel.getT(x, y, z, 0.0f);
              else
                  zVal = (float) currentModel.getProject().velocityModel.getZ(x, y, z);
          }
          float[] xy = getRotatedPoint(x, y);
          return new float[] { xy[0], xy[1], zVal };
        }
        catch(Exception ex)
        {
            StsMessageFiles.infoMessage("StsSensor:getXYZValue() Unable to compute XYZ value");
        	return null;
        }
    }

    /**
     * Get an array of the unique numbers in the cluster array. Cluster array is built when clustering, proximity
     * analysis or crossplot analysis is run. The unique members in the array are required for display purposes.
     * @return  - Array of unique numbers contained in the cluster array.
     */
    public int[] getClusterNums()
    {
		if(clusters == null) return null;
    	int num = 0;
    	boolean notUnique = false;
        int[] clusterIndexRange = StsMath.arrayMinMax(clusters);
        int[] uniqueNums = new int[clusterIndexRange[1] + 1];
    	for(int i=0; i<uniqueNums.length; i++)
    		uniqueNums[i] = -1;

    	if((clusters == null) || (!getClustering()))
    		return null;
    	else
    	{
    		for(int i=0; i<clusters.length; i++)
    		{
    			if(clusters[i] != -1)
    			{
    				for(int j=0; j<uniqueNums.length; j++)
    				{
    					if(uniqueNums[j] == clusters[i])
    						notUnique = true;
    				}
    				if(!notUnique)
    					uniqueNums[num++] = clusters[i];
    			}
    			notUnique = false;
    		}
    		uniqueNums = (int[])StsMath.trimArray(uniqueNums, num);
    		return uniqueNums;
    	}
    }

    /**
     * Compute the color for the supplied value
     * @param value - attribute value
     * @return StsColor
     */
    public StsColor defineColor(float value)
    {
        // Define the color
    	int bin = 0;
        double interval = (colorMax- colorMin)/(double)(colorscale.getNColors()-1);
        if(value > colorMax) value = colorMax;
        if(value < colorMin) value = colorMin;
        if(colorFloats != null)
        {
            bin = (int)(Math.abs((value - colorMin))/interval);
            int colorIdx = bin % colorscale.getNColors();
            return colorscale.getStsColor(colorIdx);
        }
        else
        {
            //System.out.println("Getting color named " + colorVector.getName());
            if(colorVector != null)
                return StsColor.getColorByName32(colorVector.getName());
            else
                return StsColor.RED;
        }
    }

    /**
     * Set the value of a member of the cluster array. This may be assigned via a number of analysis functions including
     * cluster analysis, proximity analysis and crossplot analysis.
     * @param idx - index of the event
     * @param val - value to place in the cluster array.
     */
    public void setClusterValue(int idx, int val)
    {
    	if(clusters != null)
    		clusters[idx] = val;
    }

    /**
     * Get the array containing the clustering data.
     * @return - array of cluster results
     */
    public boolean getClustering() { return clustering; }

    /**
     * Turn clustering on or off.
     * @param val - true (clustering on), false (clustering off)
     */
    public void setClustering(boolean val) { clustering = val; }

    /**
     * Set the values of the cluster array from an external source
     * @param vals - array of values
     */
    public void setClusters(int[] vals)
    {
        if((clusters == null) || (vals == null))
            clusters = vals;
        else
        {
            for(int i=0; i<clusters.length; i++)
            {
                if(clusters[i] > -1)
                    clusters[i] = vals[i];
            }
        }
    }

    /**
     * Reset Cluster Array and turn clusters off
     */
    public void resetClusters()
    {
        clusters = null;
        clustering = false;
    }
    /**
     * Temporarily set the name and values of an attribute to effect display changes.
     * @param aVals - Values of the attribute
     * @param name - The name of the attribute
     */
    public void setAttribute(float[] aVals, String name)
    {
    	attributeName = name;
    	attribute = aVals;
    }

    /**
     * Get the array of clusters
     * @return integer array of cluster values
     */
    public int[] getClusters() { return clusters; }

    /**
     * Clear an element in the cluster array.
     * @param idx - index of an event
     */
    public void clearClusters(int idx)
    {
    	if(clusters == null) return;
    	for(int i=0; i<clusters.length; i++)
    	{
    		if(clusters[i] == idx)
    			clusters[i] = -1;
    	}
    }

    /**
     * Use the cluster array to limit the display to a particular attribute range
     * @param crv - Attribute curve
     * @param low - minimum attribute value
     * @param high - maximum attribute value
     * @return success
     */
    public boolean setAttributeLimits(StsTimeCurve crv, float low, float high)
    {
    	StsTimeCurve curve = getTimeCurve(crv.getName());
    	if(clusters == null)
    	{
    		clusters = new int[getNumValues()];
    		for(int i=0; i<getNumValues(); i++) clusters[i] = 99;
    	}
    	float[] vals = curve.getValuesVectorFloats();
    	for(int i=0; i<vals.length; i++)
    	{
            float value = vals[i] + (float)curve.getValueVector().getOrigin();
    		if((value < low) || (value > high))
    			clusters[i] = -1;
    	}
    	return true;
    }

    /**
     * Assign a color to all the events within a defined polygon.
     * @param xcurve - The X values for all the events in the sensor
     * @param xrange - The X range of the polygon
     * @param ycurve - the Y values for all the events in the sensor
     * @param yrange - The Y range of the polygon
     * @param color - Color to assign
     * @return - success
     */
    public boolean setClusterXplotRange(StsTimeCurve xcurve, double[] xrange, StsTimeCurve ycurve, double[] yrange, int color)
    {
    	if(clusters == null)
    	{
    		clusters = new int[xcurve.getNumValues()];
    		for(int i=0; i<xcurve.getNumValues(); i++) clusters[i] = -1;
    	}
    	float[] xvals = xcurve.getValuesVectorFloats();
    	float[] yvals = ycurve.getValuesVectorFloats();
    	for(int i=0; i<xvals.length; i++)
    	{
    		if((xvals[i] >= xrange[0]) && (xvals[i] <= xrange[1]) &&
    				(yvals[i] >= yrange[0]) && (yvals[i] <= yrange[1]))
    			clusters[i] = color;
    	}
    	return true;
    }

    /**
     * Determine if an event contained in this sensor meets the time criteria to display
     * @param timeLongs - All time values
     * @param index - Index of current event
     * @param nPoints - Number of points within range
     * @param duration - Time duration of display
     * @param plotType - Plot type - single event, increment or all
     * @return success
     */
    public boolean checkTime(long[] timeLongs, int index, int nPoints, long duration, byte plotType)
    {
        // Is it time enabled, if not plot all points
        if(getSensorClass().getEnableTime())
        {
            if ((timeLongs[index] > currentModel.getProject().getProjectTime())
            		|| (timeLongs[index] == StsParameters.nullLongValue))
                return false;

            if (plotType == SHOW_INCREMENT)
            {
                long endTime = currentModel.getProject().getProjectTime() - duration;
                if(timeLongs[index] < endTime)
                    return false;
            }
            else if (plotType == SHOW_SINGLE)
            {
                if (nPoints > 0)
                    return false;
                if ((index + 1) != timeLongs.length)
                {
                    if (timeLongs[index + 1] < currentModel.getProject().getProjectTime())
                        return false;
                }
            }
        }
        return true;
    }

    /**
     * Check that current event is within the displayed time range.
     * @param time - time of current event
     * @return success
     */
    public boolean checkValue(long time)
    {
        StsSensor valueSensor = getSensorClass().getValueSensor();
        StsTimeCurve valueTimeCurve = null;
        if(valueSensor != null)
           valueTimeCurve = (StsTimeCurve)valueSensor.getTimeCurves().getElement(getSensorClass().getValueIndex());

        // See if current sensor meets value criteria. Might have been set on a different sensor.
        if(valueTimeCurve != null)
        {
            float value = valueTimeCurve.getValueAt(getSensorClass().getValueIndex(), time);
            if((value < getSensorClass().getValueMin()) || (value > getSensorClass().getValueMax()))
                return false;
        }
        return true;
    }

    /**
     * Draw the current point in the 3D view
     * @param glPanel3d - Graphics Panel
     * @param xyz - Position of the event
     * @param sizeFloat - Size of event
     * @param color - Color of event
     * @param highlighted - Is the event to be highlighted.
     */
    public void draw3DPoint(StsGLPanel3d glPanel3d, float[] xyz, float sizeFloat, StsColor color, boolean highlighted)
    {
        draw3DPoint(glPanel3d, xyz, null, sizeFloat, color, highlighted);
    }

    public void draw3DPoint(StsGLPanel3d glPanel3d, float[] xyz, float[] dsr, float sizeFloat, StsColor color, boolean highlighted)
    {
        float zscale = glPanel3d.getZScale();
    	float sizef = sizeFloat;
        if(color.getRGBA()[3] == 0.0)
            return;   // Transparent....don't draw
    	if(highlighted)
    		sizef = sizeFloat + 2.0f;

        if(pointType == SQUARE)
        {
        	int size = (int)(sizeFloat/10.0f);
        	if(highlighted) size += 2;
            StsGLDraw.drawPoint(xyz, StsColor.BLACK, glPanel3d, size + 2);
            StsGLDraw.drawPoint(xyz, color, glPanel3d, size, 2.0);
        }
        if(pointType == BEACHBALL)
        {
            StsBeachballShader beachballShader = StsBeachballShader.getShader(glPanel3d.getGL());
            if(dsr == null)
            {
                return;
            }
            beachballShader.drawBeachball(glPanel3d, xyz, dsr, color, sizef);
        }
        if(pointType == NBSQUARE)  // Needed because when there are too many points the border takes all the pixels (solid black)
        {
        	int size = (int)(sizeFloat/10.0f);
        	if(highlighted) size += 2;
            StsGLDraw.drawPoint(xyz, color, glPanel3d, size, 2.0);
        }
        else if(pointType == CUBE)
        {
            StsGLDraw.drawCube(glPanel3d, xyz, color, sizef);
        }
        else if (pointType == SPHERE)
        {
            StsSphereShader sphereShader = StsSphereShader.getShader(glPanel3d.getGL());
            sphereShader.drawSphere(glPanel3d, xyz, color, sizef);
            // StsGLDraw.drawSphere(glPanel3d, xyz, color, sizef);
        }
        else if (pointType == CYLINDER)
            StsGLDraw.drawCylinder(glPanel3d, xyz, color, 50, sizef);
        else if (pointType == DISK)
            StsGLDraw.drawCylinder(glPanel3d, xyz, color, sizef, 50);
        else if(pointType == TRIANGLE)
        {
        	if(highlighted)
    			StsGLDraw.drawTriangle(xyz, color, glPanel3d, StsBitmap.BIG);

        	else
        	{
        		if(sizeFloat > 100.0)
        			StsGLDraw.drawTriangle(xyz, color, glPanel3d, StsBitmap.BIG);
        		else
        			StsGLDraw.drawTriangle(xyz, color, glPanel3d, StsBitmap.SMALL);
        	}
        }
        else if(pointType == STAR)
        {
        	if(highlighted)
    			StsGLDraw.drawStar(xyz, color, glPanel3d, StsBitmap.BIG);

        	else
        	{
        	if(sizeFloat > 100.0)
        		StsGLDraw.drawStar(xyz, color, glPanel3d, StsBitmap.BIG);
        	else
        		StsGLDraw.drawStar(xyz, color, glPanel3d, StsBitmap.SMALL);
        	}
        }
        else if(pointType == DIAMOND)
        {
        	if(highlighted)
    			StsGLDraw.drawDiamond(xyz, color, glPanel3d, StsBitmap.BIG);

        	else
        	{
        	if(sizeFloat > 100.0)
        		StsGLDraw.drawDiamond(xyz, color, glPanel3d,StsBitmap.BIG);
        	else
        		StsGLDraw.drawDiamond(xyz, color, glPanel3d,StsBitmap.SMALL);
        	}
        }
        else if(pointType == CIRCLE)
        {
        	if(highlighted)
    			StsGLDraw.drawCircle(xyz, color, glPanel3d, StsBitmap.BIG);
        	else
        	{
        	if(sizeFloat > 100.0)
        		StsGLDraw.drawCircle(xyz, color, glPanel3d,StsBitmap.BIG);
        	else
        		StsGLDraw.drawCircle(xyz, color, glPanel3d,StsBitmap.SMALL);
        	}
        }
    }

    /**
     * Remove a sensor from the instance list and in the 3d window
     **/
    public boolean delete()
    {
        if(timeCurves != null)
            timeCurves.deleteAll();
        return super.delete();
    }

    /**
     * Delete an event from the sensor
     *
     * UNDER DEVELOPMENT
     * @param index - index of event
     * @return success
     */
    public boolean deletePoint(Integer index)
    {
        new StsMessage(currentModel.win3d, StsMessage.WARNING, "Not available yet.");
/*
        for(int i=0; i<getNTimeCurves(); i++)
        {
            StsTimeCurve timeCurve = (StsTimeCurve)timeCurves.getElement(i);
            timeCurve.deletePoint(idx);
        }
 */
        return true;
    }

    /**
     * Get the real world coordinates from relative XY
     * @param x - x position
     * @param y - y position
     * @return - real world coordinates of the supplied relative xy
     */
    public float[] getRotatedPoint(float x, float y)
    {
        StsProject project = currentModel.getProject();
        float dXOrigin = (float) (xOrigin - project.getXOrigin());
        float dYOrigin = (float) (yOrigin - project.getYOrigin());

        float[] xy = project.getRotatedRelativeXYFromUnrotatedRelativeXY(dXOrigin + x + xShift, dYOrigin + y + yShift);
        return xy;
    }

    /**
     * Version comparasion method
     */
    static public final class VersionComparator implements Comparator
    {
        VersionComparator()
        {
        }

        // order by versions and then order alphabetically
        public int compare(Object o1, Object o2)
        {
            StsLogVector v1 = (StsLogVector) o1;
            if (v1 == null)
            {
                return -1;
            }
            StsLogVector v2 = (StsLogVector) o2;
            if (v2 == null)
            {
                return 1;
            }

            int vv1 = v1.getVersion();
            int vv2 = v2.getVersion();

            // compa
            if (vv1 > vv2)
            {
                return 1;
            }
            if (vv1 < vv2)
            {
                return -1;
            }

            String s1 = v1.getName();
            String s2 = v2.getName();
            return s1.compareTo(s2);
        }
    }

    /**
     * Create required null vectors for real-time sensor
     * @param names - array of vector names
     */
    private void createNullVectors(String[] names, StsSensorFile sFile)
    {
        int ampExist = -1;
        // Create vectors
        if(isRealtimeLogging)
        {
            for(int i=0; i<names.length; i++)
                StsMessageFiles.logMessage("   Adding " + names[i] + " attribute vector to sensor(" + getName() + ")");
        }
        StsSensorKeywordIO.filename = sFile.file.getFilename();
        StsSensorKeywordIO.parseAsciiFilename(sFile.file.getFilename());
        StsTimeVector tVector = StsSensorKeywordIO.constructTimeVector("sensor", sFile.getColLocation(sFile.TIME), sFile.getColLocation(sFile.TIME));
    	tVector.setMinMaxAndNulls(StsParameters.largeLong);

        if(computeCurves)
        {
            if(isRealtimeLogging) StsMessageFiles.logMessage("   Adding computed attributes vectors to sensor(" + getName() + ")");

            // Verify that amplitude exists before adding cumulative amplitude
            for(int i=0; i<AMPLITUDE_KEYWORDS.length; i++)
            {
                ampExist = StsMath.stringArrayContainsIndex(names,AMPLITUDE_KEYWORDS[i]);
                if(ampExist != -1) break;
            }

            if(ampExist != -1)
                names = (String[])StsMath.arrayAddElement(names, new String("CumulativeAmplitude"));
            names = (String[])StsMath.arrayAddElement(names, new String("EventTotal"));
        }
    	StsLogVector[] logVectors = StsSensorKeywordIO.constructLogVectors(names, "sensor");

        // Add the first value to the vectors
        if(computeCurves)
        {
            int numCurves = logVectors.length-2;
            if(ampExist == -1)
                numCurves = logVectors.length-1;

            for(int i=0; i<numCurves; i++)
                logVectors[i].setValues(new float[] {(float)sFile.currentValues[i]});

            if(ampExist != -1)
                logVectors[logVectors.length-2].setValues(new float[] {(float)sFile.currentValues[ampExist]});
            logVectors[logVectors.length-1].setValues(new float[] {1.0f});
        }
        else
        {
            for(int i=0; i<logVectors.length; i++)
            {
                logVectors[i].setValues(new float[] {(float)sFile.currentValues[i]});
                logVectors[i].setMinMaxAndNulls(currentModel.getProject().getLogNull());
            }
        }
        StsLongVector lVector = new StsLongVector(new long[] {(long)sFile.currentTime});
        tVector.setValues(lVector);
        tVector.setMinMaxAndNulls(StsParameters.largeLong);
        setBornDate((long)sFile.currentTime);

        // Initialize the binary files
        StsSensorKeywordIO.deleteBinaryFiles(currentModel.getProject().getBinaryFullDirString(), sFile.curveNames, StsLogVector.SENSOR_PREFIX);
        StsSensorKeywordIO.checkWriteBinaryFiles(tVector, logVectors, currentModel.getProject().getBinaryFullDirString());

        // Add the time curves to the sensor object
        StsTimeCurve[] curves = StsTimeCurve.constructTimeCurves(tVector, logVectors, 0);
        addTimeCurves(curves, false);
        addToProject(StsProject.TD_DEPTH);
        realtimePoint = true;

        clearHandEdits();

        //currentModel.viewObjectChanged(this, this);
        if(isRealtimeLogging) StsMessageFiles.logMessage("   Successfully added realtime curves to sensor(" + getName() + ")");
    }

    /**
     * Clear the handEdits
     */
    public void clearHandEdits()
    {
        if(handEdits == null)
            handEdits = new byte[getNumValues()];
        for(int i = 0; i<handEdits.length; i++)
            handEdits[i] = 1;
        return;
    }
    /**
     * Add time curves to the reference lists
     **/
    public void addTimeCurves(StsTimeCurve[] timeCurves, boolean computeCurves)
    {
        if(timeCurves == null)
            return;

        for (int n = 0; n < timeCurves.length; n++)
            addTimeCurve(timeCurves[n]);

        if(computeCurves)
            computeAttributeCurves();
    	initialize();
        StsColor.getColorByName("");

        String defaultColor = StsColor.getNameByStsColor32(getSensorClass().getDefaultSensorColor());
        StsTimeCurve colorVector = (StsTimeCurve)colorByBean.fromString(defaultColor);
        StsTimeCurve propertyVector = (StsTimeCurve)propertyListBean.fromString("None");

        setColorBy(colorVector);
        setProperty(propertyVector);

        numberBins = getSensorClass().getDefaultSize();
        pointType = getSensorClass().getDefaultPointType();

        dataMin = timeCurves[firstNonPositionalCurve()].getCurveMin();
        dataMax = timeCurves[firstNonPositionalCurve()].getCurveMax();

        if(getTimeCurve(StsLogVector.X) != null)
        {
            xOrigin = getTimeCurve(StsLogVector.X).getValueVector().getOrigin();
            yOrigin = getTimeCurve(StsLogVector.Y).getValueVector().getOrigin();
        }
        else
        {
            xOrigin = 0.0f;
            yOrigin = 0.0f;
        }
        handEdits = null;
        clearHandEdits();

        numberOfElements = timeCurves[0].getNumValues();
    }

    /**
     * Set the color of the events from supplied name
     * @param color - name of supported color
     */
    public void setColor(String color)
    {
    	setColorBy((StsTimeCurve)colorByBean.fromString(color));
    	currentModel.refreshObjectPanel();
    }

    /**
     * Set the data histogram associated with the color vector
     */
    public void setDataHistogram()
    {
        if(colorVector.checkLoadVectors())
            colorscaleBean.setHistogram(colorVector.getVectorHistogram());
    }

    /**
     * Handle colorscale changes
     * @param e - event
     */
    public void actionPerformed(ActionEvent e)
    {
        if (e.getSource() instanceof StsColorscale)
        {
            currentModel.viewObjectChangedAndRepaint(this, this);
        }
    }

    /**
     * Add an attribute to the sensor
     * @param timeCurve - curve contaning time and attribute values.
     */
    public void addTimeCurve(StsTimeCurve timeCurve)
    {
        if (timeCurve == null)
        {
            return;
        }
        if (timeCurves == null)
        {
            timeCurves = StsObjectRefList.constructor(10, 1, "timeCurves", this);
        }

        if(well != null)
            timeCurve.setWell(well);
        timeCurves.add(timeCurve);
    }

    /**
     * Get a list of the attribute curves
     * @return list of attribute curves
     */
    public StsObjectRefList getTimeCurves()
    {
        if (timeCurves == null)
        {
            timeCurves = StsObjectRefList.constructor(2, 2, "timeCurves", this);
        }
        return timeCurves;
    }

    /**
     * Get the number of attributes contained in this sensor
     * @return - number of attributes
     */
    public int getNTimeCurves()
    {
        return (timeCurves == null) ? 0 : timeCurves.getSize();
    }

    /**
     * Get the number of events
     * @return - number of events
     */
    public int getNumValues()
    {
        if(timeCurves == null) return 0;
        if(timeCurves.getSize() == 0) return 0;

        if(!checkLoadVectors())
        	return 0;
        int num = ((StsTimeCurve)timeCurves.getElement(0)).getNumValues();
        return num;
    }

    /**
     * Copy all the attribute curves to a new array
     * @return - array of attribute curves containing time and values
     */
    public StsTimeCurve[] copyTimeCurveArray()
    {
        if (timeCurves == null)
        {
            return null;
        }
        return (StsTimeCurve[]) timeCurves.copyArrayList(StsTimeCurve.class);
    }

    /**
     * Get a list containing the names of all attributes in sensor
     * @return - array of attribute names
     */
    public String[] getTimeCurveList()
    {
        StsTimeCurve[] curveArray = copyTimeCurveArray();
        if (curveArray == null)
        {
            return null;
        }
        String[] curveList = new String[curveArray.length];
        for (int i = 0; i < curveArray.length; i++)
        {
            curveList[i] = curveArray[i].getName();
        }
        return curveList;
    }

    /**
     *  Find an attribute curve in the list
     *  @return time curve with attribute and time values
     **/
    public StsTimeCurve getTimeCurve(String name)
    {
        return getTimeCurve(name, true);
    }

    /**
     *  Find an attribute curve in the list
     *  @return time curve with attribute and time values
     **/
    public StsTimeCurve getAmplitudeCurve()
    {
        String ampStg = verifyAmplitude();
        if(ampStg == null)
            return null;
        return getTimeCurve(ampStg, true);
    }

    public StsTimeCurve getTimeCurve(String name, boolean fullString)
    {
        if (name == null)
        {
            return null;
        }
        if (timeCurves == null)
        {
            return null;
        }
        int nCurves = timeCurves.getSize();
        for (int i = 0; i < nCurves; i++)
        {
            StsTimeCurve curve = (StsTimeCurve) timeCurves.getElement(i);
            if(fullString)
            {
                if(name.equalsIgnoreCase(curve.getName()))
                    return curve;
            }
            else
            {
                if(curve.getName().toLowerCase().startsWith(name.toLowerCase()))
                    return curve;
            }
        }
        return null;
    }

    /**
     * Get the time curve of a specific type (X, Y, Z, ...)
     * @param type - Recognized data types such as time and positional types
     * @return curve containing requested data
     */
    public StsTimeCurve getTimeCurve(byte type)
    {
        if (timeCurves == null)
        {
            return null;
        }
        int nCurves = timeCurves.getSize();
        for (int i = 0; i < nCurves; i++)
        {
            StsTimeCurve curve = (StsTimeCurve) timeCurves.getElement(i);
            if(curve.getValueVector().getType() == type)
                return curve;
        }
        return null;
    }
    /**
     * Get the time curve for a specific index
     * @param idx - index of desired attribute curve
     * @return curve containing requested data
     */
    public StsTimeCurve getTimeCurve(int idx)
    {
        if (timeCurves == null) return null;
        if(timeCurves.getSize() < idx) return null;
        return (StsTimeCurve) timeCurves.getElement(idx);
    }

    /**
     * Set the visibility of this sensor
     * @param b - on/off
     */
    public void setIsVisible(boolean b)
    {
        if (b == isVisible)
        {
            return;
        }
        isVisible = b;
        dbFieldChanged("isVisible", isVisible);
        currentModel.win3dDisplayAll();
    }

    /**
     * Open a popup window when an event is selected
     */
    public void openOrPopWindow()
    {
        System.out.println("Delete point");
    }

    /**
     * Highlight a point
     * @param pointIdx - index of the point to highlight
     */
    public void highlight(int pointIdx)
    {
    	if(highlightedPoints == null)
    	{
    		highlightedPoints = new int[1];
    		highlightedPoints[0] = pointIdx;
    	}
    	else
    	{
    		highlightedPoints = (int[])StsMath.arrayAddElement(highlightedPoints, pointIdx);
    	}
    }

    /**
     * Clear all highlighted points
     */
    public void clearHiglightedPoints()
    {
    	highlightedPoints = null;
    }

    /**
     * Initialize the sensor on database load.
     * @param model
     * @return success
     */
    public boolean initialize(StsModel model)
    {
        if (initialized)
        {
            return true;
        }
        initialize();
        initialized = true;
 /*
        if(timeSeriesDisplayProperties != null)
        {
        	timeSeriesDisplayProperties.setParentObject(this);
            long min = ((StsTimeCurve)timeCurves.getElement(0)).getTimeVector().getMinValue();
            long max = ((StsTimeCurve)timeCurves.getElement(0)).getTimeVector().getMaxValue();
    		timeSeriesDisplayProperties.initializeTimeRange(min, max);
        }
 */
        initColorVector();
        initPropertyVector();

        colorVector.checkLoadVectors();
        propertyVector.checkLoadVectors();

        dataMin = colorVector.getCurveMin();
        dataMax = colorVector.getCurveMax();

        // Not sure why this is needed, but does not appear to persist colorscale settings
//        colorscale = new StsColorscale("Sensor", currentModel.getSpectrum(spectrumName),(float)dataMin, (float)dataMax);
//        initializeColorscale();

        if(colorVector.valueVector != null)
           setDataHistogram();

        numberOfElements = getNumValues();
        return true;
    }

    /**
     * Initialize the property vector from the persisted property name
     */
    private void initPropertyVector()
    {
    	for(int i=0; i< propertyList.length; i++)
    	{
    		if(propertyList[i].getName().equalsIgnoreCase(propertyVectorName))
    		{
    			propertyVector = propertyList[i];
    			return;
    		}
    	}
    }

    /**
     * Initialize the color vector from the persisted color vector name
     */
    private void initColorVector()
    {
    	for(int i=0; i< colorByList.length; i++)
    	{
    		if(colorByList[i].getName().equalsIgnoreCase(colorVectorName))
    		{
    			colorVector = colorByList[i];
    			return;
    		}
    	}
    }
    /**
     * Initialize the sensor object
     **/
    public boolean initialize()
    {
    	getDisplayFields();
        initPropertyList();
        initColorByList();
        initFilterList();
        return true;
    }
    private void initFilterList()
    {
       StsToolkit.runWaitOnEventThread ( new Runnable() { public void run() { doInitFilterList(); } } );
    }
    /**
     * Initialize the list of attributes for size
     */
    private void doInitFilterList()
    {
        getDisplayFields();

        StsSensorFilter filter = new StsSensorFilter(false);
        filter.setName("None");

        StsSensorFilterClass filterClass = (StsSensorFilterClass)currentModel.getCreateStsClass("com.Sts.DBTypes.StsSensorFilter");
        filterList = (StsObject[])StsMath.arrayAddElement(filterList, filter);
        StsObject[] filters = filterClass.getAllFilters();
        filterList = (StsObject[])StsMath.arrayAddArray(filterList, filters);

        if(filterListBean != null) filterListBean.setListItems(filterList);
        if(sfilterListBean != null) sfilterListBean.setListItems(filterList);
    }

    private void initPropertyList()
    {
       StsToolkit.runWaitOnEventThread ( new Runnable() { public void run() { doInitPropertyList(); } } );
    }
    /**
     * Initialize the list of attributes for size
     */
    private void doInitPropertyList()
    {
        boolean isZ;
        getDisplayFields();
        propertyList = new StsTimeCurve[getNTimeCurves()+2];
        int cnt = 0;
        propertyList[cnt] = new StsTimeCurve(false);
        propertyList[cnt++].setName("None");
        if(getNTimeCurves() != 0)
        {
        	propertyList[cnt++] = getRelativeTimeCurve();
        }

        for(int i=0; i<getNTimeCurves(); i++)
        {
            isZ = false;
            StsTimeCurve timeCurve = (StsTimeCurve)timeCurves.getElement(i);
            if((timeCurve.getValueVector().getType() != StsLogVector.X) &&
               (timeCurve.getValueVector().getType() != StsLogVector.Y))
            {
                /*
                for(int j=0; j<StsSensorKeywordIO.Z_KEYWORDS.length; j++)
                {
                    if(timeCurve.getName().equalsIgnoreCase(StsSensorKeywordIO.Z_KEYWORDS[j]))
                        isZ = true;
                }
                if(isZ)
                    continue;
                    */
                propertyList[cnt++] = timeCurve;
            }
        }
        propertyList = (StsTimeCurve[])StsMath.trimArray(propertyList,cnt);
        if(propertyListBean != null) propertyListBean.setListItems(propertyList);    // ToDo: This causes all sensor objects to be reset to sizeBy "None" on load of a new sensor
    }

    /**
     * Get a list of attribute curves excluding positional curves.
     * @return curves with attribute and time values
     */
    public StsTimeCurve[] getPropertyCurvesExcludingXYZ()
    {
    	boolean isZ = false;
    	int cnt = 0;
    	StsTimeCurve[] pList = new StsTimeCurve[getNTimeCurves()];
        for(int i=0; i<getNTimeCurves(); i++)
        {
            isZ = false;
            StsTimeCurve timeCurve = (StsTimeCurve)timeCurves.getElement(i);
            if((timeCurve.getValueVector().getType() != StsLogVector.X) &&
               (timeCurve.getValueVector().getType() != StsLogVector.Y))
            {
                for(int j=0; j<StsSensorKeywordIO.Z_KEYWORDS.length; j++)
                {
                    if(timeCurve.getName().equalsIgnoreCase(StsSensorKeywordIO.Z_KEYWORDS[j]))
                        isZ = true;
                }
                if(isZ)
                    continue;
                pList[cnt++] = timeCurve;
            }
        }
        pList = (StsTimeCurve[])StsMath.trimArray(pList,cnt);
        return pList;
    }

    /**
     * Get a list of property curves.
     * @return - array of objects
     */
    public StsObject[] getPropertyCurves()
    {
        return timeCurves.getElements();
    }

    /**
     * Determine if the supplied curves already exist in this sensor
     * @param inCurves - list of curves
     * @return false if not found
     */
    public boolean hasCurves(StsTimeCurve[] inCurves)
    {
    	boolean found;
    	StsTimeCurve[] curves = getPropertyCurvesExcludingXYZ();
    	if(inCurves.length != curves.length)
    		return false;

    	for(int i=0; i<curves.length; i++)
    	{
    		found = false;
    		for(int j=0; j<inCurves.length; j++)
    		{
    			if(curves[i].getName().equalsIgnoreCase(inCurves[j].getName()))
    			{
    				found = true;
    				break;
    			}
    		}
    		if(!found)
    			return false;
    	}
    	return true;
    }

    /**
     * Generate a cumulative amplitude time curve
     * @return
     */
    public void computeAttributeCurves()
    {
        // Set Min and Max for Cumulative Amplitude
        if(timeCurves == null)
        	return;

        // Cumulative Amplitude
        String ampName = verifyAmplitude();
        if(ampName == null)
            return;

        attributeName = "CumulativeAmplitude";

        float[] ampFloats = getTimeCurve(ampName).getValueVector().getFloats();
        attribute = new float[ampFloats.length];
        attribute[0] = ampFloats[0];
        for(int i=1; i<attribute.length; i++)
            attribute[i] = attribute[i-1] + ampFloats[i];

        saveAttribute();

        // Total Events
        attributeName = "EventTotal";
        attribute = new float[ampFloats.length];
        attribute[0] = 1;
        for(int i=1; i<attribute.length; i++)
            attribute[i] = i+1;

        saveAttribute();
        return;
    }

    /**
     * Verify that an amplitude attribute exists in this sensor.
     * @return name of amplitude curve
     */
    public String verifyAmplitude()
    {
    	// Verify that the amplitude data is available if required.
        for(int i=0; i<AMPLITUDE_KEYWORDS.length; i++)
        {
            if(getTimeCurve(AMPLITUDE_KEYWORDS[i]) != null)
                return AMPLITUDE_KEYWORDS[i];
        }
    	return null;
    }

    /**
     * Generate a relative time curve
     * @return - relative time curve
     */
    public StsTimeCurve getRelativeTimeCurve()
    {
        StsTimeCurve tc = new StsTimeCurve(false);
        tc.setName("Relative Time");
        // Set Min and Max for Relative Time Vector
        if(timeCurves == null)
        	return null;
        long min = ((StsTimeCurve)timeCurves.getElement(0)).getTimeVector().getMinValue()/1000;
        long max = ((StsTimeCurve)timeCurves.getElement(0)).getTimeVector().getMaxValue()/1000;

        tc.setCurveMax((float)(max-min));
        tc.setCurveMin(0.0f);
        return tc;
    }

    private void initColorByList()
    {
       StsToolkit.runWaitOnEventThread ( new Runnable() { public void run() { doInitColorByList(); } } );
    }

    /** Initialize the color by list which will include monochrome colors and attributes. */
    private void doInitColorByList()
    {
        boolean isZ;
        getDisplayFields();
        String[] colorNames = StsColor.colorNames32;
        colorByList = new StsTimeCurve[getNTimeCurves()+2];
        int cnt = 0;
        colorByList[cnt] = new StsTimeCurve(false);
        colorByList[cnt++].setName("None");
        if(getNTimeCurves() != 0)
        	colorByList[cnt++] = getRelativeTimeCurve();
        for(int i=0; i<getNTimeCurves(); i++)
        {
            isZ = false;
            StsTimeCurve timeCurve = (StsTimeCurve)timeCurves.getElement(i);
            if((timeCurve.getValueVector().getType() != StsLogVector.X) &&
               (timeCurve.getValueVector().getType() != StsLogVector.Y))
            {
                /*
                for(int j=0; j<StsSensorKeywordIO.Z_KEYWORDS.length; j++)
                {
                    if(timeCurve.getName().equalsIgnoreCase(StsSensorKeywordIO.Z_KEYWORDS[j]))
                        isZ = true;
                }
                if(isZ)
                    continue;
                    */
                colorByList[cnt++] = timeCurve;
            }
        }
        colorByList = (StsTimeCurve[])StsMath.trimArray(colorByList,cnt);
        colorByList = (StsTimeCurve[])StsMath.arrayGrow(colorByList,colorNames.length);
        // Add all colors. Easier to add as timeCurves
        for(int i=0; i<colorNames.length; i++)
        {
            colorByList[cnt] = new StsTimeCurve(false);
            colorByList[cnt++].setName(colorNames[i]);
        }
        colorByList = (StsTimeCurve[])StsMath.trimArray(colorByList,cnt);
        if((colorByBean != null) && (colorByList != null))
        {
            colorByBean.setListItems(colorByList);
        }
    }
    /**
     * Get the filter by list.
     * @return - List of filters
     */
    public Object getFilterList()
    {
    	if(filterList == null)
    		initFilterList();
        return filterList;
    }
    /**
     * Get the color by list.
     * @return - List of color by curves
     */
    public Object getColorByList()
    {
    	if(colorByList == null)
    		initColorByList();
        return colorByList;
    }

    /**
     * Get the property (size by) list
     * @return List of property curves
     */
    public Object getPropertyList()
    {
    	if(propertyList == null)
    		initPropertyList();
        return propertyList;
    }

    /**
     * Get the currently selected property curve
     * @return curve with time and values
     */
    public StsTimeCurve getProperty()
    {
        return propertyVector;
    }

    /**
     * Update the properties list
     */
    public void updateProperties()
    {
         initPropertyList();
         propertyListBean.setListItems(propertyList);
         if(propertyVector != null) propertyListBean.setSelectedItem(propertyVector);
         initColorByList();
         colorByBean.setListItems(colorByList);
         if(colorVector != null) colorByBean.setSelectedItem(colorVector);
    }

    /**
     * Update the filters list
     */
    public void updateFilters()
    {
        if(filterListBean == null) return;
        filterListBean.removeAll();
        sfilterListBean.removeAll();
        filterList = new StsObject[0];
        initFilterList();
        if(StsMath.arrayContains(filterList, primaryFilter))
        {
            if(primaryFilter != null)
                setPrimaryFilter(primaryFilter);
        }
        else
            setPrimaryFilter(filterList[0]);

        if(StsMath.arrayContains(filterList, secondaryFilter))
        {
            if(secondaryFilter != null)
                setSecondaryFilter(secondaryFilter);
        }
        else
            setSecondaryFilter(filterList[0]);
    }

    /**
     * Find the selected property curve by vector
     * @param vector - The vector we are looking for
     * @return - time curve of matching curve
     */
    public StsTimeCurve findPropertyVectorByName(StsTimeCurve vector)
    {
        if(vector == null) return null;
    	if(propertyList == null) return null;
    	for(int i=0; i<propertyList.length; i++)
    	{
    		if(vector.getName().equals(propertyList[i].getName()))
    		{
    			vector = propertyList[i];
    			return vector;
    		}
    	}
    	return null;
    }
    /**
     * Find the selected color by curve by vector
     * @param vector - The vector we are looking for
     * @return - time curve of matching curve
     */
    public StsTimeCurve findColorVectorByName(StsTimeCurve vector)
    {
    	if(colorByList == null) return null;
        return findColorVectorByName(vector.getName());
    }
    public StsTimeCurve findColorVectorByName(String name)
    {
    	if(colorByList == null) return null;
    	for(int i=0; i<colorByList.length; i++)
    	{
    		if(name.equals(colorByList[i].getName()))
    		{
    			return colorByList[i];
    		}
    	}
    	return null;
    }

    /**
     * Set the current color by curve
     * @param vector - selected vector
     */
    public void setColorBy(StsTimeCurve vector)
    {
        vector = findColorVectorByName(vector);
        if((vector == null) || (colorVector == vector))
            return;
        if(vector.getName() == "None")
            return;

        colorVector = vector;
        colorVector.checkLoadVectors();        
        if(!initialized)
            return;

        resetColorMinMax();

        if(isPersistent())
        {
        	colorVectorName = colorVector.getName();
            dbFieldChanged("colorVectorName", colorVectorName);
        }
        resetColorscale();
        currentModel.viewObjectRepaint(this, this);
    }

    public void resetColorscale()
    {
        if(colorscale != null)
        {
            dataMin = colorVector.getCurveMin();
            dataMax = colorVector.getCurveMax();

            resetColorMinMax();

            colorscale.setRange((float)dataMin, (float)dataMax);
            if(getColorBy().valueVector != null)
                setDataHistogram();
        }
    }
    /**
     * Get the current color by curve
     * @return - curve with time and attribute values
     */
    public StsTimeCurve getColorBy()
    {
        return colorVector;
    }

    public float getColorMin() { return colorMin; }
    public float getColorMax() { return colorMax; }
    public float getScaleMin() { return scaleMin; }
    public float getScaleMax() { return scaleMax; }

    /**
     * Set the current filter to apply to the data
     *  * @param vector - selected filter
     */
    public void setPrimaryFilter(StsObject selection)
    {
        StsSensorFilter tFilter = (StsSensorFilter)selection;
        if(primaryFilter == selection) return;
        this.primaryFilter = tFilter;
        runFilters();
        currentModel.viewObjectChanged(this, this);
    }

    /**
     * Set the current filter to apply to the data
     *  * @param vector - selected filter
     */
    public void setSecondaryFilter(StsObject selection)
    {
        StsSensorFilter tFilter = (StsSensorFilter)selection;
        if(secondaryFilter == selection) return;
        this.secondaryFilter = tFilter;
        runFilters();
        currentModel.viewObjectChanged(this, this);        
    }

    /**
     * Get current filter
     */
    public StsObject getPrimaryFilter() { return primaryFilter; }
    /**
     * Set the current property (size by) curve
     * @param vector - selected vector
     */
    public void setProperty(StsTimeCurve vector)
    {
        vector = findPropertyVectorByName(vector);
        if((vector == null) || (propertyVector == vector))
            return;
        propertyVector = vector;
        propertyVector.checkLoadVectors();

        resetPropertyMinMax();
        
        propertyVectorName = vector.getName();
        dbFieldChanged("propertyVectorName", propertyVectorName);
        currentModel.viewObjectRepaint(this, this);
    }

    // Resets the min and max scaling to the data range in the vectors
    //
    public void resetPropertyMinMax()
    {
        if(propertyVector.getName().equalsIgnoreCase("None"))
            return;
        setScaleMin(propertyVector.getCurveMin());
        scaleMinBean.setValueObject(scaleMin);
        setScaleMax(propertyVector.getCurveMax());
        scaleMaxBean.setValueObject(scaleMax);
    }
    public void resetColorMinMax()
    {
        if(colorVector.getName().equalsIgnoreCase("None"))
            return;
        setColorMin(colorVector.getCurveMin());
        colorMinBean.setValueObject(colorMin);
        setColorMax(colorVector.getCurveMax());
        colorMaxBean.setValueObject(colorMax);
    }

    public void runFilters()
    {
        if(primaryFilter != null)
        {
            resetClusters();
            if(!primaryFilter.getName().equalsIgnoreCase("none"))
                ((StsSensorFilter)primaryFilter).filter(this);
            else
            {
                if(sfilterListBean != null) sfilterListBean.setSelectedIndex(0);
                secondaryFilter = filterList[0];
                return;
            }
            if(secondaryFilter != null)
            {
                if(!secondaryFilter.getName().equalsIgnoreCase("none"))
                {
                    ((StsSensorFilter)secondaryFilter).filter(this);
                }
            }
        }
        else
        {
            resetClusters();
            if(filterListBean != null) filterListBean.setSelectedIndex(0);
            if(sfilterListBean != null) sfilterListBean.setSelectedIndex(0);
        }
    }
    /**
     * Can this object be exported
     * @return true
     */
    public boolean canExport() { return true; }

    /**
     * Export the sensor of the sensor events in view.
     * @return success
     */
    public boolean export()
    {
        progressBarDialog = StsProgressBarDialog.constructor(currentModel.win3d, "Sensor Export", false,0,0,400,60);
        Runnable runExport = new Runnable()
        {
            public void run()
            {
                exportSensor();
            }
        };
        Thread exportThread = new Thread(runExport);

        exportThread.start();
        return true;
    }

    /**
     * Export the sensor events in view
     * @param file - output file
     * @param outputAtts - selected attributes to output
     * @return number of events exported
     */
    public int exportSensorsInView(File file, StsObject[] outputAtts)
    {
    	checkLoadVectors();

    	boolean hasHeader = false;
    	if(file.exists())
    		hasHeader = true;

        int nPoints = exportPoints.length;
        PrintWriter printWriter = null;
        try
        {
            // Output file.
            printWriter = new PrintWriter(new FileWriter(file.getAbsoluteFile(), true));
            if(!hasHeader)
            {
            	if(outputAtts == null)
                {
                    if(currentModel.getProject().velocityModel != null)
                        printWriter.println("Date, Time, X, Y, Z, zTime");
                    else
                        printWriter.println("Date, Time, X, Y, Z");
                }
            	else
            		printWriter.println(buildAsciiHeader(outputAtts, file.getName().substring(0, file.getName().indexOf(".csv"))));
            }

            // Output values.
            float x = 0.0f, y = 0.0f, z = 0.0f;
            double xOrigin = 0.0f, yOrigin = 0.0f;
            String line = "";
            for (int i = 0; i < nPoints; i++)
            {
            	line = currentModel.getProject().getTimeDateFormatForOutput().format(timeLongs[exportPoints[i]]);
                for(int j=0; j<getNTimeCurves(); j++)
                {
                    StsTimeCurve timeCurve = (StsTimeCurve)timeCurves.getElement(j);
                	if(timeCurve.getValueVector().getType() == StsLogVector.X)
                	{
                		x = timeCurve.getValuesVectorFloats()[exportPoints[i]];
                		xOrigin = timeCurve.getValueVector().getOrigin();
                		if(outputAtts != null)
                			line = line + ", " + String.valueOf(x+xOrigin);
                		continue;
                	}
                	else if(timeCurve.getValueVector().getType() == StsLogVector.Y)
                	{
                		y = timeCurve.getValuesVectorFloats()[exportPoints[i]];
                		yOrigin = timeCurve.getValueVector().getOrigin();
                		if(outputAtts != null)
                			line = line + ", " + String.valueOf(y+yOrigin);
                		continue;
                	}
                	else if(timeCurve.getValueVector().getType() == StsLogVector.DEPTH)
                	{
                		z = timeCurve.getValuesVectorFloats()[exportPoints[i]];
                		if(outputAtts != null)
                			line = line + ", " + String.valueOf(z);
                		continue;
                	}
                	else
                	{
                		if(outputAtts != null)
                		{
                			boolean output = false;
                			for(int ii=0; ii<outputAtts.length; ii++)
                			{
                				if(outputAtts[ii].getName().equalsIgnoreCase(timeCurve.getName()))
                				{
                					output = true;
                					break;
                				}
                			}
                			if(output)
                				line = line + ", " + String.valueOf(timeCurve.getValuesVectorFloats()[exportPoints[i]]);
                		}
                		continue;
                	}
                }
                if(outputAtts == null)
                {
                	line = line + "," + String.valueOf(x+xOrigin) + "," + String.valueOf(y+yOrigin) + "," + String.valueOf(z);
                }
                if(currentModel.getProject().velocityModel != null)
                    line = line + "," + (float) currentModel.getProject().velocityModel.getT(x, y, z, 0.0f);

                if(clusters != null)
                	line = line + "," + clusters[exportPoints[i]];
                printWriter.println(line);
            }
            return 0;
        }
        catch (Exception e)
        {
            StsException.outputException("StsSensor.export() failed.", e, StsException.WARNING);
            return 0;
        }
        finally
        {
            if (printWriter != null)
            {
                printWriter.flush();
                printWriter.close();
                return nPoints;
            }
        }
    }

    /**
     * Construct an ASCII header for an export file.
     * @param clusterName - Name
     * @return - the header
     */
    private String buildAsciiHeader(StsObject[] curves, String clusterName)
    {
    	String header = "Date, Time,";
        //Iterator iter = timeCurves.getIterator();
        //while(iter.hasNext())
        for(int i=0; i<curves.length; i++)
        {
        	StsTimeCurve curve = (StsTimeCurve)curves[i];
        	if(i<curves.length-1)
        		header = header + curve.getName() + ", ";
        	else
        		header = header + curve.getName();
        }
        if(currentModel.getProject().velocityModel != null)
               header = header + ", zTime";

        if(clusters != null)
        	header = header + "," + clusterName;

        /*
        if(clusters != null)
        {
        	String lbl = "Cluster";
        	int i = 1;
        	while(getTimeCurve(lbl) != null)
        		lbl = "Cluster" + i++;
    		header = header + "," + lbl;
    	}
    	*/
        return header;
    }

    /**
     * Export this sensor
     */
    public void exportSensor()
    {
        checkLoadVectors();
        int nPoints = timeLongs.length;
        PrintWriter printWriter = null;

        File file = StsParameterFile.getSaveFile("Select or Input Filename", currentModel.getProject().getRootDirString(), "csv", getName());
        if(file == null)
        {
        	StsMessageFiles.logMessage("No filename provided, sensor not exported.");
        	return;
        }

        progressBarDialog.setProgressMax(nPoints);
        progressBarDialog.setSize(300, progressBarDialog.getHeight());
        try
        {
            progressBarDialog.setLabelText("Exporting Sensor to File:    " + file.getName());
            progressBarDialog.pack();
            progressBarDialog.setVisible(true);

            // Output file.
            printWriter = new PrintWriter(new FileWriter(file.getAbsoluteFile(), false));
            String header = "Time,";
            for(int i=0; i<getNTimeCurves(); i++)
            {
                StsTimeCurve timeCurve = (StsTimeCurve)timeCurves.getElement(i);
                header = header + timeCurve.getName();
                if(i < (getNTimeCurves()-1))
                	header = header + ",";
            }
            if(currentModel.getProject().velocityModel != null)
               header = header + "zTime";
            printWriter.println(header);

            // Output values.
            float x = 0.0f, y = 0.0f, z = 0.0f;
            double xOrigin = 0.0f, yOrigin = 0.0f;
            for (int i = 0; i < nPoints; i++)
            {
            	String line = "";
            	line = currentModel.getProject().getTimeDateFormat().format(timeLongs[i]) + ",";

                for(int j=0; j<getNTimeCurves(); j++)
                {
                    StsTimeCurve timeCurve = (StsTimeCurve)timeCurves.getElement(j);
                    if(timeCurve.getValueVector().getType() == StsLogVector.X)
                    {
                    	x = timeCurve.getValuesVectorFloats()[i];
                    	xOrigin = timeCurve.getValueVector().getOrigin();
                    	line = line + String.valueOf(x+xOrigin);
                    }
                    else if(timeCurve.getValueVector().getType() == StsLogVector.Y)
                    {
                    	y = timeCurve.getValuesVectorFloats()[i];
                    	yOrigin = timeCurve.getValueVector().getOrigin();
                    	line = line + String.valueOf(y+yOrigin);
                    }
                    else if(timeCurve.getValueVector().getType() == StsLogVector.DEPTH)
                    {
                        z =  timeCurve.getValuesVectorFloats()[i];
                        line = line + String.valueOf(z);
                    }
                    else
                    	line = line + String.valueOf(timeCurve.getValuesVectorFloats()[i]);
                    if(j < (getNTimeCurves()-1))
                    	line = line + ",";
                }
                if(currentModel.getProject().velocityModel != null)
                    line = line + "," + (float) currentModel.getProject().velocityModel.getT(x, y, z, 0.0f);

                printWriter.println(line);
                progressBarDialog.setProgress(i);
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsSensor.export() failed.", e, StsException.WARNING);
            return;
        }
        finally
        {
            if (printWriter != null)
            {
                printWriter.flush();
                printWriter.close();
            }
            progressBarDialog.setProgress(nPoints);
        }
    }

    /**
     * Add this sensor to the current project
     * @param zDomain - The domain of this sensor
     * @return success
     */
    public boolean addToProject(byte zDomain)
    {
        this.zDomainOriginal = zDomain;
        return currentModel.getProject().addToProject(this, zDomainOriginal);
    }

    /**
     * Set the native domain
     * @param zDomain - The domain of this sensor
     * @return success
     */
    public void setOriginalDomain(byte zDomain)
    {
        this.zDomainOriginal = zDomain;
        return;
    }
    /**
     * Get the relative size scale
     * @return size scale
     */
    public int getNumberBins() { return numberBins; }

    /**
     * Set the best fit plane on or off
     * @param val
     */
    public void setFitPlane(boolean val)
    {
        bestFitPlane = val;
        dbFieldChanged("bestFitPlane", bestFitPlane);
        currentModel.viewObjectRepaint(this, this);
    }

    /**
     * Determine if the best fit plane option is selected
     * @return
     */
    public boolean getFitPlane()
    {
        return bestFitPlane;
    }
    
    /**
     * Set the relative size scale
     * @param bins - size scale
     */
    public void setNumberBins(int bins)
    {
        numberBins = bins;
        dbFieldChanged("numberBins", numberBins);
        currentModel.viewObjectRepaint(this, this);
    }

    /**
     * Initialize the colorscale
     */
    public void initializeColorscale()
    {
        if(colorscale == null)
        {
            spectrumName = getSensorClass().getDefaultSpectrumName();
            colorscale = new StsColorscale("Sensor", currentModel.getSpectrum(spectrumName),(float)dataMin, (float)dataMax);
        }
        colorscale.addActionListener(this);
    }

    /**
     * Initialize the data histogram
     */
    private void initializeHistogram()
    {
        int dataCnt[] = new int[255];
        int ttlSamples = 0;
        for(int i=0; i<dataCnt.length; i++)
            dataCnt[i] = 0;
        colorscaleBean.setHistogram(StsToolkit.calculateHistogram(dataCnt,255));
    }

    /**
     * Get the colorscale
     * @return colorscale
     */
    public StsColorscale getColorscale() { return colorscale; }

    /**
     * Set the colorscale to be applied to the colorby curve
     * @param colorscale - selected colorscale
     */
    public void setColorscale(StsColorscale colorscale)
	{
		this.colorscale = colorscale;
		colorscale.addActionListener(this);
	}

    /**
     * Compute the size of the event
     * @param min - minimum attribute value
     * @param max - maximum attribute value
     * @param value - value to size
     * @return size of event
     */
    public int computeSize(float min, float max, float value)
    {
        if(clipRange)
        {
        	if(value > max) value = max;
        	if(value < min) value = min;
        }
        else
        {
        	value = value * ((max-min)/(getProperty().getCurveMax() - getProperty().getCurveMin()));
        }
        if(!scaleLinear)
        {
            value = (float)Math.log10(value);
            min = (float)Math.log10(min);
            max = (float)Math.log10(max);
        }
        double	interval = (max-min)/(double)numberBins;
        int	bin = (int) ((value - min) / interval);
        int	size = bin + 1;
        return size;
    }

    /**
     * Get the start time for this sensor
     * @return minimum time for this sensor
     */
    public long getStartTime()
    {
    	return ((StsTimeCurve)timeCurves.getElement(0)).getTimeVector().getMinValue();
    }

    /**
     * Set whether this sensor has relative times
     * @param val
     */
    public void setIsRelative(boolean val) { isRelative = val; }

    /**
     * Find out if the sensor is relative time
     * @return
     */
    public boolean getIsRelative() { return isRelative; }

    /**
     * Set whether sensor times are complete with date
     * @param val - true (has date)
     */
    public void setHasDate(boolean val) { hasDate = val; }

    /**
     * Does this sensor have full times with dates
     * @return
     */
    public boolean getHasDate() { return hasDate; }

    /**
     * Set the start date when it is not supplied in the data
     * @param val - Date stirng
     */
    public void setStartDate(String val) { dateString = val; }

    /**
     * Get the start date which can be concatenated with relative times
     * @return start date
     */
    public String getStartDate() { return dateString; }

    /**
     * Turn data clipping on/off
     * @param clip
     */
    public void setClipRange(boolean clip) { this.clipRange = clip; dbFieldChanged("clipRange", clipRange);}

    /**
     * Determine if clipping is on or off
     * @return
     */
    public boolean getClipRange() { return clipRange; }

    /**
     * Get current data minimum. This will change as user selects properties.
     * @return data minimum
     */
    public double getDataMin() { return dataMin; }
     /**
     * Get current data maximum. This will change as user selects properties.
     * @return data maximum
     */
    public double getDataMax() { return dataMax; }

    public float getValueAt(StsLogVector vector, long time)
    {
        int idx = getTimeIndexGE(time);
        return vector.getValue(idx);
    }

    public int getTimeIndexGE(long time)
    {
        long timeF = (long)time;

        long[] timeArray = getTimeCurve(0).getTimeLongVector().getValues();
        return StsMath.arrayIndexAbove(timeF, timeArray);
    }

     /**
     * Get time minimum.
     * @return time minimum
     */
    public long getTimeMin() { return getTimeCurve(0).getTimeLongVector().getMinValue(); }
     /**
     * Get time maximum.
     * @return time maximum
     */
    public long getTimeMax() { return getTimeCurve(0).getTimeLongVector().getMaxValue(); }
    public void setTimeMax(long time) { getTimeCurve(0).getTimeLongVector().setMaxValue(time); }

    /**
     * Set the data minimum for the currently selected attribute
     * @param min - attribute minimum
     */
    public void setDataMin(double min)
    {
    	dataMin = min;
    }
    /**
     * Set the data maximum for the currently selected attribute
     * @param max - attribute maximum
     */
    public void setDataMax(double max)
    {
    	dataMax = max;
    }

    /**
     * Get the minimum and maximum Z or T dependent on which domain is in view
     * @return  minimum and maximum Z or T
     */
    public float[] getZTRange()
    {
    	StsTimeCurve depth = getTimeCurve(StsLogVector.DEPTH);
    	if(depth == null)
    		return null;

        if(zDomainOriginal != currentModel.getProject().getZDomain())
        {
            if(zDomainOriginal == StsProject.TD_DEPTH)
            {
            	return new float[] { currentModel.getProject().velocityModel.getSliceCoor(depth.getValueVector().getMinValue()),
            	                     currentModel.getProject().velocityModel.getSliceCoor(depth.getValueVector().getMaxValue()) };
            }
            else
            {
            	return new float[] { currentModel.getProject().velocityModel.getZCoor(depth.getValueVector().getMinValue()),
	                     currentModel.getProject().velocityModel.getZCoor(depth.getValueVector().getMaxValue()) };
            }
        }
        else
        	return new float[] {depth.getValueVector().getMinValue(), depth.getValueVector().getMaxValue()};
    }

    /**
     * Re-write the currently set attribute. This will overwrite the current binary file with new data. It uses
     * the AttVals and attributeName set in the setAttribute method
     * @return success
     */
    public boolean saveAttribute()
    {
    	if(attribute == null)
    		return false;

    	StsTimeVector timeVector = ((StsTimeCurve)timeCurves.getElement(0)).getTimeVector();

    	String binaryFilename = "sensor.bin." + getName() + "." + attributeName + ".0";
        StsLogVector curveVector = new StsLogVector(null, binaryFilename, attributeName, 0, -1);

        curveVector.setValues(attribute);
        curveVector.checkWriteBinaryFile(currentModel.getProject().getBinaryFullDirString(), false, true);
    	StsTimeCurve curve = StsTimeCurve.constructTimeCurve(timeVector, curveVector, 0);
    	addTimeCurve(curve);

        StsTimeCurve tCurve = getColorBy();
    	initColorByList();
        if(tCurve != null) setColorBy(tCurve);

        tCurve = getProperty();
    	initPropertyList();
        if(tCurve != null) setProperty(tCurve);

    	return true;
    }

    /**
     * Save the current cluster values as a new attribute.
     * @param clusterName - Provided name
     * @return success
     */
    public boolean saveClusters(String clusterName)
    {
    	if(clusters == null)
    		return false;

    	StsTimeVector timeVector = ((StsTimeCurve)timeCurves.getElement(0)).getTimeVector();

    	String binaryFilename = "sensor.bin." + getName() + "." + clusterName + ".0";
        StsLogVector curveVector = new StsLogVector(null, binaryFilename, clusterName, 0, -1);

        float[] fClusters = new float[clusters.length];
        for(int i=0; i<clusters.length; i++)
        	fClusters[i] = clusters[i];
        curveVector.setValues(fClusters);
        curveVector.checkWriteBinaryFile(currentModel.getProject().getBinaryFullDirString());

    	StsTimeCurve curve = StsTimeCurve.constructTimeCurve(timeVector, curveVector, 0);
    	addTimeCurve(curve);
    	initColorByList();
    	initPropertyList();
    	return true;
    }

    /**
     *  Zoom into the provided point.
      */
    public boolean goToRelativePoint(float[] pt, float scale, byte domain)
    {
    	float distance = numberBins * scale;
    	return goToPoint(pt, distance, domain);
    }

    /**
     * Go to the selected point and back off a specified distance
     * @param xyz - XYZ position
     * @param distance - distance from point
     * @return success
     */
    public boolean goToPoint(float[] xyz, float distance, byte domain)
    {
        Iterator viewIterator;
        StsWindowFamily mainWindowFamily = currentModel.getMainWindowFamily();
        if((domain == BOTH_D) || (domain == THREE_D))
        {
            viewIterator = mainWindowFamily.getWindowViewIteratorOfType(StsView3d.class);
            while(viewIterator.hasNext())
            {
                StsView3d view = (StsView3d)viewIterator.next();
                float[] parms = view.getCenterAndViewParameters();
                parms[0] = xyz[0];
                parms[1] = xyz[1];
                if(xyz.length > 2)
                    parms[2] = xyz[2];
                parms[3] = distance;
                view.changeModelView3d(parms);
                view.viewObjectRepaint(this, this);
            }
        }

        if((domain == BOTH_D) || (domain == TWO_D))
        {
            viewIterator = mainWindowFamily.getWindowViewIteratorOfType(StsView2d.class);
            while(viewIterator.hasNext())
            {
                StsView2d view = (StsView2d)viewIterator.next();
                // Center the view on xyz
                int dirNo = currentModel.getCursor3d().getCurrentDirNo();
                float[] pt2d = new float[2];
                switch(dirNo)
                {
                    case 1:                    // Y
                        view.reCenterOnPoint(new float[] {xyz[0],xyz[2]});
                        break;
                    case 0:                    // X
                        view.reCenterOnPoint(new float[] {xyz[1],xyz[2]});
                        break;
                    case 2:
                        view.reCenterOnPoint(new float[] {xyz[0],xyz[1]});
                        break;
                }
                view.viewObjectRepaint(this, this);
            }
        }
        return true;
    }

    /**
     *     StsMonitorable Interface
      */
    public int addNewData(StsObject object)
    {
    	StsMessageFiles.errorMessage("addNewData(StsObject) is not supportted for StsSensor type.");
    	return 0;
    }
    /**
     *     StsMonitorable Interface
      */
    public int addNewData(double[] attValues, long time, String[] attNames)
    {
    	StsMessageFiles.errorMessage("addNewData(double[], long, String[]) is not supportted for StsSensor type.");
    	return 0;
    }

    public void resetTimeCurves()
    {
         if(timeCurves != null)
         {
             for(int i=0; i<getNTimeCurves(); i++)
                ((StsTimeCurve)timeCurves.getElement(i)).resetVectors();
         }
    }
    /**
     * Add new data to this sensor from real-time monitoring
     * @param _source - source of the data
     * @param stype - source type
     * @param lastPollTime - last time data was polled
     * @return true when new data is found and successfully loaded
     */
    public int addNewData(String _source, byte stype, long lastPollTime, boolean compute, boolean reload, boolean replace)
    {
        int added = 0;
        final String source = _source;

    	// If sensor was updated (reloaded) since last poll
        if((lastModified == 0) && (timeCurves.getSize() != 0))
            lastModified = getTimeMax();

    	if(lastPollTime < lastModified)
    		lastPollTime = lastModified;

        // If the files are to be reloaded from scratch each time. Reset the vectors
        if(reload)
            resetTimeCurves();

        // Read the new data - must already be in attribute list for this sensor.
        if(stype == StsMonitor.FILE)  // It is a file
        {
        	added = processNewDataFromFile(source, reload, replace);
        }
        else       // It is a directory
        {
        	StsAbstractFile[] files = getSensorFiles(source);
        	for(int i=0; i<files.length; i++)
        	{
        		// Check if file lastModified time is less than last poll time
        		if(!StsToolkit.newData(files[i].getPathname(), lastPollTime) && !reload)
        			continue;
        		else
                {
        			// Process the file to pull out new records
               		added = added + processNewDataFromFile(files[i].getPathname(), reload, replace);
                }
        	}
        }
        if((lastPollTime == 1) && (added == 0))
        {
            if(isRealtimeLogging)
                StsMessageFiles.logMessage("No new data found in file....");

            return 0;
        }
        SimpleDateFormat format = new SimpleDateFormat(currentModel.getProject().getTimeDateFormatString());
        if(isRealtimeLogging)
           StsMessageFiles.logMessage("Setting last modified time to: " + format.format(new Date(lastModified)));
        lastModified = lastPollTime;

        // Update scaleing limits if required
        if(added > 0)
        {
            resetPropertyMinMax();
            resetColorMinMax();
            runFilters();
            currentModel.viewObjectChanged(this, this);
        }
        numberOfElements = this.getTimeCurve(0).getNumValues();
        return added;
    }

    /**
     * Add new data found in a monitored file
     * @param source - file name
     * @return new data added or not
     */
    public int processNewDataFromFile(String source, boolean reload, boolean replace)
    {
        boolean valuesReplaced = false;
    	int added = 0;
        int nCurves = 0;
        long[] timeVals = null;
        double[][] values = null;
        SimpleDateFormat format = new SimpleDateFormat(currentModel.getProject().getTimeDateFormatString());

    	try
    	{
            StsFile file = StsFile.constructor(source);
    		StsSensorFile sFile = new StsSensorFile(file);
            int nLines = sFile.getNumLinesInFile();

            // On relaod of the database the lastNumLines is always 0 so we will initialize it to the number of events in the file.
            //if((lastNumLines == 0) && (timeCurves.getSize() != 0))
            //{
            //    lastNumLines = getNumValues();  // May be more with duplicate times but certainly don't need to review entire file.
            //}


            if(lastNumLines == nLines)
            {
                currentModel.viewObjectChanged(this,this);
                return 0;
            }

            // Determine the number of lines to process.
		    int nLinesToProcess = nLines;
		    if(!reload)
			    nLinesToProcess = nLines - lastNumLines;
            if(nLinesToProcess > currentModel.getProject().getRealtimeBarAt())
            {
                processFile();
                progressPanel.progressBar.initializeImmediate(nLinesToProcess-1);
            }

    		if(!sFile.analyzeFile(currentModel, StsSensorFile.VARIABLE)) return 0;
    		//String startDate = getStartDate();

            if(timeCurves.getSize() != 0)
    		{
    			Date date = new Date(((StsTimeCurve)timeCurves.getElement(0)).getTimeVector().getMinValue());
    			//startDate = currentModel.project.getDateFormat().format(date);
    		}

    		BufferedReader bufRdr = new BufferedReader(new FileReader(source));
    		String line = sFile.readLine(bufRdr); // Header
            for(int ii=0; ii<lastNumLines; ii++)
                sFile.readLine(bufRdr); // Lines previously processed

            if(isRealtimeLogging)
                StsMessageFiles.logMessage("Processing lines " + lastNumLines + " to " + (nLines-1) + " from " + file.getFilename());

            String sTime = null;
    		// Loop in case multiple lines were added.
            for(int j=lastNumLines; j<nLines; j++)
            {
    		    line = sFile.readLine(bufRdr);
                sTime = format.format(new Date(sFile.currentTime));

                //progressPanel.incrementCount();
    			if(!sFile.getAttributeValues(currentModel, line, dateString))
                    continue;
    		    if(timeCurves.getSize() == 0)
                {
                    if(debugRealtime) System.out.println("Creating vectors: 1st time is:" + sTime);
                    if(isRealtimeLogging)
                        StsMessageFiles.logMessage("Creating vectors for sensor(" + getName() + ") and added first event at " + sTime);
    			    createNullVectors(sFile.curveNames, sFile);
                    checkAlarms(sFile);
                    continue;
                }
                else
                {
    			    // Verify the event is newer than the max event in sensor
                    if(debugRealtime) System.out.println("Verifying time:" + sTime);
    			    if(!verifyEvent(sFile))
                    {
                        if(replace)
                        {
                            if(replaceValue(sFile))
                            {
                                StsMessageFiles.logMessage("New data at " + sTime + " is in the past and does not exist in sensor. It will be ignored." );
                                valuesReplaced = true;
                                checkAlarms(sFile);
                            }
                            else
                            {
                                if(isRealtimeLogging)
                                    StsMessageFiles.logMessage("New data at " + sTime + " is in the past, old value(s) at this time have been replaced." );
                            }
                        }
                        if(progressDialog != null)
                            progressPanel.progressBar.setValueImmediate(j-lastNumLines+1);
    				    continue;
                    }
                }
                if(timeVals == null)
                {
                    nCurves = sFile.currentValues.length;
                    timeVals = null;
                    values = new double[nCurves][];
                    if(nCurves != 0)
    		        {
                        timeVals = new long[nLines];
                        for(int i=0; i<nCurves; i++)
                            values[i] = new double[nLines];
                    }
                }
                if(debugRealtime)
                    System.out.println("Adding value for time:" + sTime);

    			timeVals[added] = sFile.currentTime;
                for(int i=0; i<nCurves; i++)
                    values[i][added] = sFile.currentValues[i];

                checkAlarms(sFile);
                if(progressDialog != null)
                    progressPanel.progressBar.setValueImmediate(j-lastNumLines+1);
    			added++;
                if(added%250 == 0)
                {
                    StsMessageFiles.infoMessage("Processed " + added + " real-time events into sensor(" + getName() + ")");
                    //System.out.println("Processed " + added + " realtime events into sensor(" + getName() +")");
                }
    		}
            // Keep track of number of last read so they can be skipped on next cycle.
            if(!reload)
                lastNumLines = nLines;
            else
                lastNumLines = 0;
            // Trim the arrays
            timeVals = (long[])StsMath.trimArray(timeVals,added);
            for(int i=0; i<nCurves; i++)
                values[i] = (double[])StsMath.trimArray(values[i],added);

            if(values != null)
            {
                if(isRealtimeLogging)
                    StsMessageFiles.logMessage("Appending " + values.length + " values to sensor...");
                appendValues(sFile, timeVals, values, reload);
            }

            if(isRealtimeLogging)
                StsMessageFiles.logMessage("Updating timeseries, 2D and 3D plots... ");
            currentModel.viewObjectChanged(this, this);
            if(progressDialog != null)
                progressDialog.setVisible(false);
            return added;
    	}
    	catch(Exception e)
    	{
            StsMessageFiles.errorMessage("Failed to process monitored sensor data: " + source);
            return 0;
    	}
    }
    public boolean isReloadRequired()
    {
        return reloadRequired;
    }
    public void resetReload() { reloadRequired = false; }
    public boolean replaceValue(StsSensorFile sFile)
    {
        StsTimeCurve timeCurve = null;

          // Locate and replace value in appropriate curves
        for(int i=0; i<getNTimeCurves(); i++)
        {
            timeCurve = (StsTimeCurve)timeCurves.getElement(i);
            for(int j=0; j<sFile.validCurves.length; j++)
            {
            	if(timeCurve.getValueVector().getName().equalsIgnoreCase(sFile.validCurves[j]))
            	{
            		if(!timeCurve.replaceValueInCurve(sFile.currentTime, sFile.currentValues[j]))
                        return false;
            	}
            }
            timeCurve.getValueVector().getValues().setMinMax();
        }
        realtimePoint = true;
        reloadRequired = true;
        return true;
    }
    /**
     * Add new values to sensor curves from file
     * @param sFile - file with new sensor data to add
     * @return success
      */
    public boolean appendValues(StsSensorFile sFile, long[] times, double[][] values, boolean reload)
    {
    	// Add new data to sensor
    	boolean found = false;
    	StsTimeCurve timeCurve = null;
        //System.out.println("Appending " + sFile.currentTime + "to sensor");
        // Determine if an amplitude value exist in the realtime stream
        int ampIdx = -1;
        for(int i=0; i<AMPLITUDE_KEYWORDS.length; i++)
        {
            ampIdx = StsMath.stringArrayContainsIndex(sFile.curveNames,AMPLITUDE_KEYWORDS[i]);
            if(ampIdx != -1) break;
        }
        // Locate and add value to appropriate curves
        for(int i=0; i<getNTimeCurves(); i++)
        {
            found = false;
            timeCurve = (StsTimeCurve)timeCurves.getElement(i);
            for(int j=0; j<sFile.validCurves.length; j++)
            {
            	if(timeCurve.getValueVector().getName().equalsIgnoreCase(sFile.validCurves[j]))
            	{
            		timeCurve.addValuesToCurve(values[j]);
                    if(isRealtimeLogging)
                        StsMessageFiles.logMessage("Added " + values.length + " new values to sensor(" + getName() + ") attribute " + sFile.validCurves[j]);

            		found = true;
            		break;
            	}
            }
            // If it didn't find it, it might be a computed curve
            if((computeCurves) && (!found))
            {
                if(timeCurve.getValueVector().getName().equalsIgnoreCase("CumulativeAmplitude") && (ampIdx != -1))
                {
                    double[] cumAmp = new double[times.length];
                    if(reload)
                        cumAmp[0] = values[ampIdx][0];
                    else
                        cumAmp[0] = timeCurve.getMaxValue() + values[ampIdx][0];

                    for(int ii=1; ii<cumAmp.length; ii++)
                        cumAmp[ii] = cumAmp[ii-1] + values[ampIdx][ii];
                    timeCurve.addValuesToCurve(cumAmp);
                    if(isRealtimeLogging)
                        StsMessageFiles.logMessage("Added " + cumAmp.length + " new values to sensor(" + getName() + ") attribute CumulativeAmplitude");
                    found = true;
                }
                else if(timeCurve.getValueVector().getName().equalsIgnoreCase("EventTotal"))
                {
                    double[] ttlEvt = new double[times.length];
                    if(reload)
                        ttlEvt[0] = 1;
                    else
                        ttlEvt[0] = timeCurve.getMaxValue() + 1;
                    for(int ii=1; ii<ttlEvt.length; ii++)
                        ttlEvt[ii] = ttlEvt[ii-1] + 1;
                    timeCurve.addValuesToCurve(ttlEvt);
                    if(isRealtimeLogging)
                        StsMessageFiles.logMessage("Added " + ttlEvt.length + " new values to sensor(" + getName() + ") attribute EventTotal");
                    found = true;
                }
            }

            // Curve in sensor was not supplied in realtime stream so pad with zero.
            if(!found)
            {
                double[] zeroes = new double[times.length];
                for(int ii=0; ii<zeroes.length; ii++)
                    zeroes[ii] = 0.0f;
            	timeCurve.addValuesToCurve(zeroes);
            	StsMessageFiles.logMessage("Could not find data for attribute (" + timeCurve.getValueVector().getName() +
            			") in supplied realtime data. Setting values to 0.0");
            }
        }

		if(timeCurve != null)
        {
            StsTimeCurve curve = (StsTimeCurve)timeCurves.getElement(0);
            if(debugRealtime)
                System.out.println("Adding " + times.length + " values to time vector.... ");
			curve.addTimesToCurve(times);
            if(isRealtimeLogging)
                StsMessageFiles.logMessage("Added " + times.length + " new events to sensor(" + getName() + ")");

            // check if Monotonic
            int monotonic =  curve.getTimeVector().checkMonotonic();
            StsLogVector[] logVectors = new StsLogVector[timeCurves.getSize()];
            switch(monotonic)
            {
                case StsLongVector.MONOTONIC_NOT:
                    StsMessageFiles.errorMessage("Realtime data is not in time order, sorting....");
                    StsSensorKeywordIO.constructSortIndexes(curve.getTimeVector());
                    for(int i=0; i<timeCurves.getSize(); i++)
                        logVectors[i] = ((StsTimeCurve)timeCurves.getElement(i)).getValueVector();

                    StsSensorKeywordIO.sort(curve.getTimeVector(), logVectors);
                    curve.getTimeVector().getValues().setMinMax();
                    break;
                case StsLongVector.MONOTONIC_DECR:
                    StsMessageFiles.errorMessage("Realtime data is in reverse time order, sorting....");
                    for(int i=0; i<timeCurves.getSize(); i++)
                        logVectors[i] = ((StsTimeCurve)timeCurves.getElement(i)).getValueVector();

                    StsSensorKeywordIO.reverse(curve.getTimeVector(), logVectors);
                    curve.getTimeVector().getValues().setMinMax();
                    break;
            }
        }
        realtimePoint = true;
        if(getNumValues() < handEdits.length)   // Clear all hand edits
        {
            reloadRequired = true;
            handEdits = null;
            clearHandEdits();
        }
        else
        {
            byte[] newArray = new byte[getNumValues()];
            System.arraycopy(handEdits, 0, newArray, 0, handEdits.length);
            for(int i=handEdits.length; i<getNumValues(); i++)
                newArray[i] = 1;
            handEdits = newArray;
        }
        return true;
    }
    /**
     * Add new values to sensor curves from file
     * @param sFile - file with new sensor data to add
     * @return success
      */
    public boolean appendValues(StsSensorFile sFile)
    {
    	// Add new data to sensor
    	boolean found = false;
    	StsTimeCurve timeCurve = null;
        //System.out.println("Appending " + sFile.currentTime + "to sensor");
        // Determine if an amplitude value exist in the realtime stream
        int ampIdx = -1;
        for(int i=0; i<AMPLITUDE_KEYWORDS.length; i++)
        {
            ampIdx = StsMath.stringArrayContainsIndex(sFile.curveNames,AMPLITUDE_KEYWORDS[i]);
            if(ampIdx != -1) break;
        }
        // Locate and add value to appropriate curves
        for(int i=0; i<getNTimeCurves(); i++)
        {
            found = false;
            timeCurve = (StsTimeCurve)timeCurves.getElement(i);
            for(int j=0; j<sFile.validCurves.length; j++)
            {
            	if(timeCurve.getValueVector().getName().equalsIgnoreCase(sFile.validCurves[j]))
            	{
            		timeCurve.addValueToCurve(sFile.currentValues[j]);
            		found = true;
            		break;
            	}
            }
            // If it didn't find it, it might be a computed curve
            if((computeCurves) && (!found))
            {
                if(timeCurve.getValueVector().getName().equalsIgnoreCase("CumulativeAmplitude") && (ampIdx != -1))
                {
                    timeCurve.addValueToCurve(timeCurve.getMaxValue() + sFile.currentValues[ampIdx]);
                    found = true;
                }
                else if(timeCurve.getValueVector().getName().equalsIgnoreCase("EventTotal"))
                {
                    timeCurve.addValueToCurve(timeCurve.getMaxValue() + 1);
                    found = true;
                }
            }
            // Curve in sensor was not supplied in realtime stream so pad with zero.
            if(!found)
            {
            	timeCurve.addValueToCurve(0.0);
            	StsMessageFiles.logMessage("Could not find data for attribute (" + timeCurve.getValueVector().getName() +
            			") in supplied realtime data. Setting value to 0.0");
            }
        }
		if(timeCurve != null)
        {
            if(debugRealtime) System.out.println("Adding time value: " + sFile.currentTime);
			((StsTimeCurve)timeCurves.getElement(0)).addTimeToCurve(sFile.currentTime);
            if(isRealtimeLogging) StsMessageFiles.logMessage("Added new event to sensor(" + getName() + ") at time " + sFile.currentTime);
        }
        realtimePoint = true;
        //currentModel.viewObjectChanged(this, this);
    	return true;
    }

       /**
     *  Verify that the current events are not already in the set.
     * @param sFile - Sensor file description
     * @return successfully added data
      */
    public boolean verifyEvent(StsSensorFile sFile)
    {
        long maxTime = getTimeMax();
        //System.out.println("currentTime= " + sFile.currentTime + " maxTime= " + maxTime);
        if(sFile.currentTime <= maxTime)
        {
            //if(debugRealtime) System.out.println("Realtime event already in sensor for time stamp: " +
        	//		currentModel.project.getTimeDateFormat().format(new Date(sFile.currentTime)));
        	return false;
        }
        if(debugRealtime) System.out.println("Event at " + sFile.currentTime + " is new.");
    	return true;
    }


    /**
     * Get a list of all recognized sensor files within the supplied directory
     * @param sourceDir - Sensor directory
     * @return - list of sensor files
     */
    public StsAbstractFile[] getSensorFiles(String sourceDir)
    {
		String[] filterStrings = new String[]	{"csv", "txt", "Csv", "Txt", "CSV", "TXT"};
		StsFilenameEndingFilter filter = new StsFilenameEndingFilter(filterStrings);
    	StsFileSet fileSet = StsFileSet.constructor(sourceDir, filter);
    	return fileSet.getFiles();
    }

    /**
     * Set the picked event
     * @param index - index of picked event
     */
    public void setPicked(int index) { pickedIdx = index; }

    /**
     * Get the picked event
     * @return index of picked event
     */
    public int getPicked() { return pickedIdx; }

    /**
     * Add new data to sensor
     * @param point - point containing positional and attribute data
     * @param time - time of event
     * @param names - names of curves
     * @return success
     */
    public int addNewData(StsPoint point, long time, String[] names)
    {
    	StsMessageFiles.errorMessage("addNewData(StsPoint, long, String[]) is not supportted for StsSensor type.");
    	return 0;
    }

    /**
     * Get all new data added to the sensor
     * @param attribute - name of attribute
     * @param fromTime - get data from this time forward
     * @return attribute values for time range
     */
    public float[] getNewData(String attribute, long fromTime)
    {
    	if(fromTime >= this.getTimeMax())
    		return null;
    	else
    	{
            return getData(attribute, fromTime, getTimeMax());
    	}
    }

    public float[] getDataInRange(String attribute, int[] range)
    {
    	float[] values = new float[range[1]-range[0]+1];
    	float[] allValues = getTimeCurve(attribute).getValuesVectorFloats();
    	System.arraycopy(allValues, range[0], values, 0, range[1]-range[0]+1);
    	return values;
    }

    public float[] getData(String attribute, long fromTime, long toTime)
    {
    	int[] range = getTimeCurve(0).getTimeLongVector().getIndicesInValueRange(fromTime + 1, toTime + 1);
    	float[] values = new float[range[1]-range[0]+1];
    	float[] allValues = getTimeCurve(attribute).getValuesVectorFloats();
    	System.arraycopy(allValues, range[0], values, 0, range[1]-range[0]+1);
    	return values;
    }
    /**
     * Get the number of points within the latest interval
     * @param seconds - number of seconds to compute number of points.
     * @return number of points within the last (interval) seconds.
     */
    public int getNumPoints(int seconds, long endTime)
    {
        long startTime = endTime - (long)(seconds * 1000L);
        StsTimeCurve tCurve = getTimeCurve(0);
        if(tCurve == null)
            return 0;
        int[] range = tCurve.getTimeLongVector().getIndicesInValueRange(startTime, endTime);
        if(range == null)
            return 0;
        return range[1]-range[0]+1;
    }

    /**
     * Get the number of points within the latest interval
     * @param seconds - number of seconds to compute number of points.
     * @return number of points within the last (interval) seconds.
     */
    public int getNumPoints(int seconds)
    {
        return getNumPoints(seconds, getTimeMax());
    }

    /**
     * Get all new time data added to the sensor
     * @param fromTime - get data from this time forward
     * @return time values for time range
     */
    public long[] getNewTimeData(long fromTime, long toTime)
    {
        int[] range = getTimeCurve(0).getTimeLongVector().getIndicesInValueRange(fromTime + 1, toTime + 1);
        if(range == null)
        {
            if(debugRealtime) System.out.println("No new time data found.");
            return null;
        }
        if(debugRealtime) System.out.println("New time data index range: " + range[0] + " - " + range[1]);

    	long[] values = new long[range[1]-range[0]+1];
        long[] allValues = getTimeCurve(0).getTimeLongVector().getValues();
        System.arraycopy(allValues, range[0], values, 0, range[1]-range[0]+1);
        return values;
    }

    public long[] getNewTimeData(long fromTime)
    {
    	if(fromTime >= this.getTimeMax())
    		return null;
    	else
            return getNewTimeData(fromTime, getTimeMax());
    }

    public String toString(int index)
    {
        if(this instanceof StsStaticSensor)
            return ((StsStaticSensor)this).toString(index);
        else
            return ((StsDynamicSensor)this).toString(index);
    }

    public void processFile()
    {
        progressDialog = new JDialog(currentModel.win3d, "Loading Realtime..." + getName(), false);
        progressDialog.setPreferredSize(new Dimension(300,75));
        progressDialog.add(progressPanel);
        progressDialog.pack();
        progressDialog.setLocationRelativeTo(currentModel.win3d);
        progressDialog.setVisible(true);
        return;
    }

    public StsAlarm[] getAlarms() { return alarms; }
    public boolean addAlarm(StsAlarm alarm)
    {
        if(alarms == null)
        {
            alarms = new StsAlarm[1];
            alarms[0] = alarm;
        }
        else
            alarms = (StsAlarm[])StsMath.arrayAddElement(alarms, alarm);

        return true;
    }
    
    public void checkAlarms(StsSensorFile sFile)
    {
        if(!hasAlarms()) return;
        double[] xyz = new double[4];
        for(int i=0; i<sFile.currentValues.length; i++)
        {
             if(sFile.curveNames[i] == StsSensorKeywordIO.DEPTH)
                   xyz[2] = sFile.currentValues[i];
             else if(sFile.curveNames[i] == StsSensorKeywordIO.X)
                   xyz[0] = sFile.currentValues[i];
             else if(sFile.curveNames[i] == StsSensorKeywordIO.Y)
                   xyz[1] = sFile.currentValues[i];
             else
                   xyz[3] = sFile.currentValues[i];
        }
        for(int i=0; i<alarms.length; i++)
            alarms[i].checkAlarm(xyz);
    }

    public boolean hasAlarms()
    {
        if(alarms != null)
            return true;
        else
            return false;
    }

    public void mouseSelectedEdit(StsMouse mouse)
    {
        StsPickItem items = StsJOGLPick.pickItems[0];
        if (items.names.length < 2) return;
        int sensorIdx = items.names[1];
        setPicked(sensorIdx);
        logMessage(toString(sensorIdx));
        //getDynamicSensorClass().displayGather(this, sensorIdx, currentModel.win3d);
    }

    public boolean handDelete(StsMouse mouse)
    {
        mouseSelectedEdit(mouse);
        int index = getPicked();
        handEdits[index] = 0;
        dbFieldChanged("handEdits", handEdits);        
        return true;
    }
    public boolean handClear(StsMouse mouse)
    {
        mouseSelectedEdit(mouse);
        int index = getPicked();
        handEdits[index] = 1;
        dbFieldChanged("handEdits", handEdits);
        return true;
    }
    public boolean handClearAll()
    {
        handEdits = null;
        clearHandEdits();
        dbFieldChanged("handEdits", handEdits);
        return true;
    }
}
