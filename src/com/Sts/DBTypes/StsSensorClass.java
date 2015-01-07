package com.Sts.DBTypes;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.SensorLoad.*;
import com.Sts.DB.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.UI.Sounds.StsSound;
import com.Sts.Utilities.*;

import java.io.*;
import java.text.*;
import java.util.*;

public class StsSensorClass extends StsClass implements StsSerializable, StsClassCursorDisplayable
{
    boolean displayNames = false;
    boolean enableTime = true;
    boolean enableSound = false;
    boolean enableAccent = true;
    boolean enableGoTo = false;
    boolean showDeleted = false;
    int goToOffset = 10;
    private StsColor accentColor = new StsColor(StsColor.GREY);

    float zCriteria = -1;
    int accentDuration = 10;

    static public final String EVENTS_STRING = "Microseismic";
    static public final String PUMPS_STRING = "Pump Curves";
    static public final String TRACERS_STRING = "Tracers";
    static public final String GPS_STRING = "GPS";
    static public final String OTHER_STRING = "Other";
    public static String[] typeStrings = {EVENTS_STRING, PUMPS_STRING, TRACERS_STRING, GPS_STRING, OTHER_STRING};

    static public final byte EVENTS = 0;
    static public final byte PUMPS = 1;
    static public final byte TRACERS = 2;
    static public final byte GPS = 3;
    static public final byte OTHER = 4;

    transient StsSubType[] subTypes = null;
    transient StsWin3dBase gatherWindow = null;

    /** X coordinate size that indicates relative */
    public double xRelativeCriteria = 100000.0f;
    /** Number of seconds between stages */
    public long multiStageCriteria = 9000l;   // 2.5 hours
    /** 2D Scale factor */
    protected float scale2D = 1.0f;

   StsSensor valueSensor = null;
   int valueIndex = 0;
   float valueMax;
   float valueMin;
   
   private String defaultSpectrumName = StsSpectrumClass.SPECTRUM_RWB;
   private StsColor defaultSensorColor = new StsColor(StsColor.BLUE);
   private String defaultSound = StsSound.CLANG;
   private int defaultSize = 100;
   protected byte defaultPointType = StsSensor.SPHERE;
   protected long[] lastModified = null;
   protected int[] lastNumLines = null;

    public StsSensorProperties defaultSensorProperties = null;

   public StsSensorClass()
   {
       userName = "Sensor Data";
   }

   public void initializeFields()
   {
       displayFields = new StsFieldBean[]
       {
           new StsIntFieldBean(this, "accentDuration", 1, 100, "Accent Duration (frames):", true),
           new StsDoubleFieldBean(this, "xRelativeCriteria", "Relative Coor Criteria:"),
           new StsLongFieldBean(this, "multiStageCriteria", true, "Multi-Stage Criteria(sec):"),
       };

       defaultFields = new StsFieldBean[]
           {
                   new StsIntFieldBean(this, "defaultSize", 1, 500, "Size:", true),
                   new StsComboBoxFieldBean(this, "defaultSymbolString", "Symbol:", StsSensor.SYMBOL_TYPE_STRINGS),
                   new StsComboBoxFieldBean(this, "defaultSound", "Sound:", StsSound.sounds),
                   new StsColorComboBoxFieldBean(this, "defaultSensorColor", "Color:", StsColor.colors32)
           };
   }
   
    public boolean projectInitialize(StsModel model)
    {
        initializeSubClasses();
        defaultSensorProperties = new StsSensorProperties(model, this, "defaultSensorProperties");
        return true;
    }

    public StsSubType[] getSubTypes()
    {
        return subTypes;
    }

    static public byte getTypeForName(String typename)
    {
        for (byte n = 0; n < typeStrings.length; n++)
            if (typeStrings[n] == typename) return n;
        return (byte) 0;
    }

    private void initializeSubClasses()
    {
        subTypes = new StsSubType[typeStrings.length];
        for (int n = 0; n < typeStrings.length; n++)
            subTypes[n] = new StsSubType(this, typeStrings[n], (byte) n);
    }

    public void displayTimeClass(StsGLPanel3d glPanel3d, long time)
    {
        StsObject[] sensors = getDynamicSensors();
		if(sensors.length > 0)
		{
			for(int i=0; i<sensors.length; i++)
			{
				StsDynamicSensor sensor = (StsDynamicSensor)sensors[i];
				if ((enableTime && sensor.isAlive(time)) || (!enableTime))
                {
					sensor.display(glPanel3d, displayNames);
                    if(sensor.getFitPlane())
                        sensor.displayPolygons(glPanel3d);
                }
            }
			if(StsDynamicSensor.runTimers)
			{
				StsSumTimer.printTimers("DynamicSensor");
				StsSumTimer.clear();
			}
		}
        sensors = getStaticSensors();
        for(int i=0; i<sensors.length; i++)
        {
            StsStaticSensor sensor = (StsStaticSensor)sensors[i];
            if ((enableTime && sensor.isAlive(time)) || (!enableTime))
                sensor.display(glPanel3d, displayNames);
        }
    }

    public StsObject[] getSensors()
    {
        StsObject[] dynamicSensors = ((StsDynamicSensorClass)currentModel.getStsClass("com.Sts.DBTypes.StsDynamicSensor")).getSensors();
        StsObject[] staticSensors = ((StsStaticSensorClass)currentModel.getStsClass("com.Sts.DBTypes.StsStaticSensor")).getSensors();
        return (StsObject[])StsMath.arrayAddArray(dynamicSensors, staticSensors);
    }

    public StsObject[] getStaticSensors()
    {
        return ((StsStaticSensorClass)currentModel.getStsClass("com.Sts.DBTypes.StsStaticSensor")).getSensors();
    }

    public StsObject[] getDynamicSensors()
    {
        return ((StsDynamicSensorClass)currentModel.getStsClass("com.Sts.DBTypes.StsDynamicSensor")).getSensors();
    }

    public void drawOnCursor3d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, StsPoint[] planePoints, boolean isDragging) { }

    public void drawOnCursor2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, boolean axesFlipped, boolean xAxisReversed, boolean yAxisReversed)
    {
        long time = currentModel.getProject().getProjectTime();
        StsObject[] sensors = getDynamicSensors();
        for(int i=0; i<sensors.length; i++)
        {
            StsDynamicSensor sensor = (StsDynamicSensor) sensors[i];
            if ((enableTime && sensor.isAlive(time)) || (!enableTime))
                sensor.display2d(glPanel3d, dirNo, dirCoordinate, axesFlipped, xAxisReversed, yAxisReversed);
        }
        sensors = getStaticSensors();
        for(int i=0; i<sensors.length; i++)
        {
            StsStaticSensor sensor = (StsStaticSensor) sensors[i];
            if ((enableTime && sensor.isAlive(time)) || (!enableTime))
                sensor.display2d(glPanel3d, dirNo, dirCoordinate, axesFlipped, xAxisReversed, yAxisReversed);
        }
    }

    public void clearHighlights()
    {
        Iterator iter = getVisibleObjectIterator();
        while (iter.hasNext())
        {
            StsSensor sensor = (StsSensor) iter.next();
            sensor.clearHiglightedPoints();
            currentModel.viewObjectRepaint(this, sensor);
        }
    }

    public void drawOn3dCurtain(StsGLPanel3d glPanel3d, StsSeismicCurtain seismicCurtain) { }

    public void setDisplayNames(boolean displayNames)
    {
        if(this.displayNames == displayNames) return;
        this.displayNames = displayNames;
        currentModel.win3dDisplayAll();
    }

    public boolean getDisplayNames() {	return displayNames; }

    public void setValueRangeAndSensor(StsSensor sensor, int valueIdx, float min, float max)
    {
        valueSensor = sensor;
        valueIndex = valueIdx;
        valueMax = max;
        valueMin = min;
        currentModel.win3dDisplayAll();
    }

    public StsSensor getValueSensor() { return valueSensor; }
    public float getValueMin() { return valueMin; }

    public float getValueMax() { return valueMax; }

    public int getValueIndex() { return valueIndex; }

    public void setEnableTime(boolean enable)
    {
        if (this.enableTime == enable) return;
        this.enableTime = enable;
        currentModel.win3dDisplayAll();
    }

    public boolean getEnableTime() { return enableTime; }

    public void setAccentDuration(int duration)
    {
        if (this.accentDuration == duration) return;
        this.accentDuration = duration;
        currentModel.win3dDisplayAll();
    }
    public int getAccentDuration() { return accentDuration; }
    public void setShowDeleted(boolean enable)
    {
        if (this.showDeleted == enable) return;
        this.showDeleted = enable;
        currentModel.win3dDisplayAll();
    }

    public boolean getShowDeleted() { return showDeleted; }
    public void setEnableSound(boolean enable)
    {
        if (this.enableSound == enable) return;
        this.enableSound = enable;
        currentModel.win3dDisplayAll();
    }

    public boolean getEnableSound() { return enableSound; }

    public void setEnableAccent(boolean enable)
    {
        if (this.enableAccent == enable) return;
        this.enableAccent = enable;
        currentModel.win3dDisplayAll();
    }

    public boolean getEnableAccent() { return enableAccent; }

    public void setEnableGoTo(boolean enable)
    {
        if (this.enableGoTo == enable) return;
        this.enableGoTo = enable;
        currentModel.win3dDisplayAll();
    }

    public boolean getEnableGoTo() { return enableGoTo; }

    public void setGoToOffset(int mult)
    {
        if (this.goToOffset == mult) return;
        this.goToOffset = mult;
       currentModel.win3dDisplayAll();
   }
   public int getGoToOffset() {	return goToOffset; }

   public void setScale2D(float scale)
   {
       if (this.scale2D == scale) return;
       this.scale2D = scale;
       currentModel.win3dDisplayAll();
   }
   public float getScale2D() { return scale2D; }

   public void setAccentColor(StsColor color)
   {
       if(accentColor.equals(color)) return;
       accentColor = new StsColor(color);
   }
   public StsColor getAccentColor()  { return accentColor; }

   public String getDefaultSpectrumName()  { return defaultSpectrumName; }
   public void setDefaultSpectrumName(String value)
   {
       if(this.defaultSpectrumName.equals(value)) return;
       this.defaultSpectrumName = value;
   }
   public void setDefaultSensorColor(StsColor color)
   {
       if(defaultSensorColor.equals(color)) return;
       defaultSensorColor = new StsColor(color);
   }
   public StsColor getDefaultSensorColor()  { return defaultSensorColor; }


    public StsColor getDefaultSensorStsColor() { return defaultSensorColor; }
    public void setDefaultSound(String sound)
    {
        if (this.defaultSound.equalsIgnoreCase(sound)) return;
        StsSound.play(sound);
        this.defaultSound = sound;
    }

    public String getDefaultSound() { return defaultSound; }

    public void setDefaultSize(int size)
    {
        if (this.defaultSize == size) return;
        this.defaultSize = size;
    }

    public int getDefaultSize() { return defaultSize; }

    public String getDefaultSymbolString() { return StsSensor.SYMBOL_TYPE_STRINGS[defaultPointType]; }

    public void setDefaultSymbolString(String symbolString)
    {
        for (int i = 0; i < StsSensor.SYMBOL_TYPE_STRINGS.length; i++)
        {
            if (symbolString.equals(StsSensor.SYMBOL_TYPE_STRINGS[i]))
            {
                defaultPointType = StsSensor.SYMBOL_TYPES[i];
                break;
            }
        }
        return;
    }

    public byte getDefaultPointType() { return defaultPointType; }

    public double getXRelativeCriteria() { return xRelativeCriteria; }

    public void setXRelativeCriteria(double value) { xRelativeCriteria = value; }

    public long getMultiStageCriteria() { return multiStageCriteria; }

    public void setMultiStageCriteria(long value) { multiStageCriteria = value; }

    public boolean canExportView() { return true; }

    public int exportView()
    {
    	String filename = "SensorsInView" + System.currentTimeMillis() + ".csv";    	
    	return exportView(filename);
    }            
    public StsObject[] getCommonAttributes()
    {
        boolean foundIt = false;

        Iterator iter = getVisibleObjectIterator();
        StsSensor sensor = (StsSensor) iter.next();

        StsObject[] returnCurves = sensor.getPropertyCurves();
        StsObject[] commonCurves = sensor.getPropertyCurves();
        for (int j = 0; j < commonCurves.length; j++)
        {
            iter = getVisibleObjectIterator();
            while (iter.hasNext())
            {
                StsSensor nextSensor = (StsSensor) iter.next();
                StsObject[] curves = nextSensor.getPropertyCurves();
                for (int k = 0; k < curves.length; k++)
                {
                    if (curves[k].getName().equalsIgnoreCase(commonCurves[j].getName()))
                    {
                        foundIt = true;
                        break;
                    }
                }
                if (foundIt)
                {
                    foundIt = false;
                    continue;
                }
                else
                {
                    returnCurves = (StsObject[]) StsMath.arrayDeleteElement(returnCurves, commonCurves[j]);
                    foundIt = false;                                
                    break;
                }
            }
        }
        return returnCurves;
    }

    public int exportView(String filename)
    {
        return exportView(filename, StsProgressPanel.constructorWithCancelButton(), StsSensorClass.OTHER, false);
    }

    public boolean saveClusters(String attName)
    {
        StsSensor sensor = null;    	
    	Iterator iter = getVisibleObjectIterator();
        while(iter.hasNext())
        {
        	sensor = (StsSensor)iter.next();
        	sensor.saveClusters(attName);
        }
        new StsMessage(currentModel.win3d, StsMessage.INFO, "Added a new attribute (" + attName + ") to clustered sensors.");  
    	    	
    	return true;
    }
    
    public boolean saveAttribute()
    {
        StsSensor sensor = null;    	
    	Iterator iter = getVisibleObjectIterator();
        while(iter.hasNext())
        {
        	sensor = (StsSensor)iter.next();
        	sensor.saveAttribute();
        }
        new StsMessage(currentModel.win3d, StsMessage.INFO, "Added a new attribute to clustered sensors.");  
    	    	
    	return true;
    }    
    
    public int exportView(String filename, byte type)
    {
        return exportView(filename, StsProgressPanel.constructorWithCancelButton(), type, false);
    }

    public int exportView(String filename, StsProgressPanel panel, byte type, boolean messageOut)
    {
        StsSensor sensor = null;

        File file = StsParameterFile.getSaveFile("Select or Input Filename", currentModel.getProject().getRootDirString(), "csv", filename);
        if (file == null)
        {
            StsMessageFiles.logMessage("No filename provided, sensor not exported.");
            return 0;
        }
        else
        {
            // confirm the csv suffix is present
            if (!file.getName().endsWith(".csv"))
                file = new File(file.getAbsolutePath() + ".csv");
        }
        try {file.delete();}
        catch(Exception ex) { new StsMessage(currentModel.win3d, StsMessage.ERROR, "Unable to create file, check permissions."); return 0; }

        filename = file.getName();

        int nPoints = 0;
        StsObject[] curves = (StsObject[]) getCommonAttributes();
        Iterator iter = getVisibleObjectIterator();
        while (iter.hasNext())
        {
            sensor = (StsSensor) iter.next();
            nPoints += sensor.exportSensorsInView(file, curves);
        }
        if (messageOut)
            new StsMessage(currentModel.win3d, StsMessage.INFO, "Exported " + nPoints + " to " + filename);

        // Export the values
        if (StsYesNoDialog.questionValue(currentModel.win3d, "Do you want to create a sensor with the exported data?"))
        {
            StsAbstractSensorFactory sensorFactory = new StsSensorFactory();
            StsFile stsFile = StsFile.constructor(file.getPath());
            StsSensorFile sFile = new StsSensorFile(stsFile);
            sFile.analyzeFile(currentModel, StsSensorFile.VARIABLE);

            Date date = new Date(0L);
            SimpleDateFormat tDformat = currentModel.getProject().getTimeDateFormat();
            String time = tDformat.format(date);
            sFile.numStages(currentModel);

            StsSensorKeywordIO.parseAsciiFilename(sFile.file.getFilename());
            String[] fileNames = new String[]{sFile.file.getFilename()};
            StsSensorImport.addSensorFilenameSets(fileNames, StsSensorImport.ASCIIFILES);
            try { StsSensorImport.setCurrentDirectory(StsFile.getDirectoryFromPathname(file.getPath()));  }
            catch(StsException ex) { System.err.println("File access error: " + ex); }
            byte vUnits = currentModel.getProject().getDepthUnits();    // Need to implement
            byte hUnits = currentModel.getProject().getXyUnits();       // Need to implement
            StsSensorKeywordIO.setTimeType(StsSensorKeywordIO.TIME_OR_DATE);
            StsSensor[] sensors = StsSensorImport.createSensors(currentModel, panel, sFile, false, vUnits, hUnits, sensorFactory, false);                                   
            if (sensors == null)
            {
                StsMessageFiles.logMessage("Failed to create sensor from exported data.");
                return nPoints;
            }
            if (sensors[0] == null)
            {
                StsMessageFiles.logMessage("Failed to create sensor from exported data.");
                return nPoints;
            }

            sensors[0].setType(type);
        	sensors[0].setStsColor(sensors[0].getSensorClass().getDefaultSensorColor());
        	panel.setValue(100);
        	panel.setDescription("Created sensor " + name);       	
    	}
        return nPoints;
    }

    public void close()
    {
        int nValues = getSize();
        lastModified = new long[nValues];
        lastNumLines = new int[nValues];
        for(int n = 0; n < nValues; n++)
        {
            StsSensor sensor = (StsSensor)getElement(n);
            if(sensor != null)
            {
                lastModified[n] = sensor.lastModified;
                lastNumLines[n] = sensor.lastNumLines;
            }
        }
    }

    public void finalInitialize()
    {
        if(lastModified == null) return;
        int length = lastModified.length;
        int nValues = getSize();
        if(length != nValues)
            StsException.systemError(this, "finalInitialize", "Saved lastModified.length " + length + " doesn't agree with number of StsSensor instances " + nValues);

        nValues = Math.min(nValues, length);
        for(int n = 0; n < nValues; n++)
        {
            StsSensor sensor = (StsSensor)getElement(n);
            if(sensor != null)
            {
                sensor.lastModified = lastModified[n];
                sensor.lastNumLines = lastNumLines[n];
            }
        }

        /* Initialize Colorscale realtive to data
        StsObject[] sensors = getDynamicSensors();
        for(int n = 0; n < sensors.length; n++)
        {
            StsSensor sensor = (StsSensor)sensors[n];
            sensor.initializeColorscale();
        }
              */
    }


    public void initializePropertyFields()
    {
        propertyFields = new StsFieldBean[]
        {
            new StsButtonFieldBean("Sensor Properties", "Edit Sensor Properties.", this, "displaySensorProperties")
        };
    }

   public void displaySensorProperties()
   {
       StsDialogFace[] dialogs = {defaultSensorProperties};
       new StsOkApplyCancelDialog(currentModel.win3d, dialogs , "Edit Sensor Properties", false);
   }

    public StsColor getDefaultColor(int index)
    {
        if(defaultSensorProperties != null)
            return defaultSensorProperties.getColor(index);
        else
            return StsColor.colors32[index];
    }

    public String getDefaultColorName(int index)
    {
        if(defaultSensorProperties != null)
            return defaultSensorProperties.getColorName(index);
        else
            return StsColor.colorNames32[index];
    }
}