//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Interfaces.StsTreeObjectI;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.Utilities.Shaders.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.io.*;
import java.util.Date;

public class StsStaticSensor extends StsSensor implements StsTreeObjectI
{
    /** The X location for a static sensor */
    protected double xLoc = 0.0;
    /** The Y location for a static sensor */
    protected double yLoc = 0.0;
    /** The Z location for a static sensor */
    protected double zLoc = 0.0;
    /** local coordinate location of static sensor */
    transient StsPoint point;

    /** default constructor */
    public StsStaticSensor()
    {
        super();
    }

    public StsStaticSensor(StsWell well, String name)
    {
        this(well, name, 0.0, 0.0, 0.0);
    }

    public StsStaticSensor(StsWell well, String name, double x, double y, double z)
    {
        super(well, name);
        xLoc = x;
        yLoc = y;
        zLoc = z;
        displayType = SHOW_SINGLE;
    }

    static public StsStaticSensor nullSensorConstructor(String name)
    {
        return new StsStaticSensor(null, name);
    }

    /**
     * Compute the points that are in view in preparation for export.
     */
    public void computePointsInView()
    {
    }

    /**
     * Display the sensor in 2D view
     * @param glPanel3d - Graphics handle
     * @param dirNo - X, Y or Z slice
     * @param dirCoordinate - Current slice position
     * @param axesFlipped - Is the axis flipped
     * @param xAxisReversed - Is the x axis reversed
     * @param yAxisReversed - Is the y axis reversed
     */
    public void display2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, boolean axesFlipped, boolean xAxisReversed, boolean yAxisReversed)
    {
    	StsColor aColor = getSensorClass().getAccentColor();
        if(!verifyDisplayable())
        	return;
        if(getSensorClass().getEnableTime())
        {
            if(!isAlive(currentModel.getProject().getProjectTime()))
                return;
        }
        long duration = this.displayDuration;
        byte type = this.displayType;

        if(colorscale == null) initializeColorscale();
        if(!checkLoadVectors())
        	return;
        int points = 0;

        GL gl = glPanel3d.getGL();
        glPanel3d.setViewShift(gl, StsGraphicParameters.edgeShift);
        StsSensorClass sensorClass = (StsSensorClass)currentModel.getCreateStsClass("com.Sts.DBTypes.StsSensor");
        try
        {
            float[] xyz = null;
            float[] xy = null;

            float sizeFloat, colorFloat;

            // If exporting the sensors in view, allocate array
            int nExportPoints = 0;
            if(exportView)
            	exportPoints = new int[timeLongs.length];

            // If time is disabled always display the last element only.
            int start = 0;
            if(!getStaticSensorClass().getEnableTime())
                start = timeLongs.length-1;

            for (int i = start; i < timeLongs.length; i++)
            {
               if(colorFloats == null)
                   colorFloat = 1;
               else
                   colorFloat = colorFloats[i];

               float zValue = (float)getZLoc() + zShift;
               float xValue = (float)getXLoc() + xShift;
               float yValue = (float)getYLoc() + yShift;
               xyz = getXYZValue(xValue, yValue, zValue);

               // If time is enabled, verify the the time of the event is before the current time and in range.
               if(getStaticSensorClass().getEnableTime())
               {
                   // Check if before current time
                   if(!checkTime(timeLongs, i, points, duration, type) && !clustering)
               	        continue;
                   // Check if in range
                   if(!checkValue(timeLongs[i]) && !clustering)
               	        continue;
               }
               points++;

                float[] pt2d = new float[2];
                switch (dirNo)
                {
                    case 1:                    // Y
                        pt2d[0] = xyz[0];
                        pt2d[1] = xyz[2];
                        break;
                    case 0:                    // X
                        pt2d[0] = xyz[1];
                        pt2d[1] = xyz[2];
                        break;
                    case 2:
                        pt2d[0] = xyz[0];
                        pt2d[1] = xyz[1];
                        break;
                }

                // Compute the size of the symbol
                if(sizeFloats == null)
                   sizeFloat = getNumberBins();
                else
                    sizeFloat = computeSize(scaleMin, scaleMax, sizeFloats[i]);

                final StsColor color = defineColor(colorFloat);

                if(!clustering || (clusters == null))
                {
                	draw2DPoint(pt2d, color, sizeFloat, sensorClass.getScale2D(), glPanel3d, gl);
            		if(exportView)
            			exportPoints[nExportPoints++] = i;
                }
                else
                {
                	if(clusters[i] != -1)
                	{
                		if(clusters[i] != 99)
                			draw2DPoint(pt2d, StsColor.colors32[clusters[i]%32], sizeFloat, sensorClass.getScale2D(), glPanel3d, gl );
                		else
                			draw2DPoint(pt2d, color, sizeFloat, sensorClass.getScale2D(), glPanel3d, gl );
                		if(exportView)
                			exportPoints[nExportPoints++] = i;
                	}
                }
                accentEvents(i, getStaticSensorClass(), glPanel3d, gl, xyz, sizeFloat*2, sizeFloat, aColor, TWO_D);
            }

            if(exportView)
            {
            	exportPoints = (int[])StsMath.trimArray(exportPoints, nExportPoints);
            	exportView = false;
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsSensor.display failed.", e, StsException.WARNING);
            return;
        }
        finally
        {
            glPanel3d.resetViewShift(gl);
        }
    }

    /**
     * Display sensor events in 3D
     * @param glPanel3d - Graphics panel
     * @param displayName - Name of the display
     */
    public void display(StsGLPanel3d glPanel3d, boolean displayName)
    {
        display(glPanel3d, displayName, false);
    }

    /**
     * display sensor events in 3d
     * @param gl - Graphics handle
     * @param glPanel - Graphics panel
     */
    public void pick(GL gl, StsGLPanel glPanel)
    {
        display((StsGLPanel3d)glPanel, true, true);
    }
    /**
     * Display Sensor in 3D
     * @param glPanel3d - Graphics panel
     * @param displayName - Name of the display
     * @param pick - Is picking enabled
     */
    public void display(StsGLPanel3d glPanel3d, boolean displayName, boolean pick)
    {
        if(!verifyDisplayable())
            return;
        if(getSensorClass().getEnableTime())
        {
            if(!isAlive(currentModel.getProject().getProjectTime()))
                return;
        }
        GL gl = glPanel3d.getGL();

        // Temporarily ignore settings for duration and type if user is interactively selecting a region on a time series plot.
        long duration = this.displayDuration;
        byte type = this.displayType;

        if(colorscale == null) initializeColorscale();
        if(!checkLoadVectors())
        	return;

        int points = 0;
        try
        {
            float[] xyz, dsr;
            float[] xy;
            glu = glPanel3d.getGLU();
            float sizeFloat, colorFloat;
            StsColor aColor = getSensorClass().getAccentColor();

            // If exporting the sensors in view, allocate array
            int nExportPoints = 0;
            if(exportView)
            {
            	System.out.println("Allocating points");
            	exportPoints = new int[timeLongs.length];
            }

            // If time is disabled always display the last element only.
            int start = 0;
            if(!getStaticSensorClass().getEnableTime())
                start = timeLongs.length-1;

            for (int i = start; i < timeLongs.length; i++)
            {
                if(colorFloats == null)
                    colorFloat = 1;
                else
                    colorFloat = colorFloats[i];

                // Determine the xyz values.
                float zValue = (float)getZLoc() + zShift;
                float xValue = (float)getXLoc() + xShift;
                float yValue = (float)getYLoc() + yShift;
                xyz = getXYZValue(xValue, yValue, zValue);

                if(pointType == BEACHBALL)
                    dsr = new float[] {dipFloats[i], strikeFloats[i], rakeFloats[i]};
                else
                    dsr = new float[] {0.0f, 0.0f, 0.0f};

                // If time is enabled, verify the the time of the event is before the current time and in range.
                if(getStaticSensorClass().getEnableTime())
                {
                    // Check if before current time
                    if(!checkTime(timeLongs, i, points, duration, type) && !clustering)
                	    continue;
                    // Check if in range
                    if(!checkValue(timeLongs[i]) && !clustering)
                	    continue;
                }
                points++;

                // Compute the size of the symbol
                if(sizeFloats == null)
                   sizeFloat = getNumberBins();
                else
                   sizeFloat = computeSize(scaleMin, scaleMax, sizeFloats[i]);

                if(pick) gl.glPushName(i);
                final StsColor color = defineColor(colorFloat);

                if(!clustering || (clusters == null))
                {
                	draw3DPoint(glPanel3d, xyz, dsr, sizeFloat, color, false);
                    // Accumulate the points to export.
            		if(exportView)
            			exportPoints[nExportPoints++] = i;
                }
                else
                {
                	if(clusters[i] != -1)
                	{
                		if(clusters[i] != 99)
                			draw3DPoint(glPanel3d, xyz, dsr, sizeFloat, StsColor.colors32[clusters[i]%32], false);
                		else
                			draw3DPoint(glPanel3d, xyz, dsr, sizeFloat, color, false);
                        // Accumulate the points to export.
                		if(exportView)
                			exportPoints[nExportPoints++] = i;
                	}
                }
                accentEvents(i, getStaticSensorClass(), glPanel3d, gl, xyz, numberBins*2, sizeFloat, aColor, THREE_D);
                if(pick) gl.glPopName();

            }
            drawHighlightedPoints(glPanel3d);

            if(displayName)
            {
                float zValue = (float)getZLoc() + zShift;
                float xValue = (float)getXLoc() + xShift;
                float yValue = (float)getYLoc() + yShift;
                xyz = getXYZValue(xValue, yValue, zValue);
                StsGLDraw.fontHelvetica12(gl, xyz, name);
            }

            if(exportView)
            {
            	exportPoints = (int[]) StsMath.trimArray(exportPoints, nExportPoints-1);
            	exportView = false;
            	System.out.println("# of points in view is " + nExportPoints);
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsSensor.display failed.", e, StsException.WARNING);
            return;
        }
        finally
        {
            StsShader.disableCurrentShader(gl);
            if(pointType != BEACHBALL) return;
            drawBeachballAxes(glPanel3d);
        }
    }

    /**
     * check and load the positioning, time and required attribute arrays.
     * @return - success
     */
    public boolean checkLoadVectors()
    {
        colorFloats = null;
        sizeFloats = null;
        if(timeCurves == null)
        if(((StsTimeCurve)timeCurves.getElement(0)).getNumValues() == 0) return true;  // Real-time will be zero initially.

        timeLongs = ((StsTimeCurve)timeCurves.getElement(0)).getTimeVectorLongs();
        sizeFloats = getPropertyAsFloats();
        colorFloats = getColorByAsFloats();

        return super.checkLoadVectors();
    }

    public void setType(byte type)
    {
        this.type = type;
        if(type == StsSensorClass.TRACERS)
        {
            setSymbolString(SYMBOL_TYPE_STRINGS[CYLINDER]);
            setColorBy(findColorVectorByName("Copper"));
            setBornDeathWithVector();            
            setIsVisible(true);
        }
        else if(type == StsSensorClass.PUMPS)
        {
            setDisplayTypeString(SYMBOL_TYPE_STRINGS[CYLINDER]);
            setIsVisible(false);
        }
        else
        {
            setIsVisible(false);
        }
        setDisplayTypeString(displayTypeStrings[SHOW_SINGLE]);
        setNumberBins(25);        
        getStaticSensorClass().setEnableAccent(false);
    }

    /**
     * Determine if this sensor is displayable. May not be displayable because it does not support
     * the current domain or it is turned off.
     * @return
     */
    public boolean verifyDisplayable()
    {
        if(!super.verifyDisplayable())
            return false;
        if(!getStaticSensorClass().getEnable())
            return false;
        //if(processing)
        //    return false;
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
        StsMessageFiles.infoMessage("Cannot export static sensors from view. Use export option instead.");
        return 0;
    }
     /**
     * Get the sensor class containing this object
     * @return
     */
    static public StsDynamicSensorClass getDynamicSensorClass()
    {
        return (StsDynamicSensorClass)currentModel.getCreateStsClass(StsDynamicSensor.class);
    }
    /**
     * Convert the coordinates and attributes to a string
     * @param pointIdx - index of the desired point
     * @return - string containing the point information
     */
    public String toString(int pointIdx)
    {
        timeLongs = ((StsTimeCurve)timeCurves.getElement(0)).getTimeVector().getLongs();
        sizeFloats = getProperty().getValuesVectorFloats();
        colorFloats = getColorBy().getValuesVectorFloats();
        float sizeFloat, colorFloat;
        if(sizeFloats == null)
            sizeFloat = 0;
        else
            sizeFloat = sizeFloats[pointIdx];
        if(colorFloats == null)
            colorFloat = 0;
        else
            colorFloat = colorFloats[pointIdx];
        double x = getTimeCurve(StsLogVector.X).getValueVector().getOrigin() + xShift + (float)getXLoc();
        double y = getTimeCurve(StsLogVector.Y).getValueVector().getOrigin() + yShift + (float)getYLoc();
        String timeStg = currentModel.getProject().getTimeDateFormat().format(new Date(timeLongs[pointIdx]));
        return new String("Sensor(" + this.getName() + "): Time=" + timeStg + " X=" + (xLoc + xShift) + " Y=" + (yLoc + yShift) +
                                       " Z=" + (zLoc + zShift) + " " + propertyVector.getName() + "=" + sizeFloat +
                                   " " + colorVector.getName() + "=" + colorFloat);
    }

    /**
     * Determine if an event was selected by the mouse.
     * @param gl - Graphics handle
     * @param glPanel - Graphcis panel
     */
    public void pickX(GL gl, StsGLPanel glPanel)
    {
        if(!verifyDisplayable()) return;
        StsGLPanel3d glPanel3d = (StsGLPanel3d)glPanel;

        // Temporarily ignore settings for duration and type if user is interactively selecting a region on a time series plot.
        long duration = this.displayDuration;
        byte type = this.displayType;
        if(currentModel.getProject().getProjectTimeDuration() != 0)
        {
            duration = currentModel.getProject().getProjectTimeDuration();
            type = SHOW_INCREMENT;
        }

        if(!checkLoadVectors())
        	return;

        int points = 0;
        try
        {
            float[] xyz = null;
            float sizeFloat;
            for (int i = 0; i < timeLongs.length; i++)
            {
                if(sizeFloats == null)
                   sizeFloat = getNumberBins();
                else
                   sizeFloat = sizeFloats[i];

                // Determine the xyz values.
                float zValue = (float)getZLoc() + zShift;
                float xValue = (float)getXLoc() + xShift;
                float yValue = (float)getYLoc() + yShift;
                xyz = getXYZValue(xValue, yValue, zValue);

                if(getStaticSensorClass().getEnableTime())
                {
                    // Check if before current time
                    if(!checkTime(timeLongs, i, points, duration, type) && !clustering)
                	    continue;
                    // Check if in range
                    if(!checkValue(timeLongs[i]) && !clustering)
                	    continue;
                }
                points++;

                gl.glPushName(i);
                if(sizeFloats != null)
                    sizeFloat = computeSize(getProperty().getCurveMin(), getProperty().getCurveMax(), sizeFloats[i]);
                draw3DPoint(glPanel3d, xyz, sizeFloat, StsColor.BLACK, false);
                gl.glPopName();
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsSensor.display failed.", e, StsException.WARNING);
            return;
        }
    }

    /**
     * Draw highlighted points
     * @param glPanel3d - Graphics panel
     */
    public void drawHighlightedPoints(StsGLPanel3d glPanel3d)
    {
    	float size = this.getNumberBins();
    	StsColor color = getStaticSensorClass().getAccentColor();
    	if(highlightedPoints == null) return;
    	for(int i=0; i<highlightedPoints.length; i++)
    	{
    		float[] xyz = getXYZ(highlightedPoints[i]);
    		if(sizeFloats != null)
    			size = computeSize(scaleMin, scaleMax, sizeFloats[highlightedPoints[i]]);
    		draw3DPoint(glPanel3d, xyz, size, color, true);
    	}
    }
    /**
     * Get the XYZ values
     * @param pointIdx - index of the desired point
     * @return - XYZ value for desired point
     */
    public float[] getXYZ(int pointIdx)
    {
        // Determine the xyz values.
        float zValue = (float)getZLoc() + zShift;
        float xValue = (float)getXLoc() + xShift;
        float yValue = (float)getYLoc() + yShift;
        float[] xyz = getXYZValue(xValue, yValue, zValue);
    	return xyz;
    }

    /**
     * Get event as a point
     * @param pointIdx - index of the desired point
     * @return point
     */
    public StsPoint getPoint(int pointIdx)
    {
    	float[] xyz = getXYZ(pointIdx);
    	return new StsPoint(xyz);
    }
     /**
     * Get the sensor class containing this object
     * @return
     */
    static public StsStaticSensorClass getStaticSensorClass()
    {
        return (StsStaticSensorClass)currentModel.getCreateStsClass(StsStaticSensor.class);
    }

    /**
     * Set the static position for this sensor
     * @param x - X coordinate
     * @param y - Y coordinate
     * @param z - Z coordinate
     */
    public void setPosition(double x, double y, double z)
    {
        xLoc = x;
        yLoc = y;
        zLoc = z;
    }

    /**
     * Set the static position for this sensor
     * @param xyz - XYZ position
     */
    public void setPosition(double[] xyz)
    {
        xLoc = xyz[0];
        yLoc = xyz[1];
        zLoc = xyz[2];
    }
    /**
     * Set the static position for this sensor
     * @param xyz - XYZ position
     */
    public void setPosition(float[] xyz)
    {
        xLoc = xyz[0];
        yLoc = xyz[1];
        zLoc = xyz[2];
    }
    /**
     * Get the X location of a static sensor
     * @return static X location
     */
    public double getXLoc()
    {
        return xLoc;
    }
    /**
     * Get the Y location of a static sensor
     * @return static Y location
     */
    public double getYLoc()
    {
        return yLoc;
    }
    /**
     * Get the Z location of a static sensor
     * @return static Z location
     */
    public double getZLoc()
    {
        return zLoc;
    }

    /**
     * Set the X location of a static sensor
     * @param x coordinate
     */
    public void setXLoc(double x)
    {
        xLoc = x;
    }
    /**
     * Set the Y location of a static sensor
     * @param y coordinate
     */
    public void setYLoc(double y)
    {
        yLoc = y;
    }
    /**
     * Set the Z location of a static sensor
     * @param z coordinate
     */
    public void setZLoc(double z)
    {
        zLoc = z;
    }

    /**
     * Get the relative X origin
     * @return X origin
     */
    public double getXOrigin()
    {
    	return xLoc;
    }

    /**
     * Get the relative Y origin
     * @return Y Origin
     */
    public double getYOrigin()
    {
        return yLoc;
    }


    /**
     *  Zoom into the center of all the points in the series.
     */
    public boolean goTo()
    {
        // Determine the xyz values.
        float zValue = (float)getZLoc() + zShift;
        float xValue = (float)getXLoc() + xShift;
        float yValue = (float)getYLoc() + yShift;
        float[] xyz = getXYZValue(xValue, yValue, zValue);
        float distance = numberBins * getSensorClass().getGoToOffset();
    	return goToPoint(xyz, distance, BOTH_D);
    }

    /**
     *  Zoom into the center of all the points in the series.
     */
    public boolean goTo(double x, double y, double z)
    {
        // Determine the xyz values.
        float zValue = (float)x + zShift;
        float xValue = (float)y + xShift;
        float yValue = (float)z + yShift;
        float[] xyz = getXYZValue(xValue, yValue, zValue);
        float distance = numberBins * getSensorClass().getGoToOffset();
    	return goToPoint(xyz, distance, BOTH_D);
    }

     /**
     *   Zoom into the last point added to the series.
     */
    public boolean goToLastPoint()
    {
    	return goTo();
    }
      /**
     * Add new data to this sensor from real-time monitoring
     * @param source - source of the data
     * @param stype - source type
     * @param lastPollTime - last time data was polled
     * @return true when new data is found and successfully loaded
     */
    public int addNewData(String source, byte stype, long lastPollTime, boolean compute, boolean reload, boolean replace)
    {
    	int added = super.addNewData(source, stype, lastPollTime, compute, reload, replace);
        if(getSensorClass().getEnableGoTo())
        	goToLastPoint();
    	return added;
    }
}