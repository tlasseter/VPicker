//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Actions.Import.StsSensorKeywordIO;
import com.Sts.Actions.Import.StsUTKeywordIO;
import com.Sts.Actions.Import.StsWellKeywordIO;
import com.Sts.Actions.Wizards.SensorLoad.StsSensorFile;
import com.Sts.IO.StsAbstractFile;
import com.Sts.IO.StsFile;
import com.Sts.Interfaces.StsMonitorable;
import com.Sts.Interfaces.StsTimeEnabled;
import com.Sts.Interfaces.StsTreeObjectI;
import com.Sts.Interfaces.StsViewSelectable;
import com.Sts.MVC.StsMessageFiles;
import com.Sts.MVC.StsModel;
import com.Sts.MVC.StsProject;
import com.Sts.MVC.View3d.StsGLPanel3d;
import com.Sts.MVC.View3d.StsView;
import com.Sts.MVC.View3d.StsView3d;
import com.Sts.Types.StsPoint;
import com.Sts.Types.StsWellFile;
import com.Sts.UI.Progress.StsProgressPanel;
import com.Sts.Utilities.*;
import com.magician.fonts.GLBitmapFont;
import com.magician.fonts.GLHelvetica10BitmapFont;

import javax.media.opengl.GL;
import javax.swing.*;
import java.awt.*;
import java.io.BufferedReader;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

public class StsLiveWell extends StsWell implements StsMonitorable, StsViewSelectable
{
    transient StsProgressPanel progressPanel = StsProgressPanel.constructor();
    transient JDialog progressDialog = null;
    transient StsProgressPanel progressLogPanel = StsProgressPanel.constructor();
    transient JDialog progressLogDialog = null;
    transient boolean debugRealtime = false;
    transient int lastNumLines = 0;
    transient int lastLogNumLines = 0;
    transient boolean reloadRequired = false;
    transient boolean realtimeUpdate = false;
    transient float[] previousXYZ = null;
    transient byte previousDomain;

    // instance fields
    protected StsTimeVector timeVector;
    // current time index based on project time setting
    protected int timeIndex;
    /** Last modified date when deviation is real-time accumulated */
    protected long lastModified = 0;
    /** Last modified date when log is real-time accumulated */
    protected long lastLogModified = 0;
    /** Enable verbose Real-time Logging*/
    protected boolean isRealtimeLogging = false;

    /** default constructor */
    public StsLiveWell() // throws StsException
    {
    }

    public StsLiveWell(StsWell well)
    {
        super(false);
        StsToolkit.copyAllObjectFields(well, this, true);
    }

    public StsLiveWell(boolean persistent) // throws StsException
    {
        super(persistent);
    }

    public StsLiveWell(String name, boolean persistent)
    {
        this(name, persistent, getWellClass().getDefaultWellStsColor());
    }

    public StsLiveWell(String name, boolean persistent, StsColor color)
    {
        super(name, persistent, color);
    }

    static public StsLiveWell nullWellConstructor(String name)
    {
        return new StsLiveWell(name, false);
    }

    public boolean initialize(StsModel model)
    {
        boolean val = super.initialize();
        checkLoadVectors();   // Needed for continued real-time loading.
        timeIndex = timeVector.getLongs().length -1;
        return val;
    }
    private boolean checkLoadVectors()
    {
        if(getNLogCurves() == 0) return true;
        for(int i=0; i<getNLogCurves(); i++)
        {
            ((StsLogCurve)getLogCurves().getElement(i)).getMDepthVector().checkLoadVector();
            ((StsLogCurve)getLogCurves().getElement(i)).getValueVector().checkLoadVector();            
        }
        timeVector.checkLoadVector();
        return true;
    }
    public boolean setTimeIndex(long time)
    {
        if(timeVector == null) // Realtime no data yet.
            return false;
        timeIndex = (int)timeVector.getIndexNearest(time);
        if(timeIndex > (timeVector.getLongs().length-1))
            timeIndex = timeVector.getLongs().length -1;
        return true;
    }

    public void setTimeVector(StsTimeVector tVector)
    {
        timeVector = tVector;
    }

    public int addNewData(String source, byte sourceType, long lastPollTime, boolean compute, boolean reload, boolean replace)
    {
        if(source.contains("dev"))
            return addNewDeviationData(source, lastPollTime, compute, reload, replace);
        else
            return addNewLogData(source, lastPollTime, compute, reload, replace);
    }

    public int addNewLogData(String _source, long lastPollTime, boolean compute, boolean reload, boolean replace)
    {
        int added = 0;
        final String source = _source;

    	// If well was updated (reloaded) since last poll
        if((lastLogModified == 0) && (getLogCurves().getSize() != 0))
            lastLogModified = getTimeMax();

    	if(lastPollTime < lastLogModified)
    		lastPollTime = lastLogModified;

        // If the files are to be reloaded from scratch each time. Reset the vectors
        if(reload)
            resetCurves();

        // Read the new data - must already be in attribute list for this well.
        added = processNewLogDataFromFile(source, reload, replace);
        if((lastPollTime == 1) && (added == 0))
        {
            if(isRealtimeLogging)
                StsMessageFiles.logMessage("No new data found in file....");

            return 0;
        }
        SimpleDateFormat format = new SimpleDateFormat(currentModel.getProject().getTimeDateFormatString());
        if(isRealtimeLogging)
           StsMessageFiles.logMessage("Setting last modified time to: " + format.format(new Date(lastLogModified)));
        lastLogModified = lastPollTime;
        realtimeUpdate = true;
        currentModel.viewObjectChanged(this, this);
        if(lineVertices != null)
            numberOfElements = lineVertices.getSize();        
    	return added;
    }

    public int addNewDeviationData(String _source, long lastPollTime, boolean compute, boolean reload, boolean replace)
    {
        int added = 0;
        final String source = _source;

    	// If well was updated (reloaded) since last poll
        if((lastModified == 0) && (getLogCurves().getSize() != 0))
            lastModified = getTimeMax();

    	if(lastPollTime < lastModified)
    		lastPollTime = lastModified;

        // If the files are to be reloaded from scratch each time. Reset the vectors
        if(reload)
            resetCurves();

        // Read the new data - must already be in attribute list for this well.
        added = processNewDeviationDataFromFile(source, reload, replace);
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

        computePoints();
        StsSeismicVelocityModel velocityModel = currentModel.getProject().getSeismicVelocityModel();
        if (velocityModel != null)
        {
		     adjustFromVelocityModel(velocityModel);
             computeMarkerTimesFromMDepth(velocityModel);
        }
        realtimeUpdate = true;
        currentModel.viewObjectChanged(this, this);
    	return added;
    }

    public void display(StsGLPanel3d glPanel3d, boolean displayName, boolean displayAllMarkers,
                        boolean displayPerfMarkers, boolean displayFmiMarkers, StsLiveWellClass wellClass)
    {
        if (!currentModel.getProject().canDisplayZDomain(zDomainSupported)) return;
        if (glPanel3d == null) return;
        if (rotatedPoints == null) return;

        float displayToMd = getPoint(timeIndex).getM();

        labelFormat = new DecimalFormat(getWellClass().getLabelFormatAsString());
        if (isDrawingCurtain)
        {
            glPanel3d.setViewShift(glPanel3d.getGL(), 2.0f);
            if (displayName)
                super.display(glPanel3d, true, getName(), rotatedPoints);
            else
                super.display(glPanel3d, true, null, rotatedPoints);
        }
        else
        {
            if (displayName)
                super.display(glPanel3d, getName(), rotatedPoints);
            else
                super.display(glPanel3d, null, rotatedPoints);
        }
        if (isDrawingCurtain)
        {
            glPanel3d.resetViewShift(glPanel3d.getGL());
        }
        GL gl = glPanel3d.getGL();

        if (markers != null)
        {
            int nMarkers = markers.getSize();
            int nPerfs = 0;
            for (int n = 0; n < nMarkers; n++)
            {
                StsWellMarker marker = (StsWellMarker) markers.getElement(n);

                // Ignore Marker if outside display range
                if(marker.getMDepth() > displayToMd)
                    continue;

                // Draw the marker
                boolean drawDifferent = marker.getMarker().getModelSurface() != null;
                if (marker instanceof StsPerforationMarker)
                {
                    if ((displayPerfMarkers) && (drawPerfMarkers))
                        ((StsPerforationMarker)marker).display(glPanel3d, displayName, isDrawingCurtain, StsLiveWell.getLiveWellClass().getDefaultColor(nPerfs));
                    nPerfs++;
                }
                else if (marker instanceof StsEquipmentMarker)
                {
                    if (drawEquipmentMarkers)
                        marker.display(glPanel3d, displayName, isDrawingCurtain, drawDifferent);
                }
                else if (marker instanceof StsFMIMarker)
                {
                    if (displayFmiMarkers)
                        marker.display(glPanel3d, displayName, isDrawingCurtain, drawDifferent);
                }
                else
                {
                    if ((displayAllMarkers) && (drawMarkers))
                        marker.display(glPanel3d, displayName, isDrawingCurtain, drawDifferent);
                }
            }
        }

        logDisplay3d(leftDisplayLogCurves, -1f, glPanel3d);
        logDisplay3d(rightDisplayLogCurves, 0f, glPanel3d);

        if (!drawLabelString.equalsIgnoreCase(NO_LABEL) && (labelInterval >= 0.0f))
        {
            StsPoint point = null;
            float md = 0.0f;
            String label = null;
            int nLabels = (int) (getMaxMDepth() / labelInterval);

            if (isDrawingCurtain)
            {
                StsColor.BLACK.setGLColor(gl);
                glPanel3d.setViewShift(gl, 10.0f);
            }
            else
            {
                stsColor.setGLColor(gl);
                glPanel3d.setViewShift(gl, 1.0f);
            }

            GLBitmapFont font = GLHelvetica10BitmapFont.getInstance(gl);
            int numChars = font.getNumChars();
            for (int i = 0; i < nLabels; i++, md += labelInterval)
            {
                 // No label beyond timeIndex
                if(md > displayToMd)
                    break;
                point = getPointAtMDepth((float) (i * labelInterval), true);
                float[] xyz = point.getXYZorT();
                if ((md % (5.0f * labelInterval)) != 0.0f)
                {
                    StsGLDraw.drawPoint(xyz, null, glPanel3d, 5, 1, 0.0f);
                }
                else
                {
                    StsGLDraw.drawPoint(xyz, null, glPanel3d, 10, 2, 0.0f);
                    float value = 0.0f;
                    if (drawLabelString.equals(MDEPTH_LABEL))
                    {
                        value = md;
                        //                       label = Float.toString(md);
                    }
                    else if (drawLabelString.equals(DEPTH_LABEL))
                    {
                        value = point.getZ();
                        //                        label = Float.toString(point.getZ());
                    }
                    else if (drawLabelString.equals(TIME_LABEL))
                    {
                        value = point.getT();
                        //                       label = Float.toString(point.getT());
                    }
                    label = "    " + labelFormat.format(value);
                    StsGLDraw.fontOutput(gl, xyz, label, font);
                }
            }
            glPanel3d.resetViewShift(gl);
        }
        displaySeismicCurtain(glPanel3d);


        // Realtime additions
        if((lineVertices == null) || (timeIndex < 1))
            return;

        if((isDrawingCurtain) && (realtimeUpdate))
        {
            deleteSeismicCurtain();
            createCurtain();
            realtimeUpdate = false;
        }
        GLBitmapFont font = GLHelvetica10BitmapFont.getInstance(glPanel3d.getGL());
        float rop = 0.0f;
        StsPoint point = rotatedPoints[timeIndex];
        // Compute the rate of penetration
        if((timeVector != null) && (timeIndex > 0))
        {
            timeVector.checkLoadVector();
            long[] times = timeVector.getLongs();
            if(times.length > 2)
            {
                float deltaDepth = rotatedPoints[timeIndex].getZ() - rotatedPoints[timeIndex-1].getZ();
                long deltaTime = times[timeIndex] - times[timeIndex-1];
                double deltaHrs = (float)deltaTime/3600000.f;
                rop = (float)(deltaDepth/deltaHrs);
            }
        }
        float md = (float)point.getM();
        float depth = (float)point.getZ();
        float[] xyz = point.getXYZorT();

        glPanel3d.setViewShift(glPanel3d.getGL(), 2.0f);
        glPanel3d.getGL().glDisable(GL.GL_LIGHTING);

        StsGLDraw.drawSphere(glPanel3d, xyz, StsColor.BRONZE, 12.0f);

        String label = new String("    MD: " + md + " Depth: " + depth + " ROP: " + rop);
        StsGLDraw.fontOutput(glPanel3d.getGL(), xyz, label, font);

        glPanel3d.getGL().glEnable(GL.GL_LIGHTING);
        glPanel3d.resetViewShift(glPanel3d.getGL());

        // Re-position the user relative to the well
        if(wellClass.getRealtimeTrack() && previousXYZ != null && (previousDomain == currentModel.getProject().getZDomain()))
        {
            StsView view = glPanel3d.getView();
            if(view instanceof StsView3d)
            {
                float[] viewParms = ((StsView3d)view).getCenterAndViewParameters();
                viewParms[0] = viewParms[0] + (xyz[0]-previousXYZ[0]);
                viewParms[1] = viewParms[1] + (xyz[1]-previousXYZ[1]);
                viewParms[2] = viewParms[2] + (xyz[2]-previousXYZ[2]);
                ((StsView3d)view).setCenterAndViewParameters(viewParms);
            }
        }
        previousXYZ = xyz;
        previousDomain = currentModel.getProject().getZDomain();
    }

    static public StsLiveWellClass getLiveWellClass()
    {
        return (StsLiveWellClass) currentModel.getCreateStsClass(StsLiveWell.class);
    }

    public void display2d(StsGLPanel3d glPanel3d, boolean displayName, int dirNo,
                          float dirCoordinate, boolean axesFlipped, boolean xAxisReversed, boolean yAxisReversed,
                          boolean displayIn2d, boolean displayAllMarkers,
                          boolean displayPerfMarkers, boolean displayFmiMarkers)
    {
        if(timeVector == null) // No data in well yet. Realtime stub
            return;
        if (!currentModel.getProject().canDisplayZDomain(zDomainSupported))
        {
            return;
        }
        if (glPanel3d == null)
        {
            return;
        }
        float displayToMd = getPoint(timeIndex).getM();
        super.display2d(glPanel3d, timeIndex, displayName, dirNo, dirCoordinate, axesFlipped, xAxisReversed, yAxisReversed);

        // Display BH statistics
        float rop = 0.0f;
        StsPoint point = rotatedPoints[timeIndex];
        if((timeVector != null) && (timeIndex > 0))
        {
            timeVector.checkLoadVector();
            long[] times = timeVector.getLongs();
            if(times.length > 2)
            {
                float deltaDepth = rotatedPoints[timeIndex].getZ() - rotatedPoints[timeIndex-1].getZ();
                long deltaTime = times[timeIndex] - times[timeIndex-1];
                double deltaHrs = (float)deltaTime/3600000.f;
                rop = (float)(deltaDepth/deltaHrs);
            }
        }
        float md = (float)point.getM();
        float depth = (float)point.getZ();
        float[] xyz = point.getXYZorT();
        String label = new String("    MD: " + md + " Depth: " + depth + " ROP: " + rop);
        
        displayLabel2d(glPanel3d.getGL(), xyz, label, dirNo, axesFlipped);

        if ((markers != null) && (displayIn2d))
        {
            int nMarkers = markers.getSize();
            int nPerfs = 0;

            for (int n = 0; n < nMarkers; n++)
            {
                StsWellMarker marker = (StsWellMarker) markers.getElement(n);
                // Ignore Marker if outside display range
                if(marker.getMDepth() > displayToMd)
                    continue;
                boolean drawDifferent = marker.getMarker().getModelSurface() != null;
                if (marker instanceof StsPerforationMarker)
                {
                    if ((displayPerfMarkers) && (drawPerfMarkers))
                        ((StsPerforationMarker) marker).display2d(glPanel3d, dirNo, displayName, StsLiveWell.getLiveWellClass().getDefaultColor(nPerfs) , 1.0f);
                    nPerfs++;
                }
                else if (marker instanceof StsEquipmentMarker)
                {
                    if (drawEquipmentMarkers)
                        ((StsEquipmentMarker) marker).display2d(glPanel3d, dirNo, displayName, drawDifferent);
                }
                else if (marker instanceof StsFMIMarker)
                {
                    if (displayFmiMarkers)
                        ((StsFMIMarker) marker).display2d(glPanel3d, dirNo, displayName);
                }
                else
                {
                    if ((displayAllMarkers) && (drawMarkers))
                        marker.display2d(glPanel3d, dirNo, displayName, drawDifferent);
                }
            }
        }
    }

    public void display(StsGLPanel3d glPanel3d, boolean highlighted, String name, StsPoint[] points, boolean drawDotted)
    {
        if (glPanel3d == null) return;
        if (!isVisible) return;
        if(timeIndex < 1) return;
        GL gl = glPanel3d.getGL();
        if (gl == null) return;

        if (points == null)
        {
            if (drawVertices) displayVertices(glPanel3d);
            return;
        }

        if(!currentModel.getProject().supportsZDomain(zDomainSupported)) return;

        int zIndex = getPointsZIndex();
        StsPoint topPoint = points[0];
        StsPoint botPoint = points[timeIndex];

        super.drawLine(gl, stsColor, highlighted, points, 0, timeIndex, topPoint, botPoint, zIndex, drawDotted);

        if (name != null)
        {
            displayName(glPanel3d, name, topPoint, botPoint, zIndex);
        }
        if (drawVertices)
        {
            displayVertices(glPanel3d);
        }
        return;
    }
    /**
     * Add new data found in a monitored file
     * @param source - file name
     * @return new data added or not
     */
    public int processNewDeviationDataFromFile(String source, boolean reload, boolean replace)
    {
        boolean valuesReplaced = false;
    	int added = 0;
        int nCurves = 0;
        long[] timeVals = null;
        float[][] values = null;
        SimpleDateFormat format = new SimpleDateFormat(currentModel.getProject().getTimeDateFormatString());

    	try
    	{
            StsFile file = StsFile.constructor(source);
            StsWellFile wellFile = new StsWellFile(file);
            int nLines = wellFile.getNumLinesInFile();
            wellFile.analyzeFile(currentModel);

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
            if(!wellFile.analyzeFile(currentModel))
                return 0;

            for(int ii=0; ii<lastNumLines; ii++)
                wellFile.readLine(); // Lines previously processed

            if(isRealtimeLogging)
                StsMessageFiles.logMessage("Processing lines " + lastNumLines + " to " + (nLines-1) + " from " + file.getFilename());

            String sTime = null;
            String line = null;
    		// Loop in case multiple lines were added.
            for(int j=lastNumLines; j<nLines; j++)
            {
    		    line = wellFile.readLine();

    			if(!wellFile.getAttributeValues(currentModel, line))
                    continue;
                long time = wellFile.currentTime;
                sTime = format.format(new Date(time));

    		    if(lineVertices == null)
                {
                    if(debugRealtime) System.out.println("Creating vectors: 1st time is:" + sTime);
                    if(isRealtimeLogging)
                        StsMessageFiles.logMessage("Creating vectors for well (" + getName() + ") and added first value at " + sTime);
    			    createNullVectors(wellFile.curveNames, wellFile);
                    continue;
                }
                else
                {
    			    // Verify the event is newer than the max event in well
                    if(debugRealtime) System.out.println("Verifying time:" + sTime);
    			    if(!verifyEvent(wellFile))
                    {
                        if(replace)
                        {
                            if(replaceValue(wellFile))
                            {
                                StsMessageFiles.logMessage("New data at " + sTime + " is in the past and does not exist in well. It will be ignored." );
                                valuesReplaced = true;
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
                    nCurves = wellFile.currentValues.length;
                    timeVals = null;
                    values = new float[nCurves][];
                    if(nCurves != 0)
    		        {
                        timeVals = new long[nLines];
                        for(int i=0; i<nCurves; i++)
                            values[i] = new float[nLines];
                    }
                }
                if(debugRealtime)
                    System.out.println("Adding value for time:" + sTime);

    			timeVals[added] = wellFile.currentTime;
                for(int i=0; i<nCurves; i++)
                    values[i][added] = (float)wellFile.currentValues[i];

                if(progressDialog != null)
                    progressPanel.progressBar.setValueImmediate(j-lastNumLines+1);
    			added++;
                if(added%250 == 0)
                {
                    StsMessageFiles.infoMessage("Processed " + added + " real-time points into deviation(" + getName() + ")");
                }
    		}
            // Keep track of number of last read so they can be skipped on next cycle.
            if(!reload)
                lastNumLines = nLines;
            else
                lastNumLines = 0;
            // Trim the arrays
            timeVals = (long[]) StsMath.trimArray(timeVals,added);
            for(int i=0; i<nCurves; i++)
                values[i] = (float[])StsMath.trimArray(values[i],added);

            if(values != null)
            {
                if(isRealtimeLogging)
                    StsMessageFiles.logMessage("Appending " + values.length + " values to deviation...");
                appendDevValues(wellFile, timeVals, values, reload);
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
            StsMessageFiles.errorMessage("Failed to process monitored well data: " + source);
            return 0;
    	}
    }

    /**
     * Add new data found in a monitored log file
     * @param source - file name
     * @return new data added or not
     */
    public int processNewLogDataFromFile(String source, boolean reload, boolean replace)
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
            StsWellFile wellFile = new StsWellFile(file);
            int nLines = wellFile.getNumLinesInFile();
            wellFile.analyzeFile(currentModel);

            if(lastLogNumLines == nLines)
            {
                currentModel.viewObjectChanged(this,this);
                return 0;
            }

            // Determine the number of lines to process.
		    int nLinesToProcess = nLines;
		    if(!reload)
			    nLinesToProcess = nLines - lastLogNumLines;
            if(nLinesToProcess > currentModel.getProject().getRealtimeBarAt())
            {
                processFile();
                progressLogPanel.progressBar.initializeImmediate(nLinesToProcess-1);
            }
            if(!wellFile.analyzeFile(currentModel))
                return 0;

            for(int ii=0; ii<lastLogNumLines; ii++)
                wellFile.readLine(); // Lines previously processed

            if(isRealtimeLogging)
                StsMessageFiles.logMessage("Processing lines " + lastLogNumLines + " to " + (nLines-1) + " from " + file.getFilename());

            String sTime = null;
            String line = null;
    		// Loop in case multiple lines were added.
            for(int j=lastLogNumLines; j<nLines; j++)
            {
    		    line = wellFile.readLine();

    			if(!wellFile.getAttributeValues(currentModel, line))
                    continue;

    		    if(getNLogCurves() == 0)
                {
                    if(debugRealtime) System.out.println("Creating log vectors: 1st time is:" + sTime);
                    if(isRealtimeLogging)
                        StsMessageFiles.logMessage("Creating log vectors for well (" + getName() + ") and added first value at " + sTime);
    			    createNullLogVectors(wellFile.curveNames, wellFile);
                    continue;
                }
                if(values == null)
                {
                    nCurves = wellFile.curveNames.length;
                    values = new double[nCurves][];
                    for(int i=0; i<nCurves; i++)
                        values[i] = new double[nLines];
                }
                for(int i=0; i<nCurves; i++)
                    values[i][added] = wellFile.currentValues[i];

                if(progressLogDialog != null)
                    progressLogPanel.progressBar.setValueImmediate(j-lastLogNumLines+1);
    			added++;
                if(added%250 == 0)
                {
                    StsMessageFiles.infoMessage("Processed " + added + " real-time points into log curves(" + getName() + ")");
                }
    		}
            // Keep track of number of last read so they can be kipped on next cycle.
            if(!reload)
                lastLogNumLines = nLines;
            else
                lastLogNumLines = 0;
            // Trim the arrays
            for(int i=0; i<nCurves; i++)
                values[i] = (double[])StsMath.trimArray(values[i],added);

            if(values != null)
            {
                if(isRealtimeLogging)
                    StsMessageFiles.logMessage("Appending " + values.length + " values to logs...");
                appendLogValues(wellFile, values, reload);
            }

            if(isRealtimeLogging)
                StsMessageFiles.logMessage("Updating timeseries, 2D and 3D plots... ");
            currentModel.viewObjectChanged(this, this);
            if(progressLogDialog != null)
                progressLogDialog.setVisible(false);
            return added;
    	}
    	catch(Exception e)
    	{
            StsMessageFiles.errorMessage("Failed to process monitored well data: " + source);
            return 0;
    	}
    }

     /**
     *  Verify that the current events are not already in the set.
     * @param wFile - Well file description
     * @return successfully added data
      */
    public boolean verifyEvent(StsWellFile wFile)
    {
        long maxTime = getTimeMax();
        //System.out.println("currentTime= " + sFile.currentTime + " maxTime= " + maxTime);
        if(wFile.currentTime <= maxTime)
        {
            //if(debugRealtime) System.out.println("Realtime event already in well for time stamp: " +
        	//		currentModel.project.getTimeDateFormat().format(new Date(sFile.currentTime)));
        	return false;
        }
        if(debugRealtime) System.out.println("Event at " + wFile.currentTime + " is new.");
    	return true;
    }
       /**
     * Add new values to deviation curves from file
     * @param wFile - file with new well data to add
     * @return success
      */
    public boolean appendDevValues(StsWellFile wFile, long[] times, float[][] values, boolean reload)
    {
        // Create vectors
        if(isRealtimeLogging)
        {
            for(int i=0; i<wFile.curveNames.length; i++)
                StsMessageFiles.logMessage("   Adding " + wFile.curveNames[i] + " deviation vector values to well(" + getName() + ")");
        }
        String binaryDataDir = currentModel.getProject().getBinaryFullDirString();
        timeVector.writeAppend(times, binaryDataDir);
        dbFieldChanged("timeVector", timeVector);
         
        StsWellKeywordIO.filename = wFile.file.getFilename();
        StsWellKeywordIO.parseAsciiFilename(wFile.file.getFilename());
    	StsLogVector[] logVectors = StsWellKeywordIO.constructLogVectors(wFile.curveNames, StsLogVector.WELL_DEV_PREFIX);
        StsLogVector.getVectorOfType(logVectors, StsLogVector.X).setOrigin(wFile.xOrigin);
        StsLogVector.getVectorOfType(logVectors, StsLogVector.X).setUnits(currentModel.getProject().getXyUnits());
        StsLogVector.getVectorOfType(logVectors, StsLogVector.Y).setOrigin(wFile.yOrigin);
        StsLogVector.getVectorOfType(logVectors, StsLogVector.Y).setUnits(currentModel.getProject().getXyUnits());
        for(int i=0; i<logVectors.length; i++)
        {
            logVectors[i].setValues(values[i]);
            logVectors[i].setMinMaxAndNulls(currentModel.getProject().getLogNull());
        }

        // Initialize the binary files
        StsSensorKeywordIO.deleteBinaryFiles(currentModel.getProject().getBinaryFullDirString(), wFile.curveNames, StsLogVector.WELL_DEV_PREFIX);
        StsSensorKeywordIO.deleteBinaryFiles(currentModel.getProject().getBinaryFullDirString(), wFile.curveNames, StsLogVector.WELL_LOG_PREFIX);
        StsSensorKeywordIO.checkWriteBinaryFiles(timeVector, logVectors, currentModel.getProject().getBinaryFullDirString());

        currentModel.getCreateCurrentTransaction("realTimeWellUpdate");
        constructWellDevCurves(logVectors, currentModel.getProject().getLogNull(), null);
        currentModel.commit();
        
        if(isRealtimeLogging) StsMessageFiles.logMessage("   Successfully added deviation values to well(" + getName() + ")");
        return true;
    }
      /**
     * Add new values to deviation curves from file
     * @param wFile - file with new well data to add
     * @return success
      */
    public boolean appendLogValues(StsWellFile wFile, double[][] values, boolean reload)
    {
        // Create vectors
        if(isRealtimeLogging)
        {
            for(int i=0; i<wFile.curveNames.length; i++)
                StsMessageFiles.logMessage("   Adding " + wFile.curveNames[i] + " log vector values to well(" + getName() + ")");
        }
        //
        // ToDo: Currently assumes mdepth is the first array. Must allow any column mdepth
        //
        for(int i = 0; i<getNLogCurves(); i++)
        {
            if(i == 0)
                ((StsLogCurve)logCurves.getElement(i)).addValuesToMDepthCurveWrite(values[0]);
            //else            
            //    ((StsLogCurve)logCurves.getElement(i)).addValuesToMDepthCurve(values[0]);
            ((StsLogCurve)logCurves.getElement(i)).addValuesToCurveWrite(values[i]);            
        }

        if(isRealtimeLogging) StsMessageFiles.logMessage("   Successfully added real-time curves to sensor(" + getName() + ")");
        return true;
    }

    public void processFile()
    {
        progressDialog = new JDialog(currentModel.win3d, "Loading Realtime Deviation Data..." + getName(), false);
        progressDialog.setPreferredSize(new Dimension(300,50));
        progressDialog.add(progressPanel);
        progressDialog.pack();
        progressDialog.setLocationRelativeTo(currentModel.win3d);
        progressDialog.setVisible(true);
        return;
    }

    public void processLogFile()
    {
        progressLogDialog = new JDialog(currentModel.win3d, "Loading Realtime Log Data..." + getName(), false);
        progressLogDialog.setPreferredSize(new Dimension(300,50));
        progressLogDialog.add(progressLogPanel);
        progressLogDialog.pack();
        progressLogDialog.setLocationRelativeTo(currentModel.win3d);
        progressLogDialog.setVisible(true);
        return;
    }
    public long getTimeMax() { return timeVector.getMaxValue(); }
    public void resetCurves()
    {
         if(timeVector != null)
         {
             timeVector.resetVector();
             lineVertices = null;
             for(int i=0; i<getNLogCurves(); i++)
                  ((StsLogCurve)logCurves.getElement(i)).resetVectors();
         }
    }

    /**
     * Create required null vectors for real-time sensor
     * @param names - array of vector names
     */
    private void createNullVectors(String[] names, StsWellFile wFile)
    {
        // Create vectors
        if(isRealtimeLogging)
        {
            for(int i=0; i<names.length; i++)
                StsMessageFiles.logMessage("   Adding " + names[i] + " deviation vectors to well(" + getName() + ")");
        }
    	timeVector = StsWellKeywordIO.constructTimeVector("LiveWell", wFile.getColLocation(wFile.TIME), wFile.getColLocation(wFile.TIME));
    	timeVector.setMinMaxAndNulls(StsParameters.largeLong);
        StsWellKeywordIO.filename = wFile.file.getFilename();
        StsWellKeywordIO.parseAsciiFilename(wFile.file.getFilename());
        
    	StsLogVector[] logVectors = StsWellKeywordIO.constructLogVectors(names, StsLogVector.WELL_DEV_PREFIX);
        StsLogVector.getVectorOfType(logVectors, StsLogVector.X).setOrigin(wFile.xOrigin);

        StsLogVector.getVectorOfType(logVectors, StsLogVector.X).setUnits(currentModel.getProject().getXyUnits());
        StsLogVector.getVectorOfType(logVectors, StsLogVector.Y).setOrigin(wFile.yOrigin);

        StsLogVector.getVectorOfType(logVectors, StsLogVector.Y).setUnits(currentModel.getProject().getXyUnits());
        for(int i=0; i<logVectors.length; i++)
        {
            logVectors[i].setValues(new float[] {(float)wFile.currentValues[i]});
            logVectors[i].setMinMaxAndNulls(currentModel.getProject().getLogNull());
        }
        StsLongVector lVector = new StsLongVector(new long[] {(long)wFile.currentTime});
        timeVector.setValues(lVector);
        timeVector.setMinMaxAndNulls(StsParameters.largeLong);
        dbFieldChanged("timeVector", timeVector);
        setBornDate((long)wFile.currentTime);
        dbFieldChanged("bornDate", bornDate);
        
        // Initialize the binary files
        StsSensorKeywordIO.deleteBinaryFiles(currentModel.getProject().getBinaryFullDirString(), wFile.curveNames, StsLogVector.WELL_DEV_PREFIX);
        StsSensorKeywordIO.deleteBinaryFiles(currentModel.getProject().getBinaryFullDirString(), wFile.curveNames, StsLogVector.WELL_LOG_PREFIX);
        StsSensorKeywordIO.checkWriteBinaryFiles(timeVector, logVectors, currentModel.getProject().getBinaryFullDirString());

        currentModel.getCreateCurrentTransaction("realTimeWellUpdate");
        constructWellDevCurves(logVectors, currentModel.getProject().getLogNull(), null);
        currentModel.commit();

        dbFieldChanged("xOrigin", xOrigin);
        dbFieldChanged("yOrigin", yOrigin);

        currentModel.getProject().adjustBoundingBoxes(true, false); // extend displayBoundingBox as needed and set cursor3d box accordingly
        currentModel.getProject().checkAddUnrotatedClass(StsWell.class);

        //currentModel.viewObjectChanged(this, this);
        if(isRealtimeLogging) StsMessageFiles.logMessage("   Successfully added real-time deviation curves to well(" + getName() + ")");
    }

    /**
     * Create required null vectors for real-time sensor
     * @param names - array of vector names
     */
    private void createNullLogVectors(String[] names, StsWellFile wFile)
    {
        // Create vectors
        if(isRealtimeLogging)
        {
            for(int i=0; i<names.length; i++)
                StsMessageFiles.logMessage("   Adding " + names[i] + " log vectors to well(" + getName() + ")");
        }
    	StsLogVector[] logVectors = StsWellKeywordIO.constructLogVectors(names, StsLogVector.WELL_LOG_PREFIX);
        for(int i=0; i<logVectors.length; i++)
        {
            logVectors[i].setValues(new float[] {(float)wFile.currentValues[i]});
            logVectors[i].setMinMaxAndNulls(currentModel.getProject().getLogNull());
        }
        StsLogCurve[] newLogCurves = StsLogCurve.constructLogCurves(this, logVectors, currentModel.getProject().getLogNull(), 0);

        // Initialize the binary files
        StsSensorKeywordIO.deleteBinaryFiles(currentModel.getProject().getBinaryFullDirString(), wFile.curveNames, StsLogVector.WELL_LOG_PREFIX);

        currentModel.getCreateCurrentTransaction("realTimeLogUpdate");
        addLogCurves(newLogCurves);
        currentModel.commit();

        if(isRealtimeLogging) StsMessageFiles.logMessage("   Successfully added real-time curves to sensor(" + getName() + ")");
    }

    public boolean replaceValue(StsWellFile wFile)
    {
        StsLogCurve logCurve = null;

          // Locate and replace value in appropriate curves
        for(int i=0; i<getNLogCurves(); i++)
        {
            logCurve = (StsLogCurve)logCurves.getElement(i);
            for(int j=0; j<wFile.validCurves.length; j++)
            {
            	if(logCurve.getValueVector().getName().equalsIgnoreCase(wFile.validCurves[j]))
            	{
                    int index = timeVector.getIndex(wFile.currentTime);
            		if(!logCurve.replaceValueInCurve(wFile.currentTime, wFile.currentValues[j], index))
                        return false;
            	}
            }
            logCurve.getValueVector().getValues().setMinMax();
        }
        reloadRequired = true;
        return true;
    }


    public int getNPoints()
    {
        return (timeVector == null) ? 0 : timeVector.getLongs().length;
    }
    public StsTimeVector getTimeVector() { return timeVector; }
    public float[] getMDepthFloats()
    {
        if(lineVertices == null) return null;
        if(lineVertices.getSize() < 2) return null;

        float[] flts = new float[lineVertices.getSize()];
        for(int i=0;i<lineVertices.getSize(); i++)
            flts[i] = ((StsSurfaceVertex)lineVertices.getElements()[i]).getPoint().getM();
        return flts;
    }
    public float[] getDepthFloats()
    {
        if(lineVertices == null) return null;
        if(lineVertices.getSize() < 2) return null;

        float[] flts = new float[lineVertices.getSize()];
        for(int i=0;i<lineVertices.getSize(); i++)
            flts[i] = ((StsSurfaceVertex)lineVertices.getElements()[i]).getPoint().getZ();
        return flts;
    }

    public void logDisplay3d(StsLogCurve[] logCurves, float origin, StsGLPanel3d glPanel3d)
    {
        if (logCurves == null)return;
        for(StsLogCurve logCurve : logCurves)
            logCurve.display3d(glPanel3d, this, origin, ((StsSurfaceVertex)lineVertices.getElement(timeIndex)).getPoint().getM());
    }

    public boolean addToProject()
    {
        return super.addToProject();
    }

    public long getBornDateLong()
    {
        if(timeVector != null)
            return timeVector.getMinValue();
        return bornDate;
    }

    public long getDeathDateLong()
    {
        if(timeVector != null)
            return timeVector.getMaxValue();
        return deathDate;
    }

    public boolean isAlive(long time)
    {
        boolean val = super.isAlive(time);
        if(timeVector == null)
            return false;
        else
            return val;
    }
}