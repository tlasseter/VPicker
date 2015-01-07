package com.Sts.MVC.View3d;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: Base class used to present two-dimensional data.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.1
 */

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import info.monitorenter.gui.chart.*;
import info.monitorenter.gui.chart.io.*;
import info.monitorenter.gui.chart.traces.*;
import info.monitorenter.gui.chart.traces.painters.*;

import javax.swing.*;
import java.awt.*;
import java.util.*;

public class StsViewRealtimeSeries extends StsViewTimeSeries
{
    static public final String viewClassname = "RealtimeSeries_View";

    private boolean emptyPlot;
    //private boolean labelOn = false;
    transient JLabel rtLabel;
    transient boolean debugRealtime = true;
    transient StsStaticDataCollector collectors[] = null;
    private int plotDuration;

    static public final String viewName = "RealtimeSeries_View";

    /**
     * Default constructor
     */
    public StsViewRealtimeSeries()
    {
     //   System.out.println("Default Realtime Series View Constructor.");
    }

    public StsViewRealtimeSeries(StsWin3dBase window)
    {
        super(window);
        if(selectedObjects == null) return;
        if(selectedObjects[0] == null) return;
        plotDuration = ((StsSensor)selectedObjects[0]).getTimeSeriesDuration();
        createDataCollectors();   
    }

    public boolean buildDialog()
    {
        try
        {
        	if((StsSensor)selectedObjects[0] == null)
                return false;
            this.setTitle("Realtime Series Plot - "  + ((StsSensor)selectedObjects[0]).toString());
            getAttributeData();
            setAlwaysOnTop(true);
            jbInit();
            if(emptyPlot)
            {
                rtLabel = new JLabel("Realtime Timeseries Plot - Awaiting data");
                jPanel1.add(rtLabel);
                //labelOn = true;
            }
            else
            {
                initialize();
                //setCurrentTraceToAmplitude();
            }
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
            return false;
        }
        return true;
    }

    public void buildTraces(Chart2D chart1)
    {
        boolean ampSet = false;
        ITracePainter barPainter = new TracePainterVerticalBar(2, chart1);
        for(int i=0; i<names.length; i++)
        {
            final Color color = StsColor.colors32[i%32].getColor();
            traces[i] = createDataset(names[i], attributes[i]);
            traces[i].addPropertyChangeListener(ITrace2D.PROPERTY_VISIBLE, this);
            traces[i].setColor(color);
            traces[i].setVisible(visible[i]);
            if((names[i].toLowerCase().contains("amp") || names[i].toLowerCase().contains("mag")) && !ampSet)
            {
                traces[i].setTracePainter(barPainter);
                ampSet = true;
            }
            int numPointsShown = ((StsSensor)selectedObjects[0]).getNumPoints(((StsSensor)selectedObjects[0]).getTimeSeriesDuration(), project.getProjectTime());   // time interval may be irregular so must call each time
            ((Trace2DLtd)traces[i]).setMaxSize(numPointsShown);
            chart1.addTrace(traces[i]);
            if(chart1.getCurrentTrace().getMaxY() != chart1.getCurrentTrace().getMinY())
                configureAxis(chart1);
        }
    }

    public void createDataCollectors()
    {
        if(traces == null) return;
        // Configure Data Collectors
        collectors = new StsStaticDataCollector[traces.length];
        for(int i=0; i<traces.length; i++)
        {
            collectors[i] = new StsStaticDataCollector(traces[i], ((StsSensor)selectedObjects[0]), plotDuration);
            collectors[i].changePolicy(chart, project.getProjectTime());
        }
    }

    public void changePolicy(Chart2D chart1) {}

    public boolean getAttributeData()
    {
        emptyPlot = false;
        names = ((StsSensor)selectedObjects[0]).getTimeCurveList();
        nAttributes = names.length;
        int numPoints = ((StsSensor)selectedObjects[0]).getNumPoints(((StsSensor)selectedObjects[0]).getTimeSeriesDuration(), project.getProjectTime());     // may have data just none in the latest time interval
        if((nAttributes == 0) || (numPoints < 2))    // Need at least two values to create plot.
        {
            if(debugRealtime) System.out.println("EmptyPlot: numAttributes=" + nAttributes + " & numValues=" + ((StsSensor)selectedObjects[0]).getNumValues());
            emptyPlot = true;
            return false;   // May be empty for realtime.
        }
        return super.getSensorAttributeData();
    }

    public boolean reconfigureView()
    {
        if(emptyPlot || (chart == null))
            return false;
        //return super.reconfigureView();
        return true;
    }

    public boolean viewObjectChanged(Object source, Object object)
    {
        if(!project.isRealtime())
            return true;
        if(object instanceof StsSensor)
        {
        	if(object == ((StsSensor)selectedObjects[0]))
        	{
                if(plotDuration != ((StsSensor)selectedObjects[0]).getTimeSeriesDuration())
                    StsMessageFiles.infoMessage("To reset realtime plot duration, close window and re-open.");

                if(emptyPlot || ((StsSensor)selectedObjects[0]).isReloadRequired())
                {
                   ((StsSensor)selectedObjects[0]).resetReload();
                   changePolicies(project.getProjectTime());                    
                   if(!checkLoadTraces())
                       return false;
                }
                else
                {
                    if(!collectData())
                        return false;
                }
                //super.reconfigureView();
                setCurrentTraceToAmplitudeOrPressure();
                chart.updateUI();
                chart.setRequestedRepaint(true);
        	}
        }
        return super.viewObjectChanged(source, object);
    }

    public boolean collectData()
    {
        boolean foundNewData = false;
        if(traces == null || collectors == null || !project.isRealtime())
            return false;
        for(int i=0; i<traces.length; i++)
            foundNewData = collectors[i].collectData(chart, project.getProjectTime());

        return true;
    }
    
    public boolean reloadData()
    {
        StsObjectRefList timeCurves = ((StsSensor)selectedObjects[0]).getTimeCurves();
        for(int i=0; i<timeCurves.getSize(); i++)
        {
            StsTimeCurve curve = (StsTimeCurve)timeCurves.getElement(i);
            if(i==0)
                time = curve.getTimeVectorLongs();
            attributes[i] = curve.getValuesVectorFloats();
        }

        for(int i=0; i<traces.length; i++)
        {
            chart.setCurrentTrace(i);
            traces[i].removeAllPoints();
            //System.out.println("Time length= " + time.length);
            for (int j = 0; j < time.length; j++)
            {
        	    if(maxTime < time[j]) maxTime = time[j];
        	    if(minTime > time[j]) minTime = time[j];
                if(attributes[i][j] == model.getProject().getLogNull())
                    traces[i].addPoint(time[j], traces[i].getMinY());
                else
                    traces[i].addPoint(time[j], attributes[i][j]);
            }
        }
        return true;
    }

    public boolean checkLoadTraces()
    {
        if(chart != null)
        {
            StsMessageFiles.infoMessage("Reloading time series data due to reload or replacement of events.");
            return reloadData();
        }
         if(!getAttributeData())
         {
             // Open empty time series view for realtime.
             if(debugRealtime) System.out.println("Not enough data found.");
             StsMessageFiles.infoMessage("Not enough data found in selected sensor...");
             rtLabel.setText("Not enough data found in selected sensor yet, awaiting more...");
             //labelOn = true;
             return false;
         }
         if(debugRealtime) System.out.println("Initializing realtime plot....");
         initialize();
         //super.reconfigureView();
         createDataCollectors();
         jPanel1.remove(rtLabel);
         //labelOn = false;
         setCurrentTraceToAmplitudeOrPressure();
         if(debugRealtime) System.out.println("Values in traces=" + traces[0].getSize() + " should be 2");
         return true;
    }

    public boolean viewObjectRepaint(Object source, Object object)
    {
        if(!emptyPlot && traces != null && collectors != null)
        {
            // realtime is no longer running, user clicked inside plot.
            //if(project.getProjectTime() < traces[0].getMaxX())
            //{
            if(!project.isRealtime())
            {
                //System.out.println("Setting offline time");
                setTime(project.getProjectTime());
                return false;
            }
            else
            {
                //System.out.println("Project time is: " + project.getProjectTimeAsString());
                //if(object instanceof StsProject)
                //    changePolicies(project.getProjectTime());
            }
        }
        return super.viewObjectRepaint(source, object);
    }

    public boolean changePolicies(long time)
    {
        if(traces == null || chart == null)
            return false;
        // Adjust the number of points shown in the times series. Even if no new data has been detected.
        // Only want to look at the most recent (user specified) duration;
        int numPts = 0;
        for(int i=0; i<traces.length; i++)
        {
           numPts = collectors[i].changePolicy(chart, time);
           if(numPts == 0)
               break;
        }
        //System.out.println("At (" + time + ") set the number of points in view to " + numPts);
        //setCurrentTraceToAmplitudeOrPressure();
        //reconfigureView();
        //chart.updateUI();
        //chart.setRequestedRepaint(true);
        return true;
    }

    public void setTime(long v_)
    {
        if(chart == null) return;
        setTitle("Realtime Series Plot (" + plotDuration + " seconds) of: "  + ((StsSensor)selectedObjects[0]).toString() + "     Active Attribute (" +  chart.getCurrentTrace().getName() + ")     Last update: " + project.getTimeDateFormat().format(new Date(project.getProjectTime())));
        super.setTime(v_);
    }

    public void restoreWindowState()
    {
    	super.restoreWindowState();
        createDataCollectors();
        emptyPlot = false;
        changePolicies(project.getProjectTime());
        super.reconfigureView();
        setCurrentTraceToAmplitudeOrPressure();
        chart.updateUI();
        chart.setRequestedRepaint(true);
    }
}

class StsStaticDataCollector extends AStaticDataCollector
{
    StsSensor sensor = null;
    int duration;
    long lastTime = 0;

    public StsStaticDataCollector(ITrace2D trace, StsSensor tSensor, int plotDuration)
    {
        super(trace);
        sensor = tSensor;
        duration = plotDuration;
    }

    public void collectData() {}
    public boolean collectData(Chart2D chart, long projectTime)
    {
        long processedTime = (long)m_trace.getMaxX();
    	long[] times = sensor.getNewTimeData(processedTime);

        // No new data found
        if((times == null) || (times.length == 0))
            return false;        // No new data found

        // Expand number of points shown to accomadate new data, even if no new data.
        int numPts = changePolicy(chart, projectTime);
        //if(numPts > 0)
        //    lastTime = projectTime;
        //System.out.println("3+ Collecting " + numPts + " points.");
        
        // Load new data into trace.
        chart.setCurrentTrace(m_trace);
        float[] data = sensor.getNewData(m_trace.getLabel(), processedTime);
        if(data != null)
    	{
            //System.out.println("Number of data collected=" + data.length);
    		for(int j=0; j<data.length; j++)
    			m_trace.addPoint(times[j], data[j]);
    	}
        return true;
    }


    public int changePolicy(Chart2D chart1, long projectTime)
    {
    	if(chart1 == null || sensor == null || lastTime > projectTime)
            return -1;
        lastTime = projectTime;
        chart1.setCurrentTrace(m_trace);
        int numPointsShown = sensor.getNumPoints(duration, projectTime);   // time interval may be irregular so must call each time
		if(numPointsShown == 0)
            return 0;
        //System.out.println("Increasing the number of points shown to " + numPointsShown + " at time: " + projectTime);
        ((Trace2DLtd)m_trace).setMaxSize(numPointsShown);
        return numPointsShown;
        //System.out.println("Current trace=" + chart1.getCurrentTrace().getName());
    }
}