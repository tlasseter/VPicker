
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.FractureTrack;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.event.*;
import java.util.*;

public class StsFractureTrackWizard extends StsWizard implements ActionListener
{
    public StsDynamicSensor selectedSensor = null;
    public StsSelectSensors selectSensors = null;
    public Object[] selectedSensors = new Object[0];
    public StsDefineTracker defineTracker = null;
    public StsWell selectedWell = null;
    
    transient StsDynamicSensor trackSensor = null;
    transient int sensorIdx = 0;
    transient boolean[] liveSensors;
    transient int totalSensors = 0;
    
    float distanceLimit = -1.0f;
    long timeLimit = -1l;
    boolean ignoreVertical = false;
    boolean showNotTracked = false;
    int[][] clusters = null;
    long[] startTimes = null;
    
    Thread trackFractureThread = null;
    
    static String NONE = "none";

    public StsFractureTrackWizard(StsActionManager actionManager)
    {
        super(actionManager, 600, 400);
        addSteps
        (
            new StsWizardStep[]
            {
            	selectSensors = new StsSelectSensors(this),
            	defineTracker = new StsDefineTracker(this)
            }
        );
    }
    
    public boolean start()
    {
        disableFinish();
        dialog.addActionListener(this);
        System.runFinalization();
        System.gc();
        dialog.setTitle("Track Fractures");
        return super.start();
    }
    
    public void actionPerformed(ActionEvent e)
    {
    	if(e.getID() == WindowEvent.WINDOW_CLOSING)
    		cancel();
    }
    
    public void setSelectedSensor(StsDynamicSensor sensor)
    {
        if(sensor == null) return;
        selectedSensor = sensor;
    }

    public StsDynamicSensor getSelectedSensor()
    {
        return selectedSensor;
    }
    
    public Object getSelectedSensors()
    {
        return selectedSensors;
    }
    
    public void setSelectedSensors(Object sensor)
    {
        if(sensor == null) return;
        selectedSensors = selectSensors.getSelectedSensors();
    }
    
    public void setSelectedWell(StsWell well)
    {
        if(well == null) return;
        selectedWell = well;
    }

    public StsWell getSelectedWell()
    {
        return selectedWell;
    }
    
    public boolean end()
    {
    	restoreSensors();
        // setClustering(false);
        return super.end();
    }
    
    public void previous()
    {
       gotoPreviousStep();
    }

    public void next()
    {
    	if(currentStep == selectSensors)
    	{
    		setClustering(true);
    		enableSensors();
    	}
    	gotoNextStep();
    }
    
    public float getDistanceLimit() { return distanceLimit; }
    public void setDistanceLimit(float limit) {	distanceLimit = limit; }
    public float getTimeLimit() { return timeLimit; }
    public void setTimeLimit(long time) {	timeLimit = time; }    
    
    public void finish()
    {    	
        super.finish();
    }

    public void cancel()
    {
    	restoreSensors();    	
        setClustering(false);    	
    	super.cancel();
    }
    
    // Enable only the selected Sensor
    public void enableSensors()
    {
    	StsObject[] sensors = (StsObject[])model.getObjectList(StsSensor.class);
    	totalSensors = sensors.length;
    	if(sensors == null) 
    		return;
    	liveSensors = new boolean[sensors.length];
    	for(int i=0; i<sensors.length; i++)
    	{
    		liveSensors[i] = ((StsSensor)sensors[i]).getIsVisible();
    		((StsSensor)sensors[i]).setIsVisible(false);
    	}
    	// Turn on selected sensors
    	for(int i=0; i<selectedSensors.length; i++)
    		((StsSensor)selectedSensors[i]).setIsVisible(true);
    }
    
    // Restore the visibility of the sensors before wizard
    public void restoreSensors()
    {
    	if(liveSensors == null) return;
    	
    	StsObject[] sensors = (StsObject[])model.getObjectList(StsSensor.class);
    	if(sensors == null) 
    		return;
    	for(int i=0; i<totalSensors; i++)
    		((StsSensor)sensors[i]).setIsVisible(liveSensors[i]);  	
    }
    
    public void setClustering(boolean val)
    {
    	for(int i=0; i<selectedSensors.length; i++)
    	{
    		((StsSensor)selectedSensors[i]).setClustering(val);
    		((StsSensor)selectedSensors[i]).setClusters(null);
    	}
    }
    
    private boolean verifyXYZ()
    {
    	// Verify the XYZ information exists if clustering by distance.    	
    	for(int i=0; i<selectedSensors.length; i++)
    	{
    		StsSensor sensor = (StsSensor)selectedSensors[i];
        	if((sensor instanceof StsStaticSensor) && (distanceLimit != -1))
        		return false;    		
    	}
    	return true;
    }
    
    public void setShowNotTracked(boolean val) 
    {
    	if(showNotTracked == val) return;
    	showNotTracked = val;
    	runTrackFracture();
    }
    public boolean getShowNotTracked() { return showNotTracked; }
    public boolean exportView()
    {
    	String filename = "TrackedSensors" + System.currentTimeMillis() + ".csv";
    	StsSensorClass sensorClass = ((StsSensor)selectedSensors[0]).getSensorClass();
    	int nPoints = sensorClass.exportView(filename, StsSensorClass.EVENTS);
    	defineTracker.panel.setMessage("Output " + nPoints + " points to file: " + filename);
    	return true;
    }
    
    public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
    {
        StsGLPanel3d glPanel3d = (StsGLPanel3d)glPanel;
        trackSensor = (StsDynamicSensor) StsJOGLPick.pickVisibleClass3d(glPanel3d, StsDynamicSensor.class, StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_FIRST);
        if (trackSensor != null)
        {                        	 
            StsPickItem items = StsJOGLPick.pickItems[0];
            if(items.names.length < 2)
            	return false;
            sensorIdx = items.names[1];
            trackSensor.setPicked(sensorIdx);
            trackSensor.logMessage(trackSensor.toString(sensorIdx));
        }
        else
        	return false;
        
        if( mouse.getCurrentButton() != StsMouse.LEFT ) return true;
        int buttonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

        if(buttonState == StsMouse.RELEASED)
        {
        	setClustering(true);
        	initClusters();
        	trackSensor.setClusterValue(sensorIdx, 2);
        	
            StsView currentView = glPanel3d.getView();
            if(currentView instanceof StsView3d)
            {
            	runTrackFracture();
            }
            else if(currentView instanceof StsViewCursor)
                ;
            else
                return true;
        }
        model.viewObjectRepaint(this, trackSensor);
        return true;
    }
    
    public void clearTracking()
    {
    	initClusters();
    }

    public void initClusters()
    {
    	clusters = new int[selectedSensors.length][];
    	startTimes = new long[selectedSensors.length];
    	for(int i=0; i<selectedSensors.length; i++)
    	{
    		startTimes[i] = ((StsSensor)selectedSensors[i]).getStartTime();
    		clusters[i] = new int[((StsSensor)selectedSensors[i]).getNumValues()];
    		for(int j=0; j<clusters[i].length; j++)
    			clusters[i][j] = 1;
    		((StsSensor)selectedSensors[i]).setClusters(clusters[i]);
    		model.viewObjectChanged(this, ((StsSensor)selectedSensors[i]));
    	} 
    	Arrays.sort(startTimes);
    	return;
    }
    
    public void runTrackFracture()
    {
        Runnable runTracker = new Runnable()
        {
            public void run()
            {
                trackFractureProcess();
            }
        };

        trackFractureThread = new Thread(runTracker);
        trackFractureThread.start();
    }

    private void trackFractureProcess()
    {
    	float[] nextPt = null;
    	int processIdx = 0;
    	StsDynamicSensor currentSensor = null;
    	
    	int notTrackedColorIdx = 1;
    	if(getShowNotTracked())
    		notTrackedColorIdx = -1;
    	
    	int[] distanceRange = defineTracker.panel.getDistanceRange();
    	int[] timeRange = defineTracker.panel.getTimeRange();
    	int[] azimuthRange = defineTracker.panel.getAzimuthRange();
    	
        try
        {
        	// Verify the XYZ information exists if clustering by distance.
        	if(!verifyXYZ())
        	{
        		new StsMessage(frame, StsMessage.ERROR, "Selected sensors are static so distance criteria cannot be verified. Aborting tracking.");
        		return;    		
        	}
        	if(startTimes == null)
            {
                return;
            }
        	for(int i=0; i<startTimes.length; i++)
        	{
        		for(int j=0; j<selectedSensors.length; j++)
        		{
        			StsDynamicSensor sensor = ((StsDynamicSensor)selectedSensors[j]);
        			if(startTimes[i] == sensor.getStartTime())
        			{
        				processIdx = i;
        				currentSensor = (StsDynamicSensor)selectedSensors[processIdx];
        				break;
        			}
        		}
        		
        		long[] times = currentSensor.getTimeCurve(0).getTimeVectorLongs();
        		
        		long pickTime = trackSensor.getTimeCurve(0).getTimeVectorLongs()[sensorIdx];
        		float[] pickPt = trackSensor.getXYZ(sensorIdx);

        		for(int j=0; j<clusters[processIdx].length; j++)
        		{
        			nextPt = currentSensor.getXYZ(j);
        			// Is it inside the time range?
        			if(((times[j]-pickTime) > timeRange[0]*1000) && ((times[j]-pickTime) < timeRange[1]*1000))
        			{
        				// is it inside the distance limit?
        				if((StsMath.distance(pickPt, nextPt, 3) > distanceRange[0]) &&
        						(StsMath.distance(pickPt, nextPt, 3) < distanceRange[1]))
        				{
        					clusters[processIdx][j] = 3;
        					pickTime = times[j];
        					pickPt = nextPt;
        					continue;
        				}
        				else
        					clusters[processIdx][j] = notTrackedColorIdx;
        			}
    				else
    					clusters[processIdx][j] = notTrackedColorIdx;
        		}
                currentSensor.setClusters(clusters[processIdx]);
				model.viewObjectRepaint(this, currentSensor);
                
        		trackSensor.setClusterValue(sensorIdx, 2);
        		model.viewObjectRepaint(this, trackSensor);
        	}
        }
        catch(Exception e)
        {
            StsException.outputException("StsFractureTrackWizard.trackFractureProcess() failed.", e, StsException.WARNING);
        }
        finally
        {
            ;
        }
    }
}
