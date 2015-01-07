
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.ClusterAnalysis;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.event.*;

public class StsClusterAnalysisWizard extends StsWizard implements ActionListener
{
    public StsDynamicSensor selectedSensor = null;
    public StsSelectSensors selectSensors = null;
    public Object[] selectedSensors = new Object[0];
    public StsDefineClustering defineClustering = null;
    StsProgressPanel progressPanel;
    
    float distanceLimit = -1.0f;
    int minimumClusterSize = 2;
    int timeLimit = -1; // seconds
    float amplitude = 0.0f;
    boolean relateTime = false;
    boolean relateDistance = false;
    boolean ignoreVertical = false;
    
    transient boolean[] liveSensors;
    transient int totalSensors = 0;
    
    static String NONE = "none";

    public StsClusterAnalysisWizard(StsActionManager actionManager)
    {
        super(actionManager, 400, 500);
        addSteps
        (
            new StsWizardStep[]
            {
            	selectSensors = new StsSelectSensors(this),
            	defineClustering = new StsDefineClustering(this)
            }
        );
    }
    
    public boolean start()
    {
        disableFinish();
        dialog.addActionListener(this);
        System.runFinalization();
        System.gc();
        dialog.setTitle("Cluster Analysis");
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
    
    // Enable only the selected Sensor
    public void enableSensors()
    {
    	StsObject[] sensors = (StsObject[])model.getObjectList(StsDynamicSensor.class);
    	totalSensors = sensors.length;
    	if(sensors == null) 
    		return;
    	liveSensors = new boolean[sensors.length];
    	for(int i=0; i<sensors.length; i++)
    	{
    		liveSensors[i] = ((StsDynamicSensor)sensors[i]).getIsVisible();
    		((StsDynamicSensor)sensors[i]).setIsVisible(false);
    	}
    	// Turn on selected sensors
    	for(int i=0; i<selectedSensors.length; i++)
    		((StsDynamicSensor)selectedSensors[i]).setIsVisible(true);
    }
    
    // Restore the visibility of the sensors before wizard
    public void restoreSensors()
    {
    	if(liveSensors == null) return;
    	
    	StsObject[] sensors = (StsObject[])model.getObjectList(StsDynamicSensor.class);
    	if(sensors == null) 
    		return;
    	for(int i=0; i<totalSensors; i++)
    		((StsDynamicSensor)sensors[i]).setIsVisible(liveSensors[i]);
    }
    
    public boolean end()
    {
    	restoreSensors();
        setClustering(false);    	
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
    		enableSensors();
    		setClustering(true);
    	}
    	gotoNextStep();
    }
    
    public int getTimeLimit() { return timeLimit; }
    public void setTimeLimit(int limit) {	timeLimit = limit; }
    public float getDistanceLimit() { return distanceLimit; }
    public void setDistanceLimit(float limit) {	distanceLimit = limit; }
    public int getMinimumClusterSize() { return minimumClusterSize; }
    public void setMinimumClusterSize(int min) { minimumClusterSize = min; }    
    public float getAmplitude() { return amplitude; }
    public void setAmplitude(float amp) {	amplitude = amp; }
    public boolean getRelateTime() { return relateTime; }
    public void setRelateTime(boolean relate) {	relateTime = relate; }
    public boolean getRelateDistance() { return relateDistance; }
    public void setRelateDistance(boolean relate) {	relateDistance = relate; }
    public boolean getIgnoreVertical() { return ignoreVertical; }
    public void setIgnoreVertical(boolean ignore) {	ignoreVertical = ignore; }
    
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
    
    public void setClustering(boolean val)
    {
    	for(int i=0; i<selectedSensors.length; i++)
    	{
    		((StsDynamicSensor)selectedSensors[i]).setClustering(val);
    		((StsDynamicSensor)selectedSensors[i]).setClusters(null);
    	}
    }
    
    public boolean exportView()
    {
    	String filename = "ClusteredSensors" + System.currentTimeMillis() + ".csv";
    	StsDynamicSensorClass sensorClass = ((StsDynamicSensor)selectedSensors[0]).getDynamicSensorClass();
    	int nPoints = sensorClass.exportView(filename, StsDynamicSensorClass.EVENTS);
    	defineClustering.panel.setMessage("Output " + nPoints + " points to file: " + filename);
    	
    	Main.logUsage();
    	return true;
    }
    
    public boolean saveClusters()
    {  
    	StsDynamicSensorClass sensorClass = ((StsDynamicSensor)selectedSensors[0]).getDynamicSensorClass();
    	String attName = "Clustered" + System.currentTimeMillis();
    	// Ask the user for a name....
    	StsNameInputDialog dialog = new StsNameInputDialog(frame, "Input Attribute Name", "Attribute Name:", attName, true);
        dialog.setVisible(true);
        if (dialog.wasCanceled())
        {
           return false;
        }
    	// Output the cluster as attribute.
        String userName = dialog.getUserName();
        if(userName != null)
        	attName = userName;
    	sensorClass.saveClusters(attName);
    	defineClustering.panel.setMessage("Clusters have been saved.");	   	
    	Main.logUsage();   	
    	return true;
    }
    
    public boolean checkPrerequisites()
    {
    	StsDynamicSensorClass sensorClass = (StsDynamicSensorClass)model.getStsClass("com.Sts.DBTypes.StsDynamicSensor");
    	StsObject[] sensors = sensorClass.getSensors();
    	if(sensors.length == 0)
    	{
        	reasonForFailure = "No loaded and visible (enabled) sensors. Must have at least one dynamic sensor that has been loaded and is visible in the 3D view.";
    		return false;
    	}
        if(model.getProject().isRealtime())
        {
            reasonForFailure =  "Cannot execute sensor related actions while real-time is running. Stop real-time, run this wizard and then re-start real-time.";
            return false;
        }
    	return true;
    }
    
    private boolean verifyAmplitudes()
    {
    	// Verify that the amplitude data is available if required.
    	for(int i=0; i<selectedSensors.length; i++)
    	{
    		StsDynamicSensor sensor = (StsDynamicSensor)selectedSensors[i];
    		if((getRelateDistance() || getRelateTime()) && (sensor.verifyAmplitude() == null))
    			return false;
    	}
    	return true;	
    }
    
    private boolean verifyXYZ()
    {
    	// All dynamic sensors have XYZ information
    	return true;
    }
    
    public void analyze(StsProgressPanel panel)
    {
		progressPanel = panel;
        Runnable runnable = new Runnable()
        {
            public void run()
            {
                analyze();
            }
        };
        Thread processThread = new Thread(runnable);
        processThread.start();
    }
    
    public boolean analyze()
    {
    	float[] currentPt = null, evalPt = null;
    	float dist = -1.0f;
    	int span = -1;
    	int distPass, timePass;
    	int clusterNumber = 0;
    	float[] amps = null;
    	
    	Main.logUsage();
    	// Verify that the amplitude data is available if required.
    	if(!verifyAmplitudes())
    	{
    		new StsMessage(frame, StsMessage.ERROR, "Unable to relate time or distance to amplitude. Selected sensor does not have amplitue values. Disable both checkboxs and run again.");
			return false;
    	}
    	
    	// Verify the XYZ information exists if clustering by distance.
    	if(!verifyXYZ())
    	{
    		new StsMessage(frame, StsMessage.ERROR, "Selected sensor is static so distance limit is invalid, set to -1 or pick a different sensor");
    		return false;    		
    	}
    	
    	int ttlPts = 0;
    	// If linearly relating distance or time to amplitude, load the amplitudes
    	progressPanel.setIntervalCount(selectedSensors.length);    	
    	for(int j=0; j<selectedSensors.length; j++)
    	{
        	dist = -1.0f;
        	span = -1; 
        	
    		StsDynamicSensor sensor = (StsDynamicSensor)selectedSensors[j];
    		if(getRelateDistance() || getRelateTime())
    		{
    			StsTimeCurve curve = sensor.getAmplitudeCurve();    			
    			amps = curve.getValuesVectorFloats();
    		}
    	
    		// Load time vector
    		long[] times = sensor.getTimeCurve(StsLogVector.types[StsLogVector.X]).getTimeVectorLongs();
    		int[] passed = new int[times.length];
    		int nPoints = 0;
			float dLimit = distanceLimit;
			int tLimit = timeLimit;
			
			progressPanel.resetProgressBar();
			progressPanel.setInterval(j);
			progressPanel.setMaximum(times.length);
    		
    		// Loop through points
    		for(int i=0; i<times.length; i++)
    		{
				passed[i] = -1;    			
    			currentPt = sensor.getXYZ(i);
				defineClustering.setAnalysisMessage("Found " + nPoints + " points out of " + times.length);
				progressPanel.progressBar.setDescriptionandValue("Processing event " + i + " of " + times.length, i);
				if (progressPanel.isCanceled())
				{
					progressPanel.finished();
					return false;
				}
				
				for(int kk=0; kk<selectedSensors.length; kk++)
				{
					StsDynamicSensor innerSensor = (StsDynamicSensor)selectedSensors[kk];
					long[] innerTimes = innerSensor.getTimeCurve(StsLogVector.types[StsLogVector.X]).getTimeVectorLongs();
					for(int k=0; k<innerTimes.length; k++)
					{
						if((k == i) && (j == kk))
							continue;
    				
						distPass = -1;
						timePass = -1;    				
						evalPt = innerSensor.getXYZ(k);

						// Compute Distance Limit
						if(getRelateDistance() && i>0)
							dLimit = (distanceLimit/amplitude)*amps[i];
    				
						// Compute Time Limit
						if(getRelateTime() && i>0)
							tLimit = (int)(((float)timeLimit/amplitude)*amps[i]);

						if(!ignoreVertical)
							dist = StsMath.distance(currentPt, evalPt, 3);
						else
							dist = StsMath.distance(currentPt, evalPt, 2);
						span = (int)Math.abs((innerTimes[k] - times[i]));
    		
    				//System.out.println("Distance limit is " + dLimit + " and the distance is " + dist);
    				//System.out.println("Time limit is " + tLimit + " and the time span is " + span);
						if(distPass == -1)
						{
							if(dLimit == -1) 
								distPass = 0;
							else
							{
								if(dist < dLimit)
									distPass = k;
							}
						}
						if(timePass == -1)
						{
							if(tLimit == -1) 
								timePass = 0;
							else
							{
								if(span < tLimit)
									timePass = k;
							}	
						}
						// If this point past both tests 
						if((timePass > -1) && (distPass > -1))
						{
							if(passed[i] == -1)
								passed[i] = 2;
							else
								passed[i]++;
    					
							if(passed[i] == minimumClusterSize)
								nPoints++;
						}
					}
				}
    		}
    		// Test minimum cluster size
    		
    		for(int i=0; i<passed.length; i++)
    		{
    			if(passed[i] < minimumClusterSize)
    				passed[i] = -1;
    		}
	
    		sensor.setClusters(passed);
    		model.viewObjectRepaint(this, sensor);
    		ttlPts += nPoints;
    		progressPanel.progressBar.setDescription("Completed processing");
    		progressPanel.finished();
    	}
    	defineClustering.setAnalysisMessage("Found " + ttlPts + " points from " + selectedSensors.length + " sensors.");
    	Main.logUsage();
    	return true;
    }
    
    static public void main(String[] args)
    {
        StsModel model = StsModel.constructor();
        StsActionManager actionManager = new StsActionManager(model);
        StsClusterAnalysisWizard clusterWizard = new StsClusterAnalysisWizard(actionManager);
        clusterWizard.start();
    }
}
