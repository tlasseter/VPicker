
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SensorCompare;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;

import java.awt.event.*;

public class StsSensorCompareWizard extends StsWizard implements ActionListener
{
    public StsDynamicSensor selectedSensor = null;
    public StsSelectSensors selectSensors = null;
    public Object[] selectedSensors = new Object[0];
    public StsDynamicSensor primeSensor = null;
    public boolean interactiveMode = false;
    public StsCompareSensors compareSensors = null;
    public byte thresholdType = NONE;
    public float threshold = 0.0f;
    public boolean[] liveSensors;
    public String originalSymbol;
    
    public static final byte NONE = 0;
    public static final byte TIME = 1;
    public static final byte DIST = 2;
    public static final byte XYDIST = 3;
    public static final byte AMPL = 4;
    public static final String[] attList = {"None","Time","Distance","XY Distance","Amplitude"};
    
    public StsSensorCompareWizard(StsActionManager actionManager)
    {
        super(actionManager, 300, 470);
        addSteps
        (
            new StsWizardStep[]
            {
            	selectSensors = new StsSelectSensors(this),
            	compareSensors = new StsCompareSensors(this)
            }
        );
    }
    
    public boolean start()
    {
        disableFinish();
        dialog.addActionListener(this);
        System.runFinalization();
        System.gc();
        dialog.setTitle("Sensor Comparison");
        return super.start();
    }
    
    public void actionPerformed(ActionEvent e)
    {
    	if(e.getID() == WindowEvent.WINDOW_CLOSING)
    		cancel();
    }
    
    public void setInteractiveMode(boolean interactive)
    {
        interactiveMode = interactive;
    }

    public boolean getInteractiveMode()
    {
        return interactiveMode;
    }
    public void setPrimeSensor(StsSensor sensor)
    {
        if(sensor == null) return;
        primeSensor = (StsDynamicSensor)sensor;
    }

    public StsSensor getPrimeSensor()
    {
        return primeSensor;
    }

    public void saveClusters()
    {
        ;
    }
    public void setSelectedSensor(StsDynamicSensor sensor)
    {
        if(sensor == null) return;
        selectedSensor = sensor;
    }

    public StsSensor getSelectedSensor()
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

    public boolean end()
    {
        if(primeSensor != null)
        {
    	    restoreSensors();
            setClustering(false);
        }
        return super.end();
    }
    
    public void previous()
    {
       gotoPreviousStep();
    }

    public void next()
    {
    	if(currentStep == selectSensors)
    		enableSensors();
    	gotoNextStep();
    }

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
    	StsObject[] sensors = (StsObject[])model.getObjectList(StsDynamicSensor.class);
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
        
        originalSymbol = primeSensor.getSymbolString();
        primeSensor.setSymbolString(primeSensor.SYMBOL_TYPE_STRINGS[primeSensor.CUBE]);
    	primeSensor.setIsVisible(true);
    }
    
    // Restore the visibility of the sensors before wizard
    public void restoreSensors()
    {
    	if(liveSensors == null) return;
    	
    	StsObject[] sensors = (StsObject[])model.getObjectList(StsSensor.class);
    	if(sensors == null) 
    		return;
    	for(int i=0; i<sensors.length; i++)
    		((StsSensor)sensors[i]).setIsVisible(liveSensors[i]);

        primeSensor.setSymbolString(originalSymbol);
    }
    
    public float getThreshold() { return threshold; }
    public void setThreshold(float thres) { threshold = thres; }
    public String getThresholdAttribute()
    {
        return attList[thresholdType];
    }
    public void setThresholdAttribute(String att)
    {
    	for(int i=0; i<attList.length; i++)
    	{
    		if(attList[i].equalsIgnoreCase(att))
    		{
    			thresholdType = (byte)i;
    			compareSensors.panel.reinitializeThreshold();
    			return;
    		}
    	}
    }
    
    static public StsDynamicSensorClass getSensorClass(StsModel model)
    {
        return (StsDynamicSensorClass)model.getStsClass(StsDynamicSensor.class);
    }
    
    public void setClustering(boolean val)
    {
		    primeSensor.setClustering(val);
		    primeSensor.setClusters(null);
    	    for(int i=0; i<selectedSensors.length; i++)
    	    {
    		    ((StsSensor)selectedSensors[i]).setClustering(val);
    		    ((StsSensor)selectedSensors[i]).setClusters(null);
    	    }
    }
    
    public byte getComparisonType() 
    {
    	return compareSensors.panel.getComparisonType();
    }
    
    public float[] getAmplitudeValues(StsSensor sensor)
    {
		StsTimeCurve curve = sensor.getTimeCurve("Amplitude");
		if(curve == null)
			curve = sensor.getTimeCurve("Amp");
		if(curve == null)
			curve = sensor.getTimeCurve("Ampl");
		if(curve == null)
			curve = sensor.getTimeCurve("Energy");	
		if(curve == null)
			return null;
		else
			return curve.getValuesVectorFloats();
    }
    
    public float[] getAmplitudeRange(StsSensor sensor)
    {
		StsTimeCurve curve = sensor.getTimeCurve("Amplitude");
		if(curve == null)
			curve = sensor.getTimeCurve("Amp");
		if(curve == null)
			curve = sensor.getTimeCurve("Ampl");
		if(curve == null)
			curve = sensor.getTimeCurve("Energy");
		
		if(curve == null)
			return null;
		else
			return new float[] {curve.getMinValue(), curve.getMaxValue()};
    }
    
    public boolean analyze()
    {  		
    	Main.logUsage();
    	long[] ptimes = primeSensor.getTimeCurve(0).getTimeVectorLongs();
		float[] pamps = getAmplitudeValues(primeSensor);
		
		int[] ppassed = new int[primeSensor.getNumValues()];
		for(int j=0; j<primeSensor.getNumValues(); j++)
			ppassed[j] = j%32;
		primeSensor.setClusters(ppassed);
		primeSensor.setClustering(true);
		
		int[][] passed = new int[selectedSensors.length][];
		for(int j=0; j<selectedSensors.length; j++)
		{
			StsDynamicSensor sensor = (StsDynamicSensor)selectedSensors[j];
			passed[j] = new int[sensor.getNumValues()];
			for(int k=0; k<sensor.getNumValues(); k++)
			{
				passed[j][k] = -1;
			}
		}
    	for(int i=0; i<primeSensor.getNumValues(); i++)
    	{
    		StsPoint ppt = primeSensor.getZPoint(i);
    		for(int j=0; j<selectedSensors.length; j++)
    		{
    			StsDynamicSensor sensor = (StsDynamicSensor)selectedSensors[j];
    			sensor.setIsVisible(true);
    			long tDiff = Long.MAX_VALUE;
    			float distAmpDiff = Float.MAX_VALUE;
    			int idx = -1;

    			float[] amps = getAmplitudeValues(sensor);
    			long[] times = sensor.getTimeCurve(0).getTimeVectorLongs();
    			
    			for(int k=0; k<sensor.getNumValues(); k++)
    			{
    				StsPoint pt = sensor.getZPoint(k);
    				float ampDiff = 0.0f;
    				if((pamps != null) && (amps != null))
    					ampDiff = Math.abs(pamps[i]-amps[k]);
    				
    				if(!checkThreshold(ppt, pt, ampDiff, Math.abs(ptimes[i]-times[k])))
    					continue;
    				switch(getComparisonType())
    				{
    				case TIME:
    					if(tDiff >= Math.abs(ptimes[i]-times[k]))
    					{
    						tDiff = Math.abs(ptimes[i]-times[k]);
    						idx = k;
    					}
    					break;
    				case DIST:
    					if(distAmpDiff >= pt.distance(ppt))
    					{
    						distAmpDiff = pt.distance(ppt);
    						idx = k;
    					}
    					break;
    				case XYDIST:
    					if(distAmpDiff >= pt.horizontalDistance(ppt))
    					{
    						distAmpDiff = pt.horizontalDistance(ppt);
    						idx = k;
    					}    					
    					break;
    				case AMPL:
    					if((pamps == null) || (amps == null))
    						break;
    					if(distAmpDiff >= Math.abs(pamps[i]-amps[k]))
    					{
    						distAmpDiff = Math.abs(pamps[i]-amps[k]);
    						idx = k;
    					}    					
    					break;
    				default:
    					break;
    				}
    			}
    			if(idx >= 0) passed[j][idx] = i%32;
    		}
    	}
		for(int j=0; j<selectedSensors.length; j++)
		{
			StsDynamicSensor sensor = (StsDynamicSensor)selectedSensors[j];
			sensor.setClustering(true);
			sensor.setClusters(passed[j]);
            model.viewObjectRepaint(this, sensor);
		}
    	Main.logUsage();    	
    	return true;
    }
    
    public boolean checkPrerequisites()
    {
    	StsDynamicSensorClass sensorClass = (StsDynamicSensorClass)model.getStsClass("com.Sts.DBTypes.StsDynamicSensor");
    	StsObject[] sensors = sensorClass.getSensors();
    	if(sensors.length == 0)
    	{
        	reasonForFailure = "No loaded and visible (enabled) sensors. Must have at least one sensor that has been loaded and is visible in the 3D view.";
    		return false;
    	}
    	if(sensors.length < 2)
        {
           	reasonForFailure = "Need at least 2 sensors loaded to do comparison.";
        	return false;
        }
        if(model.getProject().isRealtime())
        {
            reasonForFailure =  "Cannot execute sensor related actions while real-time is running. Stop real-time, run this wizard and then re-start real-time.";
            return false;
        }
    	return true;
    }
    
    public boolean checkThreshold(StsPoint ppt, StsPoint pt, float ampDiff, float timeDiff)
    {
    	switch(thresholdType)
    	{
    	case TIME:
    		//System.out.println("Timediff=" + timeDiff + " Threshold=" + threshold);
    		if(timeDiff < threshold)
    			return true;
    		break;
    	case AMPL:
    		//System.out.println("Ampdiff=" + ampDiff + " Threshold=" + threshold);   		
    		if(ampDiff < threshold)
    			return true;   		
    		break;
    	case DIST:
    		//System.out.println("Distance=" + ppt.distance(pt) + " Threshold=" + threshold);    		
    		if(ppt.distance(pt) < threshold)
    			return true;
    		break;
    	case XYDIST:
    		//System.out.println("Distance=" + ppt.horizontalDistance(pt) + " Threshold=" + threshold);    		
    		if(ppt.horizontalDistance(pt) < threshold)
    			return true;
    		break;
    	case NONE:
    		return true;
    	default:
    		return true;
    	}
    	return false;
    }
    
    static public void main(String[] args)
    {
        StsModel model = StsModel.constructor();
        StsActionManager actionManager = new StsActionManager(model);
        StsSensorCompareWizard clusterWizard = new StsSensorCompareWizard(actionManager);
        clusterWizard.start();
    }
}
