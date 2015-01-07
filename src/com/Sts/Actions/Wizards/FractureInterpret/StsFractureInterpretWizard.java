
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.FractureInterpret;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

import java.awt.event.*;
import java.util.*;

public class StsFractureInterpretWizard extends StsWizard implements ActionListener
{
    public StsSensor selectedSensor = null;
    public StsSelectSensors selectSensors = null;
    public Object[] selectedSensors = new Object[0];
    public StsFractureConstruction constructFracture = null;
    public StsWell selectedWell = null;
    
    transient StsSensor trackSensor = null;
    transient int sensorIdx = 0;
    transient boolean[] liveSensors;
    transient int totalSensors = 0;
    
    int[][] clusters = null;
    long[] startTimes = null;

    public StsFractureInterpretWizard(StsActionManager actionManager)
    {
        super(actionManager, 500, 400);
        addSteps
        (
            new StsWizardStep[]
            {
            	selectSensors = new StsSelectSensors(this),
            	constructFracture = new StsFractureConstruction(this)
            }
        );
    }
    
    public boolean start()
    {
        disableFinish();
        dialog.addActionListener(this);
        System.runFinalization();
        System.gc();
        dialog.setTitle("Interpret Fractures");
        return super.start();
    }
    
    public void actionPerformed(ActionEvent e)
    {
    	if(e.getID() == WindowEvent.WINDOW_CLOSING)
    		cancel();
    }
    
    public void setSelectedSensor(StsSensor sensor)
    {
        if(sensor == null) return;
        selectedSensor = sensor;
    }
    
    public boolean checkPrerequisites()
    {
    	StsDynamicSensorClass sensorClass = (StsDynamicSensorClass)model.getStsClass("com.Sts.DBTypes.StsDynamicSensor");
    	StsObject[] sensors = sensorClass.getSensors();
    	if(sensors == null)
    	{
        	reasonForFailure = "No loaded and visible (enabled) sensors. Must have at least one sensor that has been loaded and is visible in the 3D view.";
    		return false;
    	}
    	return true;
    } 

    //
    // ToDo: Need to change to allow the user to selected an attribute which is the desired cluster nums
    // They can now save clustering from any step into a new attribute, so selection of cluster attribute is needed.
    //
    public int[] getUniqueClusterNums()
    {
    	int[] clusterNums = null;
    	for(int i=0; i<selectedSensors.length; i++)
    	{
    		if(i == 0)
    			clusterNums = ((StsSensor)selectedSensors[i]).getClusterNums();
    		else
    		{
    			int[] nums = ((StsSensor)selectedSensors[i]).getClusterNums();
    			for(int j=0; j<nums.length; j++)
    			{
    				if(StsMath.arrayContains(clusterNums, nums[j]))
    						continue;
    				else
    					StsMath.arrayAddElement(clusterNums, nums[j]);
    			}
    		}
    	}
    	return clusterNums;
    }
    
    public StsSensor getSelectedSensor()
    {
        return selectedSensor;
    }
    
    public Object[] getSelectedSensors()
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
			selectedSensors = selectSensors.getSelectedSensors();
			enableFinish();
    	}
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
}
