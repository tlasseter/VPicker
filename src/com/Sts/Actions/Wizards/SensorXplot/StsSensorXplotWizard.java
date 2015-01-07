
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SensorXplot;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;

import java.awt.event.*;

public class StsSensorXplotWizard extends StsWizard implements ActionListener
{
    public StsSelectSensors selectSensors = null;
    public StsSensorPlot sensorPlot = null;

    public StsDynamicSensor primeSensor = null;
    public boolean interactiveMode = false;
    
    public transient boolean[] liveSensors;
    transient int totalSensors = 0;
    
    public StsSensorXplotWizard(StsActionManager actionManager)
    {
        super(actionManager, 800, 800);
        addSteps
        (
            new StsWizardStep[]
            {
            	selectSensors = new StsSelectSensors(this),
            	sensorPlot = new StsSensorPlot(this)
            }
        );
    }
    
    public boolean start()
    {
        disableFinish();
    	Main.logUsage();        
        dialog.addActionListener(this);
        System.runFinalization();
        System.gc();
        dialog.setTitle("Sensor Crossplot");
        return super.start();
    }

    public void setPrimeSensor(StsDynamicSensor sensor)
    {
        if(sensor == null) return;
        primeSensor = sensor;
    }

    public StsDynamicSensor getPrimeSensor()
    {
        return primeSensor;
    }

    public int[] getClusterNums()
    {
    	return primeSensor.getClusterNums();
    }
    
    public boolean end()
    {    
    	restoreSensors(true);
        setClustering(false); 
    	Main.logUsage();        
        return super.end();
    }
    
    public void finish()
    {    	
        super.finish();
    }

    public void cancel()
    {
    	restoreSensors(true);    	
    	setClustering(false);  
    	Main.logUsage();    	
    	super.cancel();
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
        if(model.getProject().isRealtime())
        {
            reasonForFailure =  "Cannot execute sensor related actions while real-time is running. Stop real-time, run this wizard and then re-start real-time.";
            return false;
        }
    	return true;
    }
    
    public void actionPerformed(ActionEvent e)
    {
    	if(e.getID() == WindowEvent.WINDOW_CLOSING)
    		cancel();
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
    	for(int i=0; i<sensors.length; i++)
    	{
    		if(((StsDynamicSensor)sensors[i]).getClustering())
    			((StsDynamicSensor)sensors[i]).setIsVisible(true);
    	}  	
    	primeSensor.setIsVisible(true);
    }
    
    // Restore the visibility of the sensors before wizard
    public void restoreSensors(boolean clearClustering)
    {
    	if(liveSensors == null) return;
    	
    	StsObject[] sensors = (StsObject[])model.getObjectList(StsDynamicSensor.class);
    	if(sensors == null) 
    		return;
    	for(int i=0; i< totalSensors; i++)
    	{
    		if(clearClustering)
    		{
    			((StsDynamicSensor)sensors[i]).setClustering(false);
    			((StsDynamicSensor)sensors[i]).setClusters(null);
    		}
    		((StsDynamicSensor)sensors[i]).setIsVisible(liveSensors[i]);
    	}
    }
    
    public boolean exportView()
    {
    	String filename = "XplotSensors" + System.currentTimeMillis() + ".csv";
    	StsDynamicSensorClass sensorClass = primeSensor.getDynamicSensorClass();
    	int nPoints = sensorClass.exportView(filename, StsDynamicSensorClass.EVENTS);
    	sensorPlot.panel.updateMessage("Output " + nPoints + " points to file: " + filename);
    	
    	Main.logUsage();
    	return true;
    }
    
    public boolean saveClusters()
    {  
    	StsDynamicSensorClass sensorClass = primeSensor.getDynamicSensorClass();
    	String attName = "Xplot" + System.currentTimeMillis();
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
    	sensorPlot.panel.updateMessage("Crossplot events have been saved.");	   	
    	Main.logUsage();   	
    	return true;
    }
    
    public void setClustering(boolean val)
    {
    	if(primeSensor != null)
    	{
    		primeSensor.setClustering(val);
    		primeSensor.setClusters(null);
    	}
    }

    public StsColor getClusterStsColor(int index)
    {
        return sensorPlot.getClusterStsColor(index);
    }
    
    public void previous()
    {
    	if(currentStep == sensorPlot)
    	{
    		restoreSensors(false);
    		sensorPlot.panel.destroyChart();
    	}
    	gotoPreviousStep();
    }

    public void next()
    {
    	if(currentStep == selectSensors)
    		enableSensors();
    	gotoNextStep();
    }
    
    static public void main(String[] args)
    {
        StsModel model = StsModel.constructor();
        StsActionManager actionManager = new StsActionManager(model);
        StsSensorXplotWizard xplotWizard = new StsSensorXplotWizard(actionManager);
        xplotWizard.start();
    }
}
