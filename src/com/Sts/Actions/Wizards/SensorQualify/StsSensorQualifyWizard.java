
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SensorQualify;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.event.*;

public class StsSensorQualifyWizard extends StsWizard implements ActionListener
{
    public StsDynamicSensor selectedSensor = null;
    public StsSelectSensors selectSensors = null;
    public Object[] selectedSensors = new Object[0];
    public StsDefineRanges defineRanges = null;
    StsProgressPanel progressPanel;
    
    transient boolean[] liveSensors;
    transient int totalSensors = 0;

    public StsSensorQualifyWizard(StsActionManager actionManager)
    {
        super(actionManager, 700, 700);
        addSteps
        (
            new StsWizardStep[]
            {
            	selectSensors = new StsSelectSensors(this),
            	defineRanges = new StsDefineRanges(this)
            }
        );
    }
    
    public boolean start()
    {
        disableFinish();
        dialog.addActionListener(this);
        System.runFinalization();
        System.gc();
        dialog.setTitle("Sensor Qualification");
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
    	liveSensors = new boolean[totalSensors];
    	for(int i=0; i<totalSensors; i++)
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
        if(selectedSensors.length > 0)
            setClustering(false);

        // Must update the static beans list of filters
        StsObject[] sensors = ((StsSensorClass)model.getCreateStsClass("com.Sts.DBTypes.StsSensor")).getSensors();
        for(int i=0; i<sensors.length; i++)
            ((StsSensor)sensors[i]).updateFilters();
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
        	// Determine if there are common attributes between the selected sensors.
    		if(getCommonAttributes() == null)
    		{
    			new StsMessage(frame, StsMessage.WARNING, "There are no common attributes to all selected sensors.\n Adjust selection.");
    			return;
    		}
    		enableSensors();
    		setClustering(true);
    	}
    	gotoNextStep();
    }
    
    public StsObject[] getCommonAttributes()
    {
    	boolean foundIt = false;
    	
    	StsObject[] returnCurves = ((StsDynamicSensor)selectedSensors[0]).getPropertyCurves();
    	StsObject[] commonCurves = ((StsDynamicSensor)selectedSensors[0]).getPropertyCurves();
		for(int j=0; j<commonCurves.length; j++)
		{
	    	for(int i=1; i<selectedSensors.length; i++)
	    	{
	    		StsObject[] curves = ((StsDynamicSensor)selectedSensors[i]).getPropertyCurves();
	    		for(int k=0; k<curves.length; k++)
	    		{
	    			if(curves[k].getName().equalsIgnoreCase(commonCurves[j].getName()))
	    			{
	    				foundIt = true;
	    				break;
	    			}
	    		}
	    		if(foundIt)
	    		{
	    			foundIt = false;
	    			continue;
	    		}
	    		else
	    		{
	    			StsMessageFiles.infoMessage("Removing " + commonCurves[j].getName() + " from attribute list because it does not exist in all selected sensors.");
	    			returnCurves = (StsObject[])StsMath.arrayDeleteElement(returnCurves, commonCurves[j]);
	    		}
	    	}
		}  
		
		// Verify that the return curves have a valid data range
		StsObject[] curves = returnCurves;
		for(int i=0; i<curves.length; i++)
		{
            ((StsTimeCurve)curves[i]).checkLoadVectors();
			if(((StsTimeCurve)curves[i]).getMinValue() == ((StsTimeCurve)curves[i]).getMaxValue())
			{
				StsMessageFiles.infoMessage("Removing " + curves[i].getName() + " from attribute list because it does not have a valid range.");
				returnCurves = (StsObject[])StsMath.arrayDeleteElement(returnCurves, curves[i]);
			}
		}
		return returnCurves;
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
    
    public void setClustering(boolean val)
    {
    	for(int i=0; i<selectedSensors.length; i++)
    	{
    		((StsDynamicSensor)selectedSensors[i]).setClustering(val);
    		((StsDynamicSensor)selectedSensors[i]).setClusters(null);
    	}
    }
    
    public boolean exportView(StsProgressPanel panel)
    {
    	// Export the values
    	String name = "QualifiedSensors" + System.currentTimeMillis();
		String filename = name + ".csv";
		StsDynamicSensorClass sensorClass = ((StsDynamicSensor)selectedSensors[0]).getDynamicSensorClass();
		int nPoints = sensorClass.exportView(filename, panel, StsDynamicSensorClass.EVENTS, false);
    	defineRanges.panel.setMessage("Output " + nPoints + " points to file: " + filename);
    	
    	Main.logUsage();
    	return true;
    }
    
    public boolean saveClusters()
    {  
    	StsDynamicSensorClass sensorClass = ((StsDynamicSensor)selectedSensors[0]).getDynamicSensorClass();
    	String attName = "Qualified" + System.currentTimeMillis();
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
    	defineRanges.panel.setMessage("Qualified events have been saved.");	   	
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
        if(model.getProject().isRealtime())
        {
            reasonForFailure =  "Cannot execute sensor related actions while real-time is running. Stop real-time, run this wizard and then re-start real-time.";
            return false;
        }
    	return true;
    }  
    
    static public void main(String[] args)
    {
        StsModel model = StsModel.constructor();
        StsActionManager actionManager = new StsActionManager(model);
        StsSensorQualifyWizard clusterWizard = new StsSensorQualifyWizard(actionManager);
        clusterWizard.start();
    }
}
