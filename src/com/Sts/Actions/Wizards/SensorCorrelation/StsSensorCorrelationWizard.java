
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SensorCorrelation;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Utilities.StsMath;
import com.Sts.Utilities.StsParameters;
import com.Sts.UI.StsMessage;
import com.Sts.UI.StsYesNoDialog;

import java.awt.event.*;

public class StsSensorCorrelationWizard extends StsWizard implements ActionListener
{
    public StsDynamicSensor selectedSensor = null;
    public StsSelectSensors selectSensors = null;
    public Object[] selectedSensors = new Object[0];
    public StsSensor primeSensor = null;
    public StsCorrelateSensors correlateSensors = null;
    public byte methodType = AVERAGE;
    public int timeRange = 1; // time range to analyze in seconds
    public boolean[] liveSensors;

    public static final byte INSTANCE = 0;
    public static final byte MIN = 1;
    public static final byte MAX = 2;
    public static final byte AVERAGE = 3;
    public static final byte SUM = 4;
    public static final String[] methodList = {"Instance", "Minimum","Maximum","Average","Sum"};

    public StsSensorCorrelationWizard(StsActionManager actionManager)
    {
        super(actionManager, 400, 470);
        addSteps
        (
            new StsWizardStep[]
            {
            	selectSensors = new StsSelectSensors(this),
            	correlateSensors = new StsCorrelateSensors(this)
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

    public void setPrimeSensor(StsSensor sensor)
    {
        if(sensor == null) return;
        primeSensor = (StsSensor)sensor;
    }

    public StsSensor getPrimeSensor()
    {
        return primeSensor;
    }

    public boolean saveClusters()
    {
        return true;
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
        return super.end();
    }

    public void previous()
    {
       gotoPreviousStep();
    }

    public void next()
    {
    	gotoNextStep();
    }

    public void finish()
    {
        super.finish();
    }

    public void cancel()
    {
    	super.cancel();
    }

    public String getCorrelationMethod() { return methodList[methodType]; }
    public void setCorrelationMethod(String att)
    {
    	for(int i=0; i<methodList.length; i++)
    	{
    		if(methodList[i].equalsIgnoreCase(att))
    		{
    			methodType = (byte)i;
    			return;
    		}
    	}
    }
    public int getTimeRange() { return timeRange; }
    public void setTimeRange(int range) { timeRange = range; }

    static public StsDynamicSensorClass getSensorClass(StsModel model)
    {
        return (StsDynamicSensorClass)model.getStsClass(StsDynamicSensor.class);
    }

    public boolean analyze()
    {
        return true;
    }

    public boolean saveAttributes()
    {
    	Main.logUsage();
    	for(int j=0; j<selectedSensors.length; j++)
    	{
            StsSensor sensor = (StsSensor)selectedSensors[j];
            StsTimeCurve curve = null;
            String[] aNames = new String[sensor.getNTimeCurves()];
            boolean foundValues = false;

            int numCurves = sensor.getNTimeCurves();
            int numPrimeValues = primeSensor.getNumValues();
            float[][] attributes = new float[numCurves][numPrimeValues];
            for(int ii=0; ii<numPrimeValues; ii++)
            {
                for(int jj=0; jj<numCurves; jj++)
                    attributes[jj][ii] = StsParameters.nullValue;
            }
            long[] stimes = sensor.getTimeCurve(0).getTimeVectorLongs();
            for(int i=0; i<stimes.length; i++)
    	    {
                long startTime = stimes[i] - (getTimeRange()*1000)/2;
                long endTime = stimes[i] + (getTimeRange()*1000)/2;
                int[] idxRange = sensor.getTimeCurve(0).getTimeLongVector().getIndicesInValueRange(startTime + 1, endTime + 1);

                int primeIndex = primeSensor.getTimeIndexGE(startTime);
                for(int k=0; k<numCurves; k++)
                {
                    foundValues = true;
                    curve = sensor.getTimeCurve(k);
                    aNames[k] = curve.getName();
                    float[] values = sensor.getDataInRange(aNames[k], idxRange);
                    switch(methodType)
                    {
                            case INSTANCE:
                                attributes[k][primeIndex] = values[0];
                                aNames[k] = aNames[k] + "_Inst";
                                break;
                            case MIN:
                                attributes[k][primeIndex] = StsMath.min(values);
                                aNames[k] = aNames[k] + "_Min";
                                break;
                            case MAX:
                                attributes[k][primeIndex] = StsMath.max(values);
                                aNames[k] = aNames[k] + "_Max";
                                break;
                            case AVERAGE:
                                attributes[k][primeIndex] = StsMath.average(values, StsParameters.nullValue);
                                aNames[k] = aNames[k] + "_Avg";
                                break;
                            case SUM:
                                for(int kk=0; kk<values.length; kk++)
                                    attributes[k][primeIndex] = attributes[k][primeIndex] + values[kk];
                                aNames[k] = aNames[k] + "_Sum";
                                break;
                            default:
                                break;
                    }
                }
    		}
            // Add attributes to the prime sensor
            if(foundValues)
            {
                correlateSensors.panel.setMessage("Data found in " + sensor.getName() + " attributes added to " + primeSensor.getName());
                for(int jj=0; jj<attributes.length; jj++)
                {
                    if(primeSensor.getTimeCurve(aNames[jj]) != null)
                    {
                        if(!StsYesNoDialog.questionValue(frame, "Attribute (" + aNames[jj] + ") already exists. Overwrite?"))
                           continue;
                    }
                    primeSensor.setAttribute(attributes[jj], aNames[jj]);
                    primeSensor.saveAttribute();
                }
            }
            else
                correlateSensors.panel.setMessage("No correlating data found in " + sensor.getName());
    	}
    	Main.logUsage();
    	return true;
    }

    public boolean checkPrerequisites()
    {

    	StsSensorClass sensorClass = (StsSensorClass)model.getStsClass("com.Sts.DBTypes.StsSensor");
    	/*
        int numSensors = sensorClass.getSensors().length;
        StsStaticSensorClass staticClass = (StsStaticSensorClass)model.getStsClass("com.Sts.DBTypes.StsStaticSensor");
    	numSensors = numSensors + staticClass.getSensors().length;
        */
        StsObject[] objects = sensorClass.getSensors();
        int numSensors = objects.length;
    	if(numSensors == 0)
    	{
        	reasonForFailure = "No loaded and visible (enabled) sensors. Must have at least two sensor that has been loaded and is visible in the 3D view.";
    		return false;
    	}
    	if(numSensors < 2)
        {
           	reasonForFailure = "Need at least 2 sensors loaded to do comparison.";
        	return false;
        }
        if(sensorClass.getDynamicSensors().length < 1)
        {
            reasonForFailure = "Need at least 1 dynamic sensor to do comparison.";
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
        StsSensorCorrelationWizard clusterWizard = new StsSensorCorrelationWizard(actionManager);
        clusterWizard.start();
    }
}