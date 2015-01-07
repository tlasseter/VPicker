
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SensorPartition;

import com.Sts.Actions.Import.StsSensorKeywordIO;
import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.DBTypes.*;
import com.Sts.MVC.StsActionManager;
import com.Sts.MVC.StsMessageFiles;
import com.Sts.MVC.StsModel;
import com.Sts.MVC.StsProject;
import com.Sts.UI.Progress.StsProgressPanel;
import com.Sts.UI.StsMessage;
import com.Sts.UI.StsOkCancelDialog;
import com.Sts.Utilities.StsFloatVector;
import com.Sts.Utilities.StsLongVector;
import com.Sts.Utilities.StsMath;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;

public class StsSensorPartitionWizard extends StsWizard implements ActionListener
{
    public StsSensor selectedSensor = null;
    public StsSelectSensors selectSensors = null;
    public Object[] selectedSensors = new Object[0];
    public StsDefineRanges defineRanges = null;

    transient boolean[] liveSensors;
    transient int totalSensors = 0;

    public StsSensorPartitionWizard(StsActionManager actionManager)
    {
        super(actionManager, 500, 700);
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
        dialog.setTitle("Sensor Partitioning");
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
        // Create the partitioned sensor objects.
        new StsMessage(this.frame, StsMessage.INFO, "Creating " + defineRanges.panel.getNumPartitions() + " new sensors...", true);

        int numPartitions = defineRanges.panel.getNumPartitions();
        long[] starts =  defineRanges.panel.getStartTimes();
        long[] ends =  defineRanges.panel.getEndTimes();
        int[] colorI =  defineRanges.panel.getColorIndices();
        StsSensor sensor = (StsSensor)selectedSensors[0];
        StsObjectRefList curves = sensor.getTimeCurves();
        for(int i=0; i<numPartitions; i++)
        {
            int[] indices = sensor.getTimeCurve(0).getTimeLongVector().getIndicesInValueRange(starts[i], ends[i]);

            if(indices == null)
                continue;
            if(indices.length == 0)
                continue;

            long[] times = new long[indices[1]-indices[0]+1];
            float[] values = null;

            String binaryFilename = "sensor.bin." + sensor.getName()+"_"+i + ".TIME.0";
            StsTimeVector timeVector = new StsTimeVector(null, binaryFilename, "TIME", 0, -1, -1);

            System.arraycopy(sensor.getTimeCurve(0).getTimeLongVector().getValues(),indices[0],times,0,indices[1]-indices[0]+1);
            timeVector.setValues(new StsLongVector(times, times.length, 1));
            timeVector.setNameAndType(sensor.getTimeCurve(0).getTimeVector().getName(), binaryFilename);
            timeVector.checkWriteBinaryFile(model.getProject().getBinaryFullDirString(), false);

            StsTimeCurve[] outCurves = new StsTimeCurve[curves.getSize()];
            if(sensor instanceof StsDynamicSensor)
            {
                StsDynamicSensor dsensor = new StsDynamicSensor(null, sensor.getName() + "_" + i);
                for(int j=0; j<curves.getSize(); j++)
                {
                    StsTimeCurve curve = (StsTimeCurve)curves.getElement(j);
                    binaryFilename = "sensor.bin." + dsensor.getName() + "." + curve.getName() + ".0";
                    StsLogVector curveVector = new StsLogVector(null, binaryFilename, curve.getName(), 0, -1);
                    curveVector.setOrigin(curve.getValueVector().getOrigin());
                    curveVector.setNameAndType(curve.getName(), binaryFilename);

                    values = new float[indices[1]-indices[0]+1];
    	            System.arraycopy(sensor.getTimeCurve(curve.getName()).getValuesVectorFloats(),indices[0],values,0,indices[1]-indices[0]+1);                    
                    curveVector.setValues(values);

                    curveVector.checkWriteBinaryFile(model.getProject().getBinaryFullDirString(), false, true);
                    outCurves[j] = StsTimeCurve.constructTimeCurve(timeVector,curveVector,0);
                }
                dsensor.addTimeCurves(outCurves, false);
                dsensor.setType(StsSensor.DYNAMIC);
                dsensor.addToProject(StsProject.TD_DEPTH);
                dsensor.setColorFromString(StsColor.colorNames32[colorI[i]]);
            }
            else
            {
                StsStaticSensor ssensor = new StsStaticSensor(null, sensor.getName() + "_" + i);
                for(int j=0; j<curves.getSize(); j++)
                {
                    StsTimeCurve curve = (StsTimeCurve)curves.getElement(j);
                    binaryFilename = "sensor.bin." + ssensor.getName() + "." + curve.getName() + ".0";
                    StsLogVector curveVector = new StsLogVector(null, binaryFilename, curve.getName(), 0, -1);
                    curveVector.setOrigin(curve.getValueVector().getOrigin());
                    curveVector.setNameAndType(curve.getName(), binaryFilename);

                    values = new float[indices[1]-indices[0]+1];
    	            System.arraycopy(sensor.getTimeCurve(curve.getName()).getValuesVectorFloats(),indices[0],values,0,indices[1]-indices[0]+1);
                    curveVector.setValues(values);

                    curveVector.checkWriteBinaryFile(model.getProject().getBinaryFullDirString(), false, true);
                    outCurves[j] = StsTimeCurve.constructTimeCurve(timeVector,curveVector,0);
                }
                ssensor.addTimeCurves(outCurves, false);
                ssensor.setType(StsSensor.STATIC);                
                ssensor.addToProject(StsProject.TD_DEPTH);
                ssensor.setColorFromString(StsColor.colorNames32[colorI[i]]);
                ssensor.setIsVisible(false);
            }
        }

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

    public boolean checkPrerequisites()
    {
    	StsSensorClass sensorClass = (StsSensorClass)model.getStsClass("com.Sts.DBTypes.StsSensor");
    	StsObject[] sensors = sensorClass.getSensors();
    	if(sensors.length == 0)
    	{
        	reasonForFailure = "No loaded and visible (enabled) sensors. Must have at least one sensor that has been loaded it must be enabled.";
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
        StsSensorPartitionWizard clusterWizard = new StsSensorPartitionWizard(actionManager);
        clusterWizard.start();
    }
}