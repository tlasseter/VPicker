
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.ProximityAnalysis;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.event.*;

public class StsProximityAnalysisWizard extends StsWizard implements ActionListener
{
    public StsDynamicSensor selectedSensor = null;
    public StsSelectSensors selectSensors = null;
    public Object[] selectedSensors = new Object[0];
    public StsDefineProximity defineProximity = null;
    public StsWell selectedWell = null;
    
    transient boolean[] liveSensors;
    transient int totalSensors = 0;
    
    float distanceLimit = -1.0f;
    boolean ignoreVertical = false;
    boolean insideLimits = true;
    float dist[][] = null;

    static String NONE = "none";

    public StsProximityAnalysisWizard(StsActionManager actionManager)
    {
        super(actionManager, 400, 470);
        addSteps
        (
            new StsWizardStep[]
            {
            	selectSensors = new StsSelectSensors(this),
            	defineProximity = new StsDefineProximity(this)
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
		if (e.getID() == WindowEvent.WINDOW_CLOSING)
			cancel();
	}

	public void setSelectedSensor(StsDynamicSensor sensor)
	{
		if (sensor == null)
			return;
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
		if (sensor == null)
			return;
		selectedSensors = selectSensors.getSelectedSensors();
	}

	public void setSelectedWell(StsWell well) 
	{
		if (well == null)
			return;
		selectedWell = well;
	}

	public StsWell getSelectedWell() 
	{
		return selectedWell;
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
		if (currentStep == selectSensors) 
		{
			setClustering(true);
			enableSensors();
		}
		gotoNextStep();
	}

	public float getDistanceLimit() 
	{
		return distanceLimit;
	}

	public void setDistanceLimit(float limit) 
	{
		distanceLimit = limit;
	}

	public boolean getIgnoreVertical() 
	{
		return ignoreVertical;
	}

	public void setIgnoreVertical(boolean ignore) 
	{
		ignoreVertical = ignore;
	}

	public boolean getInsideLimits() 
	{
		return insideLimits;
	}

	public void setInsideLimits(boolean inside) 
	{
		insideLimits = inside;
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
		StsObject[] sensors = (StsObject[]) model.getObjectList(StsDynamicSensor.class);
		totalSensors = sensors.length;
		if (sensors == null)
			return;
		liveSensors = new boolean[sensors.length];
		for (int i = 0; i < sensors.length; i++) 
		{
			liveSensors[i] = ((StsDynamicSensor) sensors[i]).getIsVisible();
			((StsDynamicSensor) sensors[i]).setIsVisible(false);
		}
		// Turn on selected sensors
		for (int i = 0; i < selectedSensors.length; i++)
			((StsDynamicSensor) selectedSensors[i]).setIsVisible(true);
	}

	// Restore the visibility of the sensors before wizard
	public void restoreSensors() 
	{
		if (liveSensors == null)
			return;

		StsObject[] sensors = (StsObject[]) model.getObjectList(StsDynamicSensor.class);
		if (sensors == null)
			return;
		for (int i = 0; i < totalSensors; i++)
			((StsDynamicSensor) sensors[i]).setIsVisible(liveSensors[i]);
	}

	public void setClustering(boolean val) 
	{
		for (int i = 0; i < selectedSensors.length; i++)
        {
            ((StsDynamicSensor) selectedSensors[i]).resetClusters();
			((StsDynamicSensor) selectedSensors[i]).setClustering(val);
		}
	}

	private boolean verifyXYZ() 
	{
		// All Dynamic sensors will have XYZ vectors.
		return true;
	}

	public boolean exportView() 
	{
		String filename = "ProximitySensors" + System.currentTimeMillis() + ".csv";
		StsDynamicSensorClass sensorClass = ((StsDynamicSensor) selectedSensors[0]).getDynamicSensorClass();
		int nPoints = sensorClass.exportView(filename, StsDynamicSensorClass.EVENTS);
		defineProximity.panel.setMessage("Output " + nPoints
				+ " points to file: " + filename);
		return true;
	}

	public boolean checkPrerequisites() 
	{
		StsDynamicSensorClass sensorClass = (StsDynamicSensorClass) model.getStsClass("com.Sts.DBTypes.StsDynamicSensor");
		StsObject[] sensors = sensorClass.getSensors();
		if (sensors.length == 0) {
			reasonForFailure = "No loaded and visible (enabled) sensors. Must have at least one sensor that has been loaded and is visible in the 3D view.";
			return false;
		}
		StsWellClass wellClass = (StsWellClass) model.getStsClass("com.Sts.DBTypes.StsWell");
        StsLiveWellClass lwellClass = (StsLiveWellClass) model.getStsClass("com.Sts.DBTypes.StsLiveWell");

		if ((wellClass.getElements().length == 0) && (lwellClass.getElements().length == 0))
        {
			reasonForFailure = "No loaded wells. Must have at least one well that has been loaded.";
			return false;
		}
        if(model.getProject().isRealtime())
        {
            reasonForFailure =  "Cannot execute sensor related actions while real-time is running. Stop real-time, run this wizard and then re-start real-time.";
            return false;
        }
		return true;
	}

	public boolean saveClusters() 
	{
		StsDynamicSensorClass sensorClass = ((StsDynamicSensor) selectedSensors[0]).getDynamicSensorClass();
		String attName = "Proximity" + System.currentTimeMillis();
		// Ask the user for a name....
		StsNameInputDialog dialog = new StsNameInputDialog(frame,"Input Attribute Name for In or Out Value (-1 or 1)",
				"Attribute Name:", attName, true);
		dialog.setVisible(true);
		if (dialog.wasCanceled()) {
			return false;
		}
		// Output the cluster as attribute.
		String userName = dialog.getUserName();
		if (userName != null)
			attName = userName;
		sensorClass.saveClusters(attName);
		sensorClass.saveAttribute();
		defineProximity.panel.setMessage("Proximity events have been saved with 2 attributes");
		Main.logUsage();
		return true;
	}

	public boolean analyze() 
	{
		StsPoint sensorPt = null;
		StsWell currentWell = null;

		Main.logUsage();
		// Verify the XYZ information exists if clustering by distance.
		if (!verifyXYZ()) 
		{
			new StsMessage(frame,StsMessage.ERROR,"Selected sensor is static so distance limit is invalid, set to -1 or pick a different sensor");
			return false;
		}

		int ttlPts = 0;
		// If linearly relating distance or time to amplitude, load the
		// amplitudes
		dist = new float[selectedSensors.length][];
		for (int j = 0; j < selectedSensors.length; j++) 
		{
			StsDynamicSensor sensor = (StsDynamicSensor) selectedSensors[j];

			// Load time vector
			long[] times = sensor.getTimeCurve(StsLogVector.types[StsLogVector.X]).getTimeVectorLongs();
			int[] passed = new int[times.length];
			int nPoints = 0;

			// Loop through the available wells
			// StsObject[] wells = model.getObjectList(StsWell.class);
			// for(int k=0; k<wells.length; k++)
			// {
			// currentWell = (StsWell)wells[k];
			currentWell = selectedWell;
			dist[j] = new float[times.length];
			// Loop through points
			for (int i = 0; i < times.length; i++) {
				passed[i] = -1;
				dist[j][i] = -1.0f;

				sensorPt = sensor.getZPoint(i);

				StsPoint pt = null;
				if (!ignoreVertical)
					pt = StsMath.getNearestPointOnLine(sensorPt, currentWell.getRotatedPoints(), 3);
				else
					pt = StsMath.getNearestPointOnLine(sensorPt, currentWell.getRotatedPoints(), 2);

				dist[j][i] = pt.distance(sensorPt);
				// System.out.println("Distance to well (" +
				// currentWell.getName() + " is " + dist);
				if (((dist[j][i] < distanceLimit) && insideLimits)
						|| ((dist[j][i] > distanceLimit) && !insideLimits)) {
					nPoints++;
					passed[i] = 1;
				}
			}
			// }
			sensor.setClusters(passed);
			sensor.setAttribute(dist[j], "DistToWell");
			model.viewObjectRepaint(this, sensor);
			ttlPts += nPoints;
		}
		if (insideLimits)
			defineProximity.setAnalysisMessage("Found " + ttlPts
					+ " points inside proximity criteria.");
		else
			defineProximity.setAnalysisMessage("Found " + ttlPts
					+ " points outside proximity criteria.");
		Main.logUsage();
		return true;
	}

    public boolean saveAsFilter()
    {
        int numRanges = 0;
        StsSensorProximityFilter filter = null;

        if(filter == null)
            filter = new StsSensorProximityFilter("ProximityFilter");

        filter.setWell(selectedWell);
        filter.setDistanceCriteria(distanceLimit, insideLimits);
        
        // Must update the static beans list of filters
        StsObject[] sensors = ((StsSensorClass)model.getCreateStsClass("com.Sts.DBTypes.StsSensor")).getSensors();
        for(int i=0; i<sensors.length; i++)
            ((StsSensor)sensors[i]).updateFilters();

       	defineProximity.setAnalysisMessage("Proximity filter successfully saved.");
        return true;
    }
    
	static public void main(String[] args) 
	{
		StsModel model = StsModel.constructor();
		StsActionManager actionManager = new StsActionManager(model);
		StsProximityAnalysisWizard clusterWizard = new StsProximityAnalysisWizard(
				actionManager);
		clusterWizard.start();
	}
}
