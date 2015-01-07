
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.PerforationAttributes;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.MVC.*;
import com.Sts.DBTypes.*;

import java.awt.event.ActionEvent;
import java.awt.event.WindowEvent;
import java.awt.event.ActionListener;

public class StsPerforationAttributesWizard extends StsWizard implements ActionListener
{
    public StsSelectWells selectWells = null;
    public StsCorrelateMarkers correlateMarkers = null;

    transient StsWell selectedWell = null;

    public StsPerforationAttributesWizard(StsActionManager actionManager)
    {
        super(actionManager, 300, 470);
        addSteps
        (
            new StsWizardStep[]
            {
                selectWells = new StsSelectWells(this),
                correlateMarkers = new StsCorrelateMarkers(this)
            }
        );
    }

    public boolean start()
    {
		disableFinish();
		dialog.addActionListener(this);
		System.runFinalization();
		System.gc();
		dialog.setTitle("Assign Perforations to Sensors");
		return super.start();
	}

	public void actionPerformed(ActionEvent e)
	{
		if (e.getID() == WindowEvent.WINDOW_CLOSING)
			cancel();
	}

	public Object getSelectedWells()
	{
		return selectedWell;
	}

    // Highlight the selected Wells
	public void setSelectedWells(Object well)
	{
		if (well == null)
			return;
        selectedWell = (StsWell)well;
	}

    // Rebuild the Perf list for the selected well
	public void setSelectedWell(StsWell well)
	{
		if (well == null)
			return;
        selectedWell = well;
        // Rebuild the perf list based on selection
        if(selectedWell != null)
            correlateMarkers.refreshLists();
	}

	public StsWell getSelectedWell()
	{
		return selectedWell;
	}

	public void previous()
	{
		gotoPreviousStep();
	}

	public void next()
	{
		gotoNextStep();
	}

	public boolean checkPrerequisites()
	{
        boolean foundPerfs = false;
		StsDynamicSensorClass sensorClass = (StsDynamicSensorClass) model.getStsClass("com.Sts.DBTypes.StsDynamicSensor");
		StsObject[] sensors = sensorClass.getSensors();
		if (sensors.length == 0) {
			reasonForFailure = "No loaded and visible (enabled) sensors. Must have at least one sensor that has been loaded and is visible in the 3D view.";
			return false;
		}
		StsWellClass wellClass = (StsWellClass) model.getStsClass("com.Sts.DBTypes.StsWell");
		if (wellClass.getElements().length == 0) {
			reasonForFailure = "No loaded wells. Must have at least one well that has been loaded.";
			return false;
		}
        else
        {
            StsObject[] wells = wellClass.getElements();
            for(int i=0; i<wells.length; i++)
            {
                if(((StsWell)wells[i]).getPerforationMarkers() != null)
                {
                    foundPerfs = true;
                    break;
                }
            }
        }
        if(!foundPerfs)
        {
            reasonForFailure = "There are no perforations in the loaded wells.";
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
		StsPerforationAttributesWizard wellZoneWizard = new StsPerforationAttributesWizard(actionManager);
		wellZoneWizard.start();
	}
}