
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.ActiveWell;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.Actions.Wizards.ProximityAnalysis.StsSelectSensors;
import com.Sts.Actions.Wizards.ProximityAnalysis.StsDefineProximity;
import com.Sts.Actions.Wizards.ProximityAnalysis.StsProximityAnalysisWizard;
import com.Sts.MVC.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;
import com.Sts.Types.StsPoint;

import java.awt.event.ActionEvent;
import java.awt.event.WindowEvent;
import java.awt.event.ActionListener;
import java.util.Iterator;

public class StsActiveWellWizard extends StsWizard implements ActionListener
{
    public StsSelectWellZone selectWellZone = null;
    public StsSelectWells selectWells = null;
    
    transient StsWell selectedWell = null;
    transient StsPerforationMarker activePerf = null;

    public StsActiveWellWizard(StsActionManager actionManager)
    {
        super(actionManager, 300, 470);
        addSteps
        (
            new StsWizardStep[]
            {
                selectWells = new StsSelectWells(this),
                selectWellZone = new StsSelectWellZone(this)
            }
        );
    }

    public boolean start()
    {
		disableFinish();
		dialog.addActionListener(this);
		System.runFinalization();
		System.gc();
		dialog.setTitle("Set Active Well Zone");
		return super.start();
	}

	public void actionPerformed(ActionEvent e)
	{
		if (e.getID() == WindowEvent.WINDOW_CLOSING)
			cancel();
	}

	public Object getSelectedPerforation()
	{
		return activePerf;
	}

    // Highlight the selected perforation
	public void setSelectedPerforation(Object perf)
	{
		if (perf == null)
			return;
		activePerf = (StsPerforationMarker)perf;
        activePerf.setHighlighted(true);
	}

    // Highlight the selected Wells
	public void setSelectedWells(Object well)
	{
		if (well == null)
			return;
        // Adjust highlights
        ((StsWell)well).setHighlighted(true);
	}

	public Object getSelectedWells()
	{
		return selectedWell;
	}

    // Rebuild the Perf list for the selected well
	public void setSelectedWell(StsWell well)
	{
		if (well == null)
			return;
        selectedWell = well;
        // Rebuild the perf list based on selection
        if((selectedWell != null) && (selectWellZone != null))
            selectWellZone.panel.rebuildPerfList(selectedWell);
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
		StsWellClass wellClass = (StsWellClass) model.getStsClass("com.Sts.DBTypes.StsWell");
		if (wellClass.getElements().length == 0) {
			reasonForFailure = "No loaded wells. Must have at least one well that has been loaded.";
			return false;
		}
		return true;
	}

	static public void main(String[] args)
	{
		StsModel model = StsModel.constructor();
		StsActionManager actionManager = new StsActionManager(model);
		StsActiveWellWizard wellZoneWizard = new StsActiveWellWizard(actionManager);
		wellZoneWizard.start();
	}
}