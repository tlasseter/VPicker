//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.LogToVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.*;

public class StsLogToVolumeWizard extends StsWizard
{
    int nWells;
    boolean resampleLog = false;
    int nsampleInterval = 1;
    public StsRotatedGridBoundingBox boundingBox;

    public StsWell selectedWell = null;
    public StsSelectWells selectWells = null;
    public Object[] selectedWells = new Object[0];

    StsSelectLog selectLog = null;
    StsBuildVolume buildVolume = null;
    StsWizardStep[] mySteps;

    public StsLogToVolumeWizard(StsActionManager actionManager)
    {
        super(actionManager, 500, 600);
        setupBoundingBox();
        addSteps();
    }

    private void addSteps()
    {
        selectWells = new StsSelectWells(this);
        buildVolume = new StsBuildVolume(this);
        selectLog = new StsSelectLog(this);
        mySteps = new StsWizardStep[] {selectWells, selectLog, buildVolume};
        addSteps(mySteps);
    }

    public boolean start()
    {

        System.runFinalization();
        System.gc();
        dialog.setTitle("Create Volume from Log");
        initialize();
        this.disableFinish();
        return super.start();
    }

    public void initialize()
    {
    }

    public StsLogCurve getSelectedLog()
    {
        return selectLog.getSelectedLog();
    }

    public int getNumberOfSampleInterval() { return nsampleInterval; }
    public boolean getResampleLog() { return resampleLog; }
    public void setNumberOfSampleInterval(int interval) { nsampleInterval = interval; }
    public void setResampleLog(boolean resample) { resampleLog = resample; }

    private void setupBoundingBox()
    {
        boundingBox = model.getProject().getRotatedBoundingBox().getGridBoundingBoxClone();
        boundingBox.adjustBoundingBox(200);
    }
    public boolean gridDefinitionExists()
    {
        // If valid depth grid already exists sample interval must be equal to that of grid.
        if((model.getProject().getRotatedBoundingBox().getNRows() != 0) && (model.getProject().supportsDepth()))
            return true;
        return false;
    }
 	public double getXOrigin() { return boundingBox.xOrigin; }
	public double getYOrigin() { return boundingBox.yOrigin; }
	public float getXInc() { return boundingBox.xInc; }
	public float getYInc() { return boundingBox.yInc; }
    public float getZInc() { return boundingBox.zInc; }
	public int getNRows() { return boundingBox.nRows; }
	public int getNCols() { return boundingBox.nCols; }
    public int getNSlices() { return boundingBox.nSlices; }

	public void setXOrigin(double value) { boundingBox.xOrigin = value; }
	public void setYOrigin(double value) { boundingBox.yOrigin = value; }
	public void setXInc(float value) { boundingBox.xInc = value; }
	public void setYInc(float value) { boundingBox.yInc = value; }
    public void setZInc(float value) { boundingBox.zInc = value; }
	public void setNRows(int value) { boundingBox.nRows = value; }
	public void setNCols(int value) { boundingBox.nCols = value; }
    public void setNSlices(int value) { boundingBox.nSlices = value; }

    public Object getSelectedWells()
    {
        return selectedWells;
    }

    public void setSelectedWells(Object sensor)
    {
        if(sensor == null) return;
        selectedWells = selectWells.getSelectedWells();
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
        if (success)
        {
            model.setActionStatus(getClass().getName(), StsModel.STARTED);
        }
        return super.end();
    }

    public void previous()
    {
        gotoPreviousStep();
    }

    public void next()
    {
        if(currentStep == selectWells)
        {
            Object[] wells = selectWells.getSelectedWells();
            if(wells.length == 0)
            {
                new StsMessage(frame, StsMessage.ERROR, "Must select at least one well with logs.");
                return;
            }
        }
         gotoNextStep();
    }

    public void finish()
    {
        success = true;
        super.finish();
    }

    static public void main(String[] args)
    {
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
        StsLogToVolumeWizard wellWizard = new StsLogToVolumeWizard(actionManager);
        wellWizard.start();
    }

}