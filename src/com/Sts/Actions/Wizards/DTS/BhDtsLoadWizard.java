package com.Sts.Actions.Wizards.DTS;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

public class BhDtsLoadWizard extends StsWizard
{
    public BhDtsSelect fileSelect = null;
    public BhAssignDtsWell assignWells = null;
    public BhDtsLoad fileLoad = null;

    byte vUnits = StsParameters.DIST_NONE;
    byte hUnits = StsParameters.TEMP_F;

    private StsFile[] files = null;
    private StsTimeLogCurve[] datasets;

    static public final String curveTypeName = "Temperature";

    private StsWizardStep[] mySteps =
    {
        fileSelect = new BhDtsSelect(this),
        assignWells = new BhAssignDtsWell(this),
        fileLoad = new BhDtsLoad(this)
    };

    public BhDtsLoadWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("Distributed Temperature Sensor Wizard");
//        super.initialize();
        disableFinish();
        if(!super.start()) return false;
        return true;
    }


    public boolean end()
    {
        if(datasets == null) return false;
        for(StsTimeLogCurve dataset : datasets)
        {
            StsWell well = dataset.getWell();
            if(well == null) continue;
            well.addLogCurve(dataset);
        //    dataset.addToModel();
        }
        model.setActionStatus(getClass().getName(), StsModel.STARTED);
        return super.end();
    }

    public void previous()
    {
        gotoPreviousStep();
    }

    public void next()
    {
        if(currentStep == fileSelect)
        {
            initializeDatasets();
            gotoStep(assignWells);
        }
        else if(currentStep == assignWells)
            gotoLoad();
        else
            gotoNextStep();
    }

    private void initializeDatasets()
    {
        StsSpectrum spectrum = getSpectrum();
        int nFiles = files.length;
        datasets = new StsTimeLogCurve[nFiles];
        for(int n = 0; n < nFiles; n++)
        {
            String curvename = StsStringUtils.trimSuffix(files[n].getFilename());
            datasets[n] = new StsTimeLogCurve(files[n], curvename, curveTypeName, spectrum, 0.0f, 250.0f);
        }
    }

    private StsSpectrum getSpectrum()
    {
        StsSpectrumClass spectrumClass = model.getSpectrumClass();
        return spectrumClass.getSpectrum(StsSpectrumClass.SPECTRUM_TEMPERATURES);
    }

    public void gotoLoad()
    {
//        super.prepareFileSets();
        fileLoad.constructPanel();
        gotoStep(fileLoad);
    }

    public void initialize()
    {
//        hUnits = model.getProject().getXyUnits();
//        vUnits = model.getProject().getDepthUnits();
    }

    public boolean checkPrerequisites()
    {
    	StsWellClass wellClass = (StsWellClass)model.getCreateStsClass("com.Sts.DBTypes.StsWell");
    	StsObject[] wells = wellClass.getObjectList();
        StsLiveWellClass lwellClass = (StsLiveWellClass)model.getCreateStsClass("com.Sts.DBTypes.StsLiveWell");
    	StsObject[] lwells = lwellClass.getObjectList();
    	if((wells.length == 0) && (lwells.length == 0))
    	{
        	reasonForFailure = "No loaded wells were found. Must have at least one well that has been loaded.";
    		return false;
    	}
    	return true;
    }

    public void addFile(StsFile file)
    {
        files = (StsFile[])StsMath.arrayAddElement(files, file);
    }

    public void removeAllFiles() { files = null; }
    
    public void removeFile(StsAbstractFile file)
    {
    	if(file == null) return;
        files = (StsFile[])StsMath.arrayDeleteElement(files, (StsFile)file);
    }

    public StsFile[] getFiles() { return files; }

    public StsTimeLogCurve[] getDatasets() { return datasets; }
}
