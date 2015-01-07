
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2006
//Author:       TJLasseter
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.VspLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;

public class StsVspLoadWizard extends StsWizard
{
    private StsVspSelect selectVolumes;
    private StsVspWellAssign assignWells;
    private StsVspLoad loadVolumes;
    private StsAbstractFile[] selectedFiles = null;
    private String[] filenameEndings = null;
    private StsVsp[] seismicVolumes = null;

    private StsWizardStep[] mySteps =
    {
        selectVolumes = new StsVspSelect(this),
        loadVolumes = new StsVspLoad(this),
        assignWells = new StsVspWellAssign(this)
    };

    public StsVspLoadWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("Load VSP Data");
        dialog.getContentPane().setSize(500, 600);
        return super.start();
    }

    public boolean end()
    {
        model.setActionStatus(getClass().getName(), StsModel.STARTED);
        return super.end();
    }

    public void previous()
    {
        if((currentStep == loadVolumes) || (currentStep == assignWells))
        {
            removeSeismicVolumes();
            gotoStep(selectVolumes);
        }
        else
            gotoPreviousStep();
    }

    public void next()
    {
        if(currentStep == selectVolumes)
        {
            selectedFiles = selectVolumes.getSelectedFiles();
            if(selectedFiles == null || selectedFiles.length == 0)
            {
                gotoStep(assignWells);
            }
            else
            {
                filenameEndings = selectVolumes.getFilenameEndings(selectedFiles);
                loadVolumes.constructPanel();
                gotoStep(loadVolumes);
            }
        }
        else
            gotoNextStep();
    }

    private void removeSeismicVolumes()
    {
        if(seismicVolumes == null) return;
        for(int n = 0; n < seismicVolumes.length; n++)
            seismicVolumes[n].delete();
    }

    public StsSpectrum getSpectrum()
    {
        return selectVolumes.panel.getSpectrum();
    }

    public StsVsp[] getSeismicVolumes() { return seismicVolumes; }

    public void finish()
    {
        super.finish();
        model.enableDisplay();
        model.win3dDisplayAll(model.win3d);
    }

    public StsAbstractFile[] getSelectedFiles() { return selectedFiles; }
    public String[] getFilenameEndings() { return filenameEndings; }

    public boolean getArchiveIt() { return selectVolumes.panel.getArchiveIt(); }

    static void main(String[] args)
    {
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
        StsVspLoadWizard detailsWizard = new StsVspLoadWizard(actionManager);
        detailsWizard.start();
    }
}
