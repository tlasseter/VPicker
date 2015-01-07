
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.PostStack2dLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;

public class StsLine2dWizard extends StsWizard
{
    private StsLine2dSelect selectLines2d;
//    private StsLine2dCrop cropLines2d;
    private StsLine2dLoad loadLines2d;
    private StsFile[] selectedFiles = null;
    private String[] filenameEndings = null;
    private StsSeismicLine2d[] seismicLines2d = null;

    private StsWizardStep[] mySteps =
    {
        selectLines2d = new StsLine2dSelect(this),
//        cropLines2d = new StsLine2dCrop(this),
        loadLines2d = new StsLine2dLoad(this)
    };

    public StsLine2dWizard(StsActionManager actionManager)
    {
        super(actionManager, 500, 768);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("Define PostStack2d");
        disableFinish();
        return super.start();
    }

    public boolean end()
    {
        model.setActionStatus(getClass().getName(), StsModel.STARTED);
        return super.end();
    }

    public void previous()
    {
        if(currentStep == loadLines2d)
        {
            removeSeismicLines2d();
            gotoStep(selectLines2d);
        }
        else

            gotoPreviousStep();
    }

    public void next()
    {

        if(currentStep == selectLines2d)
        {
            selectedFiles = selectLines2d.getSelectedFiles();
            if(selectedFiles == null || selectedFiles.length == 0)
            {
                new StsMessage(frame, StsMessage.ERROR, "No volumes selected: select or cancel.");
                return;
            }

            filenameEndings = selectLines2d.getFilenameEndings(selectedFiles);

            loadLines2d.constructPanel();
            if(success) gotoStep(loadLines2d);
        }
    /*
        else if(currentStep == cropLines2d)
        {
            seismicLines2d = cropLines2d.getSeismicLines2d();
            loadLines2d.addFilenamesToPanel(filenameEndings);
            for(int n = 0; n < seismicLines2d.length; n++)
                setCroppedValues(seismicLines2d[n]);
            gotoStep(loadLines2d);
        }
    */
      //  else
            gotoNextStep();
    }
/*
    public boolean setCroppedValues(StsSeismicLine2d seismicVol)
    {
        float[] vals = new float[4];

        vals[1] = seismicVol.getRowNumMin();
        if(cropLines2d.getRowNumMin() > vals[1])
            vals[1] = cropLines2d.getRowNumMin();
        vals[2] = seismicVol.getRowNumMax();
        if(cropLines2d.getRowNumMax() < vals[2])
            vals[2] = cropLines2d.getRowNumMax();
        vals[0] = (int)((vals[2] - vals[1])/seismicVol.getRowNumInc()) + 1;
        seismicVol.setCroppedValues(seismicVol.YDIR, vals);

        vals[1] = seismicVol.getColNumMin();
        if(cropLines2d.getColNumMin() > vals[1])
            vals[1] = cropLines2d.getColNumMin();
        vals[2] = seismicVol.getColNumMax();
        if(cropLines2d.getColNumMax() < vals[2])
            vals[2] = cropLines2d.getColNumMax();
        vals[0] = (int)((vals[2] - vals[1])/seismicVol.getColNumInc()) + 1;
        seismicVol.setCroppedValues(seismicVol.XDIR, vals);

        vals[1] = seismicVol.getZMin();
        if(cropLines2d.getZMin() > vals[1])
            vals[1] = cropLines2d.getZMin();
        vals[2] = seismicVol.getZMax();
        if(cropLines2d.getZMax() < vals[2])
            vals[2] = cropLines2d.getZMax();
        vals[0] = (int)((vals[2] - vals[1])/seismicVol.getZInc()) + 1;
        seismicVol.setCroppedValues(seismicVol.ZDIR, vals);
        return true;
    }
*/
    private void removeSeismicLines2d()
    {
        if(seismicLines2d == null) return;
        for(int n = 0; n < seismicLines2d.length; n++)
           seismicLines2d[n].delete();
    }

    public StsSpectrum getSpectrum()
    {
       return null; // return selectLines2d.panel.getSpectrum();
    }

    public StsSeismicLine2d[] getSeismicLines2d() { return null; } //seismicLines2d; }

    public void finish()
    {
//      next();
        super.finish();
        model.enableDisplay();
        model.win3dDisplayAll(model.win3d);
    }

    public StsFile[] getSelectedFiles() { return selectedFiles; }
    public String[] getFilenameEndings() { return filenameEndings; }
}
