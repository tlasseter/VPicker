
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.PostStack3dLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;

public class StsVolumeWizard extends StsWizard
{
    private StsVolumeSelect selectVolumes;
//    private StsVolumeCrop cropVolumes;
    private StsVolumeLoad loadVolumes;
    private StsFile[] selectedFiles = null;
    private String[] filenameEndings = null;
    private StsSeismicVolume[] seismicVolumes = null;

    private StsWizardStep[] mySteps =
    {
        selectVolumes = new StsVolumeSelect(this),
//        cropVolumes = new StsVolumeCrop(this),
        loadVolumes = new StsVolumeLoad(this)
    };

    public StsVolumeWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
        dialog.setPreferredSize(500, 600);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("Define PostStack3d");
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
        if(currentStep == loadVolumes)
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
                new StsMessage(frame, StsMessage.ERROR, "No volumes selected: select or cancel.");
                return;
            }

            filenameEndings = selectVolumes.getFilenameEndings(selectedFiles);
//            int nFiles = selectedFiles.length;

//            seismicVolumes = new StsSeismicVolume[nFiles];
//            success = true;
        /*
            for(int n = 0; n < nFiles; n++)
            {
                seismicVolumes[n] = new StsSeismicVolume(false);
                if(!seismicVolumes[n].classInitialize(selectedFiles[n], model))
                    success = false;
            }
         */
        loadVolumes.constructPanel();
            gotoStep(loadVolumes);
//            if(success) gotoStep(loadVolumes);
        }
    /*
        else if(currentStep == cropVolumes)
        {
            seismicVolumes = cropVolumes.getSeismicVolumes();
            loadVolumes.addFilenamesToPanel(filenameEndings);
            for(int n = 0; n < seismicVolumes.length; n++)
                setCroppedValues(seismicVolumes[n]);
            gotoStep(loadVolumes);
        }
    */
        else
            gotoNextStep();
    }
/*
    public boolean setCroppedValues(StsSeismicVolume seismicVol)
    {
        float[] vals = new float[4];

        vals[1] = seismicVol.getRowNumMin();
        if(cropVolumes.getRowNumMin() > vals[1])
            vals[1] = cropVolumes.getRowNumMin();
        vals[2] = seismicVol.getRowNumMax();
        if(cropVolumes.getRowNumMax() < vals[2])
            vals[2] = cropVolumes.getRowNumMax();
        vals[0] = (int)((vals[2] - vals[1])/seismicVol.getRowNumInc()) + 1;
        seismicVol.setCroppedValues(seismicVol.YDIR, vals);

        vals[1] = seismicVol.getColNumMin();
        if(cropVolumes.getColNumMin() > vals[1])
            vals[1] = cropVolumes.getColNumMin();
        vals[2] = seismicVol.getColNumMax();
        if(cropVolumes.getColNumMax() < vals[2])
            vals[2] = cropVolumes.getColNumMax();
        vals[0] = (int)((vals[2] - vals[1])/seismicVol.getColNumInc()) + 1;
        seismicVol.setCroppedValues(seismicVol.XDIR, vals);

        vals[1] = seismicVol.getZMin();
        if(cropVolumes.getZMin() > vals[1])
            vals[1] = cropVolumes.getZMin();
        vals[2] = seismicVol.getZMax();
        if(cropVolumes.getZMax() < vals[2])
            vals[2] = cropVolumes.getZMax();
        vals[0] = (int)((vals[2] - vals[1])/seismicVol.getZInc()) + 1;
        seismicVol.setCroppedValues(seismicVol.ZDIR, vals);
        return true;
    }
*/
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

    public StsSeismicVolume[] getSeismicVolumes() { return seismicVolumes; }

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
