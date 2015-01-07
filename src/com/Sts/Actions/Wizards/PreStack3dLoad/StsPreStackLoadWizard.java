
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.PreStack3dLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;

public class StsPreStackLoadWizard extends StsWizard
{
    private StsPreStackLoadSelect selectVolumes;
    private StsPreStackLoad loadVolumes;
//    private StsSeismicRowColRangeEdit editVolumes;
    private StsFile[] selectedFiles = null;
    private String[] filenameEndings = null;
    private StsPreStackLineSet3d volume;

	 private StsPreStackLineSet3d previousVolume = null;

    private StsWizardStep[] mySteps =
    {
        selectVolumes = new StsPreStackLoadSelect(this),
        loadVolumes = new StsPreStackLoad(this)
//        editVolumes = new StsSeismicRowColRangeEdit(this)
    };

    public StsPreStackLoadWizard(StsActionManager actionManager)
    {
        super(actionManager, 500, 768);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("Define PreStack3d");
        dialog.getContentPane().setSize(600, 600);
        return super.start();
    }

    public boolean end()
    {
        model.setActionStatus(getClass().getName(), StsModel.STARTED);
        return super.end();
    }

    public void previous()
    {
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
        }
        else if(currentStep == loadVolumes)
        {
 //           volume = loadVolumes.preStackLineSet3d;
 //           editVolumes.initializeExtend(volume);
        }
        gotoNextStep();
    }

	 //TODO: predecessor needs to be set on the lineSet ultimately built.  Perhaps move this selection to last wizard step.
	 public StsPreStackLineSet3d getPredecessor()
	 {
		 return previousVolume;
	 }

	 public void setPredecessor(StsPreStackLineSet3d volume)
	 {
		 previousVolume = volume;
	 }

    public StsSpectrum getSpectrum()
    {
        return selectVolumes.panel.getSpectrum();
    }

    // public StsSeismicVolume[] getSeismicVolumes() { return seismicVolumes; }

    public void finish()
    {
    /*
        StsEditedBoundingBox editedBox = editVolumes.getEditedBox();
        if(editedBox != null)
        {
            volume.applyEditedBox(editedBox);
        }
   */
        volume = loadVolumes.preStackLineSet3d;
        if(volume != null)
        {
            volume.completeLoading();
            model.viewObjectChanged(this, volume);
        }
        model.getProject().runCompleteLoading();
        super.finish();
        model.win3d.validateToolbars();

//		model.win3d.checkAddView(StsViewSelectToolbar.GATHERPLOT);
//        StsViewSelectToolbar vstb = model.win3d.getViewSelectToolbar();
//        vstb.setViewItemByNamed(StsViewSelectToolbar.GATHERPLOT);
//        vstb.constructViewList();

//        model.enableDisplay();
//        model.win3dDisplayAll(model.win3d);
        completeLoading(success);
    }

    public StsFile[] getSelectedFiles() { return selectedFiles; }
    public String[] getFilenameEndings() { return filenameEndings; }
    public String getVolumeName() { return selectVolumes.panel.getName(); }
}
