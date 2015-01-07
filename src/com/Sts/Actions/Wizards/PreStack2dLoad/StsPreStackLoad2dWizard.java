
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.PreStack2dLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.*;

public class StsPreStackLoad2dWizard extends StsWizard
{
    private StsPreStackLoadSelect2d selectVolumes;
    private StsPreStack2dSurveyDefinition surveyDef;
    private StsPreStackLoad2d loadVolumes;

    private StsFile[] selectedFiles = null;
    private String[] filenameEndings = null;
    private StsPreStackLineSet2d volume;
    private int nVolumes = 0;


	private StsPreStackLineSet2d previousVolume = null;

    private StsWizardStep[] mySteps =
    {
        selectVolumes = new StsPreStackLoadSelect2d(this),
		surveyDef  = new StsPreStack2dSurveyDefinition(this),
		loadVolumes = new StsPreStackLoad2d(this),
    };

    public StsPreStackLoad2dWizard(StsActionManager actionManager)
    {
        super(actionManager, 500, 600);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("Select 2D lines");
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
        if(currentStep == loadVolumes)
        {
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
        }
        else if(currentStep == loadVolumes)
        {
        }
        gotoNextStep();
    }

   public StsRotatedGridBoundingBox getBoundingBox()
   {
	   return surveyDef.getBoundingBox();
   }

	public StsPreStackLineSet2d getPredecessor()
	{
		return previousVolume;
	}

	public void setPredecessor(StsPreStackLineSet2d volume)
	{
		previousVolume = volume;
	}

    public StsSpectrum getSpectrum()
    {
        return selectVolumes.panel.getSpectrum();
    }

	/* create a dummy volume, create seismicLines & get bounding box info from them */
	public StsRotatedGridBoundingBox getLinesBoundingBoxes()
	{
		StsRotatedGridBoundingBox v = new StsRotatedGridBoundingBox();
		for(int n = 0; n < selectedFiles.length; n++)
		{
			StsPreStackLine2d seismicLine = StsPreStackLine2d.constructor(selectedFiles[n], n, model, null);
			if(seismicLine == null)continue;
			v.addLineRotatedBoundingBox(seismicLine);
			if(!v.originSet)
			{
				v.checkSetOriginAndAngle(seismicLine.xOrigin, seismicLine.yOrigin, seismicLine.angle);
			}
			seismicLine = null;
		}
		return v;
	}
    public void finish()
    {
//      next();
        super.finish();
//		model.win3d.glPanel3d.checkAddView(StsViewGather3d.class);
//		model.setViewPreferred(StsViewGather3d.class);
        model.win3d.validateToolbars();
        model.viewObjectChanged(this, volume);
//		model.win3d.checkAddView(StsViewSelectToolbar.GATHERPLOT);
//        StsViewSelectToolbar vstb = model.win3d.getViewSelectToolbar();
//        vstb.setViewItemByNamed(StsViewSelectToolbar.GATHERPLOT);
//        vstb.constructViewList();

        model.enableDisplay();
        model.win3dDisplayAll(model.win3d);
    }

	public void addVolume(StsPreStackLineSet2d volume)
	{
	    this.volume = volume;
    }

    public StsFile[] getSelectedFiles() { return selectedFiles; }
    public String[] getFilenameEndings() { return filenameEndings; }
    public String getLineName() { return selectVolumes.panel.getName(); }

    static void main(String[] args)
    {
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
        StsPreStackLoad2dWizard detailsWizard = new StsPreStackLoad2dWizard(actionManager);
        detailsWizard.start();
    }
}
