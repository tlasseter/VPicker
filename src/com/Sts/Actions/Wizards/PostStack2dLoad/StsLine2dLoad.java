package com.Sts.Actions.Wizards.PostStack2dLoad;

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System


import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;

public class StsLine2dLoad extends StsWizardStep implements Runnable
{
	private StsProgressPanel panel;
	private StsHeaderPanel header;
	StsSeismicLineSet[] seismicVolumes;
	private StsLine2dWizard wizard = null;
	private boolean isDone = false;
	private boolean canceled = false;

	public StsLine2dLoad(StsWizard wizard)
	{
		super(wizard);
		this.wizard = (StsLine2dWizard)wizard;
	}

	public void constructPanel()
	{
		panel = StsProgressPanel.constructorWithCancelButton();
		header = new StsHeaderPanel();
		setPanels(panel, header);
		header.setTitle("Post-Stack2d Selection");
		header.setSubtitle("Load Post-Stack2d Files(s)");
		header.setInfoText(wizardDialog, "(1) Once complete, press the Finish Button to dismiss the screen");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#PostStack2dLoad");
	}

	public boolean start()
	{
		run();
		return true;
	}

	public void run()
	{
		if(canceled)
		{
			success = false;
			return;
		}
		StsModel model = wizard.getModel();
		StsProject project = model.getProject();
		int nExistingVolumes = model.getObjectList(StsSeismicLineSet.class).length;
		String name = "null";

		try
		{
			StsFile[] selectedFiles = wizard.getSelectedFiles();
			String[] names = wizard.getFilenameEndings();
			disablePrevious();
			model.disableDisplay();
        /*
            StsSeismicLineSet lineSet = null;
            StsObject[] existingSeismicVolumes = model.getObjectList(StsSeismicLineSet.class);
			if (existingSeismicVolumes.length > 0)
			{
				panel.setDescriptionAndLevel("Merging Line2 with current VOLUME \n");
				StsSeismicLineSetClass vc = (StsSeismicLineSetClass)model.getStsClass(StsSeismicLineSet.class);
				lineSet = (StsSeismicLineSet)vc.getCurrentObject();
            }
        */
            StsSeismicLineSet lineSet = new StsSeismicLineSet(false);
            lineSet.add2dLines(selectedFiles, panel, this);
            lineSet.initializeColorscale();

            // Need to make excludeValue user defined.
            lineSet.initializeHistogram(0.0f);
            project.addToProject(lineSet, false);
 //           lineSet.addToProject(false);
            lineSet.addToModel();
            statusArea.setText("2D line set " + lineSet.getName() + " load successful.");
			//checkSetZDomain(line2d);

            isDone = true;
            panel.appendLine("Post-Stack2d " + name + " load successful.\n");
            panel.finished();
		}
		catch(Exception e)
		{
            StsException.outputWarningException(this, "run", e);
            String message;
			success = false;
			return;
		}
		try
		{
            project.runCompleteLoading();
            model.win3d.cursor3dPanel.setGridCheckboxState(false);
			success = true;
			return;
		}
		catch(Exception e)
		{
			StsException.outputException("StsVolumeLoad.start() failed.", e, StsException.WARNING);
			panel.setDescriptionAndLevel("StsVolumeLoad.start() failed.", StsProgressBar.ERROR);
			success = false;
			return;
		}
        finally
        {
            wizard.completeLoading(success);
        }
    }

	// thread-unsafe
/*
	public void checkSetZDomain(StsSeismicLineSet seismicVolume)
	{
		if(seismicVolume.getZDomain() != StsParameters.TD_NONE)return;
		StsComboBoxFieldBean fieldBean = new StsComboBoxFieldBean(seismicVolume, "zDomain", null, StsParameters.TD_STRINGS);
		fieldBean.setSelectedIndex(0);
		StsComboBoxDialog dialog = new StsComboBoxDialog(wizard.frame, "Select z domain", "Select time or depth domain.", fieldBean, true, true);
		seismicVolume.writeHeaderFile();
	}
*/
	public boolean end()
	{
//        wizard.finish();
//        model.win3d.getViewSelectToolbar().updateList(model.glPanel3d.getCurrentView());
		return true;
	}

	public boolean isDone()
	{
		return isDone;
	}
}
