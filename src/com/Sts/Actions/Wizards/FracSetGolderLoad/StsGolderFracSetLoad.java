
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.FracSetGolderLoad;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardDialog;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.Actions.Wizards.WizardHeaders.StsHeaderPanel;
import com.Sts.DBTypes.StsGolderFractureSet;
import com.Sts.DBTypes.StsGolderFractureSetClass;
import com.Sts.DBTypes.StsSeismicVolume;
import com.Sts.IO.StsFile;
import com.Sts.MVC.StsProject;
import com.Sts.MVC.View3d.StsView;
import com.Sts.MVC.View3d.StsView3d;
import com.Sts.UI.Progress.StsProgressPanel;
import com.Sts.Utilities.StsException;

public class StsGolderFracSetLoad extends StsWizardStep implements Runnable
{
    private StsProgressPanel panel;
    private StsHeaderPanel header;
    StsSeismicVolume[] seismicVolumes;
    private StsGolderFracSetLoadWizard wizard = null;

    private boolean isDone = false;
	private boolean canceled = false;

    public StsGolderFracSetLoad(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsGolderFracSetLoadWizard)wizard;
    }

    public void constructPanel()
    {
        panel = StsProgressPanel.constructorWithCancelButton(10, 50);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("FracSet Load");
        header.setSubtitle("Load Selected FracSet(s)");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#FracSet");
        header.setInfoText(wizardDialog,"(1) Once complete, press the Finish Button to dismiss the screen");
    }

    public boolean start()
    {
        run();
        return true;
    }


    public void run()
    {
        try
        {
            StsFile[] selectedFiles = wizard.getFracSetFiles();

            if (canceled)
            {
                success = false;
                return;
            }
            wizardDialog.enableButton(StsWizardDialog.PREVIOUS, false);

            // turn off redisplay
            model.disableDisplay();
            panel.appendLine("Starting FracSet loading...");

            StsView currentView = model.getGlPanel3d().view;
            StsProject project = model.getProject();
            boolean projectInitialized = project.isInitialized();
            panel.initialize(selectedFiles.length);
			// Determine if the project origin is set; if not, then the first fracture loaded will set the origin.
			// All subsequent coordinates will use this origin, i.e., x and y values will be offset from it.
			// If ultimately multithreaded, allow the first load operation to set the origin if need be and then
			// launch all the other threads which will use this initial origin
			boolean originSet = model.getProject().isOriginSet();
			StsGolderFracSetImport fracSetLoader = new StsGolderFracSetImport(model, originSet);
            // Process the FracSet files.
            for(int i = 0;  i< selectedFiles.length; i++)
            {
                // Process the files.
                fracSetLoader.loadFile(selectedFiles[i]);
            	panel.setValue(i+1);
            	panel.setDescription("Loaded FracSet #" + (i+1) + " of " + selectedFiles.length);
            }
			// this adds the fractureSets to the project boundingBoxes and to the model fracSetClass
			fracSetLoader.addSetsToProjectAndModel(model);

//            success = (FracSets != null);
            panel.appendLine("FracSet loading is complete. Press the Finish> button");

            project.adjustBoundingBoxes(true, true); // extend displayBoundingBox as needed and set cursor3d box accordingly
            project.checkAddUnrotatedClass(StsGolderFractureSet.class);
            project.rangeChanged();
            if(!projectInitialized) model.getGlPanel3d().setDefaultView();

            panel.setDescription("Loading Complete");
            panel.finished();
            isDone = true;

            model.win3d.getCursor3d().initialize();
            model.win3d.cursor3dPanel.setSliderValues();

            // Setup Time Toolbar
            model.win3d.checkAddTimeActionToolbar();
            model.getProject().setProjectTimeToCurrentTime(true);

            wizard.enableFinish();

            model.enableDisplay();
            if(currentView instanceof StsView3d)
                ((StsView3d)currentView).adjustView();
            model.win3dDisplay();
            model.refreshObjectPanel();
        }
        catch (Exception e)
        {
        	panel.appendLine("    Unable to load selected files. Please review format.");
            panel.appendLine("    Error message: " + e.getMessage());
            StsException.outputWarningException(this, "run", e);
            model.enableDisplay();
            wizard.enableFinish();
            success = false;
        }
    }
    public boolean end()
    {
        return true;
    }

    public boolean isDone() { return isDone; }

}
