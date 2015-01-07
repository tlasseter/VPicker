
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.FracSetMVLoad;

//import com.Sts.Actions.Import.StsFracSetImport;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;

public class StsMVFracSetLoad extends StsWizardStep implements Runnable
{
    private StsProgressPanel panel;
    private StsHeaderPanel header;
    StsSeismicVolume[] seismicVolumes;
    private StsMVFracSetLoadWizard wizard = null;

    private boolean isDone = false;
	private boolean canceled = false;

    public StsMVFracSetLoad(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsMVFracSetLoadWizard)wizard;
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
            StsMVFracSetFile[] selectedFiles = wizard.getFracSetFiles();

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
            // Process the FracSet files.
            for(int i = 0;  i< selectedFiles.length; i++)
            {
                // Process the files.
                String path = selectedFiles[i].file.getDirectory();
                StsMVFractureSetImport.setCurrentDirectory(path);
//                StsFractureSetImport.setStartDateAndType(wizard.getAsciiStartTime(),  type);
                StsMVFractureSet fracSet = StsMVFractureSetImport.createFractureSet(model, panel, selectedFiles[i]);

                fracSet.finish();
                if(fracSet.isPersistent())
                {
            	    panel.appendLine("Setting FracSet (" + fracSet.name + ") to " + fracSet.getStsColor().toLabelString() + "");
            	    panel.setValue(i+1);
            	    panel.setDescription("Loaded FracSet #" + (i+1) + " of " + selectedFiles.length);
                }
                else
                {
            	    panel.appendLine("FracSet (" + fracSet + ") NOT loaded. User aborted load.");
            	    panel.setValue(i+1);
            	    panel.setDescription("Canceled load FracSet #" + (i+1) + " of " + selectedFiles.length);
                }
            }

//            success = (FracSets != null);
            panel.appendLine("FracSet loading is complete. Press the Finish> button");

            project.adjustBoundingBoxes(true, true); // extend displayBoundingBox as needed and set cursor3d box accordingly
            project.checkAddUnrotatedClass(StsMVFractureSet.class);
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
