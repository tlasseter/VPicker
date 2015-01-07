
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.StsObjectLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;

public class StsObjectLoad extends StsWizardStep implements Runnable
{
    private StsProgressPanel panel;
    private StsHeaderPanel header;
    private StsObjectLoadWizard wizard = null;

    private boolean isDone = false;
	private boolean canceled = false;

    public StsObjectLoad(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsObjectLoadWizard)wizard;
        constructPanel();
    }

    public void constructPanel()
    {
        panel = StsProgressPanel.constructorWithCancelButton(10, 50);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("StsObject Load");
        header.setSubtitle("Load Selected Object(s)");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#StsObjectLoad");
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
            StsFile[] selectedFiles = wizard.getObjectFiles();

            if (canceled)
            {
                success = false;
                return;
            }
            wizardDialog.enableButton(wizardDialog.PREVIOUS, false);

            // turn off redisplay
            model.disableDisplay();
            panel.appendLine("Starting object loading...");

            panel.initialize(selectedFiles.length);
            // Process the sensor files.
            for(int i = 0;  i< selectedFiles.length; i++)
            {
            	// Load objects here........................
            	String filename = selectedFiles[i].getPathname();
                StsObject object = StsDBFileObjectTrader.importStsObject(filename, null);
//     			object.setIndex(-1);
//    			object.addToModel();
                panel.setValue(i+1);
                panel.setDescription("Loaded object #" + (i+1) + " of " + selectedFiles.length);
            }
            panel.appendLine("Sensor loading is complete. Press the Finish> button");

            panel.setDescription("Loading Complete");
            panel.finished();
            isDone = true;
            
			model.win3d.cursor3d.initialize();
			model.win3d.cursor3dPanel.setSliderValues();
			model.win3d.getGlPanel3d().setDefaultView();
			model.enableDisplay();
			model.win3dDisplay();
            wizard.enableFinish();
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
