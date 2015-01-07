package com.Sts.Actions.Wizards.EclipseLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;
import com.Sts.IO.*;

public class StsEclipseFileLoad extends StsWizardStep implements Runnable
{
    private StsProgressPanel panel;
    private StsHeaderPanel header;
    private StsEclipseLoadWizard wizard = null;

    private boolean isDone = false;
	private boolean canceled = false;

    public StsEclipseFileLoad(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsEclipseLoadWizard)wizard;
    }

    public void constructPanel()
    {
        panel = StsProgressPanel.constructorWithCancelButton(10, 50);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Load Eclipse Restart Files");
        header.setSubtitle("Loading Selected File(s)");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#EclipseLoad");
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
            success = true;
            StsEclipseModel dataset = wizard.getEclipseData();
            if (canceled || dataset == null)
            {
                panel.appendLine("Aborting Eclipse restart file loading...");
                success = false;
                return;
            }
            StsAbstractFile restartFile = wizard.getRestartFile();
            wizardDialog.enableButton(wizardDialog.PREVIOUS, false);
            // turn off redisplay
            model.disableDisplay();
            panel.appendLine("Starting Eclipse Restart file loading...");

            panel.initialize(1);

            // currently on a single file can be selected

            if(dataset.readFile(restartFile))
            {
                /*
                // Initialize the object
                long[] times = dataset.getTimes();
                if(times != null)
                {
                    dataset.setBornDate(times[0]);
                    dataset.setDeathDate(times[times.length -1]);
                }
                */
            }
            panel.appendLine("Built restart dataset for " + dataset.name + " from file " + restartFile.getPathname());

            if(!success)
                 panel.appendLine("     Failed to create time logs of pressure.");
            panel.appendLine("Eclipse restart data loading is complete. Press the Finish> button");

            panel.setDescription("Loading Complete");
            panel.finished();
            isDone = true;
            
            wizardDialog.enableFinish();

            // Setup Time Toolbar
            model.win3d.checkAddTimeActionToolbar();
            model.getProject().setProjectTimeToCurrentTime(true);
            model.enableDisplay();
            model.refreshObjectPanel();
        }
        catch (Exception e)
        {
        	panel.appendLine("    Unable to load selected files. Please review format.");
            panel.appendLine("    Error message: " + e.getMessage());
            StsException.outputWarningException(this, "run", e);
            success = false;
        }
    }

    public boolean end()
    {
        return true;
    }

    public boolean isDone() { return isDone; }

}