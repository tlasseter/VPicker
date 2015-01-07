package com.Sts.Actions.Wizards.DTS;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.DateTime.*;
import com.Sts.Utilities.*;

import java.util.*;

public class BhDtsLoad extends StsWizardStep implements Runnable
{
    private StsProgressPanel panel;
    private StsHeaderPanel header;
    private BhDtsLoadWizard wizard = null;

    int[] attIndices = null;
    private StsTimeLogCurve[] dtsTimeLogCurves = null;
    // StsTankLevelsFactory tankFactory = null;

    private boolean isDone = false;
	private boolean canceled = false;

    public BhDtsLoad(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (BhDtsLoadWizard)wizard;
    }

    public void constructPanel()
    {
        panel = StsProgressPanel.constructorWithCancelButton(10, 50);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Load DTS Files");
        header.setSubtitle("Load Selected File(s)");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#DTS");
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
            StsTimeLogCurve[] datasets = wizard.getDatasets();

            if (canceled)
            {
                success = false;
                return;
            }
            wizardDialog.enableButton(wizardDialog.PREVIOUS, false);

            // turn off redisplay
            model.disableDisplay();
            panel.appendLine("Starting DTS file loading...");

            StsProject project = model.getProject();
            panel.initialize(datasets.length);

            for(StsTimeLogCurve dataset : datasets)
            {
                dataset.readFile();
                panel.appendLine("Built Dts dataset for " + dataset.filename);
            }

            success = (dtsTimeLogCurves != null);
            panel.appendLine("Dts data loading is complete. Press the Finish> button");
         
            panel.setDescription("Loading Complete");
            panel.finished();
            isDone = true;

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

    private long[] constructTimes(String[] dates, String[] clockTimes, int nSkip)
    {
        int nTimes = dates.length - 1;
        long[] times = new long[nTimes];
        for(int n = 0; n < nTimes; n++)
        {
            String dateAndTime = dates[n+1] + " " + clockTimes[n+1];
            try
            {
                Calendar cal = CalendarParser.parse(dateAndTime, CalendarParser.YY_MM_DD, true);
                times[n] = cal.getTimeInMillis();
            }
            catch(Exception e)
            {
                times[n] = 0;
            }
        }
        return times;
    }

    public boolean end()
    {
        return true;
    }

    public boolean isDone() { return isDone; }

}
