
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2006
//Author:       TJLasseter
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.LogToVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.SurfacesFromMarkers.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsBuildVolume extends StsWizardStep implements Runnable
{
    public StsProgressPanel panel;
    private StsHeaderPanel header;
    StsLogCurve[] logs;
    StsLogToVolumeWizard wizard;

    private boolean isDone = false;
	private boolean canceled = false;

	static private final boolean debugInterpolate = false;

    public StsBuildVolume(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsLogToVolumeWizard)wizard;
        panel = StsProgressPanel.constructorWithCancelButton();
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Create Volume from Log");
        header.setSubtitle("Define and Build Volume");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#LogToVolume");                
        header.setInfoText(wizardDialog,"(1) Select the desired log\n"  +
                "(2) Set the desired volume increments and size\n"       +
                "(3) Press the Finish Button to dismiss the screen");
    }

    public boolean start()
    {
        run();
        disableNext();
        return true;
    }

    public void run()

    {
        if (canceled)
       {
           success = false;
           return;
       }
        StsModel model = wizard.getModel();
        StsProject project = model.getProject();
        int n = -1;
        String name = "null";

        try
        {
            StsLogCurve selectedLog = wizard.getSelectedLog();
			float xInc = wizard.getXInc();
			float yInc = wizard.getYInc();
            float zInc = wizard.getZInc();
			StsRotatedGridBoundingBox boundingBox = wizard.boundingBox;
            model.disableDisplay();
            panel.initialize(1);
            panel.appendLine("Creating volume from log: " + selectedLog.getName());

                // Does it already exist? If so, delete it.
            StsSeismicVolume volume = (StsSeismicVolume)model.getObjectWithName(StsSeismicVolume.class, selectedLog.getName());
            if(volume != null)
            {
                panel.appendLine("Deleting existing volume: " + selectedLog.getName());
                volume.delete();
            }

            statusArea.setText("Volume " + selectedLog.getName() + " created and loaded successfully.");
            panel.setDescription("Volume named " + selectedLog.getName() + " has been created.");
            panel.setValue(1);

            isDone = true;
            panel.appendLine("Volume created successfully.");
            panel.setDescription("Complete");
            panel.finished();

        }
        catch(Exception e)
        {
            panel.appendLine("Failed to create volume. Error: " + e.getMessage());
            new StsMessage(wizard.frame, StsMessage.WARNING, "Failed to create volume");
            panel.setDescriptionAndLevel("Exception thrown.", StsProgressBar.ERROR);
            success = false;
            return;
        }
        try
        {
            disableCancel();
            wizard.enableFinish();
            model.enableDisplay();
            model.setActionStatus(StsSurfacesFromMarkersWizard.class.getName(), StsModel.STARTED);
        }
        catch(Exception e)
        {
            panel.appendLine("Failed to create volume. Error: " + e.getMessage());
            panel.setDescriptionAndLevel("Exception thrown.", StsProgressBar.ERROR);
            StsException.outputException("StsLogToVolume.StsBuildVolume.run() failed.", e, StsException.WARNING);
            panel.setDescriptionAndLevel("StsLogToVolume.StsBuildVolume.run() failed.", StsProgressBar.WARNING);
            success = false;
            return;
        }
    }

    public boolean end()
    {
        return true;
    }
}