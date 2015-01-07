
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SensorLoad;

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsSensorLoad extends StsWizardStep implements Runnable
{
    private StsProgressPanel panel;
    private StsHeaderPanel header;
    StsSeismicVolume[] seismicVolumes;                         
    private StsSensorLoadWizard wizard = null;

    private float[][] attributes = null;
    int[] attIndices = null;
    private String[] names = null;
    private float[] time = null;
    private String[] timeString = null;
    private StsSensor[] sensors = null;
    StsAbstractSensorFactory sensorFactory = null;

    private boolean isDone = false;
	private boolean canceled = false;

    public StsSensorLoad(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsSensorLoadWizard)wizard;
    }

    public void constructPanel()
    {
        panel = StsProgressPanel.constructorWithCancelButton(10, 50);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Sensor Load");
        header.setSubtitle("Load Selected Sensor(s)");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#SensorLoad");
        header.setInfoText(wizardDialog,"(1) Once complete, press the Finish Button to dismiss the screen");
    }

    public boolean start()
    {
        run();
        return true;
    }

    public void setSensorFactory(StsAbstractSensorFactory sensorFactory)
    {
        this.sensorFactory = sensorFactory;
    }

    public void run()
    {
        try
        {
            StsSensorFile[] selectedFiles = wizard.getSensorFiles();

            if (canceled)
            {
                success = false;
                return;
            }
            wizardDialog.enableButton(wizardDialog.PREVIOUS, false);

            // turn off redisplay
            model.disableDisplay();
            panel.appendLine("Starting sensor loading...");

            StsView currentView = model.getGlPanel3d().view;
            StsProject project = model.getProject();
            boolean projectInitialized = project.isInitialized();
            panel.initialize(selectedFiles.length);
            boolean acceptNumStages = false;        // If unusual number of stages in file ask user -- ONCE.
            boolean askMultipleStages = false;

            StsSensorClass sensorClass = (StsSensorClass)model.getCreateStsClass(StsSensor.class);
            // Process the sensor files.
            for(int i = 0;  i< selectedFiles.length; i++)
            {
                // Verify that the project date format works for the supplied sensor data.
                if(selectedFiles[i].timeType == -1)
                {
                    panel.appendLine("Unable to load sensor from file: " + selectedFiles[i].file.getFilename() +
                    		"\n      Verify that a Time column exists in the file....");
                    continue;
                }

                String path = selectedFiles[i].file.getDirectory();
                StsSensorImport.setCurrentDirectory(path);
                //StsSensorImport.setStartDateAndType(wizard.getAsciiStartTime(),  type);

                // Is it multi-stage
                boolean multiStage = false;
                int numStages = selectedFiles[i].numStages(model);
                if((numStages > 25) && (!acceptNumStages)) // Seems unusually large
                {
                	if(!StsYesNoDialog.questionValue(model.win3d,"Based on the multi-stage criteria for the project,\n" +
                			"File: " + selectedFiles[i].file.getFilename() + " appears to have " + numStages + " stages." +
                					"\n\nThis number of stages seems unusually large.\n\n" +
        			"       Do you want to continue?"))
                	{
                		panel.appendLine("     User has elected to cancel load due to excessive number of stages detected.");
                		success = false;
                        model.enableDisplay();
                		return;
                	}
                	else
                	{
                		panel.appendLine("     User has elected to continue loading even though number of stages seems excessive.");
                        acceptNumStages = true;
                	}
                }
                if((numStages > 1)  && (!askMultipleStages))
                {
                	panel.appendLine(selectedFiles[i].file.getFilename() + " appears to contain " + numStages + " stages.");
                	if(StsYesNoDialog.questionValue(model.win3d,"File: " + selectedFiles[i].file.getFilename() + " appears to contain " + numStages + " stages.\n\n" +
                			"     Do you want to load it as separate stages?"))
                	{
                		multiStage = true;
                		panel.initialize(numStages);
                		panel.appendLine("     User has elected to load as separate stages.");
                	}
                	else
                    {
                		panel.appendLine("     User has elected to load multi-stage file as single stage.");
                    }
                    askMultipleStages = true;
                }

                // Process the files.

                byte vUnits = model.getProject().getDepthUnits();    // Need to implement
                byte hUnits = model.getProject().getXyUnits();       // Need to implement
                sensors = StsSensorImport.createSensors(model, panel, selectedFiles[i], multiStage, vUnits, hUnits, sensorFactory, wizard.getComputeAttributes());
                if(sensors == null)
                {
                    panel.appendLine("Load failed or cancelled by user");
                	panel.finished();
                	panel.setDescription("Load failed or cancelled by user.");
                    success = false;
                    return;
                }
                // Multiple stages per file
                if(sensors.length > 1)
                {
                	for(int j=0; j<sensors.length; j++)
                	{
                        sensors[j].setColorFromString(sensorClass.getDefaultColorName(j%32));
                		sensors[j].setType(wizard.getSensorType());
                        sensors[j].setComputeCurves(wizard.getComputeAttributes());
                		if(sensors[j] instanceof StsStaticSensor)
                            continue;
                		panel.appendLine("Setting sensor (" + sensors[j].getName() + ") to " + StsColor.colorNames32[j%32] + "");
                	}
                }
                // One stage per file
                else
                {
                    sensors[0].setColorFromString(sensorClass.getDefaultColorName(i%32));
                	sensors[0].setType(wizard.getSensorType());
                    sensors[0].setComputeCurves(wizard.getComputeAttributes());
                	panel.appendLine("Setting sensor (" + sensors[0].getName() + ") to " + StsColor.colorNames32[i%32] + "");
                	panel.setValue(i+1);
                	panel.setDescription("Loaded sensor #" + (i+1) + " of " + selectedFiles.length);
                }

            }

            success = (sensors != null);
            panel.appendLine("Sensor loading is complete. Press the Finish> button");

            project.adjustBoundingBoxes(true, true); // extend displayBoundingBox as needed and set cursor3d box accordingly
            project.checkAddUnrotatedClass(StsSensor.class);
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
