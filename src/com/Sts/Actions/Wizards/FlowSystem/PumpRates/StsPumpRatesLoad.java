
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.FlowSystem.PumpRates;

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.SensorLoad.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsPumpRatesLoad extends StsWizardStep implements Runnable
{
    private StsProgressPanel panel;
    private StsHeaderPanel header;
    private StsPumpRatesWizard wizard = null;

    int[] attIndices = null;
    private StsSensor[] pumps = null;
    StsPumpRatesFactory pumpFactory = null;

    private boolean isDone = false;
	private boolean canceled = false;

    public StsPumpRatesLoad(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsPumpRatesWizard)wizard;
    }

    public void constructPanel()
    {
        panel = StsProgressPanel.constructorWithCancelButton(10, 50);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Tank Levels Load");
        header.setSubtitle("Load Selected Pump Rates File(s)");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#PumpRates");
        header.setInfoText(wizardDialog,"(1) Once complete, press the Finish Button to dismiss the screen");
    }

    public boolean start()
    {
        run();
        return true;
    }

    public void setPumpRatesFactory(StsPumpRatesFactory sensorFactory)
    {
        this.pumpFactory = sensorFactory;
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
            panel.appendLine("Starting pump loading...");

            StsView currentView = model.win3d.getCurrentView();
            StsProject project = model.getProject();
            boolean projectInitialized = project.isInitialized();
            panel.initialize(selectedFiles.length);
            
            // Enable Static Sensors
    		StsPumpClass sc = (StsPumpClass)model.getCreateStsClass(StsPump.class);
    		sc.setEnable(true);
    		
            // Process the sensor files.
            for(int i = 0;  i< selectedFiles.length; i++)
            {

                // Verify that the project date format works for the supplied sensor data.
                byte type = selectedFiles[i].timeType;
                if(type == -1)
                {
                    panel.appendLine("Unable to load sensor from file: " + selectedFiles[i].file.getFilename() +
                    		"\n      Verify that a Time column exists in the file....");
                    continue;
                }
                // Is it multi-stage
                boolean multiStage = false;
                int numStages = selectedFiles[i].numStages(model);
                if(numStages > 25) // Seems unusually large
                {
                	if(!StsYesNoDialog.questionValue(model.win3d,"Based on the multi-stage criteria for the project,\n" +
                			"File: " + selectedFiles[i].file.getFilename() + " appears to have " + numStages + " stages." +
                					"\n\nThis number of stages seems unusually large." +
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
                	}
                }
                if(numStages > 1)
                {
                	panel.appendLine(selectedFiles[i].file.getFilename() + " appears to contain multiple stages.");
                	if(StsYesNoDialog.questionValue(model.win3d,"File: " + selectedFiles[i].file.getFilename() + " appears to be multi-stage.\n\n" +
                			"     Do you want to load it as seperate stages?"))
                	{
                		multiStage = true;
                		panel.initialize(numStages);
                		panel.appendLine("     User has elected to load as seperate stages.");
                	}
                	else
                		panel.appendLine("     User has elected to load multi-stage file as single stage.");
                }

                // Process the files.
                String path = selectedFiles[i].file.getDirectory();
                StsSensorImport.setCurrentDirectory(path);
 //               pumps = StsSensorImport.createSensors(model, panel, selectedFiles[i], multiStage,
 //               		                          wizard.getPosition(), wizard.getRelativePosition(), pumpFactory);

                // Multiple stages per file
                if(pumps.length > 1)
                {
                	for(int j=0; j<pumps.length; j++)
                	{
                		((StsPump)pumps[j]).setType(wizard.getSensorType());
                		((StsPump)pumps[j]).setDisplayTypeString(StsSensor.displayTypeStrings[StsSensor.SHOW_SINGLE]);
                		((StsPump)pumps[j]).setSymbolString(StsSensor.SYMBOL_TYPE_STRINGS[StsSensor.SPHERE]);
                		((StsPump)pumps[j]).setScaleMin(0.0f);
                		((StsPump)pumps[j]).setScaleMax(100.0f);
                		((StsPump)pumps[0]).setProperty(((StsPump)pumps[0]).getTimeCurve("Rate"));
                		 panel.appendLine("Setting pump (" + pumps[j].getName() + ") properties.");
                	}
                }
                // One stage per file
                else
                {
            		((StsPump)pumps[0]).setType(wizard.getSensorType());
            		((StsPump)pumps[0]).setDisplayTypeString(StsSensor.displayTypeStrings[StsSensor.SHOW_SINGLE]);
            		((StsPump)pumps[0]).setSymbolString(StsSensor.SYMBOL_TYPE_STRINGS[StsSensor.SPHERE]);
            		((StsPump)pumps[0]).setScaleMin(0.0f);
            		((StsPump)pumps[0]).setScaleMax(100.0f);
            		((StsPump)pumps[0]).setProperty(((StsPump)pumps[0]).getTimeCurve("Rate"));
                	panel.appendLine("Setting tank (" + pumps[0].getName() + ") properties.");
                	panel.setValue(i+1);
                	panel.setDescription("Loaded tank #" + (i+1) + " of " + selectedFiles.length);
                }

            }

            success = (pumps != null);
            panel.appendLine("Pump loading is complete. Press the Finish> button");

            project.adjustBoundingBoxes(true, true); // extend displayBoundingBox as needed and set cursor3d box accordingly
            project.checkAddUnrotatedClass(StsPump.class);
            project.rangeChanged();
            if(!projectInitialized) model.win3d.setDefaultView();

            panel.setDescription("Loading Complete");
            panel.finished();
            isDone = true;

            model.win3d.cursor3d.initialize();
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
