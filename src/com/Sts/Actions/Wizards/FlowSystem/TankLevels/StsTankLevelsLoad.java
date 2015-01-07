
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.FlowSystem.TankLevels;

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

public class StsTankLevelsLoad extends StsWizardStep implements Runnable
{
    private StsProgressPanel panel;
    private StsHeaderPanel header;
    StsSeismicVolume[] seismicVolumes;
    private StsTankLevelsWizard wizard = null;

    int[] attIndices = null;
    private StsSensor[] tanks = null;
    StsTankLevelsFactory tankFactory = null;

    private boolean isDone = false;
	private boolean canceled = false;

    public StsTankLevelsLoad(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsTankLevelsWizard)wizard;
    }

    public void constructPanel()
    {
        panel = StsProgressPanel.constructorWithCancelButton(10, 50);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Tank Levels Load");
        header.setSubtitle("Load Selected Tank Level File(s)");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#TankLevels");
        header.setInfoText(wizardDialog,"(1) Once complete, press the Finish Button to dismiss the screen");
    }

    public boolean start()
    {
        run();
        return true;
    }

    public void setTankLevelsFactory(StsTankLevelsFactory sensorFactory)
    {
        this.tankFactory = sensorFactory;
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

            StsView currentView = model.win3d.getCurrentView();
            StsProject project = model.getProject();
            boolean projectInitialized = project.isInitialized();
            panel.initialize(selectedFiles.length);
            
            // Enable Static Sensors
    		StsTankClass sc = (StsTankClass)model.getCreateStsClass(StsTank.class);
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
               // StsSensorImport.setStartDateAndType(wizard.getAsciiStartTime(),  type);
                //tanks = StsSensorImport.createSensors(model, panel, selectedFiles[i], multiStage,
                //		                          wizard.getPosition(), wizard.getRelativePosition(), tankFactory);
                
                // Multiple stages per file
                if(tanks.length > 1)
                {
                	for(int j=0; j<tanks.length; j++)
                	{
                		((StsTank)tanks[j]).setType(wizard.getSensorType());
                		((StsTank)tanks[j]).setDisplayTypeString(StsSensor.displayTypeStrings[StsSensor.SHOW_SINGLE]);
                		((StsTank)tanks[j]).setSymbolString(StsSensor.SYMBOL_TYPE_STRINGS[StsSensor.CYLINDER]);
                		((StsTank)tanks[j]).setScaleMin(0.0f);
                		((StsTank)tanks[j]).setScaleMax(100.0f);
                		((StsTank)tanks[0]).setProperty(((StsTank)tanks[0]).getTimeCurve("Volume"));
                		 panel.appendLine("Setting tank (" + tanks[j].getName() + ") properties.");
                	}
                }
                // One stage per file
                else
                {
            		((StsTank)tanks[0]).setType(wizard.getSensorType());
            		((StsTank)tanks[0]).setDisplayTypeString(StsSensor.displayTypeStrings[StsSensor.SHOW_SINGLE]);
            		((StsTank)tanks[0]).setSymbolString(StsSensor.SYMBOL_TYPE_STRINGS[StsSensor.CYLINDER]);
            		((StsTank)tanks[0]).setScaleMin(0.0f);
            		((StsTank)tanks[0]).setScaleMax(100.0f);
            		((StsTank)tanks[0]).setProperty(((StsTank)tanks[0]).getTimeCurve("Volume"));
                	panel.appendLine("Setting tank (" + tanks[0].getName() + ") properties.");
                	panel.setValue(i+1);
                	panel.setDescription("Loaded tank #" + (i+1) + " of " + selectedFiles.length);
                }

            }

            success = (tanks != null);
            panel.appendLine("Tank loading is complete. Press the Finish> button");

            project.adjustBoundingBoxes(true, true); // extend displayBoundingBox as needed and set cursor3d box accordingly
            project.checkAddUnrotatedClass(StsTank.class);
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
