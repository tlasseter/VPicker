
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Monitor;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;
import com.Sts.IO.*;

import javax.swing.*;
import java.io.*;
import java.util.*;

public class StsMonitorWizard extends StsWizard
{
	private StsSelectSensor selectSensor;
    private StsSelectWell selectWell;
    private StsDefineMonitor defineMonitor;
    private StsDefineSensor defineSensor;
    private StsDefineWell defineWell;
    private StsCreateMonitor createMonitor;
    private StsMonitorSource monitorSource;

    private StsSensor sensor = null;
    private StsWell well = null;

    private JFileChooser chooseFile = null;

	public StsMonitorWizard(StsActionManager actionManager)
    {
        super(actionManager, 700, 700);
        monitorSource = new StsMonitorSource(this);
        defineMonitor = new StsDefineMonitor(this);
        defineWell = new StsDefineWell(this);
        selectSensor = new StsSelectSensor(this);
        selectWell = new StsSelectWell(this);
        defineSensor = new StsDefineSensor(this);
        createMonitor = new StsCreateMonitor(this);
		addSteps(new StsWizardStep[] { monitorSource, defineMonitor, selectSensor, selectWell, defineSensor, defineWell, createMonitor } );
    }

    public boolean start()
    {
        if(model.getProject().isRealtime())
        {
            new StsMessage(frame, StsMessage.WARNING, "Cannot define a monitor object while running realtime. Time will be stopped.");
            model.getProject().stopProjectTime();
        }
    	return super.start();
    }

    public void next()
    {
        if(currentStep == monitorSource)
        {
            if(!monitorSource.panel.isManual())
            {
                if(processMonitorFile())
                    gotoStep(createMonitor);
            }
            else
                gotoNextStep();
        }
        else if(currentStep == defineMonitor)
        {
            if(getMonitorType() == StsMonitor.SENSOR)
                gotoStep(selectSensor);
            else
                gotoStep(selectWell);
        }
        else if(currentStep == selectSensor)
        {
            if(selectSensor.panel.getSelectedSensor() == null)
            {
                new StsMessage(frame, StsMessage.ERROR, "Must select an existing sensor or create a new one.");
                return;
            }
            else
            {
                sensor = selectSensor.panel.getSelectedSensor();
                buildMonitor();
                gotoStep(createMonitor);
            }
        }
        else if(currentStep == selectWell)
        {
            if(selectWell.panel.getSelectedWell() == null)
            {
                new StsMessage(frame, StsMessage.ERROR, "Must select an existing well or create a new one.");
                return;
            }
            else
            {
                well = selectWell.panel.getSelectedWell();
                buildMonitor();
                gotoStep(createMonitor);
            }
        }
        else if((currentStep == defineSensor) || (currentStep == defineWell))
        {
            if(getMonitorType() == StsMonitor.SENSOR)
                createSensor();
            else
                createWell();
            buildMonitor();
            gotoStep(createMonitor);
        }
    	else
    		gotoNextStep();
    }

    public void previous()
    {
        gotoStep(monitorSource);
    }

    public boolean processMonitorFile()
    {
        File newFile = null;
        byte SOURCE = 0;
        byte COMPUTE = 1;
        byte RELOAD = 2;
        byte REPLACE = 3;
        byte POLLBY = 4;
        byte NAME = 5;
        byte ISSTATIC = 6;
        byte ISRELATIVE = 7;
        byte REL_X = 8;
        byte REL_Y = 9;
        byte REL_Z = 10;
        byte[] colOrder = { SOURCE, COMPUTE, RELOAD, REPLACE, POLLBY, NAME, ISSTATIC, ISRELATIVE, REL_X, REL_Y, REL_Z };
        String[] colHeaders = { "Source", "Compute", "Reload", "Replace", "PollBy", "Name", "isStatic", "IsRelative", "X", "Y", "Z"};

        if(!monitorSource.panel.isManual())
        {
            // Open file selection dialog
            chooseFile = new JFileChooser(model.getProject().getRootDirectory());
            chooseFile.setFileSelectionMode(JFileChooser.FILES_ONLY);

            chooseFile.setDialogTitle("Select the file with the defined monitors.");
            chooseFile.setApproveButtonText("Open");
            while(true)
            {
               chooseFile.showOpenDialog(null);
               if(chooseFile.getSelectedFile() == null)
                   return false;
               newFile = chooseFile.getSelectedFile();
               if(newFile.isFile() && newFile.exists())
                   break;
               else
               {
                   new StsMessage(frame, StsMessage.ERROR, "Selected file must exist and be ASCII format.");
                   chooseFile.setSelectedFile(null);
               }
            }
            // Read and process the selected file
            System.out.println("Process the selected file.");
            try
            {
                BufferedReader bufRdr = new BufferedReader(new FileReader(newFile));
                String line = bufRdr.readLine();   // Header
                line = StsStringUtils.detabString(line);
                if(line == null)
                {
                    new StsMessage(frame, StsMessage.ERROR, "Unable to read selected file.");
                    return false;
                }
                StringTokenizer stok = new StringTokenizer(line, " ,");
                int cnt = 0;
                while(stok.hasMoreTokens())
                {
                    String token = stok.nextToken();
                    for(int i=0; i<colHeaders.length; i++)
                    {
                        if(colHeaders[i].equalsIgnoreCase(token))
                        {
                            colOrder[i] = (byte)cnt;
                            cnt++;
                            break;
                        }
                    }
                }
                while(line != null)
                {
                    line = bufRdr.readLine();
				    if(line == null)
					    break;
                    line = line.trim();
                    line = StsStringUtils.detabString(line);

                    stok = new StringTokenizer(line, " ,");
				    int nTokens = stok.countTokens();
                    if(nTokens < 8)
                    {
                        StsMessageFiles.errorMessage("Monitor definition must have at least 8 fields. This line has " + nTokens + " Line: " + line);
                        continue;
                    }
                    // Parse values
                    String[] tokens = new String[stok.countTokens()];
                    for(int i=0; i<nTokens; i++)
                        tokens[i] = stok.nextToken();
                    String sourceLocation = tokens[colOrder[SOURCE]];
                    String fileName = null;
                    byte sourceType;
                    File tFile = new File(sourceLocation);
                    if(!sourceLocation.endsWith(".txt") && !sourceLocation.endsWith(".csv"))
                    {
                        sourceType = StsMonitor.DIRECTORY;
                    }
                    else
                    {

                	    sourceType = StsMonitor.FILE;
                        sourceLocation = StsFile.getDirectoryFromPathname(tFile.getPath());
                        fileName = tFile.getName();
                    }
                    StsFile.checkCreateDirectory(sourceLocation);
                    boolean computeAttributes = false;
                    boolean reloadFile = false;
                    boolean replaceEvents = false;
                    if(tokens[colOrder[COMPUTE]].equalsIgnoreCase("true"))
                        computeAttributes = true;
                    if(tokens[colOrder[RELOAD]].equalsIgnoreCase("true"))
                        reloadFile = true;
                    if(tokens[colOrder[REPLACE]].equalsIgnoreCase("true"))
                        replaceEvents = true;

                    byte pollBy = StsMonitor.TIME;
                    if(tokens[colOrder[POLLBY]].equalsIgnoreCase("Size"))
                        pollBy = StsMonitor.SIZE;
                    if(pollBy == StsMonitor.SIZE)
                    {
                        if(fileName == null)
                        {
                            new StsMessage(frame, StsMessage.ERROR, "Can only monitor files for size change, not directories.\n" +
                                    " Ignoring line: " + line + "\n\n");
                            continue;
                        }
                    }

                    String sensorName = tokens[colOrder[NAME]];
                    boolean isStatic = false;
                    boolean isRelative = false;
                    if(tokens[colOrder[ISSTATIC]].equalsIgnoreCase("true"))
                        isStatic = true;
                    if(tokens[colOrder[ISRELATIVE]].equalsIgnoreCase("true"))
                        isRelative = true;                    

                    double xRel = model.getProject().getXOrigin();
                    double yRel = model.getProject().getYOrigin();
                    double zRel = model.getProject().getDepthMin();
                    if(((isRelative) || (isStatic))  && (tokens.length > 8))
                    {
                        xRel = Double.parseDouble(tokens[colOrder[REL_X]]);
                        yRel = Double.parseDouble(tokens[colOrder[REL_Y]]);
                        zRel = Double.parseDouble(tokens[colOrder[REL_Z]]);
                        continue;
                    }
                    
                    // Build Sensor
                    if(isStatic)
                    {
                        sensor = new StsStaticSensor(null, sensorName);
                        ((StsStaticSensor)sensor).setXLoc(xRel);
    		            ((StsStaticSensor)sensor).setYLoc(yRel);
    		            ((StsStaticSensor)sensor).setZLoc(zRel);
                        sensor.setIsVisible(false);

                    }
    		        else
                        sensor = new StsDynamicSensor(null, sensorName);
                    sensor.setComputeCurves(computeAttributes);
    		        sensor.setIsRelative(isRelative);
    		        sensor.setHasDate(true);
    		        sensor.setStartDate("");
                    sensor.getTimeCurves();    // Initializes the timeCurves
                    sensor.setOriginalDomain(StsProject.TD_DEPTH);
                    //sensor.setBornDate(System.currentTimeMillis());
                    //sensor.addToProject(StsProject.TD_DEPTH);
    		        setSensor(sensor);

                    // Build the monitor
                    StsMonitor monitor = new StsMonitor(sensorName + "_Monitor", sensor, sourceLocation,
    			        sourceType, fileName, computeAttributes, reloadFile, replaceEvents, pollBy);
                }
            }
            catch(Exception ex)
            {
                ;
            }
        }
        return true;
    }

    public boolean createWell()
    {
    	try
    	{
    		// Create the empty well
    		StsLiveWell well = null;
            well = new StsLiveWell(defineWell.panel.getName(), true, StsColor.BLUE);
    		setWell(well);
    		return true;
    	}
    	catch(Exception e)
    	{
    		new StsMessage(this.frame,StsMessage.ERROR,"Failed to create empty well. Back up and try again or select an existing well.");
    		return false;
    	}
    }

    public boolean createSensor()
    {
    	try
    	{
    		// Create the empty sensor
    		StsSensor sensor = null;

    		if(defineSensor.panel.getIsStatic())
            {
                sensor = new StsStaticSensor(null, defineSensor.panel.getName());
                ((StsStaticSensor)sensor).setXLoc(defineSensor.panel.getXOrigin());
    		    ((StsStaticSensor)sensor).setYLoc(defineSensor.panel.getYOrigin());
    		    ((StsStaticSensor)sensor).setZLoc(defineSensor.panel.getZOrigin());
                sensor.setIsVisible(false);

            }
    		else
    			sensor = new StsDynamicSensor(null, defineSensor.panel.getName());

            sensor.setComputeCurves(defineMonitor.panel.getComputeAttributes());
    		sensor.setIsRelative(defineSensor.panel.getIsRelative());
    		sensor.setHasDate(defineSensor.panel.getHasDate());
    		sensor.setStartDate(defineSensor.panel.getStartDate());
            sensor.addToProject(StsProject.TD_DEPTH);
    		setSensor(sensor); 
    		return true;
    	}
    	catch(Exception e)
    	{
    		new StsMessage(this.frame,StsMessage.ERROR,"Failed to create empty sensor. Back up and try again or select an existing sensor.");
    		return false;
    	}
    }
    public void newSensor()
    {
    	gotoStep(defineSensor);
    }
    public void newWell()
    {
    	gotoStep(defineWell);
    }
    public boolean end()
    {
    	if(!super.end()) 
    		return false;
    	
    	// Build the monitor object
    	StsSensor sensor = this.getSensor();
        StsWell well = this.getWell();
    	if((sensor == null) && (well == null))
    	{
    		new StsMessage(model.win3d, StsMessage.ERROR, "Failed to find or build object to monitor.");
    		return false;
    	}
		return true;
    }

    public boolean buildMonitor()
    {
        StsMonitor monitor = null;
        if(getMonitorType() == StsMonitor.SENSOR)
    	    monitor = new StsMonitor(sensor.getName() + "_Monitor", sensor, defineMonitor.panel.getSourceLocation(),
    			defineMonitor.panel.getSourceType(), defineMonitor.panel.getFilename(), defineMonitor.panel.getComputeAttributes(),
                defineMonitor.panel.getReloadFile(), defineMonitor.panel.getReplaceEvents(),
                defineMonitor.panel.getPollBy());
        else
        {
    	    monitor = new StsMonitor(well.getName() + "Dev_Monitor", well, StsMonitor.WELL, defineMonitor.panel.getSourceLocation(),
    			defineMonitor.panel.getFilename(), defineMonitor.panel.getComputeAttributes(),
                defineMonitor.panel.getReloadFile(), defineMonitor.panel.getReplaceEvents(),
                defineMonitor.panel.getPollBy());
            StsMonitor logMonitor = new StsMonitor(well.getName() + "Log_Monitor", well, StsMonitor.LOG, defineMonitor.panel.getSourceLocation2(),
    			defineMonitor.panel.getFilename2(), defineMonitor.panel.getComputeAttributes(),
                defineMonitor.panel.getReloadFile(), defineMonitor.panel.getReplaceEvents(),
                defineMonitor.panel.getPollBy());
        }
    	if(monitor != null)
    		return true;
    	else
    		return false;
    }

    public void setWell(StsWell well)
    {
    	this.well = well;
    }
    public StsWell getWell()
    {
    	return well;
    }
    public void setSensor(StsSensor sensor)
    {
    	this.sensor = sensor;
    }
    public StsSensor getSensor()
    {
    	return sensor;
    }
    public byte getMonitorType()
    {
        return monitorSource.panel.getMonitorType();
    }
}