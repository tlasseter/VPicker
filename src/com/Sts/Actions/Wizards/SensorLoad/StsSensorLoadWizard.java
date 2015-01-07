
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SensorLoad;

import com.Sts.Actions.Import.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsSensorLoadWizard extends StsWizard // implements Runnable
{
    public StsSensorSelect sensorSelect = null;
    public StsSensorDefineRows sensorDefineRows = null;
    public StsSensorDefineColumns sensorDefineCols = null;
    public StsSensorDefineStages defineStages = null;
    public StsSensorTime sensorTime = null;
    public StsSensorStatic sensorStatic = null;
    public StsSensorRelative sensorRelative = null;
    public StsSensorLoad sensorLoad = null;

    byte vUnits = StsParameters.DIST_NONE;
    byte hUnits = StsParameters.DIST_NONE;
    byte binaryHorzUnits = StsParameters.DIST_NONE;
    byte binaryVertUnits = StsParameters.DIST_NONE;
    
    private StsSensorFile[] sensorFiles = null;

    private StsWizardStep[] mySteps =
    {
        sensorSelect = new StsSensorSelect(this),
        sensorDefineRows = new StsSensorDefineRows(this),
        sensorDefineCols = new StsSensorDefineColumns(this),
        sensorTime = new StsSensorTime(this),
        defineStages = new StsSensorDefineStages(this),
        sensorStatic = new StsSensorStatic(this),
        sensorRelative = new StsSensorRelative(this),

        sensorLoad = new StsSensorLoad(this)
    };

    public StsSensorLoadWizard(StsActionManager actionManager)
    {
        super(actionManager, 800, 600);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("Sensor Load Wizard");
        initialize();
        disableFinish();
        if(!super.start()) return false;
        return true;
    }

    public void initialize()
    {
        sensorLoad.setSensorFactory(new StsSensorFactory());
        hUnits = model.getProject().getXyUnits();
        vUnits = model.getProject().getDepthUnits();
    }

    public boolean end()
    {
        model.setActionStatus(getClass().getName(), StsModel.STARTED);    	
        return super.end();
    }

    public void previous()
    {
        gotoPreviousStep();
    }

    public void next()
    {
        if(currentStep == sensorSelect)
        {
            gotoNextStep();
        }
        else if(currentStep == sensorDefineRows)
        {
            // Re-evaluate the time type after the user has provided the number of header rows.
            StsSensorFile[] files = getSensorFiles();
            for(int i=0; i<files.length; i++)
            {
                if(!files[i].determineTimeType(model))
                {
                    new StsMessage(frame, StsMessage.WARNING, "Unable to determine valid time & date column(s). Please verify \n" +
                                                              "that the column(s) exist and that the number of header rows are \n" +
                                                              "properly set, or remove the file from the load list.\n\n" +
                                                              files[i].getName());
                    return;
                }
            }
            gotoNextStep();
        }
        else if(currentStep == sensorDefineCols)
        {
            if(!sensorSelect.panel.getAutoConfigure())
                gotoStep(defineStages);
            else if(defineTimeRequired() == true)
                gotoStep(sensorTime);
            else if(defineStaticLocationRequired())
                gotoStep(sensorStatic);
            else if(defineRelativeLocationRequired())
            	gotoStep(sensorRelative);
            else
               gotoLoad();
        }
        else if(currentStep == defineStages)
        {
            if(defineTimeRequired() == true)
                gotoStep(sensorTime);
            else if(defineStaticLocationRequired())
                gotoStep(sensorStatic);
            else if(defineRelativeLocationRequired())
            	gotoStep(sensorRelative);
            else
               gotoLoad();
        }
        else if(currentStep == sensorTime)
        {
            if(defineStaticLocationRequired())
               gotoStep(sensorStatic);
            else if(defineRelativeLocationRequired())
            	gotoStep(sensorRelative);            
            else
               gotoLoad();
        }
        else if(currentStep == sensorStatic)
        {
        	if(defineRelativeLocationRequired())
            	gotoStep(sensorRelative); 
        	else
        		gotoLoad();
        }
        else
            gotoLoad();

    }

    public void gotoLoad()
    {
        prepareFileSets();
        sensorLoad.constructPanel();
        gotoStep(sensorLoad);
    }
    
    public byte getSensorType() { return sensorSelect.panel.getSensorType(); }
    public boolean getComputeAttributes() { return sensorSelect.panel.getComputeAttributes(); }
    private boolean defineTimeRequired()
    {
        StsSensorFile[] selectedFiles = getSensorFiles();
        for(int i=0; i<selectedFiles.length; i++)
        {
            if((selectedFiles[i].timeType != StsSensorKeywordIO.TIME_AND_DATE) &&
            	(selectedFiles[i].timeType != StsSensorKeywordIO.TIME_OR_DATE))
            {
               return true;
            }
        }
        return false;
    }
    
    private boolean defineStaticLocationRequired()
    {
        StsSensorFile[] selectedFiles = getSensorFiles();
        for(int i=0; i<selectedFiles.length; i++)
        {
            if(selectedFiles[i].positionType == StsSensor.STATIC)
            {
               return true;
            }
        }
        return false;
    }
    
    private boolean defineRelativeLocationRequired()
    {
        StsSensorFile[] selectedFiles = getSensorFiles();
        for(int i=0; i<selectedFiles.length; i++)
        {
            if(selectedFiles[i].relativeCoordinates)
            {
               return true;
            }
        }
        return false;
    }    
    
    public void prepareFileSets()
    {
       // Loop through each sensor file and create its own sensor set.
       StsSensorFile[] selectedFiles = getSensorFiles();
       for(int i=0; i<selectedFiles.length; i++)
       {
           StsSensorKeywordIO.parseAsciiFilename(sensorFiles[i].file.getFilename());
           String[] fileNames = new String[] { selectedFiles[i].file.getFilename() };
           StsSensorImport.addSensorFilenameSets(fileNames, StsSensorImport.ASCIIFILES);
       }
    }

    public void finish()
    {
        super.finish();
    }

    public String getSensorName()
    {
        return "SensorName";
    }
    
    public boolean getArchiveIt()
    {
        return false;
    }
/*
    public String getAsciiStartTime()
    {
        Date date = new Date(0L);
        SimpleDateFormat format = model.getProject().getTimeDateFormat();
        String time = format.format(date);
        try
        {
            time = sensorTime.panel.getStartTime();
            if(time == null)
            {
                new StsMessage(frame, StsMessage.ERROR,
                    "Invalid date & time value (" + format.format(date) + ") in time series file selection step.\n" +
                    "\n   Solution: Return to step and re-enter valid start time.");
                return null;
            }
            Calendar cal = CalendarParser.parse(time, model.getProject().getDateOrder(), true);
            long lvalue = cal.getTimeInMillis();
        }
        catch (Exception e)
        {
            StsMessageFiles.logMessage("Failed to create date, setting to 01-01-71 00:00:00.0");
            return "01-01-71 00:00:00.0";
        }
        return time;
    }
*/
    public boolean addFile(StsFile file)
    {
    	StsSensorFile sFile = new StsSensorFile(file);
        byte type = StsSensorFile.VARIABLE;
        if(file.filename.endsWith("src"))
           type = StsSensorFile.FIXED;
        sFile.setType(type);
    	if(!sFile.analyzeFile(model, type))
        {
            new StsMessage(frame, StsMessage.ERROR, "Sensor file (" + file.getFilename() + ") will not be processed.\n Required data is not available in file.\n\n File must have header describing columns, and a column of time values.");
            return false;
        }
        if(sFile.startTime > System.currentTimeMillis())
        {
            new StsMessage(frame,StsMessage.ERROR, "Time in the selected file (" + file.getFilename() + ") is in the future.\n File will not be processed, check the Date Order field in the project folder and try again.");
            return false;
        }
        sensorFiles = (StsSensorFile[]) StsMath.arrayAddElement(sensorFiles, sFile);
        return true;
    }

    public void removeFile(StsAbstractFile file)
    {
    	StsSensorFile sFile = getSensorFile((StsFile)file);
    	if(sFile != null)
    		sensorFiles = (StsSensorFile[])StsMath.arrayDeleteElement(sensorFiles, sFile);
    }

    public StsSensorFile getSensorFile(StsFile inFile)
    {
    	for(int i=0; i<sensorFiles.length; i++)
    		if(sensorFiles[i].file == inFile)
    			return sensorFiles[i];
    	return null;
    }
    
    public StsSensorFile[] getSensorFiles() { return sensorFiles; }
    public StsSensorFile[] getRelativeXYZSensorFiles()
    {
        StsSensorFile[] relativeFiles = new StsSensorFile[sensorFiles.length];
        int cnt = 0;
        for(int i=0; i< sensorFiles.length; i++)
        {
            if(sensorFiles[i].relativeCoordinates)
                relativeFiles[cnt++] = sensorFiles[i];
        }
        return (StsSensorFile[])StsMath.trimArray(relativeFiles, cnt);
    }
    
    public StsSensorFile[] getRelativeTimeSensorFiles()
    {
        StsSensorFile[] relativeFiles = new StsSensorFile[sensorFiles.length];
        int cnt = 0;
        for(int i=0; i< sensorFiles.length; i++)
        {
            if((sensorFiles[i].timeType != StsSensorKeywordIO.TIME_AND_DATE) &&
            	(sensorFiles[i].timeType != StsSensorKeywordIO.TIME_OR_DATE))
                relativeFiles[cnt++] = sensorFiles[i];
        }
        return (StsSensorFile[])StsMath.trimArray(relativeFiles, cnt);
    }

    public StsSensorFile[] getStaticSensorFiles()
    {
        StsSensorFile[] staticFiles = new StsSensorFile[sensorFiles.length];
        int cnt = 0;
        for(int i=0; i< sensorFiles.length; i++)
        {
            if(sensorFiles[i].positionType == StsSensor.STATIC)
                staticFiles[cnt++] = sensorFiles[i];
        }
        return (StsSensorFile[])StsMath.trimArray(staticFiles, cnt);
    }

    static String cleanTimeString(String ts)
    {
        if(ts.length() == 19)
            return ts;
        if(ts.indexOf(".") == -1)
            ts = ts + ".0";
        if(ts.indexOf(".") == 7)
            ts = "0" + ts;
        if((ts.indexOf(".") < 7) || (ts.indexOf(".") > 8))
            return null;
        ts = "2000-01-01 " + ts;
        return ts;
    }

}
