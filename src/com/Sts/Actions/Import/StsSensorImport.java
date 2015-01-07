//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Import;

import com.Sts.Actions.Wizards.SensorLoad.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.util.*;

public class StsSensorImport
{
    static StsModel model;
    static String currentDirectory = null;
    static String binaryDataDir;
    static float sensorNull;
    public static int nLoaded = 0;
    static boolean success = false;
    static String fullBinaryDirName;
    static TreeMap sensorFilenameSets;
    static StsAbstractFileSet[] fileSets = new StsAbstractFileSet[2];

    static String[] asciiFiles = null;
    static private String[] selectedWells = null;
    static private StsFile[] selectedFiles = null;
    static private boolean reloadAscii = false;
    static private boolean multiStage = false;
    //static private double falseX = 0.0f;
    //static private double falseY = 0.0f;
    static private StsPoint[] sensorLocations = null;
    static private byte[] sensorTypes = null;

    static final String SENSOR_PREFIX = "sensor";
    static public final int SENSOR = 0;
    static String[] filePrefixs = new String[] {SENSOR_PREFIX};

    static public final byte BINARYFILES = 2;
    static public final byte ASCIIFILES = 3;

    static byte vUnits = StsParameters.DIST_NONE;
    static byte hUnits = StsParameters.DIST_NONE;
    static byte binaryHorzUnits = StsParameters.DIST_NONE;
    static byte binaryVertUnits = StsParameters.DIST_NONE;

    static public final boolean debug = false;
    static StsTimer timer = new StsTimer();
    static final String[] typeNames = new String[] {"BINARYFILES", "ASCIIFILES"};

    static double modelXOrigin, modelYOrigin;

    static public void initialize(StsModel model_)
    {
        model = model_;
        StsProject project = model.getProject();
        currentDirectory = project.getRootDirString();
        binaryDataDir = project.getBinaryFullDirString();
        sensorNull = project.getLogNull();
        nLoaded = 0;
    }

    static public String getCurrentDirectory()
    {
        return currentDirectory;
    }

    static public void setCurrentDirectory(String dirPath)
    {
        currentDirectory = dirPath;
    }

    /** Create sensors
     *
     */
    static public StsSensor[] createSensors(StsModel model, StsProgressPanel progressPanel, StsSensorFile selectedFile, byte vertUnits, byte horzUnits,
                                            StsAbstractSensorFactory sensorFactory, boolean computeAttributes)
    {
        binaryDataDir = model.getProject().getBinaryFullDirString();

        try
        {
            int nSensors = 1;
            StsSensor[] sensors = new StsSensor[nSensors];
            nLoaded = 0;

            StsSensorKeywordIO.parseAsciiFilename(selectedFile.file.getFilename());
            String sensorname = StsSensorKeywordIO.name;
            SensorFilenameSet sensorFilenameSet = (SensorFilenameSet) sensorFilenameSets.get(StsSensorKeywordIO.name);

            if ((model.getObjectWithName(StsSensor.class, sensorname) == null) || reloadAscii)
            {
                StsSensor sensorTemp = (StsSensor) model.getObjectWithName(StsSensor.class, sensorname);
                if ( (sensorTemp != null) && reloadAscii)
                {
                    sensorTemp.delete();
                }

                byte formatType = sensorFilenameSet.sensorType;
                vUnits = vertUnits;
                hUnits = horzUnits;
                if(formatType == ASCIIFILES)
                {
                    sensors = createSensorFromAsciiFiles(model.getProject().getRotatedBoundingBox(),sensorFilenameSet, currentDirectory,
                            binaryDataDir, selectedFile, progressPanel, sensorFactory,
                            model.getProject().getDateOrder(), computeAttributes, model.getProject().getTimeString(selectedFile.startTime));
                }
            }
            else
            {
                if (progressPanel != null)
                {
                    progressPanel.appendLine("Sensor: " + sensorname + " already loaded...\n");
                }
            }
            return sensors;

        }
        catch (Exception e)
        {
            StsException.outputException("StsSensorImport.createSensors() failed.", e, StsException.WARNING);
            progressPanel.setDescriptionAndLevel("StsSensorImport.createWells() failed.\n", StsProgressBar.WARNING);
            return null;
        }
    }

    /** read in multiple well deviations from a list of Ascii files */
    static public StsSensor[] createSensorFromAsciiFiles(StsRotatedBoundingBox rBox, SensorFilenameSet sensorFilenameSet, String dataDir,
        String binDir, StsSensorFile sensorFile, StsProgressPanel progressPanel, StsAbstractSensorFactory sensorFactory, int dateOrder, boolean computeAttributes, String startTime)
    {
        StsSensor sensor = null;
        StsSensor[] sensors = null;
        if(multiStage)
        	sensors = new StsSensor[sensorFile.numStages];
        else
        	sensors = new StsSensor[1];

        String sensorname = sensorFilenameSet.sensorName;

        if (dataDir == null)
        {
            return null;
        }

        try
        {
            StsTimeCurve[] timeCurves = null;
            boolean deleteBinaries = reloadAscii;

            // Process time curve
            if(multiStage)
            {
            	for(int i=0; i<sensorFile.numStages; i++)
            	{
            		timeCurves = StsSensorKeywordIO.readTimeCurves(rBox, sensorFilenameSet.sensorFilenames, sensorname, sensorFile, dataDir,
            						binDir, deleteBinaries, vUnits, dateOrder, multiStage, i+1, startTime);
            		sensor = buildSensor(sensorname+(i+1), progressPanel, sensorFactory, sensorFile.positionType);
            		if (sensor == null)
            			return null;
            		else
            			sensors[i] = sensor;

            		// Add the TD Curve
            		if (timeCurves != null)
            			sensors[i].addTimeCurves(timeCurves, computeAttributes);
            		progressPanel.setValue(i+1);          		
                    progressPanel.appendLine("Loaded sensor for stage #" + (i+1) + " of " + sensorFile.numStages);
            	}
            }
            else
            {
        		timeCurves = StsSensorKeywordIO.readTimeCurves(rBox, sensorFilenameSet.sensorFilenames, sensorname, sensorFile, dataDir,
						binDir, deleteBinaries, vUnits, dateOrder, startTime);
        		sensor = buildSensor(sensorname, progressPanel, sensorFactory, sensorFile.positionType);
        		if (sensor == null)
        			return null;
        		else
        			sensors[0] = sensor;

        		// Add the TD Curve
        		if (timeCurves != null)
        			sensors[0].addTimeCurves(timeCurves, computeAttributes);
            }
            if (progressPanel != null)
            {
                progressPanel.appendLine("Processing ASCII Sensor: " + sensorname + "...");
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsSensorImport failed.", e, StsException.WARNING);
            StsMessageFiles.logMessage("Unable to read Ascii sensor files.");
            progressPanel.appendLine("Unable to read Ascii sensor files.");

            return null;
        }
        return sensors;
    }

    static public boolean archiveSensors(String[] sensorNames, String archiveDir)
    {
        try
        {
            for(int i=0; i<sensorNames.length; i++)
            {
                SensorFilenameSet sensorFilenameSet = (SensorFilenameSet) sensorFilenameSets.get(sensorNames[i]);
                for(int j=0; j<sensorFilenameSet.sensorFilenames.length; j++)
                {
                    if(sensorFilenameSet.sensorFilenames[j] != null)
                    {
                        if (sensorFilenameSet.sensorFilenames[j].indexOf("bin") == -1)
                            StsFile.copy(currentDirectory + File.separator + sensorFilenameSet.sensorFilenames[j],
                                         archiveDir + File.separator + sensorFilenameSet.sensorFilenames[j]);
                    }
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsSensorImport.archiveWells Error", e, StsException.WARNING);
            return false;
        }
        return true;
    }

    static public String[] getSensorFilenames()
    {
        String[] filenames = new String[0];
        SensorFilenameSet sensorFilenameSet = null;
        Iterator iter;

        int nSensors = sensorFilenameSets.size();
        if (nSensors > 0)
        {
            filenames = new String[nSensors];
            int i = 0;
            iter = sensorFilenameSets.entrySet().iterator();
            while (iter.hasNext())
            {
                Map.Entry entry = (Map.Entry) iter.next();
                sensorFilenameSet = (SensorFilenameSet) entry.getValue();
                filenames[i] = sensorFilenameSet.sensorName;
                i++;
            }
        }
        if (debug)
        {
            timer.stopPrint("Time to find " + filenames.length + " files.");
        }
        return filenames;
    }

    static public String getSensorFilename(int idx)
    {
        SensorFilenameSet sensorFilenameSet = null;
        Iterator iter;

        int nSensors = sensorFilenameSets.size();
        if (nSensors > 0)
        {
            int i = 0;
            iter = sensorFilenameSets.entrySet().iterator();
            while (iter.hasNext())
            {
                Map.Entry entry = (Map.Entry) iter.next();
                sensorFilenameSet = (SensorFilenameSet) entry.getValue();
                if (i == idx)
                {
                    return sensorFilenameSet.sensorFilenames[0];
                }
                i++;
            }
        }
        return null;
    }

    static public String[] getSelectedSensorFilenames(String[] selectedSensornames, int type, boolean includeBinaries)
    {
        String sensorname = null;
        String[] filenames = null;
        SensorFilenameSet sensorFilenameSet = null;
        String fmt = null;

        for (int n = 0; n < selectedSensornames.length; n++)
        {
            sensorname = selectedSensornames[n];
            sensorFilenameSet = (SensorFilenameSet) sensorFilenameSets.get(sensorname);
            if (sensorFilenameSet.sensorType == type)
            {
                filenames = (String[]) StsMath.arrayAddElement(filenames, sensorFilenameSet.sensorFilenames[0]);
            }
        }
        return filenames;
    }

    static public SensorFilenameSet getSensorFilename(String sensorname)
    {
        return (SensorFilenameSet) sensorFilenameSets.get(sensorname);
    }

    static public String getFilePrefix(int type) { return filePrefixs[type]; }
    static public void setFilePrefixs(String[] prefixs)
    {
        filePrefixs = prefixs;
    }
    static public boolean constructSensorFilenames(byte type)
    {
        String dirname = null;
        try
        {
            switch (type)
            {
                case BINARYFILES:
                    if (currentDirectory == null)
                    {
                        fullBinaryDirName = model.getProject().getBinaryFullDirString();
                    }
                    else
                    {
                        fullBinaryDirName = currentDirectory + File.separator;
                    }
                    fileSets[type] = StsFileSet.constructor(fullBinaryDirName, new SensorNameFilter("bin"));
                    break;

                case ASCIIFILES:
                    if (currentDirectory == null)
                    {
                        dirname = model.getProject().getRootDirectory() + "/";
                    }
                    else
                    {
                        dirname = currentDirectory + "/";
                    }
                    fileSets[type] = StsFileSet.constructor(dirname, new SensorNameFilter("txt"));
                    break;
                }
            if (debug)
            {
                timer.stopPrint("Time to build well filenameSets of type " + typeNames[type]);
            }

            if (fileSets[type] == null)
            {
                return false;
            }
            String[] filenames = fileSets[type].getFilenames();
            if (filenames == null)
            {
                return false;
            }
            return addSensorFilenameSets(filenames, type);
        }
        catch (Exception e)
        {
            StsException.outputException("StsSensorImport.constructSensorFilenameSets() failed.", e, StsException.WARNING);
            return false;
        }
    }

    static public boolean selectedOfType(int type, String[] selectedSensorNames)
    {
        String sensorname = null;
        String[] filenames = null;
        SensorFilenameSet sensorFilenameSet = null;
        for (int n = 0; n < selectedSensorNames.length; n++)
        {
            sensorname = selectedSensorNames[n];
            sensorFilenameSet = (SensorFilenameSet) sensorFilenameSets.get(sensorname);
            if (sensorFilenameSet.sensorType == type)
            {
                return true;
            }
        }
        return false;
    }

    static public void setReloadAscii(boolean reload)
    {
        reloadAscii = reload;
    }

    static public boolean getReloadAscii()
    {
        return reloadAscii;
    }

    static public void compressSensorSet()
    {
        Iterator iter;
        SensorFilenameSet sensorFilenameSet;
        if (debug)
        {
            timer.start();
        }

        iter = sensorFilenameSets.entrySet().iterator();
        while (iter.hasNext())
        {
            Map.Entry entry = (Map.Entry) iter.next();
            sensorFilenameSet = (SensorFilenameSet) entry.getValue();
            if (!sensorFilenameSet.isOK())
            {
                iter.remove();
                continue;
            }

            // If we don't want to reload it, delete it from the list of available wells.
            StsSensor existingSensor = (StsSensor) model.getObjectWithName(StsSensor.class, sensorFilenameSet.sensorName);
            if (existingSensor != null)
            {
                if(!reloadAscii)
                {
                    iter.remove();
                }
            }
        }
        if (debug)
        {
            timer.stopPrint("Time for compressSensorSet.");
        }
        return;
    }

    static public boolean addSensorFilenameSets(String[] filenames, byte type)
    {
        String filename = null;
        SensorFilenameSet sensorFilenameSet;

        try
        {
            if (debug)
            {
                timer.start();
            }

            if (filenames == null)
            {
                return false;
            }

            if (sensorFilenameSets == null)
            {
                sensorFilenameSets = new TreeMap();
            }

            int nFilenames = filenames.length;
            for (int n = 0; n < nFilenames; n++)
            {
                filename = filenames[n];
                if (type == ASCIIFILES)
                {
                    StsSensorKeywordIO.parseAsciiFilename(filename);
                }
                else
                {
                    StsSensorKeywordIO.parseBinaryFilename(filename);
                }
                String sensorname = StsSensorKeywordIO.getSensorName();
                if (sensorFilenameSets.containsKey(sensorname))
                {
                    sensorFilenameSet = (SensorFilenameSet) sensorFilenameSets.get(sensorname);
                    if(type == ASCIIFILES)
                    {
                        if (!sensorFilenameSet.asciiAlreadyInList(filename))
                        {
                            sensorFilenameSet.addAsciiFilename(filename);
                        }
                    }
                    else
                    {
                        if (!sensorFilenameSet.binaryAlreadyInList(filename))
                        {
                            sensorFilenameSet.addBinaryFilename(filename);
                        }
                    }
                }
                else
                {
                    sensorFilenameSet = new SensorFilenameSet(sensorname, type);
                    sensorFilenameSets.put(sensorname, sensorFilenameSet);
                    if (type == ASCIIFILES)
                    {
                        sensorFilenameSet.addAsciiFilename(filename);
                    }
                    else
                    {
                        sensorFilenameSet.addBinaryFilename(filename);
                    }
                }
            }
            if (debug)
            {
                timer.stopPrint("Time to initializeSensorFilenameSets for " + nFilenames + "files");
            }

            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsSensorImport.initializeSensorFilenameSets() failed for file: " + filename + ".", e,
                                         StsException.WARNING);
            return false;
        }
    }

    static public class SensorFilenameSet
    {
        public String sensorName;
        public byte sensorType = -1;
        public int nFilenames = 0;
        public String[] sensorFilenames = null;
        public int increment = 10;

        public boolean hasDepth = false;

        public SensorFilenameSet(String sensorName, byte sensorType)
        {
            this.sensorName = sensorName;
            this.sensorType = sensorType;
            sensorFilenames = null;
        }

        public void addBinaryFilename(String filename)
        {
            sensorFilenames = (String[]) StsMath.arrayAddElement(sensorFilenames, filename, nFilenames, increment);
            nFilenames++;
        }

        public void addAsciiFilename(String filename)
        {
            sensorFilenames = (String[]) StsMath.arrayAddElement(sensorFilenames, filename, nFilenames, increment);
            nFilenames++;
        }

        public boolean binaryAlreadyInList(String filename)
        {
            for (int i = 0; i < nFilenames; i++)
            {
                if (sensorFilenames[i].equals(filename))
                {
                    return true;
                }
            }
            return false;
        }

        public boolean asciiAlreadyInList(String filename)
        {
            for (int i = 0; i < nFilenames; i++)
            {
                if (sensorFilenames[i].equals(filename))
                {
                    return true;
                }
            }
            return false;
        }

        public boolean isOK()
        {
                trimLogFilenames();
                StsMath.qsort(sensorFilenames);
                return true;
        }

        void trimLogFilenames()
        {
            sensorFilenames = (String[]) StsMath.trimArray(sensorFilenames, nFilenames);
        }
    }

    static public void setLocations(StsPoint[] locations, byte[] types)
    {
        sensorLocations = locations;
        sensorTypes = types;
    }

    static public StsSensor[] createSensors(StsModel model, StsProgressPanel progressPanel, StsSensorFile selectedFile, boolean loadMultiStage, byte vUnits, byte hUnits,
                           StsAbstractSensorFactory sensorFactory, boolean computeAttributes)
    {
        StsSensor sensor = null;
        try
        {

            if (selectedFile == null)
            {
                return null;
            }

            progressPanel.appendLine("Preparing to load sensors from " + selectedFile.file.getFilename() + "...");
            StsMessageFiles.logMessage("Preparing to load sensor from " + selectedFile.file.getFilename() + "...");

            StsProject project = model.getProject();

            //StsSensorKeywordIO.setFalseX(model.getXOrigin());
            //StsSensorKeywordIO.setFalseY(model.getYOrigin());

            // Reload from ASCII always until switch is provide for user to decide
            reloadAscii = true;
            multiStage = loadMultiStage;

            StsSensor[] sensors = StsSensorImport.createSensors(model, progressPanel, selectedFile, vUnits, hUnits, sensorFactory, computeAttributes);
            success = (sensors != null);
            if(!success)
            {
            	progressPanel.setDescriptionAndLevel("Failed to load sensor file " + selectedFile.file.getFilename(), StsProgressBar.ERROR);
            	return null;
            }

            for(int i=0; i< sensors.length; i++)
            {
               if(sensors[i] instanceof StsStaticSensor)
               {
                   ((StsStaticSensor)sensors[i]).setPosition(new double[] {selectedFile.getStaticX(), selectedFile.getStaticY(), selectedFile.getStaticZ()});
                   sensors[i].setIsVisible(false);
               }
               if(!sensors[i].addToProject(StsProject.TD_DEPTH))
                   return null;
            }
            progressPanel.appendLine("Loaded sensor file " + selectedFile.file.getFilename() + " ...");
            return sensors;
        }
        catch (Exception e)
        {
            progressPanel.setDescriptionAndLevel("Exception thrown.", StsProgressBar.ERROR);
            progressPanel.appendLine("Exception thrown: " + e.getMessage());
            StsException.outputException("StsSensorImport.createSensors() failed.", e, StsException.WARNING);
            return null;
        }
    }

    /** build the well entity with its well line displayed in the 3d window */
    static private StsSensor buildSensor(String sensorName, StsProgressPanel progressPanel, StsAbstractSensorFactory sensorFactory, byte type)
    {
        StsSensor sensor = sensorFactory.createSensorInstance(sensorName, type);

        // Add xyz, time and other vectors to sensor
        try
        {
            return sensor;
        }
        catch (Exception e)
        {
            StsException.outputException("StsSensorImport.buildSensor() failed.", e, StsException.WARNING);
            progressPanel.setDescriptionAndLevel("StsSensorImport.buildSensor() failed.\n", StsProgressBar.ERROR);
            return null;
        }
    }

    static final class SensorNameFilter implements FilenameFilter
    {
        String format = null;
        public SensorNameFilter(String fmt)
        {
            this.format = fmt;
        }

        public boolean accept(File dir, String name)
        {
            if (!name.startsWith("sensor"))
            {
                return false;
            }
            if(format.equals("bin") && name.startsWith("sensor"))
            {
                StsWellKeywordIO.parseBinaryFilename(name);
            }
            else
            {
                StsWellKeywordIO.parseAsciiFilename(name);
            }
            return StsWellKeywordIO.format.equals(format);
        }
    }
}