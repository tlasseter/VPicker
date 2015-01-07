//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Import;

import com.Sts.Actions.Wizards.OSWell.*;
import com.Sts.DBTypes.OpenSpirit.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;
import com.Sts.Utilities.DateTime.CalendarParser;

import java.io.*;
import java.util.*;

public class StsWellImport
{
    static StsModel model;
    static String currentDirectory = null;
    static String binaryDataDir;
    static float logNull;
    public static int nLoaded = 0;
    static boolean success = false;
    static String fullBinaryDirName;
    static TreeMap wellFilenameSets;
    static StsAbstractFileSet[] fileSets = new StsAbstractFileSet[7];
    static StsAbstractFileSet fileSet;
    static String jarFilename = "wells.jar";

    private transient double progressValue = 0.0;
    private transient String progressDescription = "";

    static String[] jarFiles = null;
    static String[] asciiFiles = null;
    static private String[] selectedWells = null;
    static private StsFile[] selectedFiles = null;
    static private float datumShift = 0.0f;
    static private int[] selectedIdxs = null;
    static private int dateOrder = CalendarParser.DD_MM_YY;
    static private StsPoint[] topHoleLocations = null;
    static private byte[] wellTypes = null;
    static private boolean reloadAscii = false;
    static private boolean overrideFilter = false;

    static final String WELL_DEV_PREFIX = StsLogVector.WELL_DEV_PREFIX;
    static final String WELL_LOG_PREFIX = StsLogVector.WELL_LOG_PREFIX;
    static final String WELL_TD_PREFIX = StsLogVector.WELL_TD_PREFIX;
    static final String WELL_REF_PREFIX = "well-ref";
    static final String WELL_PERF_PREFIX = "well-perf";
    static final String WELL_FMI_PREFIX = "well-fmi";
    static final String WELL_EQUIP_PREFIX = "well-equipment";
    static public final int DEV = 0;
    static public final int LOG = 1;
    static public final int TD = 2;
    static public final int REF = 3;
    static public final int PERF = 4;
    static public final int FMI = 5;
    static public final int EQUIPMENT = 6;
    static public final int HDR = 5;
    static String[] filePrefixs = new String[] {WELL_DEV_PREFIX, WELL_LOG_PREFIX, WELL_TD_PREFIX, WELL_REF_PREFIX, WELL_PERF_PREFIX, WELL_FMI_PREFIX, WELL_EQUIP_PREFIX};
    static boolean[] applyShift = new boolean[] {false, false, false, false, false, false, false };

    static public final byte WEBJARFILE = 0;
    static public final byte JARFILE = 1;
    static public final byte BINARYFILES = 2;
    static public final byte ASCIIFILES = 3;
    static public final byte LASFILES = 4;
    static public final byte UTFILES = 5;
    static public final byte GEOGRAPHIX = 6;

    static byte vUnits = StsParameters.DIST_NONE;
    static byte hUnits = StsParameters.DIST_NONE;
    static byte binaryHorzUnits = StsParameters.DIST_NONE;
    static byte binaryVertUnits = StsParameters.DIST_NONE;

    static public final boolean debug = false;
    static StsTimer timer = new StsTimer();
    static final String[] typeNames = new String[]
        {"WEBJARFILE", "JARFILE", "BINARYFILES", "ASCIIFILES", "LASFILES", "UTFILES"};

    static public void initialize(StsModel model_)
    {
        model = model_;
        StsProject project = model.getProject();
        currentDirectory = project.getRootDirString();
        binaryDataDir = project.getBinaryFullDirString();
        logNull = project.getLogNull();
        nLoaded = 0;
        overrideFilter = false;
        resetFilePrefixs();
    }

    static public String getCurrentDirectory()
    {
        return currentDirectory;
    }

    static public void setCurrentDirectory(String dirPath)
    {
        currentDirectory = dirPath;
    }

    /** Create wells
     *
     */
    static public StsWell[] createWells(StsProgressPanel progressPanel, String[] selectedWellnames, double[][] wellXYs, boolean[] applyKbs, byte vUnits, byte hUnits,
                                        StsPoint[] topHolePoints, byte[] utTypes, StsAbstractWellFactory wellFactory)
    {
        StsWell well = null;

        String[] filenames = null;
        StsAbstractFile[] logFiles = null;

        try
        {
            int nWells = selectedWellnames.length;
            int utCounter = 0;
            StsWell[] wells = new StsWell[0];
            progressPanel.initialize(nWells);
            for (int n = 0; n < nWells; n++)
            {
                String wellname = selectedWellnames[n];
                StsPoint topHolePoint = null;
                byte utType = StsUTKeywordIO.MINCURVE;

                // Set the right fileSet
                WellFilenameSet wellFilenameSet = (WellFilenameSet) wellFilenameSets.get(wellname);
                filenames = wellFilenameSet.devFilenames;

                // read and build a well (if we don't already have it)
                if ( (model.getObjectWithName(StsWell.class, wellname) == null) || reloadAscii)
                {
                    StsWell wellTemp = (StsWell) model.getObjectWithName(StsWell.class, wellname);
                    if ((wellTemp != null) && reloadAscii)
                    {
                        wellTemp.delete();
                    }

                    byte formatType = wellFilenameSet.welltype;
                    if(formatType == ASCIIFILES || formatType == LASFILES || formatType == UTFILES)
                    {
                        if (topHolePoints != null)
                        {
                            if (formatType == (int) UTFILES)
                            {
                                topHolePoint = topHolePoints[utCounter];
                                utType = utTypes[utCounter];
                                utCounter++;
                            }
                        }
                        well = createWellFromAsciiFiles(wellFilenameSet, currentDirectory, binaryDataDir, topHolePoint, utType, logNull,
                            progressPanel, wellFactory,  vUnits, hUnits);
                        if (well != null)
                        {
                            wells = (StsWell[])StsMath.arrayAddElement(wells, well);
                            nLoaded++;
                        }
                    }
                    else if(formatType == GEOGRAPHIX)
                    {
                        StsWell[] geographixWells = createWellFromGeographixFiles(wellFilenameSet, currentDirectory, binaryDataDir, wellXYs, applyKbs, logNull, progressPanel, wellFactory, vUnits, hUnits);
                        if(geographixWells != null)
                        {
                             wells = (StsWell[])StsMath.arrayAddArray(wells, geographixWells);
                             nLoaded = nLoaded + geographixWells.length;
                        }
                    }
                }
                else
                {
                    if (progressPanel != null)
                    {
                        progressPanel.appendLine("Well: " + wellname + " already loaded...");
                        progressPanel.setDescriptionAndLevel("Well: " + wellname + " already loaded...", StsProgressBar.WARNING);
                    }
                }
                progressPanel.setValue(n+1);
                progressPanel.setDescription("Loaded well #" + nLoaded + " of " + nWells);
            }
            return wells;

        }
        catch (Exception e)
        {
            StsException.outputException("StsWellWizard.createWells() failed.", e, StsException.WARNING);
            progressPanel.setDescriptionAndLevel("StsWellWizard.createWells() failed.\n", StsProgressBar.WARNING);
            return null;
        }
    }

    /** read in multiple well deviations from a list of Ascii files */
    static public StsWell createWellFromAsciiFiles(WellFilenameSet wellFilenameSet, String dataDir,
        String binDir, StsPoint point, byte utType, float logNull,
        StsProgressPanel progressPanel, StsAbstractWellFactory wellFactory, byte vUnits, byte hUnits)
    {
        StsWell well = null;
        String wellname = wellFilenameSet.wellname;

        if (dataDir == null)
        {
            return null;
        }

        try
        {
            StsLogVector[] xyzVectors = null;
            StsLogVector[] logVectors = null;
            StsLogCurve[] logCurves = null;
            StsLogCurve tdCurve = null;
            StsLogCurve[] tdCurves = null;

            boolean deleteBinaries = reloadAscii;
            
            if (progressPanel != null) progressPanel.appendLine("Processing S2S formatted well: " + wellname + "...");
            
            // Process TD Curves --- Currently only supports S2S Format.
            float shift = 0.0f;
            if(applyShift[TD])
                shift = datumShift;

            StsWellKeywordIO.parseAsciiFilename(wellFilenameSet.devFilenames[0]);
            tdCurves = StsWellKeywordIO.readLogCurves(well, wellname, dataDir, binDir, filePrefixs[TD], logNull, deleteBinaries, vUnits, shift, dateOrder);
            if (tdCurves.length > 0)
            {
            	if (progressPanel != null) progressPanel.appendLine("   Successfully processed txt formatted time-depth curve for well " + wellname + "...");
                tdCurve = tdCurves[0];
            }
            else
            {
                // Try LAS formated TD Curves
                tdCurves = StsLasKeywordIO.readLogCurves(well, wellname, dataDir, binDir, filePrefixs[TD],
                    logNull, deleteBinaries, vUnits, shift);
                if (tdCurves.length > 0)
                {
                	if (progressPanel != null) progressPanel.appendLine("   Successfully processed las formatted time-depth curve for well " + wellname + "...");              	
                    tdCurve = tdCurves[0];
                }
            }

            // Process Deviation Curves
            shift = 0.0f;
            if(applyShift[DEV])
                shift = datumShift;
            StsWellKeywordIO.parseAsciiFilename(wellFilenameSet.devFilenames[0]);
            Object fileHeader = null;
            String format = StsWellKeywordIO.format;
            if (format.equals("txt"))
            {
                xyzVectors = StsWellKeywordIO.readDeviation(dataDir, binDir, wellFilenameSet.devFilenames[0],
                    logNull, deleteBinaries, vUnits, hUnits, filePrefixs[DEV], shift);
                fileHeader = StsWellKeywordIO.getFileHeader();
                wellFactory.setIsLive(StsWellKeywordIO.fileHeader.liveWell);
            }
            else if (format.equals("ut"))
            {
                xyzVectors = StsUTKeywordIO.readDeviation(dataDir, binDir, wellFilenameSet.devFilenames[0], point,
                    utType, logNull, deleteBinaries, vUnits, hUnits, filePrefixs[DEV], shift);
                fileHeader = StsUTKeywordIO.getFileHeader();
                wellFactory.setIsLive(false);
            }
            else if (format.equals("las"))
            {
                xyzVectors = StsLasKeywordIO.readDeviation(dataDir, binDir, wellFilenameSet.devFilenames[0],
                    logNull, deleteBinaries, vUnits, hUnits, filePrefixs[DEV], shift);
                fileHeader = StsLasKeywordIO.getFileHeader();
                wellFactory.setIsLive(false);
            }
            if (xyzVectors == null)
            {
            	if (progressPanel != null) progressPanel.appendLine("Failed to read deviation vectors (" + wellname + "), check file header.");
            	if (progressPanel != null) progressPanel.setDescriptionAndLevel("Failed to load " + wellname + ".", StsProgressBar.ERROR);
                return null;
            }
            else
            	if (progressPanel != null) progressPanel.appendLine("   Successfully processed " + format + " formatted deviation file for well " + wellname + "...");

            // Build the Well
            well = buildWell(wellname, xyzVectors, tdCurve, progressPanel, wellFactory);

            if (well == null)
            {
            	if (progressPanel != null) progressPanel.appendLine("Failed to build well from supplied information (" + wellname + "), check file contents.");
            	if (progressPanel != null) progressPanel.setDescriptionAndLevel("Failed to build" + wellname + ".", StsProgressBar.ERROR);
                return null;
            }

            if(StsWellKeywordIO.fileHeader != null)
            {
                if(StsWellKeywordIO.fileHeader.liveWell)
                    ((StsLiveWell)well).setTimeVector(StsWellKeywordIO.getTimeVector());
            }
            // Add the TD Curve
            if (tdCurve != null)
            {
                well.addLogCurve(tdCurve);
            }

            // Add Drift and Azimuth as Logs if UT Format
            StsWellKeywordIO.parseAsciiFilename(wellFilenameSet.devFilenames[0]);
            format = StsWellKeywordIO.format;
            if (format.equals("ut"))
            {
                logCurves = StsUTKeywordIO.addLogCurves(well, xyzVectors, well.getName(), logNull);
                fileHeader = StsUTKeywordIO.fileHeader;
                if (logCurves != null)
                {
                    well.addLogCurves(logCurves);
                    if (progressPanel != null) progressPanel.appendLine("   Successfully processed drift & azimuth from UT deviation file for well " + wellname + "...");
                }
            }
            // Set units to the ones from the binaries, whether just created or old.
            vUnits = StsLogVector.getVectorOfType(xyzVectors, StsLogVector.DEPTH).getUnits();
            hUnits = StsLogVector.getVectorOfType(xyzVectors, StsLogVector.X).getUnits();

            // Process Log Curves
            shift = 0.0f;
            if(applyShift[LOG])
                shift = datumShift;
			if(wellFilenameSet.logFilenames != null && wellFilenameSet.logFilenames.length > 0)
			{
                if(wellFilenameSet.logFilenames[0] != null)
                {
                    StsWellKeywordIO.parseAsciiFilename(wellFilenameSet.logFilenames[0]);
                    format = StsWellKeywordIO.format;
                    if (format.equals("txt"))
                    {
                        logCurves = StsWellKeywordIO.readLogCurves(well, wellname, dataDir, binDir, filePrefixs[LOG], logNull, deleteBinaries, vUnits, shift, dateOrder);
                        fileHeader = StsWellKeywordIO.fileHeader;
                    }
                    else if (format.equals("las"))
                    {
                        logCurves = StsLasKeywordIO.readLogCurves(well,  wellname, dataDir, binDir,
                                    filePrefixs[LOG], logNull, deleteBinaries, vUnits, shift);
                        fileHeader = StsLasKeywordIO.fileHeader;
                    }
                    if (logCurves.length == 0)
                    {
                        return well;
                    }
                    well.addLogCurves(logCurves);
                    if (progressPanel != null) progressPanel.appendLine("   Successfully processed " + logCurves.length + " logs from file(s) for well " + wellname + "...");                    
                }
			}
            // If a file header exists (usually in the log file since it is the only standard format, store the values

            if(fileHeader != null)
            {
                if(fileHeader instanceof StsLasKeywordIO.FileHeader)
                {
                    StsLasKeywordIO.FileHeader header = (StsLasKeywordIO.FileHeader)fileHeader;
                    well.setCompany(header.company);
                    well.setOperator(header.operator);
                    well.setDate(header.spudDate);
                    well.setApi(header.uwiNumber);
                    well.setUwi(header.wellIdNumber);
                    well.setField(header.field);
                    well.setWellLabel(header.wellLabel);
                    well.setWellNumber(header.wellNumber);
                }
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsWellWizard failed.", e, StsException.WARNING);
            StsMessageFiles.logMessage("Unable to read well deviation files.");
            if (progressPanel != null) progressPanel.setDescriptionAndLevel("Unable to read well deviation files.\n", StsProgressBar.ERROR);

            return null;
        }
        return well;
    }

    /** read in multiple well deviations from a list of Ascii files */
    static public StsWell[] createWellFromGeographixFiles(WellFilenameSet wellFilenameSet, String dataDir, String binDir, double[][] wellXYs, boolean[] applyKb, float logNull,
                                StsProgressPanel progressPanel, StsAbstractWellFactory wellFactory, byte vUnits, byte hUnits)
    {
        StsWell well = null;
        StsWell[] wells = new StsWell[0];
        StsGeographixKeywordIO.WellHeader wellHeader = null;
        BufferedReader bufRdr = null;

        if (dataDir == null)
        {
            return null;
        }

        try
        {
            StsLogVector[] xyzVectors = null;
            StsLogVector[] logVectors = null;
            StsLogCurve[] logCurves = null;
            StsLogCurve tdCurve = null;
            StsLogCurve[] tdCurves = null;

            boolean deleteBinaries = reloadAscii;

            if (progressPanel != null) progressPanel.appendLine("Processing Geographix formatted well: " + wellFilenameSet.wellname + "...");
            int numberOfWells = StsGeographixKeywordIO.getFileHeader().getNumWellsInFile();
            progressPanel.initialize(numberOfWells);
            if(selectedIdxs == null)
            {
                selectedIdxs = new int[numberOfWells];
                for(int i=0; i<numberOfWells; i++)
                     selectedIdxs[i] = i;
            }
            float shift = 0.0f;
            if(applyShift[DEV])
                shift = datumShift;

            //int nWells = 0;


            for(int i=0; i<selectedIdxs.length; i++)
            {
                bufRdr = StsGeographixKeywordIO.readWellHeader(dataDir, wellFilenameSet.devFilenames[0], selectedIdxs[i]);
                wellHeader = StsGeographixKeywordIO.getWellHeader();
                if(wellHeader == null)
                    return wells;

                if(wellHeader != null)
                {
                    wellHeader.xOrigin = wellXYs[selectedIdxs[i]][0];
                    wellHeader.yOrigin = wellXYs[selectedIdxs[i]][1];
                    wellHeader.datumElevation = (float)wellXYs[selectedIdxs[i]][2];
                    wellHeader.applyDatum = applyKb[selectedIdxs[i]];
                }
                wellFilenameSet.wellname = StsGeographixKeywordIO.getWellName();
                String wellname = wellFilenameSet.wellname;

                // Read the TD curves
                //tdCurves = StsGeographixKeywordIO.readLogCurves(well, wellname, StsGeographixKeywordIO.TD, dataDir, binDir, logNull, deleteBinaries, vUnits, shift);
                //if (tdCurves.length > 0)
                //{
                //	if (progressPanel != null) progressPanel.appendLine("   Successfully processed time-depth curve for well " + wellname + "...");
                //    tdCurve = tdCurves[0];
                //}
                // Read XYZ Vectors (O records)
                xyzVectors = StsGeographixKeywordIO.readDeviation(bufRdr, dataDir, binDir, wellname, wellFilenameSet.devFilenames[0], logNull, deleteBinaries, vUnits, hUnits, shift);
                if (xyzVectors == null)
                {
                    if(wells.length == 0)
                    {
            	        if (progressPanel != null)
                        {
                            progressPanel.appendLine("Failed to read deviation vectors (" + wellname + "), check file header.");
                            progressPanel.setDescriptionAndLevel("Failed to load " + wellname + ".", StsProgressBar.ERROR);
                        }
                    }
                    if(bufRdr != null)
                        bufRdr.close();
                    return wells;
                }
                else
            	    if (progressPanel != null) progressPanel.appendLine("   Successfully processed Geographix formatted deviation file for well " + wellname + "...");

                vUnits = StsLogVector.getVectorOfType(xyzVectors, StsLogVector.DEPTH).getUnits();
                hUnits = StsLogVector.getVectorOfType(xyzVectors, StsLogVector.X).getUnits();

                // Build Well
                well = buildWell(wellname, xyzVectors, tdCurve, progressPanel, wellFactory);
                if (well == null)
                {
                    if(wells.length == 0)
                    {
            	        if (progressPanel != null) progressPanel.appendLine("Failed to build well from supplied information (" + wellname + "), check file contents.");
            	        if (progressPanel != null) progressPanel.setDescriptionAndLevel("Failed to build" + wellname + ".", StsProgressBar.ERROR);
                    }
                    if(bufRdr != null)
                        bufRdr.close();
                    return wells;
                }
                progressPanel.setValue(i);
                progressPanel.setDescription("Processed deviations for well #" + i + " of " + selectedIdxs.length);
                // Add the TD Curve
                if (tdCurve != null)
                    well.addLogCurve(tdCurve);

                // Read and add Log Curves
                //logCurves = StsGeographixKeywordIO.readLogCurves(well,  wellname, StsGeographixKeywordIO.LOGS, dataDir, binDir, logNull, deleteBinaries, vUnits, shift);
                //if (logCurves.length == 0)
                //    return well;
                //else
                //    well.addLogCurves(logCurves);

                // Process Header
                if(wellHeader != null)
                {
                    well.setCompany(wellHeader.operator);
                    well.setField(wellHeader.field);
                    well.setOperator(wellHeader.operator);
                    well.setSpudDate(wellHeader.spudDate.getTime());
                    well.setCompletionDate(wellHeader.completionDate.getTime());
                    well.setApi(wellHeader.altIdNumber);
                    well.setUwi(wellHeader.uwiNumber);
                    well.setState(wellHeader.state);
                    well.setCounty(wellHeader.county);
                    well.setField(wellHeader.field);
                    well.setWellLabel(wellHeader.wellName);
                    well.setWellNumber(wellHeader.wellNumber);
                    well.setElev(wellHeader.datumElevation);
                    well.setElevDatum(wellHeader.datumDescription);
                }
                wells = (StsWell[])StsMath.arrayAddElement(wells, well);
            }
            if(bufRdr != null)
                bufRdr.close();
        }
        catch (Exception e)
        {
            StsException.outputException("StsWellWizard failed.", e, StsException.WARNING);
            StsMessageFiles.logMessage("Unable to read well deviation files.");
            if (progressPanel != null) progressPanel.setDescriptionAndLevel("Unable to read well deviation files.\n", StsProgressBar.ERROR);

            return wells;
        }
        return wells;
    }

    static public boolean archiveWells(String[] wellNames, String archiveDir)
    {
        try
        {
            for(int i=0; i<wellNames.length; i++)
            {
                WellFilenameSet wellFilenameSet = (WellFilenameSet) wellFilenameSets.get(wellNames[i]);
                for(int j=0; j<wellFilenameSet.devFilenames.length; j++)
                {
                    if(wellFilenameSet.devFilenames[j] != null)
                    {
                        if (wellFilenameSet.devFilenames[j].indexOf("bin") == -1)
                            StsFile.copy(currentDirectory + File.separator + wellFilenameSet.devFilenames[j],
                                         archiveDir + File.separator + wellFilenameSet.devFilenames[j]);
//                       System.out.println("Archiving deviation file: " + wellFilenameSet.devFilenames[j]);
                    }
                }
                for(int j=0; j<wellFilenameSet.nLogFilenames; j++)
                {
                    if(wellFilenameSet.logFilenames[j] != null)
                    {
                        if(wellFilenameSet.logFilenames[j].indexOf("bin") == -1)
                            StsFile.copy(currentDirectory + File.separator + wellFilenameSet.logFilenames[j],
                                         archiveDir + File.separator + wellFilenameSet.logFilenames[j]);
//                            System.out.println("Archiving log file: " + wellFilenameSet.logFilenames[j]);
                    }
                }
                for(int j=0; j<wellFilenameSet.nRefFilenames; j++)
                {
                    if(wellFilenameSet.refFilenames[j] != null)
                    {
                        if(wellFilenameSet.refFilenames[j].indexOf("bin") == -1)
                            StsFile.copy(currentDirectory + File.separator + wellFilenameSet.refFilenames[j],
                                         archiveDir + File.separator + wellFilenameSet.refFilenames[j]);
 //                           System.out.println("Archiving ref file: " + wellFilenameSet.refFilenames[j]);
                    }
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsWellImport.archiveWells Error", e, StsException.WARNING);
            return false;
        }
        return true;
    }

    static public String[] getWellFilenames()
    {
        String[] filenames = new String[0];
        WellFilenameSet wellFilenameSet = null;
        Iterator iter;

        if (debug)
        {
            timer.start();
        }

        int nWells = wellFilenameSets.size();
        if (nWells > 0)
        {
            filenames = new String[nWells];
            int i = 0;
            iter = wellFilenameSets.entrySet().iterator();
            while (iter.hasNext())
            {
                Map.Entry entry = (Map.Entry) iter.next();
                wellFilenameSet = (WellFilenameSet) entry.getValue();
                filenames[i] = wellFilenameSet.wellname;
                i++;
            }
        }
        if (debug)
        {
            timer.stopPrint("Time to find " + filenames.length + " files.");
        }
        return filenames;
    }

    static public String getWellFilename(int idx)
    {
        String filename = null;
        WellFilenameSet wellFilenameSet = null;
        Iterator iter;

        int nWells = wellFilenameSets.size();
        if (nWells > 0)
        {
            int i = 0;
            iter = wellFilenameSets.entrySet().iterator();
            while (iter.hasNext())
            {
                Map.Entry entry = (Map.Entry) iter.next();
                wellFilenameSet = (WellFilenameSet) entry.getValue();
                if (i == idx)
                {
                    return wellFilenameSet.devFilenames[0];
                }
                i++;
            }
        }
        return null;
    }

    static public String[] getSelectedWellFilenames(String[] selectedWellnames, int type, boolean includeBinaries)
    {
        String wellname = null;
        String[] filenames = null;
        StsWellImport.WellFilenameSet wellFilenameSet = null;
        String fmt = null;

//        String binDir = model.getProject().getBinaryFullDirString();
//        String asciiDir = currentDirectory;
        for (int n = 0; n < selectedWellnames.length; n++)
        {
            wellname = selectedWellnames[n];
            wellFilenameSet = (StsWellImport.WellFilenameSet) wellFilenameSets.get(wellname);
            if (wellFilenameSet.welltype == type)
            {
//                if(!StsUTKeywordIO.verifyDeviationBinaries(asciiDir, binDir,  wellFilenameSet.devFilenames[0])
//                   || selectWells.getReloadASCII())
                filenames = (String[]) StsMath.arrayAddElement(filenames, wellFilenameSet.devFilenames[0]);
            }
        }
        return filenames;
    }

    static public WellFilenameSet getWellFilenameSet(String wellname)
    {
        return (WellFilenameSet) wellFilenameSets.get(wellname);
    }

    /*
        static public void initFileSets()
        {
            fileSets = new StsAbstractFileSet[6];
        }
     */

    static public void initializeWellFilenameSets()
    {
        if (wellFilenameSets == null)
        {
            wellFilenameSets = new TreeMap();
        }
        else
        {
            wellFilenameSets.clear();
        }
    }
    static public void resetFilePrefixs()
    {
        filePrefixs = new String[] {WELL_DEV_PREFIX, WELL_LOG_PREFIX, WELL_TD_PREFIX, WELL_REF_PREFIX, WELL_PERF_PREFIX, WELL_FMI_PREFIX, WELL_EQUIP_PREFIX};
    }
    static public String getFilePrefix(int type) { return filePrefixs[type]; }
    static public void setFilePrefixs(String[] prefixs)
    {
        filePrefixs = prefixs;
    }
    static public void setApplyDatumShift(boolean[] apply)
    {
        applyShift = apply;
    }
    static public boolean constructWellFilenameSets(byte type)
    {
        return constructWellFilenameSets(null, type);
    }

    static public boolean constructWellFilenameSets(StsFilenameFilter filter, byte type)
    {
        String dirname = null;
        try
        {
            if (debug)
            {
                timer.start();
            }
            if(filter != null)
            {
                if (currentDirectory == null)
                {
                    dirname = model.getProject().getRootDirectory() + File.separator;
                }
                else
                {
                    dirname = currentDirectory + File.separator;
                }
                type = LASFILES;
                fileSets[type] = StsFileSet.constructor(dirname, new AsciiFilter("las", filePrefixs));
//                fileSets[type] = StsFileSet.constructor(dirname, filter);
            }
            else
            {
                switch (type)
                {
                    case WEBJARFILE:
                        fileSets[type] = StsWebStartJar.constructor(jarFilename);
                        break;

                    case JARFILE:
                        if (currentDirectory == null)
                        {
                            dirname = new String(model.getProject().getRootDirString() + "/");
                        }
                        else
                        {
                            dirname = currentDirectory + "/";
                        }
                        fileSets[type] = StsJar.constructor(dirname, jarFilename);
                        break;

                    case BINARYFILES:
                        if (currentDirectory == null)
                        {
                            fullBinaryDirName = model.getProject().getBinaryFullDirString();
                        }
                        else
                        {
                            fullBinaryDirName = currentDirectory + File.separator;
                        }
                        fileSets[type] = StsFileSet.constructor(fullBinaryDirName, new WellnameFilter("bin"));
                        break;

                    case ASCIIFILES:
                        if (currentDirectory == null)
                        {
                            dirname = model.getProject().getRootDirectory() + File.separator;
                        }
                        else
                        {
                            dirname = currentDirectory + File.separator;
                        }
                        fileSets[type] = StsFileSet.constructor(dirname, new WellnameFilter("txt"));
                        break;
                    case LASFILES:
                        if (currentDirectory == null)
                        {
                            dirname = model.getProject().getRootDirectory() + File.separator;
                        }
                        else
                        {
                            dirname = currentDirectory + File.separator;
                        }
                        fileSets[type] = StsFileSet.constructor(dirname, new AsciiFilter("las", filePrefixs));
                        break;
                    case UTFILES:
                        if (currentDirectory == null)
                        {
                            dirname = model.getProject().getRootDirectory() + File.separator;
                        }
                        else
                        {
                            dirname = currentDirectory + File.separator;
                        }
                        fileSets[type] = StsFileSet.constructor(dirname, new AsciiFilter("ut", filePrefixs));
                        break;
                    case GEOGRAPHIX:
                        if (currentDirectory == null)
                        {
                            dirname = model.getProject().getRootDirectory() + File.separator;
                        }
                        else
                        {
                            dirname = currentDirectory + File.separator;
                        }
                        fileSets[type] = StsFileSet.constructor(dirname, new SuffixAsciiFilter(new String[] {"wba", "wls"}));
                        break;

                }
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
            return addWellFilenameSets(filenames, type);
        }
        catch (Exception e)
        {
            StsException.outputException("StsWellWizard.constructWellFilenameSets() failed.", e, StsException.WARNING);
            return false;
        }
    }

    static public boolean selectedOfType(int type, String[] selectedWellnames)
    {
        String wellname = null;
        String[] filenames = null;
        StsWellImport.WellFilenameSet wellFilenameSet = null;
        for (int n = 0; n < selectedWellnames.length; n++)
        {
            wellname = selectedWellnames[n];
            wellFilenameSet = (StsWellImport.WellFilenameSet) wellFilenameSets.get(wellname);
            if (wellFilenameSet.welltype == type)
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

    static public void setOverrideFilter(boolean override)
    {
        overrideFilter = override;
    }

    static public boolean getOverrideFilter()
    {
        return overrideFilter;
    }

    static public boolean getReloadAscii()
    {
        return reloadAscii;
    }

    static final class WellnameFilter implements FilenameFilter
    {
        String format = null;

        public WellnameFilter(String fmt)
        {
            this.format = fmt;
        }

        public boolean accept(File dir, String name)
        {
            if (!name.startsWith("well"))
            {
                return false;
            }
            if (format.equals("txt"))
            {
                StsWellKeywordIO.parseAsciiFilename(name);
            }
            else
            {
                StsWellKeywordIO.parseBinaryFilename(name);
            }
            String group = StsWellKeywordIO.group;
            boolean value = StsWellKeywordIO.format.equals(format) && (group.equals(filePrefixs[DEV])
                || group.equals(filePrefixs[LOG]) || group.equals(filePrefixs[TD]) || group.equals(filePrefixs[REF])
                || group.equals(filePrefixs[PERF]) || group.equals(filePrefixs[FMI]) || group.equals(filePrefixs[EQUIPMENT]));
            return value;
        }
    }

    static final class AsciiFilter implements FilenameFilter
    {
        String format = null;
        String[] prefixs = null;

        public AsciiFilter(String fmt, String[] prefixs)
        {
            this.format = fmt;
            this.prefixs = prefixs;
        }

        public boolean accept(File dir, String name)
        {
            StsWellKeywordIO.parseAsciiFilename(name);
            String group = StsWellKeywordIO.group;
            return StsWellKeywordIO.format.equals(format) && (group.equals(prefixs[DEV])
                || group.equals(prefixs[LOG]) || group.equals(prefixs[TD]));
        }
    }

    static final class SuffixAsciiFilter implements FilenameFilter
    {
        String[] suffixs = null;

        public SuffixAsciiFilter(String[] suffixs)
        {
            this.suffixs = suffixs;
        }

        public boolean accept(File dir, String name)
        {
            for(int i=0; i<suffixs.length; i++)
            {
                if(name.endsWith(suffixs[i]))
                    return true;
            }
            return false;
        }
    }

    static final class binaryWellnameFilter implements FilenameFilter
    {
        public binaryWellnameFilter()
        {}

        public boolean accept(File dir, String name)
        {
            StsWellKeywordIO.parseBinaryFilename(name);
            String group = StsWellKeywordIO.group;
            return StsWellKeywordIO.format.equals("bin") && (group.equals(WELL_DEV_PREFIX) || group.equals(WELL_LOG_PREFIX));
        }
    }

    static public void compressWellSet()
    {
        Iterator iter;
        WellFilenameSet wellFilenameSet;
        if (debug)
        {
            timer.start();
        }

        iter = wellFilenameSets.entrySet().iterator();
        while (iter.hasNext())
        {
            Map.Entry entry = (Map.Entry) iter.next();
            wellFilenameSet = (WellFilenameSet) entry.getValue();
            if (!wellFilenameSet.isOK())
            {
                iter.remove();
                continue;
            }

            // If we don't want to reload it, delete it from the list of available wells.
            StsWell existingWell = (StsWell) model.getObjectWithName(StsWell.class, wellFilenameSet.wellname);
            if (existingWell != null)
            {
                if(!reloadAscii)
                {
                    iter.remove();
                }
            }
        }
        if (debug)
        {
            timer.stopPrint("Time for compressWellSet.");
        }
        return;
    }

    static public boolean addWellFilenameSets(String[] filenames, byte type)
    {
        String filename = null;
        WellFilenameSet wellFilenameSet;

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

            if (wellFilenameSets == null)
            {
                wellFilenameSets = new TreeMap();
            }

            int nFilenames = filenames.length;
            for (int n = 0; n < nFilenames; n++)
            {
                filename = filenames[n];
                if ( (type == ASCIIFILES) || (type == LASFILES) || (type == UTFILES))
                {
                    StsWellKeywordIO.parseAsciiFilename(filename);
                }
                else if(type == GEOGRAPHIX)
                {
                    StsWellKeywordIO.name = filename.substring(0,filename.indexOf("."));
                }
                else
                {
                    StsWellKeywordIO.parseBinaryFilename(filename);
                }
                String wellname = StsWellKeywordIO.getWellName();
                if (wellFilenameSets.containsKey(wellname))
                {
                    wellFilenameSet = (WellFilenameSet) wellFilenameSets.get(wellname);
                    if ( (type == ASCIIFILES) || (type == LASFILES) || (type == UTFILES))
                    {
                        if (!wellFilenameSet.asciiAlreadyInList(filename, StsWellKeywordIO.group))
                        {
                            wellFilenameSet.addAsciiFilename(filename, StsWellKeywordIO.group, type, filePrefixs);
                        }
                    }

                    else
                    {
                        if (!wellFilenameSet.binaryAlreadyInList(filename, StsWellKeywordIO.group, StsWellKeywordIO.getCurveName()))
                        {
                            wellFilenameSet.addBinaryFilename(filename, StsWellKeywordIO.group, StsWellKeywordIO.getCurveName());
                        }
                    }
                }
                else
                {
                    wellFilenameSet = new WellFilenameSet(wellname, type);
                    wellFilenameSets.put(wellname, wellFilenameSet);
                    if ( (type == ASCIIFILES) || (type == LASFILES) || (type == UTFILES))
                    {
                        wellFilenameSet.addAsciiFilename(filename, StsWellKeywordIO.group, type, filePrefixs);
                    }
                    else if(type == GEOGRAPHIX)
                    {
                        wellFilenameSet.devFilenames[0] = filename;
                        wellFilenameSet.welltype = type;
                    }
                    else
                    {
                        wellFilenameSet.addBinaryFilename(filename, StsWellKeywordIO.group, StsWellKeywordIO.getCurveName());
                    }
                }
            }
            if (debug)
            {
                timer.stopPrint("Time to initializeWellFielnameSets for " + nFilenames + "files");
            }

            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsWellWizard.initializeWellFilenameSets() failed for file: " + filename + ".", e,
                                         StsException.WARNING);
            return false;
        }
    }

    static public class WellFilenameSet
    {
        public String wellname;
        public byte welltype = -1;
        public int nLogFilenames = 0;
        public String[] devFilenames;
        public String[] logFilenames = null;
        public int nRefFilenames = 0;
        public String[] refFilenames = null;
        public int nPerfFilenames = 0;
        public String[] perfFilenames = null;
        public int nFmiFilenames = 0;
        public String[] fmiFilenames = null;
        public int nEquipFilenames = 0;
        public String[] equipFilenames = null;
        public int increment = 10;

        public boolean hasDepth = false;

        public WellFilenameSet(String wellname, byte welltype)
        {
            this.wellname = wellname;
            this.welltype = welltype;
            devFilenames = new String[4];
            logFilenames = new String[4];
			refFilenames = new String[4];
			perfFilenames = new String[4];
			equipFilenames = new String[4];
			fmiFilenames = new String[4];
        }

        // add entries according to type; see StsLogVector for type definitions
        // types: 0 = X, 1 = Y, 2 = DEPTH, 3 = MDEPTH
        public void addBinaryFilename(String filename, String group, String curvename)
        {
            byte type = StsLogVector.getTypeFromString(StsWellKeywordIO.getCurveName());

            if (group.equals(WELL_DEV_PREFIX))
            {
                devFilenames[type] = filename;
            }
            else if (group.equals(WELL_LOG_PREFIX))
            {
                logFilenames = (String[]) StsMath.arrayAddElement(logFilenames, filename, nLogFilenames, increment);
                nLogFilenames++;
                if (type == StsLogVector.DEPTH || type == StsLogVector.MDEPTH)
                {
                    hasDepth = true;
                }
            }
        }

        public void addAsciiFilename(String filename, String group, byte type, String[] filePrefixs)
        {
//			System.out.println("add Ascii "+filename+" "+prefix);
            if (group.equals(filePrefixs[DEV]))
            {
                devFilenames[0] = filename;
                welltype = type;
            }
            else if (group.equals(filePrefixs[LOG]) || group.equals(filePrefixs[TD]))
            {
                logFilenames = (String[]) StsMath.arrayAddElement(logFilenames, filename, nLogFilenames, increment);
                nLogFilenames++;
            }
            else if (group.equals(filePrefixs[REF]))
            {
                refFilenames = (String[]) StsMath.arrayAddElement(refFilenames, filename, nRefFilenames, increment);
                nRefFilenames++;
            }
            else if (group.equals(filePrefixs[PERF]))
            {
                perfFilenames = (String[]) StsMath.arrayAddElement(perfFilenames, filename, nPerfFilenames, increment);
                nPerfFilenames++;
            }
            else if (group.equals(filePrefixs[FMI]))
            {
                fmiFilenames = (String[]) StsMath.arrayAddElement(fmiFilenames, filename, nFmiFilenames, increment);
                nFmiFilenames++;
            }
            else if (group.equals(filePrefixs[EQUIPMENT]))
            {
                equipFilenames = (String[]) StsMath.arrayAddElement(equipFilenames, filename, nEquipFilenames, increment);
                nEquipFilenames++;
            }
        }

        public boolean binaryAlreadyInList(String filename, String group, String curvename)
        {
            byte type = StsLogVector.getTypeFromString(StsWellKeywordIO.getCurveName());
            if (group.equals(WELL_DEV_PREFIX))
            {
                if (devFilenames[type] != null)
                {
                    return true;
                }
            }
            else if (group.equals(WELL_LOG_PREFIX))
            {
                for (int i = 0; i < nLogFilenames; i++)
                {
                    if (logFilenames[i].equals(filename))
                    {
                        return true;
                    }
                }
            }
            return false;
        }

        public boolean asciiAlreadyInList(String filename, String group)
        {
            if (group.equals(WELL_DEV_PREFIX))
            {
                if (devFilenames[0] != null)
                {
                    return true;
                }
                else if (group.equals(WELL_LOG_PREFIX))
                {
                    if (logFilenames[0].equals(filename))
                    {
                        return true;
                    }
                }
            }
            return false;
        }

        public boolean isOK()
        {
            if ( (welltype == ASCIIFILES) || (welltype == LASFILES) || (welltype == UTFILES) || (welltype == GEOGRAPHIX))
            {
                if (devFilenames[0] == null)
                {
                    return false;
                }
                return true;
            }
            else
            {
                // check deviation vectors
                if (devFilenames[0] == null || devFilenames[1] == null || (devFilenames[2] == null && devFilenames[3] == null))
                {
                    return false;
                }

                if (!hasDepth)
                {
                    logFilenames = new String[0];
                    nLogFilenames = 0;
                    return false;
                }
                trimLogFilenames();
                StsMath.qsort(logFilenames);
                return true;
            }
        }

        void trimLogFilenames()
        {
            logFilenames = (String[]) StsMath.trimArray(logFilenames, nLogFilenames);
        }
    }

    static public final class SuffixFilenameFilter implements FilenameFilter
    {
        String suffix = null;
        public SuffixFilenameFilter(String suffix)
        {
            this.suffix = suffix;
        }

        public boolean accept(File dir, String name)
        {
            return name.endsWith(suffix);
        }
    }

    /*
     public StsWell createWellFromBinaryFile(String wellname, String dataDir, String binDir)
     {
      try
      {
       StsLogVector[] xyzmtVectors = null;

       xyzmtVectors = getBinaryWellDeviationVectors(wellname);

                StsLogVector vector = StsLogVector.getVectorOfType(xyzmtVectors, StsLogVector.X);
                vUnits = vector.getVerticalUnits();
                vector = StsLogVector.getVectorOfType(xyzmtVectors, StsLogVector.MDEPTH);
                hUnits = vector.getHorizontalUnits();

                if (xyzmtVectors == null)
                {
                    return null;
                }
                StsLogVector tVector = StsLogVector.getVectorOfType(xyzmtVectors, StsLogVector.TIME);
                StsWell well = buildWell(wellname, xyzmtVectors, null);
                if (well == null) return null;
                well.addToProject();

                StsLogVector[] logCurveVectors = getBinaryLogCurveVectors(well);
                vUnits = logCurveVectors[0].getVerticalUnits();
                hUnits = logCurveVectors[0].getHorizontalUnits();

                if (logCurveVectors != null)
                {
                    well.constructLogCurvesCheckVersions(logCurveVectors, StsParameters.nullValue);
                }
                else if (!Main.isJarFile)
                {
//             StsWellKeywordIO.classInitialize(model);
     StsLogCurve[] logCurves = StsWellKeywordIO.readLogCurves(well, wellname, dataDir, binDir, StsLogVector.WELL_LOG_PREFIX, StsParameters.nullValue);
                    well.addLogCurves(logCurves);
                }
                loadWells.panel.setText("Processing binary well: " + wellname + "...", false);

                return well;
            }
            catch (Exception e)
            {
                StsException.outputException("StsWellWizard.createWellsFromBinaryJar() failed.", e, StsException.WARNING);
                return null;
            }
        }
     */
    // read the 4 logVectors: X, Y, DEPTH, MDEPTH.  First 3 are required.
    /*	private StsLogVector[] getBinaryWellDeviationVectors(String wellname)
     {
      try
      {
       WellFilenameSet wellFilenameSet = (WellFilenameSet) wellFilenameSets.get(wellname);
       if (wellFilenameSet == null)
       {
        return null;
       }
       String[] filenames = wellFilenameSet.devFilenames;
       StsLogVector[] logVectors = getBinaryWellVectors(wellname, filenames);
       if (!StsLogVector.deviationVectorsOK(logVectors))
       {
        return null;
       }
       return logVectors;
      }
      catch (Exception e)
      {
       StsException.outputException("StsWellWizard.getBinaryWellDeviationVectors() failed.", e, StsException.WARNING);
       return null;
      }
     }
     */
    /*	private StsLogVector[] getBinaryWellVectors(String wellname, String[] filenames)
     {
      InputStream is;
      URL url;
      URLConnection urlConnection;

      try
      {
       boolean ok = true;

       if (filenames == null)
       {
        return null;
       }

       int nFilenames = filenames.length;
       StsLogVector[] logVectors = new StsLogVector[nFilenames];

       for (int n = 0; n < nFilenames; n++)
       {
        String filename = filenames[n];
        if (filename == null)
        {
         continue;
        }
        StsAbstractFile file = fileSet.getFile(filename);
        logVectors[n] = new StsLogVector(filename);
        if (!logVectors[n].readBinaryFile(file, true))
        {
         ok = false;
        }
       }
       if (!ok)
       {
        return null;
       }
       return logVectors;
      }
      catch (Exception e)
      {
       StsException.outputException("StsWellWizard.getWellDeviationVectors() failed.", e, StsException.WARNING);
       return null;
      }
     }

     private StsLogVector[] getBinaryLogCurveVectors(StsWell well)
     {
      try
      {
       if (well == null)
       {
        return null;
       }
       String wellname = well.getName();
       WellFilenameSet wellFilenameSet = (WellFilenameSet) wellFilenameSets.get(wellname);
       if (wellFilenameSet == null)
       {
        return null;
       }
       String[] filenames = wellFilenameSet.logFilenames;
       return getBinaryWellVectors(wellname, filenames);
      }
      catch (Exception e)
      {
       StsException.outputException("StsWellWizard.getBinaryLogCurveVectors() failed.", e, StsException.WARNING);
       return null;
      }
     }
     */

    static public void setTopHoleLocations(StsPoint[] topPoints, byte[] types)
    {
        topHoleLocations = topPoints;
        wellTypes = types;
    }

    static public StsWell[] createWells(StsProgressPanel progressPanel, String[] selectedWellnames,  int[] selectedIndices, double[][] wellXYs, boolean[] applyKbs, StsAbstractWellFactory wellFactory, float shift, byte hUnits, byte vUnits, int order)
    {
        datumShift = shift;
        selectedIdxs = selectedIndices;
        dateOrder = order;
        return createWells(progressPanel, selectedWellnames, wellXYs, applyKbs, wellFactory, hUnits, vUnits);
    }

    static public StsWell[] createWells(StsProgressPanel progressPanel, String[] selectedWellnames, double[][] wellXYs, boolean[] applyKbs,
                                        StsAbstractWellFactory wellFactory, byte hUnits, byte vUnits)
    {
        try
        {
            if (selectedWellnames == null)
            {
                return null;
            }

            int nSelected = selectedWellnames.length;

            progressPanel.appendLine("Preparing to load " + nSelected + " well deviations ...");
            StsMessageFiles.logMessage("Preparing to load " + nSelected + " well deviations ...");

            StsProject project = model.getProject();
            String binaryDataDir = project.getBinaryFullDirString();

            // If the user has selected both binary files and ASCII files the units may be mixed, so reload all from ASCII
            String wellname = selectedWellnames[0];
            WellFilenameSet wellFilenameSet = (WellFilenameSet) wellFilenameSets.get(wellname);
            byte formatType = wellFilenameSet.welltype;
            if(formatType != GEOGRAPHIX)    // Number of binaries is not related to number of files in Geographix format
            {
                int nWithBinaries = 0;
                for (int n = 0; n < nSelected; n++)
                {
                    // If binaries exist for some and not for others. Load all from ASCII
                    if (StsWellKeywordIO.binaryFileExist(binaryDataDir, (String) selectedWellnames[n], StsLogVector.WELL_DEV_PREFIX,StsLogVector.X))
                    {
                        nWithBinaries++;
                    }
                }
                if (nWithBinaries != nSelected)
                {
                    progressPanel.appendLine("All binary wells will be re-loaded from ASCII to apply unit selections\n");
                    StsMessageFiles.infoMessage("All binary wells will be re-loaded from ASCII to apply unit selections");
                    reloadAscii = true;
                }
                else if (nWithBinaries == nSelected)
                {
            	    if (progressPanel != null) progressPanel.appendLine("   Ignoring unit selection since all wells are in binary format. Override with ReLoad from ASCII option");
                    StsMessageFiles.infoMessage("Ignoring unit selection since all wells are in binary format. Override with ReLoad from ASCII option"); ;
                }
            }
            StsWell[] wells = StsWellImport.createWells(progressPanel, selectedWellnames, wellXYs, applyKbs, vUnits, hUnits, topHoleLocations, wellTypes, wellFactory);
            success = (wells != null);

            // add lineVertex points to project to generate boundingBoxes
            int nLoaded = StsWellImport.nLoaded;
            for (int n = 0; n < nLoaded; n++)
            {
                boolean added = false;
                if(wells[n] instanceof StsLiveWell)
                    added = ((StsLiveWell)wells[n]).addToProject();
                else
                    added = wells[n].addToProject();
                if(!added) // this adds lineVertices to projectBoundingBoxes
                {
                    progressPanel.appendLine("Well loading aborted for well (" + wells[n].getName() + ")....");
                    break;
                }
            }
            project.adjustBoundingBoxes(true, false); // extend displayBoundingBox as needed and set cursor3d box accordingly
            project.checkAddUnrotatedClass(StsWell.class);
            project.rangeChanged();

            // create splined points to define lines
            // note that we don't computePoints until wells are added to project above and origin has been set
            for (int n = 0; n < nLoaded; n++)
            {
                wells[n].computePoints(); // generate splined points between vertices
            }
            progressPanel.appendLine("Loaded " + nLoaded + " well deviation files ...");
            
            progressPanel.appendLine("Preparing to load markers associated with loaded wells ...");
            progressPanel.initialize(nLoaded);
            for (int n = 0; n < nLoaded; n++)
            {
                progressPanel.appendLine("Loading markers for well " + wells[n].getName() + "...");
                if(formatType == GEOGRAPHIX)
                {
                    StsGeographixKeywordIO.readWellMarkers(wells[n], n, currentDirectory, wellFilenameSet.devFilenames[0], progressPanel);
                    StsGeographixKeywordIO.readPerforationMarkers(wells[n], n, currentDirectory, wellFilenameSet.devFilenames[0], progressPanel);
                    StsGeographixKeywordIO.readFMIMarkers(wells[n], n, currentDirectory, wellFilenameSet.devFilenames[0], progressPanel);
                    StsGeographixKeywordIO.readEquipmentMarkers(wells[n], n, currentDirectory, wellFilenameSet.devFilenames[0], progressPanel);
                }
                else
                {
                    StsWellImport.loadWellMarkers(wells[n], progressPanel);
                    StsWellImport.loadPerforationMarkers(wells[n], progressPanel);
                    StsWellImport.loadFMIMarkers(wells[n], progressPanel);
                    StsWellImport.loadEquipmentMarkers(wells[n], progressPanel);
                }
                progressPanel.appendLine("Completed load of markers for well " + wells[n].getName() + "...");
                progressPanel.setValue(n+1);
                progressPanel.setDescription("Processed markers for well #" + (n+1) + " of " + nLoaded);
            }
            StsSeismicVelocityModel velocityModel = project.getSeismicVelocityModel();
            if (velocityModel != null)
            {
                for (int n = 0; n < nLoaded; n++)
                {
					wells[n].adjustFromVelocityModel(velocityModel);
                    wells[n].computeMarkerTimesFromMDepth(velocityModel);
                }
            }
            //        model.win3d.resetCursorPanel();
            return wells;

        }
        catch (Exception e)
        {
            StsException.outputException("StsWellWizard.createWells() failed.", e, StsException.WARNING);
            return null;
        }
    }

    static public StsWell[] createFromOSWells(StsProgressPanel progressPanel,
    						StsOSWellDatastore datastore, Object[] osWells)
    {
        try
        {
            if (osWells == null)
            {
                return null;
            }

            int nSelected = osWells.length;

            progressPanel.appendLine("Preparing to load " + nSelected + " well from OpenSpirit ...");
            StsMessageFiles.logMessage("Preparing to load " + nSelected + " well from OpenSpirit ...");

            StsWell[] wells = null;
            StsOSWell[] ospWells = (StsOSWell[])osWells;

            ospWells = datastore.getOpenSpiritImport().createStsOSWells(ospWells, progressPanel);
        /*
            StsProject project = model.getProject();
            String binaryDataDir = project.getBinaryFullDirString();

            // If the user has selected both binary files and ASCII files the units may be mixed, so reload all from ASCII
            int nWithBinaries = 0;
            for (int n = 0; n < nSelected; n++)
            {
                // If binaries exist for some and not for others. Load all from ASCII
                if (StsWellKeywordIO.binaryFileExist(binaryDataDir, (String) selectedWellnames[n], StsLogVector.WELL_DEV_PREFIX,
                    StsLogVector.X))
                {
                    nWithBinaries++;
                }
            }
            if (nWithBinaries != nSelected)
            {
                progressPanel.appendLine("All binary wells will be re-loaded from ASCII to apply unit selections\n");
                StsMessageFiles.infoMessage("All binary wells will be re-loaded from ASCII to apply unit selections");
                reloadAscii = true;
            }
            else if (nWithBinaries == nSelected)
            {
                progressPanel.appendLine("Ignoring unit selection since all wells are in binary format. Override with ReLoad from ASCII option");
                StsMessageFiles.infoMessage("Ignoring unit selection since all wells are in binary format. Override with ReLoad from ASCII option"); ;
            }

            StsWell[] wells = StsWellImport.createWells(progressPanel, selectedWellnames, vUnits, hUnits, topHoleLocations, wellTypes, wellFactory);
            success = (wells != null);

            // add lineVertex points to project to generate boundingBoxes

            int nLoaded = StsWellImport.nLoaded;
            for (int n = 0; n < nLoaded; n++)
            {
                wells[n].addToProject(); // this adds lineVertices to projectBoundingBoxes
            }

            project.adjustBoundingBoxes(true, false); // extend displayBoundingBox as needed and set cursor3d box accordingly
            project.checkAddUnrotatedClass(StsWell.class);
            project.rangeChanged();

            // create splined points to define lines
            // note that we don't computePoints until wells are added to project above and origin has been set
            for (int n = 0; n < nLoaded; n++)
            {
                wells[n].computeXYZPoints(); // generate splined points between vertices
            }
            progressPanel.appendLine("Loaded " + nLoaded + " of " + nSelected + " well deviation files ...");

            for (int n = 0; n < nLoaded; n++)
            {
                StsWellImport.loadWellMarkers(wells[n], progressPanel);
                StsWellImport.loadPerforationMarkers(wells[n]);
                StsWellImport.loadFMIMarkers(wells[n]);
                StsWellImport.loadEquipmentMarkers(wells[n]);
            }
            StsSeismicVelocityModel velocityModel = project.getSeismicVelocityModel();
            if (velocityModel != null)
            {
                for (int n = 0; n < nLoaded; n++)
                {
					wells[n].adjustWellTimes(velocityModel);
                    wells[n].convertMarkersToTime(velocityModel);
                }
            }
            //        model.win3d.resetCursorPanel();
        */
            return wells;

        }
        catch (Exception e)
        {
            StsException.outputException("StsWellWizard.createWells() failed.", e, StsException.WARNING);
            return null;
        }
    }


    /** build the well entity with its well line displayed in the 3d window */
    static private StsWell buildWell(String wellname, StsLogVector[] xyzmtVectors, StsLogCurve tdCurve, StsProgressPanel progressPanel,
                                     StsAbstractWellFactory wellFactory)
    {
        StsWell well = wellFactory.createWellInstance(wellname);
        try
        {
            // build the well & line vertices

            if (!well.constructWellDevCurves(xyzmtVectors, StsParameters.nullValue, tdCurve))
            {
                if (progressPanel != null)
                {
                    progressPanel.appendLine("Insufficient data to build well vertices for " + wellname);
                }
                StsMessageFiles.errorMessage("Insufficient data to build well vertices for " + wellname);
                return null;
            }
            well.addModelSurfaceMarkers();
            return well;
        }
        catch (Exception e)
        {
            StsException.outputException("StsWellWizard.buildWell() failed.", e, StsException.WARNING);
            progressPanel.setDescriptionAndLevel("StsWellWizard.buildWell() failed.\n", StsProgressBar.ERROR);
            return null;
        }
    }

    static public void loadWellMarkers(StsWell well, StsProgressPanel progressPanel)
    {
        StsWellImport.WellFilenameSet wellFilenameSet = getWellFilenameSet(well.getName());
        String[] markerFilenames = wellFilenameSet.refFilenames;
        if (markerFilenames == null)
        {       	
            return;
        }
        int nFilenames = markerFilenames.length;
        for (int n = 0; n < nFilenames; n++)
        {
            if(markerFilenames[n] == null) continue;
            StsWellKeywordIO.constructWellMarkers(model, well, currentDirectory, markerFilenames[n], model.getProject().getDepthScalar(vUnits), progressPanel);
        }      
    }

    static public void loadPerforationMarkers(StsWell well, StsProgressPanel progressPanel)
    {
        StsWellImport.WellFilenameSet wellFilenameSet = getWellFilenameSet(well.getName());
        String[] markerFilenames = wellFilenameSet.perfFilenames;
        if (markerFilenames == null)
        {               	        	
            return;
        }
        int nFilenames = markerFilenames.length;
        for (int n = 0; n < nFilenames; n++)
        {
            StsWellKeywordIO.constructPerforationMarkers(model, well, currentDirectory, markerFilenames[n], model.getProject().getDepthScalar(vUnits), progressPanel);
        }       
    }

    static public void loadEquipmentMarkers(StsWell well, StsProgressPanel progressPanel)
    {
        StsWellImport.WellFilenameSet wellFilenameSet = getWellFilenameSet(well.getName());
        String[] markerFilenames = wellFilenameSet.equipFilenames;
        if (markerFilenames == null)
        {        	
            return;
        }
        int nFilenames = markerFilenames.length;
        for (int n = 0; n < nFilenames; n++)
        {
            StsWellKeywordIO.constructEquipmentMarkers(model, well, currentDirectory, markerFilenames[n], model.getProject().getDepthScalar(vUnits), progressPanel);
        }               
    }

    static public void loadFMIMarkers(StsWell well, StsProgressPanel progressPanel)
    {
        StsWellImport.WellFilenameSet wellFilenameSet = getWellFilenameSet(well.getName());
        String[] markerFilenames = wellFilenameSet.fmiFilenames;
        if (markerFilenames == null)
        {              	
            return;
        }
        int nFilenames = markerFilenames.length;
        for (int n = 0; n < nFilenames; n++)
        {
            StsWellKeywordIO.constructFmiMarkers(model, well, currentDirectory, markerFilenames[n], model.getProject().getDepthScalar(vUnits), progressPanel);
        }        
    }
}
