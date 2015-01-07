//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Import;

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.util.*;

public class StsMultiWellImport
{
    static StsModel model;
    static String currentDirectory = null;
    static String binaryDataDir;
    static float logNull;
    public static int nLoaded = 0;
    static boolean success = false;
    static String fullBinaryDirName;
    static TreeMap wellFilenameSets;
    static StsAbstractFileSet[] fileSets = new StsAbstractFileSet[6];
    static StsAbstractFileSet fileSet;

    private transient double progressValue = 0.0;
    private transient String progressDescription = "";

    static String[] asciiFiles = null;
    static private String[] selectedWells = null;
    static private StsFile[] selectedFiles = null;
    static private float datumShift = 0.0f;
    static private StsPoint[] topHoleLocations = null;
    static private byte[] wellTypes = null;
    static private boolean reloadAscii = false;

    static final String WELL_DEV_PREFIX = StsLogVector.WELL_DEV_PREFIX;
    static final String WELL_LOG_PREFIX = StsLogVector.WELL_LOG_PREFIX;
    static final String WELL_TD_PREFIX = StsLogVector.WELL_TD_PREFIX;
    static final String WELL_REF_PREFIX = "well-ref";
    static public final int DEV = 0;
    static public final int LOG = 1;
    static public final int TD = 2;
    static public final int REF = 3;
    static public final int HDR = 4;
    static String[] filePrefixs = new String[] {WELL_DEV_PREFIX, WELL_LOG_PREFIX, WELL_TD_PREFIX, WELL_REF_PREFIX };
    static boolean[] applyShift = new boolean[] {false, false, false, false };

    static byte vUnits = StsParameters.DIST_NONE;
    static byte hUnits = StsParameters.DIST_NONE;
    static byte binaryHorzUnits = StsParameters.DIST_NONE;
    static byte binaryVertUnits = StsParameters.DIST_NONE;

    static public final boolean debug = false;
    static StsTimer timer = new StsTimer();

    static public void initialize(StsModel model_)
    {
        model = model_;
        StsProject project = model.getProject();
        currentDirectory = project.getRootDirString();
        binaryDataDir = project.getBinaryFullDirString();
        logNull = project.getLogNull();
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

    static public boolean skipFileHeader(BufferedReader bufRdr, int numRows)
    {
        try
        {
            bufRdr.skip(0l);
            for(int i=0; i<numRows; i++)
                bufRdr.readLine();
            return true;
        }
        catch(Exception ex)
        {
            StsException.outputException("StsMultiWellImport.skipFileHeader Error skipping header rows.", ex, StsException.WARNING);
            return false;
        }
    }

    static public String[] getHeaderRow(String filename, int numRowsToSkip)
    {
        String[] headerTokens = new String[0];
        String line;
        try
        {
            if(numRowsToSkip == 0) return null;

            BufferedReader bufRdr = new BufferedReader(new FileReader(filename));
            bufRdr.skip(0l);
            for(int i=0; i<numRowsToSkip-1; i++)
                bufRdr.readLine();

            line = bufRdr.readLine();
            if(line == null) return null;
            line = StsStringUtils.detabString(line).trim();
			StringTokenizer stok = new StringTokenizer(line,", ");
            while(stok.hasMoreTokens())
            {
                String token = StsStringUtils.deQuoteString(stok.nextToken());
                headerTokens = (String[])StsMath.arrayAddElement(headerTokens, token);
            }
            return headerTokens;
        }
        catch(Exception ex)
        {
            StsException.outputException("StsMultiWellImport.getHeaderRow Error skipping header rows.", ex, StsException.WARNING);
            return null;
        }
    }

    static public String[] getWellNamesFromHeaderFile(String filename, int numRowsToSkip, String tokens)
    {
        String[] names = new String[0];
        String line = null;
        try
        {
            BufferedReader bufRdr = new BufferedReader(new FileReader(filename));
            skipFileHeader(bufRdr, numRowsToSkip);
            line = bufRdr.readLine().trim();
            while(line != null)
            {
                line = StsStringUtils.detabString(line).trim();
			    StringTokenizer stok = new StringTokenizer(line,tokens);
                String name = stok.nextToken();
                name = StsStringUtils.deQuoteString(name);
                names = (String[])StsMath.arrayAddElement(names, name);
                line = bufRdr.readLine();
            }
        }
        catch(Exception ex)
        {
            StsException.outputException("StsMultiWellImport.getWellNamesFromHeaderFile reading " + filename, ex, StsException.WARNING);
            return null;
        }
        return names;
    }

    static public boolean isWellInFile(String filename, String wellname)
    {
        String[] names = new String[0];
        String line;
        try
        {
            BufferedReader bufRdr = new BufferedReader(new FileReader(filename));
            line = bufRdr.readLine().trim();
            while(line != null)
            {
                line = StsStringUtils.detabString(line).trim();
//                if(line.contains(wellname))
                    return true;
//                line = bufRdr.readLine();
            }
        }
        catch(Exception ex)
        {
            StsException.outputException("StsMultiWellImport.isWellInFile Error reading " + filename, ex, StsException.WARNING);
            return false;
        }
        return false;
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

    static public String[] getSelectedWellFilenames(String[] selectedWellnames, boolean includeBinaries)
    {
        String wellname = null;
        String[] filenames = null;
        WellFilenameSet wellFilenameSet = null;

        for (int n = 0; n < selectedWellnames.length; n++)
        {
            wellname = selectedWellnames[n];
            wellFilenameSet = (WellFilenameSet) wellFilenameSets.get(wellname);
            filenames = (String[]) StsMath.arrayAddElement(filenames, wellFilenameSet.devFilenames[0]);
        }
        return filenames;
    }

    static public WellFilenameSet getWellFilenameSet(String wellname)
    {
        return (WellFilenameSet) wellFilenameSets.get(wellname);
    }

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

    static public void setApplyDatumShift(boolean[] apply)
    {
        applyShift = apply;
    }

    static public boolean selectedOfType(int type, String[] selectedWellnames)
    {
        String wellname = null;
        WellFilenameSet wellFilenameSet = null;
        for (int n = 0; n < selectedWellnames.length; n++)
        {
            wellname = selectedWellnames[n];
            wellFilenameSet = (WellFilenameSet) wellFilenameSets.get(wellname);
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

    static public boolean addWellFilenameSets(String headerFile, int hdrSkip, int[] hdrCols,
                                              String surveyFile, int svySkip, int[] surveyCols,
                                              String topsFile, int topSkip, int[] topCols,
                                              String tdFile, int tdSkip, int[] tdCols,
                                              String[] logFiles)
    {
        WellFilenameSet wellFilenameSet;
        try
        {
            if((surveyFile == null) || (headerFile == null))
            {
                return false;
            }

            if (wellFilenameSets == null)
            {
                wellFilenameSets = new TreeMap();
            }
            // Parse the headerFile to determine the names of all the wells.
            StsWellKeywordIO.loadMultipleHeaderFile(headerFile, hdrSkip, hdrCols);
            String[] names = StsWellKeywordIO.getMultipleFileHeader().wellNames;

            // Create a wellfilenameSet for each well name detected
            for(int i=0; i<names.length; i++)
            {
                wellFilenameSet = new WellFilenameSet(names[i]);
                // ToDo: Verify that each file has a reference to this UWI before adding it to the wellnameSet

                wellFilenameSet.addDeviationFilename(new String[] {surveyFile});
                if(logFiles != null) wellFilenameSet.addLogFilenames(logFiles);
                if(topsFile != null) wellFilenameSet.addReferenceFilenames(new String[] {topsFile});
                if(tdFile != null) wellFilenameSet.addTdFilenames(new String[] {tdFile});
                wellFilenameSets.put(names[i], wellFilenameSet);
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsWellImport.addWellFilenameSets() failed for file", e, StsException.WARNING);
            return false;
        }
        return true;
    }

    static public class WellFilenameSet
    {
        public String wellname;
        public int nLogFilenames = 0;
        public String[] devFilenames;
        public String[] logFilenames = null;
        public int nRefFilenames = 0;
        public String[] refFilenames = null;
        public int increment = 10;

        public boolean hasDepth = false;

        public WellFilenameSet(String wellname)
        {
            this.wellname = wellname;
            devFilenames = new String[4];
            logFilenames = new String[4];
			refFilenames = new String[4];
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

        public void addDeviationFilename(String[] filename)
        {
            devFilenames = filename;
        }
        public void addLogFilenames(String[] filenames)
        {
            logFilenames = (String[]) StsMath.arrayAddArray(logFilenames, filenames);
            nLogFilenames = logFilenames.length;
        }
        public void addTdFilenames(String[] filenames)
        {
            logFilenames = (String[]) StsMath.arrayAddArray(logFilenames, filenames);
            nLogFilenames = logFilenames.length;
        }
        public void addReferenceFilenames(String[] filenames)
        {
            refFilenames = filenames;
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

        public boolean isOK()
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

    static public StsWell[] createWells(StsProgressPanel progressPanel, String headerFile, int headerSkip, int[] headerOrder,
                               String surveyFile, int surveySkip, int[] surveyOrder,
                               String[] logFiles,
                               String tdFile, int tdSkip, int[] tdOrder,
                               String topsFile, int topsSkip, int[] topsOrder,
                               StsAbstractWellFactory wellFactory, float shift)
    {
        StsWell well = null;
        String[] filenames = null;

        // Create wells filename set
        if(!addWellFilenameSets(headerFile, headerSkip, headerOrder,
                surveyFile, surveySkip, surveyOrder,
                topsFile, topsSkip, topsOrder,
                tdFile, tdSkip, tdOrder, logFiles))
        {
            progressPanel.setDescriptionAndLevel("Failed to create filename set: Missing survey or header file", StsProgressBar.WARNING);
            return null;
        }

        // Process the wells in the well sets
        try
        {
            StsWellKeywordIO.loadMultipleHeaderFile(headerFile, headerSkip, headerOrder);
            StsWellKeywordIO.MultipleFileHeader header = StsWellKeywordIO.getMultipleFileHeader();
            String[] names = header.wellNames;

            int nWells = names.length;

            progressPanel.initialize(nWells);
            StsWell[] wells = new StsWell[nWells];
            for (int n = 0; n < nWells; n++)
            {
                String wellname = names[n];
                double x = header.xOrigins[n];
                double y = header.yOrigins[n];
                double z = header.elevations[n];
                StsPoint topHolePoint = new StsPoint(x,y,z);
                byte utType = StsUTKeywordIO.MINCURVE;

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
                    well = createWellFromAsciiFiles(wellFilenameSet, currentDirectory, binaryDataDir, topHolePoint, utType, logNull, progressPanel, wellFactory);
                    if (well != null)
                        wells[nLoaded++] = well;
                }
                else
                {
                    if (progressPanel != null)
                    {
                        progressPanel.appendLine("Well: " + wellname + " already loaded...\n");
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
        StsProgressPanel progressPanel, StsAbstractWellFactory wellFactory)
    {
    	return null;
    }

    static public void loadWellMarkers(StsWell well)
    {
    	return;
    }
}
