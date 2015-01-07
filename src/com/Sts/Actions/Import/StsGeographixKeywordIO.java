
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2000
//Company:      4D Systems LLC
//Description:  Web-enabled interpretation system

/**
    This class is used to read LAS ASCII well files that use keywords
    that are followed by values
*/

package com.Sts.Actions.Import;

import com.Sts.DBTypes.*;
import com.Sts.IO.StsFilenameEndingFilter;
import com.Sts.IO.StsFilenameGroupFilter;
import com.Sts.MVC.StsMessageFiles;
import com.Sts.MVC.StsModel;
import com.Sts.Types.StsPoint;
import com.Sts.UI.Progress.StsProgressBar;
import com.Sts.UI.Progress.StsProgressPanel;
import com.Sts.Utilities.*;
import com.Sts.Utilities.Coordinates.StsCoordinateConversion;
import com.Sts.Utilities.DateTime.CalendarParser;

import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Calendar;
import java.util.Date;
import java.util.StringTokenizer;

public class StsGeographixKeywordIO extends StsKeywordIO
{
    static public final byte TD = 0;
    static public final byte LOGS = 1;

    static private StsLogVector xVector, yVector, zVector, mVector, tVector;
	static private long asciiFileDate = 0;
    static private float datumShift = 0.0f;
    static private boolean deleteBinaries = false;
    static private float vScalar = 1.0f;
    static private byte verticalUnits = StsParameters.DIST_NONE;
    static private byte horizontalUnits = StsParameters.DIST_NONE;
    static WellHeader wellHeader = null;
    static FileHeader fileHeader = null;

    static String[] suffixs = new String[] {"wls"};

    static public String getWellName()
    {
        if(wellHeader != null)
            return wellHeader.wellName;
        else
            return name; 
    }
    static public String getCurveName() { return subname; }
    static public WellHeader getWellHeader() { return wellHeader; }
    static public FileHeader getFileHeader() { return fileHeader; }

	static public void initialize(StsModel model)
	{
        ;
    }

    static public StsLogVector[] readDeviation(BufferedReader bufRdr, String dataDir, String binaryDataDir, String wellname, String filename, float logCurveNull, boolean deleteBins, byte vUnits, byte hUnits, float shift)
    {
        name = wellname;
        deleteBinaries = deleteBins;
        verticalUnits = vUnits;
        horizontalUnits = hUnits;
        datumShift = shift;
        return readDeviation(bufRdr, dataDir, binaryDataDir, filename, logCurveNull);
    }

    /** read a deviation and return values as log vectors */
    static public StsLogVector[] readDeviation(BufferedReader bufRdr, String dataDir, String binaryDataDir, String fname, float logCurveNull)
    {
        StsLogVector[] deviationVectors = null;
        filename = fname;
        try
        {
            // open the file and position to the correct location
			String[] curveNames = new String[] {"MDepth","Inclination","Azimuth","X","Y","Depth"};

			// if ascii file is newer than any of the binaries, delete all binaries as they
			// are potentially out of date
			if(!binaryFileDatesOK(dataDir, binaryDataDir, filename, curveNames, StsLogVector.WELL_DEV_PREFIX))
				deleteBinaryFiles(binaryDataDir, curveNames, StsLogVector.WELL_DEV_PREFIX);

            deviationVectors = constructLogVectors(curveNames, StsLogVector.WELL_DEV_PREFIX);
            if(!StsLogVector.deviationVectorsOK(deviationVectors))
            {
                StsException.systemError("StsGeographixKeywordIO.readDeviation() failed. Didn't find  X or Y or DEPTH|MDEPTH vectors in file: " + filename);
                return null;
            }

            xVector = StsLogVector.getVectorOfType(deviationVectors, StsLogVector.X);
            yVector = StsLogVector.getVectorOfType(deviationVectors, StsLogVector.Y);
            zVector = StsLogVector.getVectorOfType(deviationVectors, StsLogVector.DEPTH);
            zVector.setUnits(verticalUnits);
            mVector = StsLogVector.getVectorOfType(deviationVectors, StsLogVector.MDEPTH);
            if(mVector != null)
                mVector.setUnits(verticalUnits);
            tVector = StsLogVector.getVectorOfType(deviationVectors, StsLogVector.TIME);

            // classInitialize origin to well-dev file header origin
            xVector.setOrigin(wellHeader.xOrigin);
            yVector.setOrigin(wellHeader.yOrigin);
            xVector.setUnits(horizontalUnits);
            yVector.setUnits(horizontalUnits);

            // read the curve values and set the log vectors; may be from ascii or binary files
            if (!readCurveValues(bufRdr, wellHeader, binaryDataDir, curveNames, deviationVectors, logCurveNull, true))
            {
                StsMessageFiles.logMessage("Unable to read deviation vector values for"
                        + " well " + name + " from file: " + filename);
                return null;
            }

            if(zVector != null) zVector.checkMonotonic();
            if(mVector != null) mVector.checkMonotonic();

            return deviationVectors;
        }
    	catch(Exception e)
        {
            StsException.outputException("StsGeographixKeywordIO.readDeviation() failed.", e, StsException.WARNING);
            StsMessageFiles.logMessage("Unable to read deviation for well" + name + ".");
            return null;
    	}
    }

    static private boolean checkWriteBinaryFiles(StsLogVector[] vectors, String binaryDataDir)
    {
        int nVectors = vectors.length;
        boolean writeOK = true;
        for(int n = 0; n < nVectors; n++)
            if(!vectors[n].checkWriteBinaryFile(binaryDataDir)) writeOK = false;

        return writeOK;
    }

	/** check that binary file dates are all newer than ascii file date */
	 static private boolean binaryFileDatesOK(String dir, String binaryDir, String asciiFilename, String[] curveNames, String filePrefix)
	 {
		 File file;

		 try
		 {
			 file = new File(dir + File.separator + asciiFilename);
			 long asciiFileDate = file.lastModified();
			 boolean binaryDatesOK = true;
			 int nNames = curveNames.length;
			 for(int n = 0; n < nNames; n++)
			 {
				 String curveName = curveNames[n];
				 if(curveName.equals("") || curveName.equalsIgnoreCase("ignore")) continue;
				 String binaryFilename = filePrefix + ".bin." + name + "." + curveName + "." + version;
                 binaryFilename = StsStringUtils.cleanString(binaryFilename);
				 file = new File(binaryDir + File.separator + binaryFilename);
				 long binaryFileDate = file.lastModified();
				 if(binaryFileDate < asciiFileDate) binaryDatesOK = false;
			 }
			 return binaryDatesOK;
		 }
		 catch (Exception e)
		 {
			 StsException.outputException("StsGeographixKeywordIO.logVectorDatesOK() failed.",
				 e, StsException.WARNING);
			 return false;
		 }
	 }
	 /** remove all binary log curves associated with this ascii file */
	  static private void deleteBinaryFiles(String binaryDir, String[] curveNames, String filePrefix)
	  {
		  File file;

		  try
		  {
			  int nNames = curveNames.length;
			  for(int n = 0; n < nNames; n++)
			  {
				  String curveName = curveNames[n];
				  if(curveName.equals("") || curveName.equalsIgnoreCase("ignore")) continue;
				  String binaryFilename = filePrefix + ".bin." + name + "." + curveName + "." + version;
				  file = new File(binaryDir + File.separator + binaryFilename);
				  file.delete();
			  }
		  }
		  catch (Exception e)
		  {
			  StsException.outputException("StsGeographixKeywordIO.logVectorDatesOK() failed.",
				  e, StsException.WARNING);
		  }
	  }

      /* read the curve names */
      static private StsLogVector[] addLogVector(StsLogVector[] vectors, String curveName, String filePrefix)
      {
          try
          {
              String binaryFilename = filePrefix + ".bin." + name + "." + curveName + "." + version;
              StsLogVector curveVector = new StsLogVector(filename, binaryFilename, curveName, version, vectors.length);
              vectors = (StsLogVector[])StsMath.arrayAddElement(vectors, curveVector);
              return vectors;
          }
          catch (Exception e)
          {
              StsException.outputException("StsGeographixKeywordIO.addLogVector() failed.",e, StsException.WARNING);
              return null;
          }
      }

    /* read the curve names */
    static private StsLogVector[] constructLogVectors(String[] curveNames, String filePrefix)
    {
        try
        {
            int nNames = curveNames.length;

            StsLogVector[] curveVectors = new StsLogVector[0];
            for(int n = 0; n < nNames; n++)
            {
                String curveName = curveNames[n];
				if(curveName.equals("") || curveName.equalsIgnoreCase("ignore")) continue;
                String binaryFilename = filePrefix + ".bin." + name + "." + curveName + "." + version;
                StsLogVector curveVector = new StsLogVector(filename, binaryFilename, curveName, version, n);
				curveVectors = (StsLogVector[])StsMath.arrayAddElement(curveVectors, curveVector);
            }
            return curveVectors;
        }
        catch (Exception e)
        {
            StsException.outputException("StsWellKeywordIO.constructLogVectors() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

     static public StsLogCurve[] readLogCurves(StsWell well, String wellname, byte logType, String dataDir, String binaryDataDir,
                                          float logCurveNull, boolean deleteBins, byte vUnits, float shift)
     {
         datumShift = shift;
         deleteBinaries = deleteBins;
         verticalUnits = vUnits;

         String[] filenames = null;
         StsLogCurve[] logCurves = new StsLogCurve[0];

         try
         {
             File directory = new File(dataDir);
             StsFilenameEndingFilter filter = new StsFilenameEndingFilter(suffixs);
             filenames = directory.list(filter);
             if(filenames == null) return logCurves;
             int nFilenames = filenames.length;
             if(nFilenames == 0) return logCurves;

             StsLogVector[] logVectors;
             for(int n = 0; n < nFilenames; n++)
             {
                 logVectors = readLogVectors(wellname, dataDir, binaryDataDir, filenames[n], logType, logCurveNull);
                 if(logVectors == null)
                 {
                     StsMessageFiles.logMessage("StsGeographixKeywordIO.ReadLogVectors:Failed to construct log curves for file: " + filenames[n]);
                 }
                 else
                 {
                     StsLogCurve[] newLogCurves = StsLogCurve.constructLogCurves(well, logVectors, logCurveNull, version);
                     if (newLogCurves.length == 0) continue;
                     logCurves = (StsLogCurve[])StsMath.arrayAddArray(logCurves, newLogCurves);
                 }
                 if(group == StsLogVector.WELL_LOG_PREFIX) well.checkAddMDepthToDev(logVectors);
             }
             return logCurves;
         }
         catch (Exception e)
         {
             StsMessageFiles.logMessage("StsGeographixKeywordIO.ReadLogVectors:Log curve read failed for well " + wellname);
             return logCurves;
         }
    }

     static public StsLogVector[] readLogVectors(String wellname, String dataDir, String binaryDataDir, String filename,
                                                 byte logType, float logCurveNull, boolean deleteBins, byte vUnits)
	 {
         deleteBinaries = deleteBins;
         verticalUnits = vUnits;
         return readLogVectors(wellname, dataDir, binaryDataDir, filename, logType, logCurveNull);

     }

	 static public StsLogVector[] readLogVectors(String wellname, String dataDir, String binaryDataDir, String filename, byte logType, float logCurveNull)
	 {
		 BufferedReader bufRdr = null;
		 StsLogVector[] logVectors;
		 StsLogCurve[] logCurves = null;

		 try
		 {
			 bufRdr = new BufferedReader(new FileReader(dataDir + File.separator + filename));
			 if(wellHeader == null)
                 return null;

			 String[] curveNames = null; //wellHeader.curveNames;

			 logVectors = constructLogVectors(curveNames, group);
			 if (logVectors == null) return null;

			 if(!readCurveValues(bufRdr, wellHeader, binaryDataDir, curveNames, logVectors, logCurveNull, false))
                 return null;
			 return logVectors;
		 }
		 catch (Exception e)
		 {
			 StsMessageFiles.logMessage("StsGeographixKeywordIO.ReadLogVectors: Log curve read failed for well " + wellname);
			 return null;
			}
		 finally
		 {
			 closeBufRdr(bufRdr);
		 }
	 }

    static public void readFileHeader(String dataDir, String filename)
    {
        fileHeader = new FileHeader(dataDir, filename);
    }
    /* read the A and B records */
    static public BufferedReader readWellHeader(String dataDir, String filename, int wellIdx)
    {
        int nWells = 0;
        BufferedReader bufRdr = null;
        boolean foundA = false, foundB = false;

		wellHeader = new WellHeader();
        try
        {
            bufRdr = new BufferedReader(new FileReader(dataDir + File.separator + filename));
			while(nWells<=wellIdx)
			{
				String line = bufRdr.readLine();
                if(line == null)
                {
                    wellHeader = null;
                    return null;
                }
                if(line.startsWith("A"))
                {
                    foundA = true;
                    wellHeader.readRecordA(line);
                }
                else if(line.startsWith("B"))
                {
                    foundB = true;
                    wellHeader.readRecordB(line);
                }
                if(foundA && foundB)
                {
                    nWells++;
                    foundA = false;
                    foundB = false;
                }
            }
            return bufRdr;
        }
        catch(Exception e)
        {
            StsException.outputException("StsGeographixKeywordIO.readWellHeader() failed.",e, StsException.WARNING);
            wellHeader = null;
			return null;
        }
    }
	// used for wls files to capture file header info
	static public class FileHeader
	{
		int numWellsInFile = 0;
        double[][] wellheadXYZs = null;
        String[] wellNames = null;
        boolean[] applyKbCorrection = null;
        String directory = null;
        String fileName = null;

        public FileHeader(String dataDir, String filename)
        {
            boolean foundA = false, foundB = false;
            directory = dataDir;
            fileName = filename;
            numWellsInFile = 0;
            wellNames = null;
            wellheadXYZs = new double[1000][3];
            applyKbCorrection = new boolean[1000];
            BufferedReader bufRdr = null;
            try
            {
                bufRdr = new BufferedReader(new FileReader(dataDir + File.separator + filename));
                String line = bufRdr.readLine();
                WellHeader wellHeader = new WellHeader();
			    while(line != null)
			    {
                    if(line.startsWith("A"))
                    {
                        foundA = true;
                        wellHeader.readRecordA(line);
                    }
                    else if(line.startsWith("B"))
                    {
                        foundB = true;
                        wellHeader.readRecordB(line);
                    }
                    if(foundA && foundB)
                    {
                        wellheadXYZs[numWellsInFile][0] = wellHeader.xOrigin;
                        wellheadXYZs[numWellsInFile][1] = wellHeader.yOrigin;
                        wellheadXYZs[numWellsInFile][2] = wellHeader.datumElevation;
                        applyKbCorrection[numWellsInFile] = true;
                        wellNames = (String[]) StsMath.arrayAddElement(wellNames, wellHeader.wellName);
                        numWellsInFile++;
                        foundA = false;
                        foundB = false;
                    }
                    line = bufRdr.readLine();
                }
                return;
            }
            catch(Exception e)
            {
                StsException.outputException("StsGeographixKeywordIO.FileHeader() failed.",e, StsException.WARNING);
                fileHeader = null;
			    return;
            }
        }

        public String[] getWellnames() { return wellNames; }
        public double[][] getWellheads() { return wellheadXYZs; }
        public int getNumWellsInFile() { return numWellsInFile; }
        public boolean[] getApplyKbCorrection() { return applyKbCorrection; }
    }

	// used for wls files to capture well header info
	static public class WellHeader
	{
		String formatVersion = "wls";
        boolean lineWrap = false;
        float nullValue = StsParameters.nullValue;
        private int dateOrder = CalendarParser.MM_DD_YY;
        byte depthType = StsLogVector.MDEPTH;

        // A Records
        String uwiNumber = "00000000000000000000";
        String field = "unknown";
        String operator = "operator";
        String area = null;
        String wellName = "unknown";
        String wellNumber = "unknown";
        String state = "unknown";
        String county = "unknown";
        Date spudDate = new Date(0);
        String datumDescription = "unknown";
        float datumElevation = 0.0f;
        boolean applyDatum = true;
        float totalDepth = 0.0f;
        String formationAtTd = "unknown";
        String classification = "unknown";
        String status = "unknown";
        Date completionDate = new Date(0);
        String altIdNumber = "00000000000000000000";
        String wellIdNumber = "00000000000000000000";
        String permitNumber = "0000000";
        Date permitDate = new Date(0);
        static final byte ENGLISH = 0;
        static final byte METRIC = 1;
        byte units=ENGLISH;

        // B Records
        static final byte TOWNSHIP = 0;
        static final byte INTERNATIONAL = 1;
        static final byte NORTHEAST = 2;
        int locationFormat = NORTHEAST;
        float latitude = 0.0f;
        float longitude = 0.0f;
        String locationDescription = "unknown";
        float refLatitude = 0.0f;
        float refLongitude = 0.0f;
        double xOrigin = 0.0f;  // UTM
        double yOrigin = 0.0f;  // UTM

        public WellHeader()
        {
            ;
        }

        public boolean readRecordA(String line)
        {
            try
            {
            if(line.length() < 244) return false;

            uwiNumber = line.substring(1,21).trim();
            field = line.substring(21,36).trim();
            operator = line.substring(36,61).trim();
            area = line.substring(61,66).trim();
            wellName = line.substring(66,91).trim();
            wellName = uwiNumber + "-" + wellName;
            if(wellName.length() == 0)
                wellName = uwiNumber;
            wellNumber = line.substring(91,97).trim();
            wellName = wellName + " " + wellNumber;
            state = line.substring(97,101).trim();
            county = line.substring(101,105).trim();
            String dateStg = line.substring(105,115).trim();
            spudDate = parseDate(dateStg);
            datumDescription = line.substring(115,117).trim();
            try { datumElevation = Float.valueOf(line.substring(117,125).trim()).floatValue(); }
            catch(Exception ex) { datumElevation = 0.0f; }
            try { totalDepth = Float.valueOf(line.substring(125,130).trim()).floatValue(); }
            catch(Exception ex) { totalDepth = 0.0f; }
            formationAtTd = line.substring(130,150).trim();
            classification = line.substring(150,170).trim();
            status = line.substring(170,176).trim();
            dateStg = line.substring(176,186).trim();
            completionDate = parseDate(dateStg);
            altIdNumber = line.substring(186,206).trim();
            wellIdNumber = line.substring(206,226).trim();
            permitNumber = line.substring(226,233).trim();
            dateStg = line.substring(233,243).trim();
            permitDate = parseDate(dateStg);
            String unitChar = line.substring(243,line.length()-1).trim();
            if(unitChar.equalsIgnoreCase("M"))
                units = METRIC;
            }
            catch(Exception ex)
            {
                StsException.outputException("StsGeographixKeywordIO.WellHeader.ReadRecordA:Failed to read A record", ex, StsException.WARNING);
                return false;
            }
            return true;
        }

        public Date parseDate(String dateStg)
        {
            try
            {
                Calendar cal = CalendarParser.parse(dateStg, dateOrder , true);
                return cal.getTime();
            }
            catch(Exception ex)
            {
                return new Date(0);
            }
        }

        public boolean readRecordB(String line)
        {
            try
            {
                String format = line.substring(1,3).trim();
                if(format.equalsIgnoreCase("0"))
                {
                    locationFormat = TOWNSHIP;
                    if(line.length() < 72) return false;

                    try { latitude = Float.valueOf(line.substring(2,17).trim()).floatValue(); }
                    catch(Exception ex) { latitude = 0.0f; }
                    try { longitude = Float.valueOf(line.substring(17,32).trim()).floatValue(); }
                    catch(Exception ex) { longitude = 0.0f; }
                    locationDescription = line.substring(32,72).trim();
                    refLatitude = 0.0f;
                    refLongitude = 0.0f;
                }
                else if(format.equalsIgnoreCase("1"))
                {
                    locationFormat = INTERNATIONAL;
                    if(line.length() < 72) return false;

                    try { latitude = Float.valueOf(line.substring(2,17).trim()).floatValue(); }
                    catch(Exception ex) { latitude = 0.0f; }
                    try { longitude = Float.valueOf(line.substring(17,32).trim()).floatValue(); }
                    catch(Exception ex) { longitude = 0.0f; }
                    locationDescription = line.substring(32,72).trim();
                    refLatitude = 0.0f;
                    refLongitude = 0.0f;
                }
                else if(format.equalsIgnoreCase("2"))
                {
                    locationFormat = NORTHEAST;
                    if(line.length() < 88) return false;

                    try { latitude = Float.valueOf(line.substring(2,17).trim()).floatValue(); }
                    catch(Exception ex) { latitude = 0.0f; }
                    try { longitude = Float.valueOf(line.substring(17,32).trim()).floatValue(); }
                    catch(Exception ex) { longitude = 0.0f; };
                    locationDescription = line.substring(32,72).trim();
                    try { refLatitude = Float.valueOf(line.substring(72,88).trim()).floatValue(); }
                    catch(Exception ex) { refLatitude = 0.0f; }
                    try { refLongitude = Float.valueOf(line.substring(88,line.length()-1).trim()).floatValue(); }
                    catch(Exception ex) { refLongitude = 0.0f; }
                }
                else
                {
                    StsMessageFiles.errorMessage("StsGeographixKeywordIO.WellHeader.ReadRecordB:Failed to read B record - Invalid format type (0,1, or 2)");
                    return false;
                }
                yOrigin = latitude;
                xOrigin = longitude;
            }
            catch(Exception ex)
            {
                StsException.outputException("StsGeographixKeywordIO.WellHeader.ReadRecordB:Failed to read B record", ex, StsException.WARNING);
                return false;
            }
            return true;
        }
        /*
        // Read latitude and longitude and convert to UTM
        public boolean readRecordB(String line)
        {
            try
            {
                String format = line.substring(1,3).trim();
                if(format.equalsIgnoreCase("0"))
                {
                    locationFormat = TOWNSHIP;
                    try { latitude = Float.valueOf(line.substring(2,17).trim()).floatValue(); }
                    catch(Exception ex) { latitude = 0.0f; }
                    try { longitude = Float.valueOf(line.substring(17,32).trim()).floatValue(); }
                    catch(Exception ex) { longitude = 0.0f; }
                    locationDescription = line.substring(32,72).trim();
                    refLatitude = 0.0f;
                    refLongitude = 0.0f;
                }
                else if(format.equalsIgnoreCase("1"))
                {
                    locationFormat = INTERNATIONAL;
                    try { latitude = Float.valueOf(line.substring(2,17).trim()).floatValue(); }
                    catch(Exception ex) { latitude = 0.0f; }
                    try { longitude = Float.valueOf(line.substring(17,32).trim()).floatValue(); }
                    catch(Exception ex) { longitude = 0.0f; }
                    locationDescription = line.substring(32,72).trim();
                    refLatitude = 0.0f;
                    refLongitude = 0.0f;
                }
                else if(format.equalsIgnoreCase("2"))
                {
                    locationFormat = NORTHEAST;
                    try { latitude = Float.valueOf(line.substring(2,17).trim()).floatValue(); }
                    catch(Exception ex) { latitude = 0.0f; }
                    try { longitude = Float.valueOf(line.substring(17,32).trim()).floatValue(); }
                    catch(Exception ex) { longitude = 0.0f; };
                    locationDescription = line.substring(32,72).trim();
                    try { refLatitude = Float.valueOf(line.substring(72,88).trim()).floatValue(); }
                    catch(Exception ex) { refLatitude = 0.0f; }
                    try { refLongitude = Float.valueOf(line.substring(88,line.length()-1).trim()).floatValue(); }
                    catch(Exception ex) { refLongitude = 0.0f; }
                }
                else
                {
                    StsMessageFiles.errorMessage("StsGeographixKeywordIO.WellHeader.ReadRecordB:Failed to read B record - Invalid format type (0,1, or 2)");
                    return false;
                }
                // Convert to UTM
                double[] xy = latLong2Utm();
                xOrigin = xy[0];  // UTM
                yOrigin = xy[1];  // UTM
            }
            catch(Exception ex)
            {
                StsException.outputException("StsGeographixKeywordIO.WellHeader.ReadRecordB:Failed to read A record", ex, StsException.WARNING);
                return false;
            }
            return true;
        }
       */
        public  double[] latLong2Utm()
        {
            StsCoordinateConversion c = new StsCoordinateConversion();
            return c.latLon2UTM(latitude, longitude);
        }
    }

    /* read the curve values and set the log vectors */
    static private boolean readCurveValues(BufferedReader bufRdr, WellHeader fileHdr, String binaryDataDir, String[] curveNames,
										   StsLogVector[] curveVectors, float logCurveNull, boolean loadValues)
    {
        try
        {
            int nCurveVectors = curveVectors.length;

			boolean binaryFilesOK = true;
		    for(int n = 0; n < nCurveVectors; n++)
				if(!curveVectors[n].hasBinaryFile(binaryDataDir)|| deleteBinaries)
                    binaryFilesOK = false;

			if(binaryFilesOK)
			{
				for (int n = 0; n < nCurveVectors; n++)
					curveVectors[n].readBinaryFile(binaryDataDir, loadValues);
				return true;
			}
            else
                deleteBinaryFiles(binaryDataDir, curveNames, StsLogVector.WELL_LOG_PREFIX);

			if(!readAsciiCurveValues(bufRdr, fileHdr, curveNames, curveVectors, logCurveNull, loadValues))
                return false;
            return checkWriteBinaryFiles(curveVectors, binaryDataDir);
        }
        catch (Exception e)
        {
            StsException.outputException("StsGeographix.KeywordIO.readCurveValues() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    static private boolean readAsciiCurveValues(BufferedReader bufRdr, WellHeader fileHdr, String[] curveNames,
										   StsLogVector[] curveVectors, float logCurveNull, boolean loadValues)
	{
		try
		{
            // read lines until we hit the eof
            int nCurveVectors = curveVectors.length;
            int depthVectorIdx = -1, owtVectorIdx = -1;

            StsFloatVector[] curveValues = new StsFloatVector[nCurveVectors];

            for (int n = 0; n < nCurveVectors; n++)
            {
                curveValues[n] = new StsFloatVector(5000, 2000);
                curveVectors[n].setValues(curveValues[n]);
				curveVectors[n].setNullValue(logCurveNull);
                curveVectors[n].setUnits(verticalUnits);
            }

			int nCurveNames = curveNames.length;
			int[] colAssignments = new int[nCurveNames];
			for(int n = 0; n < nCurveNames; n++)
				colAssignments[n] = -1;

			for(int n = 0; n < nCurveNames; n++)
			{
				int col = curveVectors[n].getAsciiFileColumn();
				if(col < 0 || col >= nCurveNames) continue;
                if(curveVectors[n].getType() == StsLogVector.DEPTH)
                    depthVectorIdx = n;

                if(curveVectors[n].getType() == StsLogVector.OWT)
                    owtVectorIdx = n;

                if(curveVectors[n].getType() == StsLogVector.TWT)
                {
                    String binaryFilename = StsLogVector.WELL_TD_PREFIX + ".bin." + name + ".TIME." + version;
                    curveVectors[n].setNameAndType(StsLogVector.types[StsLogVector.TIME], binaryFilename);
                }
				colAssignments[col] = n;
			}

            // Locate the N Record
            String line;
            line = getRecordType(bufRdr, "N");

            // Loop through the O Records
			int nLines = 0;
            while ((line = bufRdr.readLine()) != null)
            {
                line = line.trim();
                if(line.equals("")) continue;  // blank line
                if(!line.startsWith("O"))
                    break;
                
                if(nLines == 0)
                {
                    curveVectors[0].append(0.0f);
                    curveVectors[1].append(0.0f);
                    curveVectors[2].append(0.0f);
                    curveVectors[3].append(0.0f);
                    curveVectors[4].append(0.0f);
                    if(fileHdr.applyDatum)
                        curveVectors[5].append(datumShift - fileHdr.datumElevation);
                    else
                        curveVectors[5].append(datumShift);

                }

                int start;
                int end;
                for(int j=0; j<nCurveNames; j++)
                {
                    start = (j+1) * 10;
                    if(j == (nCurveNames-1))
                        end = line.length()-1;
                    else
                        end = (j+2) * 10;
                    double value = 0.0f;
                    try { value = Double.valueOf(line.substring(start,end)); }
                    catch (Exception ex) { }

				    int nVector = colAssignments[j];
				    if(nVector == -1) continue;

                    if(nVector == depthVectorIdx)
                    {
                        if(fileHdr.applyDatum)
                            curveVectors[nVector].append(value + datumShift - fileHdr.datumElevation);
                        else
                            curveVectors[nVector].append(value + datumShift);

                    }
                    else if(nVector == owtVectorIdx)
                        curveVectors[nVector].append(value * 2.0f);
                    else
                        curveVectors[nVector].append(value);
                }
                if(++nLines%100 == 0)
                    StsMessageFiles.logMessage("File: " + filename + ": " + nLines + " values read.");
            }
            // Vertical Well or Only Top and Bottom supplied.
            if(nLines < 2)
            {
                // If N Record has coordinates for bottom hole.

                // else vertical well
                curveVectors[0].append(0.0f); curveVectors[0].append(fileHdr.totalDepth  + datumShift);
                curveVectors[1].append(0.0f); curveVectors[1].append(0.0f);
                curveVectors[2].append(0.0f); curveVectors[2].append(0.0f);
                curveVectors[3].append(0.0f); curveVectors[3].append(0.0f);
                curveVectors[4].append(0.0f); curveVectors[4].append(0.0f);
                if(fileHdr.applyDatum)
                {
                    curveVectors[5].append(datumShift - fileHdr.datumElevation);
                    curveVectors[5].append(fileHdr.totalDepth  + datumShift - fileHdr.datumElevation);
                }
                else
                {
                    curveVectors[5].append(datumShift);
                    curveVectors[5].append(fileHdr.totalDepth  + datumShift);
                }
            }

            if(owtVectorIdx != -1)
            {
                String binaryFilename = StsLogVector.WELL_TD_PREFIX + ".bin." + name + ".TIME." + version;
                curveVectors[owtVectorIdx].setNameAndType(StsLogVector.types[StsLogVector.TIME], binaryFilename);
            }

            // finish up values: trim, min, max, null
            for(int n = 0; n < nCurveVectors; n++)
			{
				curveValues[n].trimToSize();
                curveValues[n].setMinMax();
                if(curveVectors[n].isNull())  curveVectors[n].setValues((StsFloatVector)null);
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsGeographixKehywordIO.readCurveValues() failed.", e, StsException.WARNING);
            return false;
        }
    }

    static private String getRecordType(BufferedReader bufRdr, String recordType)
    {
        try
        {
            String line = bufRdr.readLine();
            while(!line.startsWith(recordType) && !line.startsWith("A"))
            {
                line = bufRdr.readLine();
                if(line == null)
                    return null;
            }
            return line;
        }
        catch (Exception e)
        {
            StsException.outputException("StsGeographixKeywordIO.getRecordType(" + recordType + ") failed.", e, StsException.WARNING);
            return null;
        }

    }
    static public int readWellMarkers(StsWell well, int wellIdx, String dataDir, String fname, StsProgressPanel progressPanel)
    {
        BufferedReader bufRdr = null;
        String line = null;
        StsPoint location = null;
        int nLoaded = 0;

        try
        {
            bufRdr = readWellHeader(dataDir, fname, wellIdx);
            line = getRecordType(bufRdr, "C");
            if(line == null)
            {
                progressPanel.appendLine("      No geologic marker records found.");
                progressPanel.appendLine("      Loaded " + nLoaded + " geologic markers for well " + well.getName() + "...");
                return 0;
            }
            while(line.startsWith("C"))
            {
                // Read the marker name
                String markerName = line.substring(5,25);
                if (well.getMarker(markerName) != null)
                {
                    line = bufRdr.readLine();
                    continue;  // already have it
                }
                // Read the value
                float value = 0.0f;
                try
                {
                    value = Float.valueOf(line.substring(37,45)).floatValue();
                }
                catch(Exception ex)
                {
                    value = 0.0f;
                    line = bufRdr.readLine();
                    continue;
                }
                // Marker outside well extent
                if((value > well.getMaxMDepth()) || (value < well.getMinMDepth()))
                {
                    value = 0.0f;
                    line = bufRdr.readLine();
                    continue;
                }
                // compute and validate the marker location
                location = StsWellKeywordIO.getWellMarkerLocation(well, value, wellHeader.depthType);
                if (location != null)
                {
                    StsWellMarker.constructor(markerName, well, StsMarker.GENERAL, location);
                    nLoaded++;
                }
                else
                {
                    String depthTypeString = StsLogVector.depth_types[StsLogVector.MDEPTH];
                    progressPanel.appendLine("Failed to load marker " + markerName + " at " + depthTypeString + ": " + value);
                    progressPanel.setDescriptionAndLevel("Errors", StsProgressBar.WARNING);
                }
                line = bufRdr.readLine();
            }
        }
        catch(Exception ex)
        {
            StsException.outputException("StsGeographixKeywordIO:readWellMarkers: failed to read well marker records.", ex, StsException.WARNING);
            progressPanel.appendLine("      Failed to load any geologic markers for well " + well.getName() + "...");
            return 0;
        }
        finally
        {
            closeBufRdr(bufRdr);
        }
        progressPanel.appendLine("      Loaded " + nLoaded + " geologic markers for well " + well.getName() + "...");
        return nLoaded;
    }
    static public int readPerforationMarkers(StsWell well, int wellIdx, String dataDir, String fname, StsProgressPanel progressPanel)
    {
        BufferedReader bufRdr = null;
        String line = null;
        StsPoint location = null;
        int nLoaded = 0;

        try
        {
            bufRdr = readWellHeader(dataDir, fname, wellIdx);
            line = getRecordType(bufRdr, "K");
            if(line == null)
            {
                progressPanel.appendLine("      No perforation records found.");
                progressPanel.appendLine("      Loaded " + nLoaded + " perforation markers for well " + well.getName() + "...");
                return 0;
            }
            while(line.startsWith("K"))
            {
                String markerName = line.substring(5,25);
                if(well.getMarker(markerName) != null)
                {
                    line = bufRdr.readLine();
                    continue;  // already have it
                }
                float top = 0.0f, btm = 0.0f;
                try
                {
                    top = Float.valueOf(line.substring(25,33)).floatValue() * vScalar;
                    btm = Float.valueOf(line.substring(33,41)).floatValue() * vScalar;
                }
                catch(Exception ex)
                {
                    line = bufRdr.readLine();                    
                    continue;
                }
                float len = btm - top;

                location = well.getPointAtMDepth(top+((top-btm)/2.0f), false);

				if(location != null)
				{
					StsPerforationMarker.constructor(markerName, well, StsMarker.GENERAL, location, len);
					nLoaded++;
				}
                line = bufRdr.readLine();
            }
        }
        catch(Exception ex)
        {
            StsException.outputException("StsGeographixKeywordIO:readPerforationMarkers: failed to read well marker records.", ex, StsException.WARNING);
            progressPanel.appendLine("      Failed to load any perforation markers for well " + well.getName() + "...");
            return 0;
        }
        finally
        {
            closeBufRdr(bufRdr);
        }
        progressPanel.appendLine("      Loaded " + nLoaded + " perforation markers for well " + well.getName() + "...");
        return nLoaded;
    }

    static public int readFMIMarkers(StsWell well, int wellIdx, String dataDir, String fname, StsProgressPanel progressPanel)
    {
        BufferedReader bufRdr = null;
        String line = null;
        StsPoint location = null;
        int nLoaded = 0;
        float[] atts = new float[1];
        StsPoint locWithAtts = null;

        try
        {
            bufRdr = readWellHeader(dataDir, fname, wellIdx);
            line = getRecordType(bufRdr, "L");
            if(line == null)
            {
                progressPanel.appendLine("      No Fault records found.");
                progressPanel.appendLine("      Loaded " + nLoaded + " fault markers for well " + well.getName() + "...");
                return 0;
            }
            while(line.startsWith("L"))
            {
                String markerName = line.substring(5,35).trim();
                if(well.getMarker(markerName) != null)
                {
                    line = bufRdr.readLine();
                    continue;  // already have it
                }
                float md = 0.0f, dip = 0.0f, azimuth = 0.0f;
                try
                {
                    md = Float.valueOf(line.substring(35,43).trim()).floatValue();
                    dip = Float.valueOf(line.substring(51,59).trim()).floatValue();
                    azimuth = Float.valueOf(line.substring(59,67).trim()).floatValue();
                }
                catch(Exception ex)
                {
                    line = bufRdr.readLine();
                    continue;
                }
                atts[0] = Float.valueOf(line.substring(43,51)).floatValue();

                location = well.getPointAtMDepth(md, false);
				if(location == null)
                {
                    line = bufRdr.readLine();
                    continue;
                }

                locWithAtts = new StsPoint(location.getLength() + atts.length, location);
                locWithAtts.v[location.getLength()] = atts[0];
				if(locWithAtts != null)
				{
					StsFMIMarker.constructor(well, locWithAtts, dip, azimuth);
					nLoaded++;
				}
                line = bufRdr.readLine();
            }
        }
        catch(Exception ex)
        {
            StsException.outputException("StsGeographixKeywordIO:readFMIMarkers: failed to read well marker records.", ex, StsException.WARNING);
            progressPanel.appendLine("      Failed to load any FMI markers for well " + well.getName() + "...");
            return 0;
        }
        finally
        {
            closeBufRdr(bufRdr);
        }
        progressPanel.appendLine("      Loaded " + nLoaded + " FMI markers for well " + well.getName() + "...");
        return nLoaded;
    }

    static public int readEquipmentMarkers(StsWell well, int wellIdx, String dataDir, String fname, StsProgressPanel progressPanel)
    {
        return 0;
    }
}