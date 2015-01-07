
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2000
//Company:      4D Systems LLC
//Description:  Web-enabled interpretation system

/**
    This class is used to read ASCII well files that use keywords
    that are followed by values
*/

package com.Sts.Actions.Import;

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;
import com.Sts.Utilities.DateTime.CalendarParser;

import java.awt.*;
import java.io.*;
import java.util.*;

public class StsWellKeywordIO extends StsKeywordIO
{
    /** Multiple Well Files Definistions */
    static public final byte UWI = 0;
    static public final byte NAME = 1;
    static public final byte MX = 2;
    static public final byte MY = 3;
    static public final byte SYMBOL = 4;
    static public final byte KB = 5;
    static public final byte GRD = 6;
    static public final byte TVD = 7;
    static public final byte DATUM = 8;
    static public final byte[] headerCols = new byte[] {UWI, NAME, MX, MY, SYMBOL, KB, GRD, TVD, DATUM};
    static public final byte MD = 0;
    static public final byte AZIM = 1;
    static public final byte DIP = 2;
    static public final byte[] surveyCols = new byte[] {MD,AZIM,DIP};
    static public final byte TD_TIME = 1;
    static public final byte[] tdCols = new byte[] {MD,TD_TIME};
    static public final byte TOP = 1;
    static public final byte[] topsCols = new byte[] {MD,TOP};

    /** keyword constants */
    static public final String WELL_NAME = "WELLNAME";
    static public final String ORIGIN = "ORIGIN";
    static public final String XY = "XY";
    static public final String YX = "YX";
    static public final String CURVE = "CURVE";
    static public final String VALUE = "VALUE";
	static public final String NULL_VALUE = "NULL";
    static public final String WELL_REF = "WELLREF";
    static public final String X = "X";
    static public final String Y = "Y";
    static public final String DEPTH = "DEPTH";
    static public final String MDEPTH = "MDEPTH";
    static public final String SSDEPTH = "SSDEPTH";
    static public final String TIME = "TIME";
    static public final String DUMMY_NAME = "IGNORE";
    static public final String[] Z_KEYWORDS = { DEPTH, MDEPTH, SSDEPTH }; //TIME };
    static public final String[] X_KEYWORDS = { "X", "EASTING", "DX", "XEVT"};
    static public final String[] Y_KEYWORDS = { "Y", "NORTHING", "DY", "YEVT" };
    static public final String[] TIME_KEYWORDS = { "CTIME", "HOUR", "HMS", "TIMESTAMP", "TIMEVT" };
    static public final String[] DATE_KEYWORDS = { "DATE", "DAY", "DATEVT" };

    static private StsLogVector xVector, yVector, zVector, mVector, tVector;
    static private StsTimeVector timeVector = null;
	static private long asciiFileDate = 0;
    static private float datumShift = 0.0f;
    static private int dateOrder = CalendarParser.DD_MM_YY;
    static private boolean deleteBinaries = false;
    static public float vScalar = 1.0f;
    static public float hScalar = 1.0f;
    static public byte verticalUnits = StsParameters.DIST_NONE;
    static public byte horizontalUnits = StsParameters.DIST_NONE;
    static public FileHeader fileHeader = null;
    static public MultipleFileHeader multipleFileHeader = null;

    static public String getWellName() { return name; }
    static public String getCurveName() { return subname; }
    static public FileHeader getFileHeader() { return fileHeader; }
    static public MultipleFileHeader getMultipleFileHeader() { return multipleFileHeader; }

	static public void initialize(StsModel model)
	{
//        StsKeywordIO.classInitialize(model);
//		logCurveNull = project.getLogNull();
    }

    static public boolean loadMultipleHeaderFile(String filename, int skippedRows, int[] colOrder)
    {
        BufferedReader bufRdr = null;
        try
        {
            bufRdr = new BufferedReader(new FileReader(filename));
            multipleFileHeader = StsWellKeywordIO.readMultipleFileHeader(bufRdr, skippedRows, colOrder);
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsWellKeywordIO.loadMultipleHeaderFile() failed.", e, StsException.WARNING);
            StsMessageFiles.logMessage("Unable to load multiple well header file" + filename + ".");
            return false;
        }
    }
    /** Verify that the deviation binaries exist and are up-to-date for this ASCII file */
    static public boolean verifyDeviationBinaries(String dataDir, String binaryDataDir, String filename)
    {
        BufferedReader bufRdr = null;
        StsLogVector[] vectors = null;
        try
        {
            // parse the filename to get the well name and file type & class
            if(!setParseFilename(filename, "well-dev", "txt"))
                return false;

            // open the file
            bufRdr = new BufferedReader(new FileReader(dataDir + File.separator + filename));

            // read the well name (not used)
            fileHeader = readFileHeader(bufRdr);
            String[] curveNames = fileHeader.curveNames;

            // if ascii file is newer than any of the binaries, delete all binaries as they
            // are potentially out of date
            return binaryFileDatesOK(dataDir, binaryDataDir, filename, curveNames, StsLogVector.WELL_DEV_PREFIX);
        }
        catch(Exception e)
        {
            StsException.outputException("StsWellKeywordIO.readDeviation() failed.",
                e, StsException.WARNING);
            StsMessageFiles.logMessage("Unable to read deviation for well" + name + ".");
            return false;
        }

    }
    static public StsLogVector[] readDeviation(String dataDir, String binaryDataDir, String filename,
                             float logCurveNull, boolean deleteBins, byte vUnits, byte hUnits, String group, float shift)
    {
        deleteBinaries = deleteBins;
        verticalUnits = vUnits;
        horizontalUnits = hUnits;
        datumShift = shift;
        return readDeviation(dataDir, binaryDataDir, filename, logCurveNull, group);
    }

    /** read a deviation and return values as log vectors */
    static public StsLogVector[] readDeviation(String dataDir, String binaryDataDir, String filename, float logCurveNull, String group)
    {
        BufferedReader bufRdr = null;

        StsLogVector[] deviationVectors = null;
        StsTimeVector clockTimeVector = null;
        try
        {
            // parse the filename to get the well name and file type & class
            if(!setParseFilename(filename, group, "txt")) return null;

            // open the file
            bufRdr = new BufferedReader(new FileReader(dataDir + File.separator + filename));

            // read the well name (not used)
            fileHeader = readFileHeader(bufRdr);
			String[] curveNames = fileHeader.curveNames;

			// if ascii file is newer than any of the binaries, delete all binaries as they
			// are potentially out of date
			if(!binaryFileDatesOK(dataDir, binaryDataDir, filename, curveNames, StsLogVector.WELL_DEV_PREFIX) || deleteBinaries)
				deleteBinaryFiles(binaryDataDir, curveNames, StsLogVector.WELL_DEV_PREFIX);

            deviationVectors = constructLogVectors(curveNames, StsLogVector.WELL_DEV_PREFIX);

            if(!StsLogVector.deviationVectorsOK(deviationVectors))
            {
                StsException.systemError("StsWellKeywordIO.readDeviation() failed." +
                    " Didn't find  X or Y or DEPTH|MDEPTH vectors in file: " + filename);
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
            if(fileHeader.liveWell)
            {
                timeVector = StsWellKeywordIO.constructTimeVector("LiveWell", fileHeader.timeIdx, fileHeader.dateIdx);
                timeVector.setValues(new StsLongVector(new long[0], 0, 1));
                timeVector.setMinMaxAndNulls(StsParameters.largeLong);
            }
            // classInitialize origin to well-dev file header origin
            xVector.setOrigin(fileHeader.xOrigin);
            xVector.setUnits(horizontalUnits);
            yVector.setOrigin(fileHeader.yOrigin);
            yVector.setUnits(horizontalUnits);

            // read the curve values and set the log vectors; may be from ascii or binary files
            if (!readCurveValues(bufRdr, binaryDataDir, curveNames, deviationVectors, logCurveNull, true))
            {
                StsMessageFiles.logMessage("Unable to read deviation vector values for"
                        + " well " + name + " from file: " + filename);
                return null;
            }

//            xVector.checkAdjustOrigin(modelXOrigin);
//            yVector.checkAdjustOrigin(modelYOrigin);

            if(zVector != null) zVector.checkMonotonic();
            if(mVector != null) mVector.checkMonotonic();
            if(timeVector != null) timeVector.checkMonotonic();
//            checkWriteBinaryFiles(deviationVectors, binaryDataDir);

            return deviationVectors;
        }
    	catch(Exception e)
        {
            StsException.outputException("StsWellKeywordIO.readDeviation() failed.",
                e, StsException.WARNING);
            if(xVector == null)
            	StsMessageFiles.logMessage("Unable to read deviation for well" + name + ". Verify that X column exists in file.");
            else if(yVector == null)
            	StsMessageFiles.logMessage("Unable to read deviation for well" + name + ". Verify that Y column exists in file.");
            else if(zVector == null)
            	StsMessageFiles.logMessage("Unable to read deviation for well" + name + ". Verify that DEPTH column exists in file.");
            else if(mVector == null)
            	StsMessageFiles.logMessage("Unable to read deviation for well" + name + ". Verify that MDEPTH column exists in file.");
            else
            	StsMessageFiles.logMessage("Unable to read deviation for well" + name + ".");
            return null;
    	}
        finally
        {
            closeBufRdr(bufRdr);
        }
    }

    static public boolean checkWriteBinaryFiles(StsLogVector[] vectors, String binaryDataDir)
    {
        int nVectors = vectors.length;
        boolean writeOK = true;
        for(int n = 0; n < nVectors; n++)
            if(!vectors[n].checkWriteBinaryFile(binaryDataDir, false)) writeOK = false;
        if(fileHeader.liveWell)
            if(!timeVector.checkWriteBinaryFile(binaryDataDir, false)) writeOK = false;
        return writeOK;
    }

    static public StsTimeVector getTimeVector() { return timeVector; }
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
				 file = new File(binaryDir + File.separator + binaryFilename);
				 long binaryFileDate = file.lastModified();
				 if(binaryFileDate < asciiFileDate) binaryDatesOK = false;
			 }
			 return binaryDatesOK;
		 }
		 catch (Exception e)
		 {
			 StsException.outputException("StsWellKeywordIO.logVectorDatesOK() failed.",
				 e, StsException.WARNING);
			 return false;
		 }
	 }

     /** check that binary files exist */
      static public boolean binaryFileExist(String binaryDir, String name, String filePrefix, byte type)
      {
          File file;

          try
          {
              String binaryFilename = filePrefix + ".bin." + name + "." + StsLogVector.getStringFromType(type) + "." + version;
              file = new File(binaryDir + File.separator + binaryFilename);
              return file.exists();
          }
          catch (Exception e)
          {
              StsException.outputException("StsWellKeywordIO.logVectorDatesOK() failed.",
                  e, StsException.WARNING);
              return false;
          }
      }

	  static private void deleteBinaryFiles(String binaryDir, String[] curveNames, String filePrefix)
      {
          	  deleteBinaryFiles(name, binaryDir, curveNames, filePrefix);
      }

	 /** remove all binary log curves associated with this ascii file */
	  static public void deleteBinaryFiles(String name, String binaryDir, String[] curveNames, String filePrefix)
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
              if(timeVector != null)
              {
                  String binaryFilename = filePrefix + ".bin." + name + ".CTIME." + version;
				  file = new File(binaryDir + File.separator + binaryFilename);
				  file.delete();
              }
		  }
		  catch (Exception e)
		  {
			  StsException.outputException("StsWellKeywordIO.logVectorDatesOK() failed.",
				  e, StsException.WARNING);
		  }
	  }

      /* read the curve names */
      static public StsTimeVector constructTimeVector(String filePrefix, int timeIdx, int dateIdx)
      {
          try
          {
              String binaryFilename = filePrefix + ".bin." + name + ".CTIME." + version;
              StsTimeVector timeVector = new StsTimeVector(filename, binaryFilename, "CTIME", version, timeIdx, dateIdx);
              return timeVector;
          }
          catch (Exception e)
          {
              StsException.outputException("StsWellKeywordIO.constructTimeVector() failed.",
                                           e, StsException.WARNING);
              return null;
          }
    }

    /* read the curve names */
    static public StsLogVector[] constructLogVectors(String[] curveNames, String filePrefix)
    {
        return constructLogVectors(curveNames, name, filePrefix);
    }
    /* read the curve names */
    static public StsLogVector[] constructLogVectors(String[] curveNames, String wellname, String filePrefix)
    {
        try
        {
            int nNames = curveNames.length;

            StsLogVector[] curveVectors = new StsLogVector[0];
            for(int n = 0; n < nNames; n++)
            {
                String curveName = curveNames[n];
				if(curveName.equals("") || curveName.equalsIgnoreCase("ignore")) continue;
                String binaryFilename = filePrefix + ".bin." + wellname + "." + curveName + "." + version;
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
    static public StsLogCurve[] readLogCurves(StsWell well, String wellname, String dataDir, String binaryDataDir, String group,
                                              float logCurveNull, boolean deleteBins, byte vUnits, float shift, int order)
    {
        datumShift = shift;
        dateOrder = order;
        return readLogCurves(well, wellname, dataDir, binaryDataDir, group, logCurveNull, deleteBins, vUnits);
    }
    static public StsLogCurve[] readLogCurves(StsWell well, String wellname, String dataDir, String binaryDataDir, String group,
                                              float logCurveNull, boolean deleteBins, byte vUnits)
    {
        deleteBinaries = deleteBins;
        verticalUnits = vUnits;
        return readLogCurves(well, wellname, dataDir, binaryDataDir, group, logCurveNull);
    }

    static public StsLogCurve[] readLogCurves(StsWell well, String wellname, String dataDir, String binaryDataDir, String group, float logCurveNull)
    {
        String[] filenames = null;
		StsLogCurve[] logCurves = new StsLogCurve[0];

        try
        {
            File directory = new File(dataDir);
            StsFilenameGroupFilter filter = new StsFilenameGroupFilter(group, "txt", wellname);
            filenames = directory.list(filter);
            if(filenames == null) return logCurves;
            int nFilenames = filenames.length;
            if(nFilenames == 0) return logCurves;

			StsLogVector[] logVectors;
           for(int n = 0; n < nFilenames; n++)
			{
                logVectors = readLogVectors(wellname, dataDir, binaryDataDir, filenames[n], group, logCurveNull);
				if(logVectors == null)
			    {
				    StsMessageFiles.logMessage("Failed to construct log curves for file: " + filenames[n]);
			    }
			    else
				{
					StsLogCurve[] newLogCurves = StsLogCurve.constructLogCurves(well, logVectors, logCurveNull, 0);
					if (newLogCurves.length == 0) continue;
				    logCurves = (StsLogCurve[])StsMath.arrayAddArray(logCurves, newLogCurves);
				}
                if(group == StsLogVector.WELL_LOG_PREFIX) well.checkAddMDepthToDev(logVectors);
 			}
			return logCurves;
        }
    	catch (Exception e)
        {
            StsMessageFiles.logMessage("Log curve read failed for well " + wellname);
			return logCurves;
       	}
    }
    static public String[] getCurveNames()
    {
        if(fileHeader == null) return null;
        else return fileHeader.curveNames;
    }
    static public double getXOrigin()
    {
        if(fileHeader == null) return 0.0f;
        else return fileHeader.xOrigin;
    }
    static public double getYOrigin()
    {
        if(fileHeader == null) return 0.0f;
        else return fileHeader.yOrigin;
    }
    static public void readHeader(String filename)
    {
        BufferedReader bufRdr = null;
        try
        {
            bufRdr = new BufferedReader(new FileReader(filename));
			fileHeader = readFileHeader(bufRdr);
            return;
        }
        catch(Exception ex)
        {
            StsMessageFiles.logMessage("Unable to read header for: " + filename);
            return;
        }
    }

    static public StsLogVector[] readLogVectors(String wellname, String dataDir, String binaryDataDir, String filename, String group, float logCurveNull)
    {
        BufferedReader bufRdr = null;
        StsLogVector[] logVectors;
        StsLogCurve[] logCurves = null;

        try
        {
            if(!setParseFilename(filename, group, "txt")) return null;
            bufRdr = new BufferedReader(new FileReader(dataDir + File.separator + filename));

			fileHeader = readFileHeader(bufRdr);
		    if(fileHeader == null) return null;

			String[] curveNames = fileHeader.curveNames;
			float curveNullValue = fileHeader.nullValue;

            logVectors = constructLogVectors(curveNames, group);
            if (logVectors == null) return null;

			if(!readCurveValues(bufRdr, binaryDataDir, curveNames, logVectors, logCurveNull, false)) return null;
			return logVectors;
//            checkWriteBinaryFiles(logVectors, binaryDataDir);
        }
    	catch (Exception e)
        {
            StsMessageFiles.logMessage("Log curve read failed for well " + wellname);
            return null;
       	}
        finally
        {
            closeBufRdr(bufRdr);
        }
    }
/*
    static private boolean checkAddMDepthToDev(StsWell well, StsLogVector[] curveLogVectors)
    {
        StsLogVector mDepthVector = StsLogVector.getVectorOfType(curveLogVectors, StsLogVector.MDEPTH);
        StsLogVector depthVector = StsLogVector.getVectorOfType(curveLogVectors, StsLogVector.DEPTH);
        if(mDepthVector == null || depthVector == null) return false;
        return well.checkAddMDepthVector(mDepthVector, depthVector);
    }
 */
/*
	public StsLogCurve[] constructLogCurvesCheckVersions(StsLogVector[] logVectors, float curveNullValue)
	{
		StsLogCurve[] logCurves = new StsLogCurve[0];

		try
		{
			if (logVectors == null)return logCurves;
			int nLogVectors = logVectors.length;
			if (nLogVectors == 0) return logCurves;

			// sort the logVectors by version
			ArrayList objects = new ArrayList(nLogVectors);
			for (int n = 0; n < nLogVectors; n++)
				objects.add(logVectors[n]);

			Comparator comparator = new VersionComparator();
			Collections.sort(objects, comparator);

			Iterator iter = objects.iterator();
			int currentVersion = StsParameters.nullInteger;
			int nVersionVectors = 0;
			StsLogCurve[] newLogCurves = new StsLogCurve[0];
			while (iter.hasNext())
			{
				StsLogVector logVector = (StsLogVector) iter.next();
				int version = logVector.getVersion();
				if (version == currentVersion)
					logVectors[nVersionVectors++] = logVector;
				else
				{
					if (nVersionVectors > 0)
					{
						logVectors = (StsLogVector[]) StsMath.trimArray(logVectors, nVersionVectors);
						newLogCurves = constructLogCurves(logVectors, curveNullValue);
						logCurves = (StsLogCurve[])StsMath.arrayAddArray(logCurves, newLogCurves);
						nVersionVectors = 0;
					}
					currentVersion = version;
					logVectors = new StsLogVector[nLogVectors];
					logVectors[nVersionVectors++] = logVector;
				}
			}
			if (nVersionVectors > 0)
			{
				logVectors = (StsLogVector[]) StsMath.trimArray(logVectors, nVersionVectors);
				newLogCurves = constructLogCurves(logVectors, curveNullValue);
				logCurves = (StsLogCurve[])StsMath.arrayAddArray(logCurves, newLogCurves);
			}
			return logCurves;
		}
		catch (Exception e)
		{
			StsException.outputException("StsLogCurve.constructLogCurvesCheckVersions() failed.",
										 e, StsException.WARNING);
			return logCurves;
		}
	}
*/

/*
	static private boolean writeLogVectors(StsLogVector[] logVectors)
	{
		boolean writeOK = true;

		if(logVectors == null) return false;
		int nVectors = logVectors.length;
		for(int n = 0; n < nVectors; n++)
		    if(!logVectors[n].checkWriteBinaryFile(binaryDataDir)) writeOK = false;
	    return writeOK;
	}
*/
    /** read a strat zone file and return values as log vectors */
    static public int constructWellStratZones(StsWell well, String dataDir, String filename)
    {
        if(filename == null) return 0;

        BufferedReader bufRdr = null;
        int nZonesLoaded = 0;

        try
        {
            // parse the filename to get the well name and file type & class
            if(!setParseFilename(filename, "zone", "strat")) return nZonesLoaded;
            // open the file
            bufRdr = new BufferedReader(new FileReader(dataDir + File.separator + filename));

            // read the well name (not used)
            String readWellName = readWellName(bufRdr);

            // read the curve names & verify we have the required ones
            String[] curveNames = readCurveNames(bufRdr);
            if(curveNames == null || curveNames.length != 1)
            {
                StsMessageFiles.logMessage("Unable to find a single depth keyword "
                        + "for well strat zones in file: " + filename);
                return nZonesLoaded;
            }
            byte depthType = StsLogVector.getTypeFromString(curveNames[0]);
            nZonesLoaded = constructWellStratZones(well, bufRdr, depthType);
            StsMessageFiles.logMessage("Constructed " + nZonesLoaded + " zones for well: " + well.getName());
        }
    	catch (Exception e)
        {
            StsMessageFiles.logMessage("Well strat zone read failed for well " + well.getName());
       	}
        finally
        {
            closeBufRdr(bufRdr);
            return nZonesLoaded;
        }
    }

    /* read the zone values */
    static private int constructWellStratZones(StsWell well, BufferedReader bufRdr, byte depthType)
    {
        int nLoaded = 0;

        try
        {
            if(well == null || bufRdr == null) return nLoaded;

            String line = bufRdr.readLine().trim();
            if (line.indexOf(VALUE) < 0)  return nLoaded; // is keyword not there?

            while ((line = bufRdr.readLine()) != null)
            {
                line = line.trim();
                if (line.equals("")) continue;  // blank line

                StringTokenizer stok = new StringTokenizer(line);
                String zoneName = stok.nextToken();
                if (well.getZone(zoneName, StsWellZoneSet.STRAT)!=null) continue;  // already have it

                float topDepth = Float.valueOf(stok.nextToken()).floatValue();
                float botDepth = Float.valueOf(stok.nextToken()).floatValue();
                if(botDepth < topDepth) continue;

                if(depthType == StsLogVector.MDEPTH)
                {
                    topDepth = well.getDepthFromMDepth(topDepth);
                    botDepth = well.getDepthFromMDepth(botDepth);
                }
                StsSurfaceVertex topVertex = new StsSurfaceVertex(well.getLinePointAtZ(topDepth, false), null);
                StsSurfaceVertex botVertex = new StsSurfaceVertex(well.getLinePointAtZ(botDepth, false), null);
                StsWellZone wellZone = new StsWellZone(well, StsWellZoneSet.STRAT, zoneName, topVertex, botVertex);
                well.addZone(wellZone);
                nLoaded++;
            }
            return nLoaded;
        }
        catch (Exception e)
        {
            StsException.outputException("StsWellKeywordIO.consructWellStratZones() failed.",
                e, StsException.WARNING);
            return nLoaded;
        }
        finally
        {
            closeBufRdr(bufRdr);
        }
    }
    /** read a marker file and construct markers */
    static public int constructFmiMarkers(StsModel model, StsWell well, String dataDir, String filename, float vScale, StsProgressPanel progressPanel)
    {
        vScalar = vScale;
        return StsWellKeywordIO.constructFmiMarkers(model, well, dataDir, filename, progressPanel);
    }

    /** read a marker file and construct markers */
    static public int constructFmiMarkers(StsModel model, StsWell well, String dataDir, String filename, StsProgressPanel progressPanel)
    {
        if(filename == null) return 0;

        BufferedReader bufRdr = null;
        int nMarkersLoaded = 0;

        try
        {
            // parse the filename to get the well name and file type & class
            if(!setParseFilename(filename, "well-fmi", "txt")) return nMarkersLoaded;
            // open the file
            bufRdr = new BufferedReader(new FileReader(dataDir + File.separator + filename));

            // read the well name (not used)
            String readWellName = readWellName(bufRdr);

            // read the curve names & verify we have the required ones
            String[] curveNames = readCurveNames(bufRdr);
            if(curveNames == null || curveNames.length < 3)
            {
                StsMessageFiles.logMessage("Unable to read FMI markers, need MDEPTH, DIP and AZIMUTH values in file.");
                return nMarkersLoaded;
            }

            byte depthType = StsLogVector.getTypeFromString(curveNames[0]);
            nMarkersLoaded = constructFmiMarkers(model, well, bufRdr, depthType, curveNames, progressPanel);
            StsMessageFiles.logMessage("Constructed " + nMarkersLoaded + " FMI markers for well: " + well.getName());
            progressPanel.appendLine("   Successfully processed FMI markers for well " + well.getName() + "...");
        }
        catch (Exception e)
        {
            progressPanel.appendLine("   Failed to process perforation markers for well " + well.getName() + "from file " + filename + "...");        	
            StsMessageFiles.logMessage("Well FMI read failed for well " + well.getName());
        }
        finally
        {
            closeBufRdr(bufRdr);
            return nMarkersLoaded;
        }
    }

    /** read a marker file and construct markers */
    static public int constructEquipmentMarkers(StsModel model, StsWell well, String dataDir, String filename, float vScale, StsProgressPanel progressPanel)
    {
        vScalar = vScale;
        return StsWellKeywordIO.constructEquipmentMarkers(model, well, dataDir, filename, progressPanel);
    }

    /** read a marker file and construct markers */
    static public int constructEquipmentMarkers(StsModel model, StsWell well, String dataDir, String filename, StsProgressPanel progressPanel)
    {
        if(filename == null) return 0;

        BufferedReader bufRdr = null;
        int nMarkersLoaded = 0;

        try
        {
            // parse the filename to get the well name and file type & class
            if(!setParseFilename(filename, "well-equipment", "txt")) return nMarkersLoaded;
            // open the file
            bufRdr = new BufferedReader(new FileReader(dataDir + File.separator + filename));
            fileHeader = readFileHeader(bufRdr);
            nMarkersLoaded = constructEquipmentMarkers(model, well, bufRdr, progressPanel);
            StsMessageFiles.logMessage("Constructed " + nMarkersLoaded + " equipment markers for well: " + well.getName());
            progressPanel.appendLine("   Successfully processed equipment markers for well " + well.getName() + "...");
        }
        catch (Exception e)
        {
            StsMessageFiles.logMessage("Failed to process equipment markers for well " + well.getName());
            progressPanel.appendLine("   Failed to process equipment markers for well " + well.getName() + "from file " + filename + "...");
        }
        finally
        {
            closeBufRdr(bufRdr);
            return nMarkersLoaded;
        }
    }

    /* read the zone values */
    static private int constructEquipmentMarkers(StsModel model, StsWell well, BufferedReader bufRdr, StsProgressPanel progressPanel)
    {
		String line = null;
		StsPoint location = null;
        int nLoaded = 0;

        try
        {
            if(well == null || bufRdr == null) return nLoaded;
            while ((line = bufRdr.readLine()) != null)
            {
                line = line.trim();
                if (line.equals("")) continue;  // blank line

                StringTokenizer stok = new StringTokenizer(line);
                String markerName = stok.nextToken();
                if(well.getMarker(markerName) != null) continue;  // already have it

                String subType = "NONE";
                int tIdx = 0;
                float x = -1.0f, y = -1.0f, z = -1.0f, len = 10.0f;
                while(stok.hasMoreTokens())
                {
                    if(tIdx == fileHeader.zIdx)
                    {
                        z = Float.valueOf(stok.nextToken()).floatValue() * vScalar;
                        if(fileHeader.depthType == StsLogVector.MDEPTH)
                            location = well.getPointAtMDepth(z, false);
                        else if(fileHeader.depthType == StsLogVector.DEPTH)  // This is a potential problem on non-vertical wells.
                            location = well.getPointFromDepth(z);
                    }
                    else if(tIdx == fileHeader.xIdx)
                    {
                        x = Float.valueOf(stok.nextToken()).floatValue();
                    }
                    else if (tIdx == fileHeader.yIdx)
                    {
                        y = Float.valueOf(stok.nextToken()).floatValue();
                    }
                    else
                    {
                        subType = stok.nextToken();
                    }
                    tIdx++;
                }
                // If x, y, and depth are supplied.
                if(x != -1.0f && y != -1.0f && fileHeader.depthType == StsLogVector.DEPTH)
                {
                    location = well.getPointFromLocation(x,y,z);
                }
				if(location != null)
				{
					StsEquipmentMarker equipmentMarker = StsEquipmentMarker.constructor(markerName, well, location,
							StsEquipmentMarker.getSubTypeFromString(subType));
					nLoaded++;
				}
            }
            return nLoaded;
        }
        catch (Exception e)
        {
            StsException.outputException("StsWellKeywordIO.consructPerforationMarkers() failed reading line: " + line + location.toString(),
                e, StsException.WARNING);
            return nLoaded;
        }
        finally
        {
            closeBufRdr(bufRdr);
        }
    }

    /** read a marker file and construct markers */
    static public int constructPerforationMarkers(StsModel model, StsWell well, String dataDir, String filename, float vScale, StsProgressPanel progressPanel)
    {
        vScalar = vScale;
        return StsWellKeywordIO.constructPerforationMarkers(model, well, dataDir, filename, progressPanel);
    }

    /** read a marker file and construct markers */
    static public int constructPerforationMarkers(StsModel model, StsWell well, String dataDir, String filename, StsProgressPanel progressPanel)
    {
        if(filename == null) return 0;

        BufferedReader bufRdr = null;
        int nMarkersLoaded = 0;

        try
        {
            // parse the filename to get the well name and file type & class
            if(!setParseFilename(filename, "well-perf", "txt")) return nMarkersLoaded;
            // open the file
            bufRdr = new BufferedReader(new FileReader(dataDir + File.separator + filename));
            fileHeader = readFileHeader(bufRdr);
            nMarkersLoaded = constructPerforationMarkers(model, well, bufRdr, progressPanel);
            StsMessageFiles.logMessage("Constructed " + nMarkersLoaded + " perforation markers for well: " + well.getName());
            progressPanel.appendLine("   Successfully processed perforation markers for well " + well.getName() + "...");
        }
        catch (Exception e)
        {
            progressPanel.appendLine("   Failed to process perforation markers for well " + well.getName() + "from file " + filename + "...");        	
            StsMessageFiles.logMessage("Failed to process perforation markers for well " + well.getName());
        }
        finally
        {
            closeBufRdr(bufRdr);
            return nMarkersLoaded;
        }
    }
    /** read a marker file and construct markers */
    static public int constructWellMarkers(StsModel model, StsWell well, String dataDir, String filename, float vScale, StsProgressPanel progressPanel)
    {
        BufferedReader bufRdr = null;

        vScalar = vScale;
        // parse the filename to get the well name and file type & class
        if(!setParseFilename(filename, "well-ref", "txt")) return 0;
        // open the file
        try
        {
            bufRdr = new BufferedReader(new FileReader(dataDir + File.separator + filename));
            fileHeader = readFileHeader(bufRdr);
            return constructWellMarkers(model, well, bufRdr, progressPanel);
        }
        catch (Exception e)
        {
            StsMessageFiles.logMessage("Well strat zone read failed for well " + well.getName());
            progressPanel.appendLine("   Failed to process geologic markers for well " + well.getName() + "from file " + filename + "...");
            return 0;
        }
        finally
        {
            closeBufRdr(bufRdr);
        }
    }

    static int constructWellMarkers(StsModel model, StsWell well, BufferedReader bufRdr, StsProgressPanel progressPanel)
    {
        String line = null;
        StsPoint location = null;
        int nLoaded = 0;

        try
        {
            if (well == null || bufRdr == null) return nLoaded;

            while ((line = bufRdr.readLine()) != null)
            {
                line = line.trim();
                if (line.equals("")) continue;  // blank line

                StringTokenizer stok = new StringTokenizer(line);
                String markerName = stok.nextToken();
                if (well.getMarker(markerName) != null) continue;  // already have it

                float value = Float.valueOf(stok.nextToken()).floatValue();
                location = getWellMarkerLocation(well, value, fileHeader.depthType);
                if (location != null)
                {
                    StsWellMarker.constructor(markerName, well, StsMarker.GENERAL, location);
                    nLoaded++;
                }
                else
                {
                    String depthTypeString = StsLogVector.depth_types[fileHeader.depthType];
                    progressPanel.appendLine("Failed to load marker " + markerName + " at " + depthTypeString + ": " + value);
                    progressPanel.setDescriptionAndLevel("Errors", StsProgressBar.WARNING);
                }
            }
            return nLoaded;
        }
        catch (Exception e)
        {
            StsException.outputException("StsWellKeywordIO.consructWellStratZones() failed reading line: " + line + location.toString(),
                e, StsException.WARNING);
            return nLoaded;
        }
        finally
        {
            closeBufRdr(bufRdr);
        }
    }
    /* read the zone values */
    static private int constructFmiMarkers(StsModel model, StsWell well, BufferedReader bufRdr, byte depthType, String[] names, StsProgressPanel progressPanel)
    {
		String line = null;
		StsPoint location = null;
        int nLoaded = 0, nAtts = 0;;
        float dip, azimuth, mdepth;
        float[] atts = new float[names.length - 3];
        StsPoint locWithAtts = null;
        try
        {
            if(well == null || bufRdr == null)
            	return nLoaded;

            line = bufRdr.readLine().trim();
            if (line.indexOf(VALUE) < 0)  return nLoaded; // is keyword not there?


            while ((line = bufRdr.readLine()) != null)
            {
                line = line.trim();
                if (line.equals("")) continue;  // blank line

                dip = -1.0f;
                azimuth = -1.0f;
                mdepth = -1.0f;
                nAtts = 0;
                int tokenNum = 0;

                StringTokenizer stok = new StringTokenizer(line);
                while(stok.hasMoreTokens())
                {
                	if(names[tokenNum].equalsIgnoreCase("MDEPTH") || names[tokenNum].equalsIgnoreCase("DEPTH") ||
                			names[tokenNum].equalsIgnoreCase("TVD") || names[tokenNum].equalsIgnoreCase("Z"))
                	{
                		mdepth = Float.valueOf(stok.nextToken()).floatValue();
                		if(depthType == StsLogVector.MDEPTH)
                			location = well.getPointAtMDepth(mdepth * vScalar, false);
                		else if(depthType == StsLogVector.DEPTH)
                			location = well.getPointFromDepth(mdepth * vScalar);
                	}
                	else if(names[tokenNum].equalsIgnoreCase("DIP"))
                	{
                		dip = Float.valueOf(stok.nextToken()).floatValue();
                	}
                	else if(names[tokenNum].equalsIgnoreCase("AZIMUTH"))
                	{
                		azimuth = Float.valueOf(stok.nextToken()).floatValue();
                	}
                	/*
                	else if(names[tokenNum].equalsIgnoreCase("HEIGHT"))
                	{
                		height = Float.valueOf(stok.nextToken()).floatValue();
                	}
                	*/
                	else
                	{
                		atts[nAtts++] = Float.valueOf(stok.nextToken()).floatValue();
                	}
                	tokenNum++;
                }
                if(location == null)
                	continue;
                locWithAtts = new StsPoint(location.getLength() + atts.length, location);
                for(int i=0; i<nAtts; i++)
                	locWithAtts.v[location.getLength() + i] = atts[i];

				if(locWithAtts != null)
				{
					StsFMIMarker fmiMarker = StsFMIMarker.constructor(well, locWithAtts, dip, azimuth);
					nLoaded++;
				}
            }
            return nLoaded;
        }
        catch (Exception e)
        {
            StsException.outputException("StsWellKeywordIO.consructFmiMarkers() failed reading line: " + line + location.toString(),
                e, StsException.WARNING);
            return nLoaded;
        }
        finally
        {
            closeBufRdr(bufRdr);
        }
    }

    /* read the zone values */
    static private int constructPerforationMarkers(StsModel model, StsWell well, BufferedReader bufRdr, StsProgressPanel progressPanel)
    {
		String line = null;
		StsPoint location = null;
        int nLoaded = 0;

        try
        { 
            for(int i=0; i<fileHeader.curveNames.length; i++)
            {
                System.out.println("Attribute= " + fileHeader.curveNames[i]);
            }
            String[] names = fileHeader.curveNames;
            if(well == null || bufRdr == null) return nLoaded;
            while ((line = bufRdr.readLine()) != null)
            {
                line = line.trim();
                if (line.equals("")) continue;  // blank line

                StringTokenizer stok = new StringTokenizer(line);
                String markerName = stok.nextToken();
                if(well.getMarker(markerName) != null) continue;  // already have it

                int tIdx = 0, nShots = 1;
                float x = -1.0f, y = -1.0f, z = -1.0f, len = 10.0f;
                String shotTime = null, shotDate = null;
                while(stok.hasMoreTokens())
                {
                    if(tIdx == fileHeader.zIdx)
                    {
                        z = Float.valueOf(stok.nextToken()).floatValue() * vScalar;
                        if(fileHeader.depthType == StsLogVector.MDEPTH)
                            location = well.getPointAtMDepth(z, false);
                        else if(fileHeader.depthType == StsLogVector.DEPTH)  // This is a potential problem on non-vertical wells.
                            location = well.getPointFromDepth(z);
                    }
                    else if(tIdx == fileHeader.xIdx)
                    {
                        x = Float.valueOf(stok.nextToken()).floatValue();
                    }
                    else if (tIdx == fileHeader.yIdx)
                    {
                        y = Float.valueOf(stok.nextToken()).floatValue();
                    }
                    else if(names[tIdx].equalsIgnoreCase("NSHOTS") || names[tIdx].equalsIgnoreCase("NUMSHOTS"))
                    {
                        nShots = Integer.valueOf(stok.nextToken()).intValue();
                    }
                    else if(names[tIdx].equalsIgnoreCase("LEN") || names[tIdx].equalsIgnoreCase("LENGTH"))
                    {
                        len = Float.valueOf(stok.nextToken()).floatValue();
                    }
                    else if(names[tIdx].equalsIgnoreCase("DATE"))
                    {
                        shotDate = stok.nextToken();
                    }
                    else if(names[tIdx].equalsIgnoreCase("TIME"))
                    {
                        shotTime = stok.nextToken();
                    }
                    tIdx++;
                }
                // If x, y, and depth are supplied.
                if(x != -1.0f && y != -1.0f && fileHeader.depthType == StsLogVector.DEPTH)
                {
                    location = well.getPointFromLocation(x,y,z);
                }
                Long perfTime = 0l;
                if(shotTime != null)
                {
                    if(shotDate != null)
                        shotTime = shotDate.trim() + " " + shotTime.trim();
                    Calendar cal = CalendarParser.parse(shotTime, model.getProject().getDateOrder(), true);
                    perfTime = cal.getTimeInMillis();
                }
				if(location != null)
				{
					StsPerforationMarker perfMarker = StsPerforationMarker.constructor(markerName, well, StsMarker.GENERAL, location, len, nShots, perfTime);
					nLoaded++;
				}
            }
            return nLoaded;
        }
        catch (Exception e)
        {
            StsException.outputException("StsWellKeywordIO.consructPerforationMarkers() failed reading line: " + line + location.toString(),
                e, StsException.WARNING);
            return nLoaded;
        }
        finally
        {
            closeBufRdr(bufRdr);
        }
    }

    static public StsPoint getWellMarkerLocation(StsWell well, float value, byte depthType)
    {
        StsPoint location;
        float z;

        if (depthType == StsLogVector.MDEPTH)
        {
            float m = value*vScalar; // z is actually mdepth
            location = well.getPointAtMDepth(m, true);
            z = location.getZ();
        }
        else // depthType == StsLogVector.DEPTH
        {
            z = value*vScalar;
            location = well.getPointAtZ(z, true);
        }

        StsLogCurve tdCurve = well.getTdCurve();
        if(tdCurve != null)
        {
            StsLogVector timeVector = tdCurve.getValueVector();
            StsLogVector depthVector = tdCurve.getDepthVector();
            float t = StsLogVector.interpolateValue(z, depthVector, timeVector);
            if(t != StsParameters.nullValue)
                location.setT(t);
        }
        return location;
    }

    /* read the well name */
    static public FileHeader readFileHeader(BufferedReader bufRdr)
    {
		FileHeader fileHeader = new FileHeader();
        try
        {
			while(true)
			{
				String line = bufRdr.readLine().trim();
                line = StsStringUtils.detabString(line);
				if(line.endsWith(WELL_NAME))
					fileHeader.wellName =  new String(bufRdr.readLine().trim());
				else if(line.indexOf(ORIGIN) >= 0)  // is origin keyword there?
				{
					boolean xyOrder = true;
					if (line.indexOf(YX) >= 0) xyOrder = false;  // determine x-y order
					line = bufRdr.readLine().trim();  // get the next line

					// tokenize the x-y values and convert to a point object
					StringTokenizer stok = new StringTokenizer(line);
					fileHeader.xOrigin = Double.valueOf(stok.nextToken()).doubleValue();
					fileHeader.yOrigin = Double.valueOf(stok.nextToken()).doubleValue();

					if(!xyOrder)
					{
						double temp = fileHeader.xOrigin;
						fileHeader.xOrigin = fileHeader.yOrigin;
						fileHeader.yOrigin = temp;
					}
				}
				else if(line.endsWith(NULL_VALUE))
				{
					line = bufRdr.readLine().trim();  // get the next line
					StringTokenizer stok = new StringTokenizer(line);
					fileHeader.nullValue = Float.valueOf(stok.nextToken()).floatValue();
				}
				else if(line.endsWith(CURVE))
				{
					String[] curveNames = new String[0];
					line = bufRdr.readLine().trim();
                    int currentIdx = 0;
					while(!lineHasKeyword(line, VALUE))
					{
                        fileHeader.xIdx = checkIndex(line, currentIdx, fileHeader.xIdx, StsWellKeywordIO.X_KEYWORDS);
                        fileHeader.yIdx = checkIndex(line, currentIdx, fileHeader.yIdx, StsWellKeywordIO.Y_KEYWORDS);
                        fileHeader.zIdx = checkIndex(line, currentIdx, fileHeader.zIdx, StsWellKeywordIO.Z_KEYWORDS);
                        fileHeader.timeIdx = checkIndex(line, currentIdx, fileHeader.timeIdx, StsWellKeywordIO.TIME_KEYWORDS);
                        fileHeader.dateIdx = checkIndex(line, currentIdx, fileHeader.dateIdx, StsWellKeywordIO.DATE_KEYWORDS);
						if((currentIdx != fileHeader.timeIdx) && (currentIdx != fileHeader.dateIdx))
                            curveNames = (String[])StsMath.arrayAddElement(curveNames, line);
                        currentIdx++;
						line = bufRdr.readLine().trim();
					}
					fileHeader.curveNames = curveNames;
                    fileHeader.depthType = StsLogVector.getTypeFromString(curveNames[fileHeader.zIdx]);
                    if((fileHeader.timeIdx != -1) && (fileHeader.dateIdx != -1))
					    fileHeader.liveWell = true;
                    return fileHeader;
				}
				else if(line.endsWith(VALUE))
					return null;
				else
					return null;
			}
        }
        catch(Exception e)
        {
            StsException.outputException("StsWellKeywordIO.readFileHeader() failed.",
				e, StsException.WARNING);
			return null;
        }
    }

	// used for both well-dev.* and well-logs.* files to capture header info
	static public class FileHeader
	{
		public String wellName = "none";
		public double xOrigin = 0.0;
		public double yOrigin = 0.0;
		public float nullValue = StsParameters.nullValue;
		public String[] curveNames;
        byte depthType = StsLogVector.MDEPTH;
        byte timeType = StsLogVector.TWT;
        boolean liveWell = false;
        public int xIdx = -1;
        public int yIdx = -1;
        public int zIdx = -1;
        public int timeIdx = -1;
        public int dateIdx = -1;
	}

    /* read the well name */
    static private MultipleFileHeader readMultipleFileHeader(BufferedReader bufRdr, int numSkipRows, int[] colOrder)
    {
		MultipleFileHeader fileHeader = new MultipleFileHeader();
        int cnt = 0;
        try
        {
			while(true)
			{
				String line = bufRdr.readLine().trim();
                if(line == null) break;
                cnt++;
                if(cnt < numSkipRows) continue;
                line = StsStringUtils.detabString(line);
                StringTokenizer stok = new StringTokenizer(line,",");
                String[] tokens = new String[stok.countTokens()];
                for(int i=0; i<tokens.length; i++)
                    tokens[i] = stok.nextToken();

                String uwi = tokens[colOrder[UWI]];
                String name = tokens[colOrder[NAME]];
                double x = Double.parseDouble(tokens[colOrder[MX]]);
                double y = Double.parseDouble(tokens[colOrder[MY]]);
                String symbol = tokens[colOrder[SYMBOL]];
                float td = Float.parseFloat(tokens[colOrder[TVD]]);
                float kb = Float.parseFloat(tokens[colOrder[KB]]);
                float elev = Float.parseFloat(tokens[colOrder[GRD]]);
                float datum = Float.parseFloat(tokens[colOrder[DATUM]]);

                // TODO: Add checks for missing and valid values and substitute defaults where not valid or missing

                fileHeader.addWell(uwi,name,x,y,symbol,td,kb,elev, datum);
            }
            return fileHeader;
        }
        catch(Exception e)
        {
            StsException.outputException("StsWellKeywordIO.readMultipleFileHeader() failed.",
				e, StsException.WARNING);
			return null;
        }
    }

    // used for multiple wells in single header file.
	static class MultipleFileHeader
	{
        int numWells;
        String[] wellNames = null;
        String[] UWIs = null;
        double[] xOrigins = new double[1];
		double[] yOrigins = new double[1];
        String[] symbols = null;
        float[] kbElevations = new float[1];
        float[] elevations = new float[1];
        float[] totalDepths = new float[1];
        float[] datums = new float[1];
        float nullValue = StsParameters.nullValue;

        public void addWell(String uwi, String name, double x, double y, String symbol, float td, float kb, float elev, float datum)
        {
            numWells++;
            UWIs = (String[]) StsMath.arrayAddElement(UWIs, uwi);
            wellNames = (String[]) StsMath.arrayAddElement(wellNames, name);
            symbols = (String[]) StsMath.arrayAddElement(symbols, symbol);

            double[] tempD = new double[numWells];
            System.arraycopy(xOrigins, 0, tempD, 0, numWells);
            xOrigins = tempD;
            xOrigins[numWells-1] = x;

            tempD = new double[numWells];
            System.arraycopy(yOrigins, 0, tempD, 0, numWells);
            yOrigins = tempD;
            yOrigins[numWells-1] = y;

            float[] tempF = new float[numWells];
            System.arraycopy(totalDepths, 0, tempF, 0, numWells);
            totalDepths = tempF;
            totalDepths[numWells-1] = td;

            tempF = new float[numWells];
            System.arraycopy(kbElevations, 0, tempF, 0, numWells);
            kbElevations = tempF;
            kbElevations[numWells-1] = kb;

            tempF = new float[numWells];
            System.arraycopy(elevations, 0, tempF, 0, numWells);
            elevations = tempF;
            elevations[numWells-1] = elev;

            tempF = new float[numWells];
            System.arraycopy(datums, 0, tempF, 0, numWells);
            datums = tempF;
            datums[numWells-1] = datum;
        }
    }

    /* read the well name */
    static private String readWellName(BufferedReader bufRdr)
    {
        try
        {
            String line = bufRdr.readLine().trim();
            if (line.endsWith(WELL_NAME))
            {
                return new String(bufRdr.readLine().trim());
            }
            return null;
        }
        catch(Exception e)
        {
            return null;
        }
    }

    /* read the deviation origin */
/*
    static private boolean readXYOrigin(BufferedReader bufRdr)
    {
        try
        {
            String line = bufRdr.readLine().trim();
            if (line.indexOf(ORIGIN) >= 0)  // is origin keyword there?
            {
                boolean xyOrder = true;
                if (line.indexOf(YX) >= 0) xyOrder = false;  // determine x-y order
                line = bufRdr.readLine().trim();  // get the next line

                // tokenize the x-y values and convert to a point object
                StringTokenizer stok = new StringTokenizer(line);
                xOrigin = Double.valueOf(stok.nextToken()).doubleValue();
                yOrigin = Double.valueOf(stok.nextToken()).doubleValue();
                if (StsTrace.getTrace())System.out.println("xOrigin = "
                        + Double.toString(xOrigin) +
                        ", yOrigin = " + Double.toString(yOrigin));
                if (!xyOrder)
                {
                    double temp = xOrigin;
                    xOrigin = yOrigin;
                    yOrigin = temp;
                }
                return true;
            }
            return false;
        }
        catch (Exception e)
        {
            return false;
        }
    }
*/
    static private String[] readCurveNames(BufferedReader bufRdr)
    {
		String line;
        String[] curveNames = new String[0];

        try
        {
			if(!findKeyword(bufRdr, CURVE)) return new String[0];
			bufRdr.mark(1024);

			line = bufRdr.readLine().trim();
			while(!lineHasKeyword(line, VALUE))
			{
				curveNames = (String[])StsMath.arrayAddElement(curveNames, line);
				bufRdr.mark(1024);
				line = bufRdr.readLine().trim();
			}
			bufRdr.reset();
			return curveNames;
		}
        catch (Exception e)
        {
            StsException.outputException("StsWellKeywordIO.readCurveNames() failed.",
                e, StsException.WARNING);
			return curveNames;
        }
    }

    /* read the curve values and set the log vectors */
    static private boolean readCurveValues(BufferedReader bufRdr, String binaryDataDir, String[] curveNames,
                                           StsLogVector[] curveVectors, float logCurveNull, boolean loadValues)
    {
        try
        {
            int nCurveVectors = curveVectors.length;

			boolean binaryFilesOK = true;
		    for(int n = 0; n < nCurveVectors; n++)
            {
				if(!curveVectors[n].hasBinaryFile(binaryDataDir) || deleteBinaries)
                    binaryFilesOK = false;
                if(fileHeader.liveWell)
                {
                    if(!timeVector.hasBinaryFile(binaryDataDir))
                        binaryFilesOK = false;
                }
            }
			if(binaryFilesOK)
			{
				for (int n = 0; n < nCurveVectors; n++)
					curveVectors[n].readBinaryFile(binaryDataDir, loadValues);
                if(fileHeader.liveWell)
                    timeVector.readBinaryFile(binaryDataDir, loadValues);
				return true;
			}
            else
                deleteBinaryFiles(binaryDataDir, curveNames, StsLogVector.WELL_LOG_PREFIX);

			if(!readAsciiCurveValues(bufRdr, curveNames, curveVectors, logCurveNull, loadValues)) return false;
            return checkWriteBinaryFiles(curveVectors, binaryDataDir);
        }
        catch (Exception e)
        {
            StsException.outputException("StsWellKeywordIO.readCurveValues() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    static private boolean readAsciiCurveValues(BufferedReader bufRdr, String[] curveNames, StsLogVector[] curveVectors, float logCurveNull, boolean loadValues)
	{
		try
		{
            int depthVectorIdx = -1, owtVectorIdx = -1;
            // read lines until we hit the eof
            int nCurveVectors = curveVectors.length;
            StsFloatVector[] curveValues = new StsFloatVector[nCurveVectors];

            for (int n = 0; n < nCurveVectors; n++)
            {
                curveValues[n] = new StsFloatVector(5000, 2000);
                curveVectors[n].setValues(curveValues[n]);
				curveVectors[n].setNullValue(logCurveNull);
            }

			int nCurveNames = curveNames.length;
			int[] colAssignments = new int[nCurveNames];
			for(int n = 0; n < nCurveNames; n++)
				colAssignments[n] = -1;

			for(int n = 0; n < nCurveVectors; n++)
			{
				int col = curveVectors[n].getAsciiFileColumn();
                if(curveVectors[n].getType() == StsLogVector.DEPTH)
                   depthVectorIdx = n;

               if(curveVectors[n].getType() == StsLogVector.OWT)
                   owtVectorIdx = n;

               if(curveVectors[n].getType() == StsLogVector.TWT)
               {
                   String binaryFilename = StsLogVector.WELL_TD_PREFIX + ".bin." + name + ".TIME." + version;
                   curveVectors[n].setNameAndType(StsLogVector.types[StsLogVector.TIME], binaryFilename);
                }
				if(col < 0 || col >= nCurveNames) continue;
				colAssignments[col] = n;
			}

            curveVectors[0].getType();
		    String line;
			int nLines = 0;
            long lvalue = 0l;
            while ((line = bufRdr.readLine()) != null)
            {
                line = line.trim();
                if (line.equals("")) continue;  // blank line
                line = StsStringUtils.detabString(line);

                StringTokenizer stok = new StringTokenizer(line," \t");
				int nTokens = stok.countTokens();
                String dateToken = null, timeToken = null;
                int subtract = 0;
                int nCols = nCurveVectors;
                if(fileHeader.liveWell)
                    nCols = nCols + 2;
                if(nTokens < nCols) 
                    continue;
                for(int col = 0; col < nCols; col++)
                {
                    double value;
                    String token = StsStringUtils.deWhiteSpaceString(stok.nextToken().trim());
                    if(col == fileHeader.timeIdx)
                    {
                        timeToken = token;
                        subtract++;
                    }
                    else if(col == fileHeader.dateIdx)
                    {
                        dateToken = token;
                        subtract++;
                    }
                    else
                    {
                        try { value = Double.parseDouble(token); }
                        catch(Exception ex) { value = logCurveNull; }
					    int nVector = colAssignments[col-subtract];
					    if(nVector == -1) continue;
                        if(nVector == depthVectorIdx)
                            curveVectors[nVector].append(value + datumShift);
                        else if(nVector == owtVectorIdx)
                            curveVectors[nVector].append(value * 2.0f);
                        else
                            curveVectors[nVector].append(value);
                    }
                }
                if(fileHeader.liveWell)
                {
                    Calendar cal = null;
                    if((timeToken != null) && (dateToken != null))
                    {
                        String dateTime = "undefined";
    	                try
    	                {
                            dateTime = dateToken + " " + timeToken;
                            cal = CalendarParser.parse(dateTime, dateOrder, true);
                            lvalue = cal.getTimeInMillis();
                            timeVector.append(lvalue);
                        }
                        catch(Exception ex)
                        {
                            StsMessageFiles.infoMessage("Unable to parse the time value (" + dateTime + "), setting to 12-31-68 16:00:00");
                            timeVector.append(0l);
                        }
                    }
                }
				if(++nLines%1000 == 0) StsMessageFiles.logMessage("File: " + filename + ": " + nLines + " values read.");
            }

            if(fileHeader.liveWell)
                timeVector.setMinMaxAndNulls(0l);

            if(owtVectorIdx != -1)
            {
                String binaryFilename = StsLogVector.WELL_TD_PREFIX + ".bin." + name + ".TIME." + version;
                curveVectors[owtVectorIdx].setNameAndType(StsLogVector.types[StsLogVector.TIME], binaryFilename);
            }
            // finish up values: trim, min, max, null
//            StsLogVector[] goodCurveVectors = new StsLogVector[nCurveVectors];
//            int nGood = 0;
            for(int n = 0; n < nCurveVectors; n++)
			{
               curveValues[n].trimToSize();
               curveValues[n].setMinMax();
               if(curveVectors[n].isNull())  curveVectors[n].setValues((StsFloatVector)null);
//				if(curveValues[n].setMinMaxAndNulls(logCurveNull)) curveVectors[nGood++] = curveVectors[n];
            }
//            if(nGood < nCurveVectors) curveVectors = (StsLogVector[])StsMath.trimArray(curveVectors, nGood);
            return true;
/*
			boolean writeOK = true;
            for(int n = 0; n < nCurveVectors; n++)
			{
				curveValues[n].trimToSize();
				curveValues[n].setMinMaxAndNulls(logCurveNull);
				if(!curveVectors[n].writeBinaryFile(binaryDataDir)) writeOK = false;
				if(!loadValues) curveVectors[n].getValues().deleteValues();
			}
            return writeOK;
*/
        }
        catch (Exception e)
        {
            StsException.outputException("StsWellKehywordIO.readCurveValues() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    // test program
  	public static void main(String[] args)
  	{
        StsLogVector[] curves = null;
    	try
        {
            // Create a file dialog to query the user for a filename.
    	    Frame frame = new Frame();
   	 	    FileDialog f = new FileDialog(frame, "choose a well deviation file",
                                FileDialog.LOAD);
            f.setVisible(true);
    	    String path = f.getDirectory();
            String filename = f.getFile();

            // make a database
        	StsModel model = new StsModel();
            // read the file
            curves =  StsWellKeywordIO.readDeviation(path, path, filename, StsParameters.nullValue, StsWellImport.WELL_TD_PREFIX);
        }
        catch (Exception e)
        {
            StsException.outputException("StsWellKeywordIO.main: ",
                                e, StsException.WARNING);
        }
  	}

	static public int checkIndex(String token, int currentIndex, int index, String[] keywords)
	{
        if(index == -1)
        {
            for(int j=0; j<keywords.length; j++)
            {
                if(token.equalsIgnoreCase(keywords[j]))
                    return currentIndex;
            }
            return index;
        }
        else
        	return index;
	}
}
