
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
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

import java.awt.*;
import java.io.*;
import java.util.*;

public class StsUTKeywordIO extends StsKeywordIO
{
    static private StsLogVector xVector, yVector, zVector, mVector;
	static private long asciiFileDate = 0;

    static private final int X = 3;          // Delta Y
    static private final int Y = 4;          // Delta X
    static private final int DEPTH = 5;      // Total Vertical Depth
    static private final int DOGLEG = 6;     // Dog Leg Severity
    static private final int LRDL = 7;       // Left Right Dog Leg
    static private final int BDDL = 8;       // Build Drop Deg Leg
    static private final int CL_DIST = 9;    // Closure Distance
    static private final int CL_DIR = 10;    // Closure Direction

    static public final byte MINCURVE = 0;
    static public final byte TANGENTIAL = 1;
    static public final byte BLNCDTANG = 2;
    static public final byte AVGANGLE = 3;
    static public final byte CURVE = 4;
    static public final String[] algorithms = {"Minimum Radius of Curvature","Tangential","Balanced Tangential","Average Angle","Radius of Curvature"};

    static boolean debug = false;
    static public boolean deleteBinaries = false;
    static public byte verticalUnits = StsParameters.DIST_NONE;
    static public byte horizontalUnits = StsParameters.DIST_NONE;
    static public double startZ = 0.0f;
    static private float datumShift = 0.0f;
    static public FileHeader fileHeader = null;

    static public FileHeader getFileHeader() { return fileHeader; }
    static public String getWellName() { return name; }
    static public String getCurveName() { return subname; }

	static public void initialize(StsModel model)
	{
        ;
    }

    /** Verify that the deviation binaries exist and are up-to-date for this ASCII file */
    static public boolean verifyDeviationBinaries(String dataDir, String binaryDataDir, String filename)
    {
        BufferedReader bufRdr = null;
        StsLogVector[] vectors = null;
        try
        {
            // parse the filename to get the well name and file type & class
            if(!setParseFilename(filename, "well-dev", "ut"))
                return false;

            // open the file
            bufRdr = new BufferedReader(new FileReader(dataDir + File.separator + filename));

            // read the well name (not used)
            fileHeader = readFileHeader(bufRdr);
			String[] curveNames = { "MDEPTH" , "DRIFT", "AZIMUTH", "X", "Y", "DEPTH", "DOGLEG", "VERTICALDL", "LATERALDL", "CL_DIST", "CL_DIR" };

            // if ascii file is newer than any of the binaries, delete all binaries as they
            // are potentially out of date
            return binaryFileDatesOK(dataDir, binaryDataDir, filename, curveNames, StsLogVector.WELL_DEV_PREFIX);
        }
        catch(Exception e)
        {
            StsException.outputException("StsUtKeywordIO.readDeviation() failed.",
                e, StsException.WARNING);
            StsMessageFiles.logMessage("Unable to read deviation for well" + name + ".");
            return false;
        }

    }
    static public StsLogVector[] readDeviation(String dataDir, String binaryDataDir, String filename, StsPoint origin,
                                               byte utType, float logCurveNull, boolean deleteBins, byte vUnits,
                                               byte hUnits, String group, float shift)
    {
        deleteBinaries = deleteBins;
        verticalUnits = vUnits;
        horizontalUnits = hUnits;
        datumShift = shift;
        return readDeviation(dataDir, binaryDataDir, filename, origin, utType, logCurveNull, group);
    }

    /** read a deviation and return values as log vectors */
    static public StsLogVector[] readDeviation(String dataDir, String binaryDataDir, String filename, StsPoint origin, byte utType, float logCurveNull, String group)
    {
        BufferedReader bufRdr = null;

        StsLogVector[] deviationVectors = null;
        try
        {
            // parse the filename to get the well name and file type & class
            if(!setParseFilename(filename, group, "ut")) return null;

            // open the file
            bufRdr = new BufferedReader(new FileReader(dataDir + File.separator + filename));

            // read the well name (not used)
            fileHeader = readFileHeader(bufRdr);
			String[] curveNames = { "MDEPTH" , "DRIFT", "AZIMUTH", "X", "Y", "DEPTH", "DOGLEG", "VERTICALDL", "LATERALDL", "CL_DIST", "CL_DIR" };

			// if ascii file is newer than any of the binaries, delete all binaries as they
			// are potentially out of date
			if(!binaryFileDatesOK(dataDir, binaryDataDir, filename, curveNames, StsLogVector.WELL_DEV_PREFIX) || deleteBinaries)
				deleteBinaryFiles(binaryDataDir, curveNames, StsLogVector.WELL_DEV_PREFIX);

            deviationVectors = constructLogVectors(curveNames, StsLogVector.WELL_DEV_PREFIX);
            if(!StsLogVector.deviationVectorsOK(deviationVectors))
            {
                StsException.systemError("StsUTKeywordIO.readDeviation() failed." +
                    " Didn't find  X or Y or DEPTH|MDEPTH vectors in file: " + filename);
                return null;
            }

            xVector = StsLogVector.getVectorOfType(deviationVectors, StsLogVector.X);
            xVector.setUnits(horizontalUnits);
            yVector = StsLogVector.getVectorOfType(deviationVectors, StsLogVector.Y);
            yVector.setUnits(horizontalUnits);
            zVector = StsLogVector.getVectorOfType(deviationVectors, StsLogVector.DEPTH);
            zVector.setUnits(verticalUnits);
            mVector = StsLogVector.getVectorOfType(deviationVectors, StsLogVector.MDEPTH);
            mVector.setUnits(verticalUnits);

            // classInitialize origin to well-dev file header origin
            xVector.setOrigin(origin.getX());
            yVector.setOrigin(origin.getY());

            startZ = origin.getZ();
            // read the curve values and set the log vectors; may be from ascii or binary files
            if (!readCurveValues(bufRdr, fileHeader, binaryDataDir, curveNames, deviationVectors, logCurveNull, true, utType))
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
            StsException.outputException("StsUTKeywordIO.readDeviation() failed.",
                e, StsException.WARNING);
            StsMessageFiles.logMessage("Unable to read deviation for well" + name + ".");
            return null;
    	}
        finally
        {
            closeBufRdr(bufRdr);
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
				 file = new File(binaryDir + File.separator + binaryFilename);
				 long binaryFileDate = file.lastModified();
				 if(binaryFileDate < asciiFileDate) binaryDatesOK = false;
			 }
			 return binaryDatesOK;
		 }
		 catch (Exception e)
		 {
			 StsException.outputException("StsUTKeywordIO.logVectorDatesOK() failed.",
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
			  StsException.outputException("StsUTKeywordIO.logVectorDatesOK() failed.",
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
              StsException.outputException("StsUTKeywordIO.addLogVector() failed.",e, StsException.WARNING);
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
            StsException.outputException("StsUTKeywordIO.constructLogVectors() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

	/** read a "multiplexed" log curve file and return values as log curves */
	 static public StsLogCurve[] readLogCurves(StsWell well, String wellname, String dataDir, String binaryDataDir, String group, float logCurveNull, byte type)
	 {
		 String[] filenames = null;
		 StsLogCurve[] logCurves = new StsLogCurve[0];

		 try
		 {
			 File directory = new File(dataDir);
			 StsFilenameGroupFilter filter = new StsFilenameGroupFilter(group, "ut", wellname);
			 filenames = directory.list(filter);
			 if(filenames == null) return logCurves;
			 int nFilenames = filenames.length;
			 if(nFilenames == 0) return logCurves;

			 boolean ok = true;
			 StsLogVector[] logVectors;
			for(int n = 0; n < nFilenames; n++)
			 {
				 logVectors = readLogVectors(wellname, dataDir, binaryDataDir, filenames[n], group, logCurveNull, type);
				 if(logVectors == null)
				 {
					 StsMessageFiles.logMessage("Failed to construct log curves for file: " + filenames[n]);
					 ok = false;
				 }
				 else
				 {
					 StsLogCurve[] newLogCurves = StsLogCurve.constructLogCurves(well, logVectors, logCurveNull, 0);
					 if (newLogCurves.length == 0) ok = false;
					 logCurves = (StsLogCurve[])StsMath.arrayAddArray(logCurves, newLogCurves);
				 }
			  }
			 return logCurves;
		 }
		 catch (Exception e)
		 {
			 StsMessageFiles.logMessage("Log curve read failed for well " + wellname);
			 return logCurves;
			}
	 }

	 static public StsLogVector[] readLogVectors(String wellname, String dataDir, String binaryDataDir, String filename, String group, float logCurveNull, byte type)
	 {
		 BufferedReader bufRdr = null;
		 StsLogVector[] logVectors;
		 StsLogCurve[] logCurves = null;

		 try
		 {
			 if(!setParseFilename(filename, group, "ut")) return null;
			 bufRdr = new BufferedReader(new FileReader(dataDir + File.separator + filename));

			 fileHeader = readFileHeader(bufRdr);
			 if(fileHeader == null) return null;

			 String[] curveNames = fileHeader.curveNames;
			 float curveNullValue = fileHeader.nullValue;

			 logVectors = constructLogVectors(curveNames, group);
			 if (logVectors == null) return null;

			 if(!readCurveValues(bufRdr, fileHeader, binaryDataDir, curveNames, logVectors, logCurveNull, false, type)) return null;
			 return logVectors;
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

    /* read the well name */
    static private FileHeader readFileHeader(BufferedReader bufRdr)
    {
		FileHeader fileHeader = new FileHeader();

        // Position to the first data record
        try
        {
            while(true)
            {
                bufRdr.mark(80);
                String line = bufRdr.readLine();
                line = StsStringUtils.detabString(line);

                StringTokenizer stok = new StringTokenizer(line.trim()," ");
                String keyWord = null;
                try
                {
                    float temp = new Float(stok.nextToken().trim()).floatValue();
                    bufRdr.reset();
                    break;
                }
                catch (Exception e)
                {
                    continue;
                }
            }
            return fileHeader;
        }
        catch(Exception e)
        {
            StsException.outputException("StsLasKeywordIO.readFileHeader() failed.",e, StsException.WARNING);
            return null;
        }
    }

	// used for las files to capture header info
	static class FileHeader
	{
		String formatVersion = "unknown";
        boolean lineWrap = false;
        String wellName = "none";
        float startZ = 0.0f;
        float stopZ = 0.0f;
        float stepZ = 0.0f;
        String company = null;
        String wellId = null;
		double xOrigin = 0.0;
		double yOrigin = 0.0;
        float elevation = 0.0f;
		float nullValue = StsParameters.nullValue;
		String[] curveNames;
    }

    /* read the curve values and set the log vectors */
    static private boolean readCurveValues(BufferedReader bufRdr, FileHeader fileHdr, String binaryDataDir, String[] curveNames,
                                            StsLogVector[] curveVectors, float logCurveNull, boolean loadValues, byte utType)
    {
        try
        {
            int nCurveVectors = curveVectors.length;

			boolean binaryFilesOK = true;
		    for(int n = 0; n < nCurveVectors; n++)
				if(!curveVectors[n].hasBinaryFile(binaryDataDir))
                    binaryFilesOK = false;

			if(binaryFilesOK)
			{
				for (int n = 0; n < nCurveVectors; n++)
					curveVectors[n].readBinaryFile(binaryDataDir, loadValues);
				return true;
			}
			if(!readAsciiCurveValues(bufRdr, fileHdr, curveNames, curveVectors, logCurveNull, loadValues, utType)) return false;
            return checkWriteBinaryFiles(curveVectors, binaryDataDir);
        }
        catch (Exception e)
        {
            StsException.outputException("StsUTKeywordIO.readCurveValues() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    static private boolean readAsciiCurveValues(BufferedReader bufRdr, FileHeader fileHdr, String[] curveNames,
            StsLogVector[] curveVectors, float logCurveNull, boolean loadValues, byte type)
	{
        double prev_md = 0.0, prev_azimuth = 0.0f, prev_drift = 0.0f;
		try
		{
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
		    String line;
			int nLines = 0;
            while ((line = bufRdr.readLine()) != null)
            {
                line = line.trim();
                if(line.equals("")) continue;  // blank line
                line = StsStringUtils.detabString(line);
                if(line.substring(0).equals("H")) continue; // header line

                StringTokenizer stok = new StringTokenizer(line);
				int nTokens = stok.countTokens();
                if(nTokens < 3) continue; // Line too short

                try
                {
                    // Assumes first three columns are MD, Drift and Azimuth, with remaining three calculated
                    double md = Double.parseDouble(stok.nextToken()) + datumShift;
                    if (md == fileHdr.nullValue)
                        md = logCurveNull;
                    curveVectors[0].append(md);

                    double drift = Double.parseDouble(stok.nextToken());
                    curveVectors[1].append(drift);

                    double azimuth = Double.parseDouble(stok.nextToken());
                    curveVectors[2].append(azimuth);

                    switch (type) {
                        case MINCURVE:
                            calculateVectorsWithMinCurve(curveVectors, md,
                                prev_md, drift, prev_drift, azimuth,
                                prev_azimuth, nLines);
                            break;
                        case TANGENTIAL:
                            calculateVectorsWithTangential(curveVectors, md,
                                prev_md, drift, prev_drift, azimuth,
                                prev_azimuth, nLines);
                            break;
                        case BLNCDTANG:
                            calculateVectorsWithBalTang(curveVectors, md,
                                prev_md, drift, prev_drift, azimuth,
                                prev_azimuth, nLines);
                            break;
                        case AVGANGLE:
                            calculateVectorsWithAvgAngle(curveVectors, md,
                                prev_md, drift, prev_drift, azimuth,
                                prev_azimuth, nLines);
                            break;
                        case CURVE:
                            calculateVectorsWithCurvature(curveVectors, md,
                                prev_md, drift, prev_drift, azimuth,
                                prev_azimuth, nLines);
                            break;
                        default:
                            break;
                    }

                    prev_md = md;
                    prev_drift = drift;
                    prev_azimuth = azimuth;

                    if (++nLines % 1000 == 0)
                        StsMessageFiles.logMessage("File: " + filename + ": " +  nLines + " values read.");
                }
                catch (Exception e)
                {
                    continue;
                }
            }

            // finish up values: trim, min, max, null
            for(int n = 0; n < nCurveVectors; n++)
			{
				curveValues[n].trimToSize();
                curveValues[n].setMinMax();
                if(curveVectors[n].isNull())
                    curveVectors[n].setValues((StsFloatVector)null);
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsUTKeywordIO.readCurveValues() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    public static void calculateVectorsWithMinCurve(StsLogVector[] vectors, double md, double pmd, double da, double pda, double az, double paz, int idx)
    {
        double pi180 = Math.PI/180.0f;
        double dl, rf, dTVD, dLAT, dDep, TNS = 0.0f, TEW = 0.0f, TVD = 0.0f;
        double courseLen, closureDir, closureDist, VertSec;

        if(idx == 0)
        {
            TNS = 0.0f;
            TEW = 0.0f;
            TVD = startZ;
        }
        else
        {
            TNS = vectors[Y].getFloats()[idx-1];
            TEW = vectors[X].getFloats()[idx-1];
            TVD = vectors[DEPTH].getFloats()[idx-1];
        }
        dl = Math.acos(Math.cos((da-pda)*pi180) - Math.sin(da*pi180) * Math.sin(pda*pi180)
                       * (1-Math.cos((az-paz) * pi180))) * 180.0f/Math.PI;
        if(dl == 0.0f)
            rf = 1.0f;
        else
            rf = 360 * (1- Math.cos(dl * pi180))/dl/Math.PI/Math.sin(dl*pi180);

        courseLen = md - pmd;
        dTVD = (float)(courseLen * rf/2.0f * (Math.cos(da * pi180) + Math.cos(pda*pi180)));
        dLAT = (float)(courseLen * rf/2.0f * (Math.sin(da*pi180) * Math.cos(az*pi180) +
                                              Math.sin(pda*pi180) * Math.cos(paz*pi180)));
        dDep = (float)(courseLen * rf/2.0f * (Math.sin(da*pi180) * Math.sin(az*pi180) +
                                             Math.sin(pda*pi180) * Math.sin(paz*pi180)));
        TNS = TNS + dLAT;
        TEW = TEW + dDep;
        TVD = TVD + dTVD;

        double[] dogLegs = null;
        dogLegs = calcDogLegs(da, pda, az, paz, courseLen);
        closureDist = calcClosureDistance(TNS, TEW);
        closureDir = calcClosureDirection(dDep, closureDist, TNS);

        vectors[X].append(TEW);
        vectors[Y].append(TNS);
        vectors[DEPTH].append(TVD);
        vectors[DOGLEG].append(dogLegs[0]);
        vectors[LRDL].append(dogLegs[1]);
        vectors[BDDL].append(dogLegs[2]);
        vectors[CL_DIST].append(closureDist);
        vectors[CL_DIR].append(closureDir);
    }

    public static void calculateVectorsWithCurvature(StsLogVector[] vectors, double md, double pmd, double da, double pda, double az, double paz, int idx)
    {
        double pi180 = Math.PI/180.0f;
        double dl, rf, dTVD, dLAT, dDep, TNS = 0.0f, TEW = 0.0f, TVD = 0.0f;
        double courseLen, closureDir, closureDist, VertSec;

        if(idx == 0)
        {
            TNS = 0.0f;
            TEW = 0.0f;
            TVD = 0.0f;
        }
        else
        {
            TNS = vectors[Y].getFloats()[idx-1];
            TEW = vectors[X].getFloats()[idx-1];
            TVD = vectors[DEPTH].getFloats()[idx-1];
        }
        courseLen = md - pmd;
        dTVD = (float)(courseLen/(da-pda)) * (180.0f/Math.PI) * (Math.sin(da*pi180) - Math.sin(pda*pi180));
        dLAT = (float)(courseLen/(da-pda)) * (180.0f/Math.PI) * (180.0f/Math.PI) *
            (((Math.cos(pda*pi180) - Math.cos(da*pi180)) * (Math.sin(paz*pi180) - Math.sin(az*pi180))) /
            (az*pi180 - paz*pi180));
        dDep = (float)(courseLen/(da-pda)) * (180.0f/Math.PI) * (180.0f/Math.PI) *
            (((Math.cos(pda*pi180) - Math.cos(da*pi180)) * (Math.cos(paz*pi180) - Math.cos(az*pi180))) /
            (az*pi180 - paz*pi180));

        TNS = TNS + dLAT;
        TEW = TEW + dDep;
        TVD = TVD + dTVD;

        double[] dogLegs = null;
        dogLegs = calcDogLegs(da, pda, az, paz, courseLen);
        closureDist = calcClosureDistance(TNS, TEW);
        closureDir = calcClosureDirection(dDep, closureDist, TNS);

        vectors[X].append(TEW);
        vectors[Y].append(TNS);
        vectors[DEPTH].append(TVD);
        vectors[DOGLEG].append(dogLegs[0]);
        vectors[LRDL].append(dogLegs[1]);
        vectors[BDDL].append(dogLegs[2]);
        vectors[CL_DIST].append(closureDist);
        vectors[CL_DIR].append(closureDir);
    }

    public static void calculateVectorsWithBalTang(StsLogVector[] vectors, double md, double pmd, double da, double pda, double az, double paz, int idx)
    {
        double pi180 = Math.PI/180.0f;
        double dl, rf, dTVD, dLAT, dDep, TNS = 0.0f, TEW = 0.0f, TVD = 0.0f;
        double courseLen, closureDir, closureDist, VertSec;

        if(idx == 0)
        {
            TNS = 0.0f;
            TEW = 0.0f;
            TVD = 0.0f;
        }
        else
        {
            TNS = vectors[Y].getFloats()[idx-1];
            TEW = vectors[X].getFloats()[idx-1];
            TVD = vectors[DEPTH].getFloats()[idx-1];
        }
        dl = Math.acos(Math.cos((da-pda)*pi180) - Math.sin(da*pi180) * Math.sin(pda*pi180)
                       * (1-Math.cos((az-paz) * pi180))) * 180.0f/Math.PI;
        if(dl == 0.0f)
            rf = 1.0f;
        else
            rf = 360 * (1- Math.cos(dl * pi180))/dl/Math.PI/Math.sin(dl*pi180);

        courseLen = md - pmd;
        dTVD = (float)(courseLen/2.0f * (Math.cos(da * pi180) + Math.cos(pda*pi180)));
        dLAT = (float)(courseLen/2.0f * (Math.sin(da*pi180) * Math.cos(az*pi180) +
                                              Math.sin(pda*pi180) * Math.cos(paz*pi180)));
        dDep = (float)(courseLen/2.0f * (Math.sin(da*pi180) * Math.sin(az*pi180) +
                                             Math.sin(pda*pi180) * Math.sin(paz*pi180)));
        TNS = TNS + dLAT;
        TEW = TEW + dDep;
        TVD = TVD + dTVD;

        double[] dogLegs = null;
        dogLegs = calcDogLegs(da, pda, az, paz, courseLen);
        closureDist = calcClosureDistance(TNS, TEW);
        closureDir = calcClosureDirection(dDep, closureDist, TNS);

        vectors[X].append(TEW);
        vectors[Y].append(TNS);
        vectors[DEPTH].append(TVD);
        vectors[DOGLEG].append(dogLegs[0]);
        vectors[LRDL].append(dogLegs[1]);
        vectors[BDDL].append(dogLegs[2]);
        vectors[CL_DIST].append(closureDist);
        vectors[CL_DIR].append(closureDir);
    }

    public static void calculateVectorsWithTangential(StsLogVector[] vectors, double md, double pmd, double da, double pda, double az, double paz, int idx)
    {
        double pi180 = Math.PI/180.0f;
        double dl, rf, dTVD, dLAT, dDep, TNS = 0.0f, TEW = 0.0f, TVD = 0.0f;
        double courseLen, closureDir, closureDist, VertSec;

        if(idx == 0)
        {
            TNS = 0.0f;
            TEW = 0.0f;
            TVD = 0.0f;
        }
        else
        {
            TNS = vectors[Y].getFloats()[idx-1];
            TEW = vectors[X].getFloats()[idx-1];
            TVD = vectors[DEPTH].getFloats()[idx-1];
        }
        courseLen = md - pmd;
        dTVD = (float)(courseLen * (Math.cos(da * pi180)));
        dLAT = (float)(courseLen * (Math.sin(da*pi180) * Math.cos(az*pi180)));
        dDep = (float)(courseLen * (Math.sin(da*pi180) * Math.sin(az*pi180)));
        TNS = TNS + dLAT;
        TEW = TEW + dDep;
        TVD = TVD + dTVD;

        double[] dogLegs = null;
        dogLegs = calcDogLegs(da, pda, az, paz, courseLen);
        closureDist = calcClosureDistance(TNS, TEW);
        closureDir = calcClosureDirection(dDep, closureDist, TNS);

        vectors[X].append(TEW);
        vectors[Y].append(TNS);
        vectors[DEPTH].append(TVD);
        vectors[DOGLEG].append(dogLegs[0]);
        vectors[LRDL].append(dogLegs[1]);
        vectors[BDDL].append(dogLegs[2]);
        vectors[CL_DIST].append(closureDist);
        vectors[CL_DIR].append(closureDir);
    }

    public static void calculateVectorsWithAvgAngle(StsLogVector[] vectors, double md, double pmd, double da, double pda, double az, double paz, int idx)
    {
        double pi180 = Math.PI/180.0f;
        double dl, rf, dTVD, dLAT, dDep, TNS = 0.0f, TEW = 0.0f, TVD = 0.0f;
        double courseLen, closureDir, closureDist, VertSec;

        if(idx == 0)
        {
            TNS = 0.0f;
            TEW = 0.0f;
            TVD = 0.0f;
        }
        else
        {
            TNS = vectors[Y].getFloats()[idx-1];
            TEW = vectors[X].getFloats()[idx-1];
            TVD = vectors[DEPTH].getFloats()[idx-1];
        }
        courseLen = md - pmd;
        dTVD = (float)(courseLen * Math.cos(((pda*pi180)+(da*pi180))/2.0f));
        dLAT = (float)(courseLen * Math.sin(((pda*pi180)+(da*pi180))/2.0f) * Math.cos(((paz*pi180)+(az*pi180))/2.0f));
        dDep = (float)(courseLen * Math.sin(((pda*pi180)+(da*pi180))/2.0f) * Math.sin(((paz*pi180)+(az*pi180))/2.0f));
        TNS = TNS + dLAT;
        TEW = TEW + dDep;
        TVD = TVD + dTVD;

        double[] dogLegs = null;
        dogLegs = calcDogLegs(da, pda, az, paz, courseLen);
        closureDist = calcClosureDistance(TNS, TEW);
        closureDir = calcClosureDirection(dDep, closureDist, TNS);

        vectors[X].append(TEW);
        vectors[Y].append(TNS);
        vectors[DEPTH].append(TVD);
        vectors[DOGLEG].append(dogLegs[0]);
        vectors[LRDL].append(dogLegs[1]);
        vectors[BDDL].append(dogLegs[2]);
        vectors[CL_DIST].append(closureDist);
        vectors[CL_DIR].append(closureDir);
    }

    static private double calcClosureDistance(double deltaNS, double deltaEW)
    {
        return Math.sqrt((Math.abs(deltaNS) * Math.abs(deltaNS)) + (Math.abs(deltaEW) * Math.abs(deltaEW)));
    }

    static private double calcClosureDirection(double deltaDepth, double closureDist, double deltaNS)
    {
        if(deltaDepth < 0.0f)
            return 360.0f - Math.acos(deltaNS/closureDist) * 180.0f/Math.PI;
        else
        {
            if(closureDist > 0.0f)
                return Math.acos(deltaNS/closureDist) * 180.0f/Math.PI;
            else
                return 0.0f;
        }
    }

    static private double[] calcDogLegs(double drift, double pdrift, double azimuth, double pazimuth, double courseLen)
    {
        double[] dogLegs = new double[3];
        double pi180 = Math.PI/180.0f;

        dogLegs[0] = (float)((Math.acos((Math.sin(pdrift*pi180) * Math.sin(drift*pi180)) * (Math.sin(pazimuth*pi180) *
                                    Math.sin(azimuth*pi180) + Math.cos(pazimuth*pi180) * Math.cos(azimuth*pi180)) +
                                    (Math.cos(pdrift*pi180) * Math.cos(drift*pi180)))) * 100.0f/courseLen) *180.0f/Math.PI;
        if(drift == 0.0f)
            dogLegs[2] = 0.0f;
        else
            dogLegs[2] = (drift - pdrift)/courseLen * 100.0f;

        if(azimuth == 0.0f)
            dogLegs[1] = 0.0f;
        else
            dogLegs[1] = (azimuth - pazimuth)/courseLen * 100.0f;

        return dogLegs;
    }

    static private boolean checkAddMDepthToDev(StsWell well, StsLogVector[] curveLogVectors)
    {
        StsLogVector mDepthVector = StsLogVector.getVectorOfType(curveLogVectors, StsLogVector.MDEPTH);
        StsLogVector depthVector = StsLogVector.getVectorOfType(curveLogVectors, StsLogVector.DEPTH);
        return mDepthVector != null && depthVector != null;
        // return well.checkAddMDepthVector(mDepthVector, depthVector);
    }

    /** Add logs to well */
    static public StsLogCurve[] addLogCurves(StsWell well, StsLogVector[] vectors, String wellname, float logCurveNull)
    {
        StsLogCurve[] logCurves = new StsLogCurve[0];
        try
        {
            StsLogCurve[] newLogCurves = StsLogCurve.constructLogCurves(well, vectors, logCurveNull, 0);
            logCurves = (StsLogCurve[])StsMath.arrayAddArray(logCurves, newLogCurves);
            if(well != null)
                checkAddMDepthToDev(well, vectors);
            return logCurves;
        }
        catch (Exception e)
        {
            StsMessageFiles.logMessage("Log curve read failed for well " + wellname);
            return logCurves;
        }
    }

    public static byte getAlgorithmTypeFromName(String name)
    {
        if(name.equalsIgnoreCase(StsUTKeywordIO.algorithms[0]))
            return StsUTKeywordIO.MINCURVE;
        else if(name.equalsIgnoreCase(StsUTKeywordIO.algorithms[1]))
            return StsUTKeywordIO.TANGENTIAL;
        else if(name.equalsIgnoreCase(StsUTKeywordIO.algorithms[2]))
            return StsUTKeywordIO.BLNCDTANG;
        else if(name.equalsIgnoreCase(StsUTKeywordIO.algorithms[3]))
            return StsUTKeywordIO.AVGANGLE;
        else if(name.equalsIgnoreCase(StsUTKeywordIO.algorithms[4]))
            return StsUTKeywordIO.CURVE;
        else
            return -1;
    }

    public static String getAlgorithmNameFromType(byte type)
    {
        if(type == MINCURVE)
            return StsUTKeywordIO.algorithms[0];
        else if(type == TANGENTIAL)
            return StsUTKeywordIO.algorithms[1];
        else if(type == BLNCDTANG)
            return StsUTKeywordIO.algorithms[2];
        else if(type == AVGANGLE)
            return StsUTKeywordIO.algorithms[3];
        else if(type == CURVE)
            return StsUTKeywordIO.algorithms[4];
        else
            return "Unknown";
    }

    // test program
  	public static void main(String[] args)
  	{
        StsLogVector[] curves = null;
        byte type = StsUTKeywordIO.MINCURVE;
    	try
        {
            // Create a file dialog to query the user for a filename.
    	    Frame frame = new Frame();
   	 	    FileDialog f = new FileDialog(frame, "choose a well deviation file", FileDialog.LOAD);
            f.setVisible(true);
    	    String path = f.getDirectory();
            String filename = f.getFile();


            // make a database
        	StsModel model = new StsModel();
            // read the file
            curves =  StsUTKeywordIO.readDeviation(path, path, filename, new StsPoint(0.0f,0.0f,0.0f), type, StsParameters.nullValue, StsWellImport.WELL_DEV_PREFIX);
        }
        catch (Exception e)
        {
            StsException.outputException("StsWellKeywordIO.main: ",
                                e, StsException.WARNING);
        }
  	}
}
