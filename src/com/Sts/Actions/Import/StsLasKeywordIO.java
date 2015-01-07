
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
import com.Sts.Utilities.*;

import java.awt.*;
import java.io.*;
import java.util.*;

public class StsLasKeywordIO extends StsKeywordIO
{
    /** string constants for well file type/class */

    /** keyword constants */
    static public final String WELL_HDR = "~WELL";
    static public final String WELL_NAME = "WELL";
    static public final String WELL_LABEL = "LABL";
    static public final String WELL_NUMBER = "NUMB";
    static public final String X_ORIGIN = "XCOORD";
    static public final String Y_ORIGIN = "YCOORD";
	static public final String NULL_VALUE = "NULL";
	static public final String Z_START = "STRT";
	static public final String Z_STOP = "STOP";
	static public final String Z_STEP = "STEP";
    static public final String API = "API";
    static public final String OPERATOR = "OPER";
    static public final String COMPANY = "COMP";
    static public final String WID = "UWI";
    static public final String ELEV = "ELEV ";
    static public final String DATE = "DATE ";

    static public final String VER_HDR = "~VER";
    static public final String FMT_VERSION = "VERS";
    static public final String LINEWRAP = "WRAP";

    static public final String CRV_HDR = "~CUR";
    static public final String PARAM_HDR = "~PAR";

    static public final int SPACE = 0;
    static public final int TAB = 1;

    static public final String DEPTH = "DEPTH";
    static public final String MDEPTH = "MDEPTH";
    static public final String SSDEPTH = "SSDEPTH";
    static public final String TIME = "TIME";
    static public final String[] Z_KEYWORDS = { DEPTH, MDEPTH, SSDEPTH, TIME };

    static public final String DATA_HDR = "~A";
    static public final String OTHER_HDR = "~OTH";

    static private StsLogVector xVector, yVector, zVector, mVector, tVector;
	static private long asciiFileDate = 0;
    static private float datumShift = 0.0f;
    static private boolean deleteBinaries = false;
    static private float vScalar = 1.0f;
    static private byte verticalUnits = StsParameters.DIST_NONE;
    static private byte horizontalUnits = StsParameters.DIST_NONE;
    static FileHeader fileHeader = null;

    static public String getWellName() { return name; }
    static public String getCurveName() { return subname; }
    static public FileHeader getFileHeader() { return fileHeader; }

	static public void initialize(StsModel model)
	{
//        StsKeywordIO.classInitialize(model);
//		logCurveNull = project.getLogNull();
    }

    static public StsLogVector[] readDeviation(String dataDir, String binaryDataDir, String filename,
                                               float logCurveNull, boolean deleteBins, byte vUnits, byte hUnits,
                                               String group, float shift)
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
        try
        {
            // parse the filename to get the well name and file type & class
            if(!setParseFilename(filename, group, "las")) return null;

            // open the file
            bufRdr = new BufferedReader(new FileReader(dataDir + File.separator + filename));

            // read the well name (not used)
            fileHeader = readFileHeader(bufRdr);
			String[] curveNames = fileHeader.curveNames;

			// if ascii file is newer than any of the binaries, delete all binaries as they
			// are potentially out of date
			if(!binaryFileDatesOK(dataDir, binaryDataDir, filename, curveNames, StsLogVector.WELL_DEV_PREFIX))
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

            // classInitialize origin to well-dev file header origin
            xVector.setOrigin(fileHeader.xOrigin);
            yVector.setOrigin(fileHeader.yOrigin);
            xVector.setUnits(horizontalUnits);
            yVector.setUnits(horizontalUnits);

            // read the curve values and set the log vectors; may be from ascii or binary files
            if (!readCurveValues(bufRdr, fileHeader, binaryDataDir, curveNames, deviationVectors, logCurveNull, true))
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
            StsException.outputException("StsWellKeywordIO.readDeviation() failed.",
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
			 StsException.outputException("StsWellKeywordIO.logVectorDatesOK() failed.",
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
			  StsException.outputException("StsWellKeywordIO.logVectorDatesOK() failed.",
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
              StsException.outputException("StsWellKeywordIO.addLogVector() failed.",e, StsException.WARNING);
              return null;
          }
      }

//    static private StsLogVector[] constructLogVectors(StsLogVector[] vectors, String[] curveNames, String filePrefix, boolean inDev)
//    {
//        StsLogVector[] vecs = constructLogVectors(curveNames, filePrefix, inDev);
//        for(int i=0; i<vecs.length; i++)
//            vectors = (StsLogVector[])StsMath.arrayAddElement(vectors, vecs[i]);
//        return vectors;
//    }

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

//    static private boolean inDevTypes(String name)
//    {
//        for(int i=0; i< DEV_KEYWORDS.length; i++)
//            if(name.equals(DEV_KEYWORDS[i])) return true;
//        return false;
//    }

	/** read a "multiplexed" log curve file and return values as log curves */
/*
	 static public StsLogCurve[] readLogCurves(String wellname, String dataDir, String binaryDataDir, String prefix, float logCurveNull)
	 {
		 String[] filenames = null;
		 StsLogCurve[] logCurves = new StsLogCurve[0];

		 try
		 {
			 File directory = new File(dataDir);
			 String filterString = prefix + ".las." + wellname;
			 StsFilenameFilter filter = new StsFilenameFilter(filterString, true);
			 filenames = directory.list(filter);
			 if(filenames == null) return logCurves;
			 int nFilenames = filenames.length;
			 if(nFilenames == 0) return logCurves;

			 boolean ok = true;
			 StsLogVector[] logVectors;
			for(int n = 0; n < nFilenames; n++)
			 {
				 logVectors = readLogVectors(null, wellname, dataDir, binaryDataDir, filenames[n], prefix, logCurveNull);
				 if(logVectors == null)
				 {
					 StsMessageFiles.logMessage("Failed to construct log curves for file: " + filenames[n]);
					 ok = false;
				 }
				 else
				 {
					 StsLogCurve[] newLogCurves = StsLogCurve.constructLogCurves(logVectors, logCurveNull);
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
*/
     /** Add logs to well */
/*
     static public StsLogCurve[] addLogCurves(StsWell well, StsLogVector[] vectors, String wellname, float logCurveNull)
     {
         StsLogCurve[] logCurves = new StsLogCurve[0];
         try
         {
             StsLogCurve[] newLogCurves = StsLogCurve.constructLogCurves(vectors, logCurveNull);
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

     static private boolean checkAddMDepthToDev(StsWell well, StsLogVector[] curveLogVectors)
     {
         StsLogVector mDepthVector = StsLogVector.getVectorOfType(curveLogVectors, StsLogVector.MDEPTH);
         StsLogVector depthVector = StsLogVector.getVectorOfType(curveLogVectors, StsLogVector.DEPTH);
         if(mDepthVector == null || depthVector == null) return false;
         return well.checkAddMDepthVector(mDepthVector, depthVector);
     }
*/
     static public StsLogCurve[] readLogCurves(StsWell well, String wellname, String dataDir, String binaryDataDir, String group,
                                          float logCurveNull, boolean deleteBins, byte vUnits, float shift)
     {
         datumShift = shift;
         return readLogCurves(well, wellname, dataDir, binaryDataDir, group, logCurveNull, deleteBins, vUnits);
     }
     static public StsLogCurve[] readLogCurves(StsWell well, String wellname, String dataDir, String binaryDataDir, String group,
                                               float logCurveNull, boolean deleteBins, byte vUnits)
     {
         deleteBinaries = deleteBins;
         verticalUnits = vUnits;
         return readLogCurves(well, wellname, dataDir, binaryDataDir, group, logCurveNull);
     }

     /** read a "multiplexed" log curve file and return values as log curves */
     static public StsLogCurve[] readLogCurves(StsWell well, String wellname, String dataDir, String binaryDataDir, String group, float logCurveNull)
     {
         String[] filenames = null;
         StsLogCurve[] logCurves = new StsLogCurve[0];

         try
         {
             File directory = new File(dataDir);
             StsFilenameGroupFilter filter = new StsFilenameGroupFilter(group, "las", wellname);
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
             StsMessageFiles.logMessage("Log curve read failed for well " + wellname);
             return logCurves;
         }
    }

     static public StsLogVector[] readLogVectors(String wellname, String dataDir, String binaryDataDir, String filename,
                                                 String group, float logCurveNull, boolean deleteBins, byte vUnits)
	 {
         deleteBinaries = deleteBins;
         verticalUnits = vUnits;
         return readLogVectors(wellname, dataDir, binaryDataDir, filename, group, logCurveNull);

     }

	 static public StsLogVector[] readLogVectors(String wellname, String dataDir, String binaryDataDir, String filename, String group, float logCurveNull)
	 {
		 BufferedReader bufRdr = null;
		 StsLogVector[] logVectors;
		 StsLogCurve[] logCurves = null;

		 try
		 {
			 if(!setParseFilename(filename, group, "las")) return null;
			 bufRdr = new BufferedReader(new FileReader(dataDir + File.separator + filename));

			 fileHeader = readFileHeader(bufRdr);
			 if(fileHeader == null)
                 return null;

			 String[] curveNames = fileHeader.curveNames;

			 logVectors = constructLogVectors(curveNames, group);
			 if (logVectors == null) return null;

			 if(!readCurveValues(bufRdr, fileHeader, binaryDataDir, curveNames, logVectors, logCurveNull, false))
                 return null;
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

    /* read the well name */
    static public FileHeader readFileHeader(BufferedReader bufRdr)
    {
		FileHeader fileHeader = new FileHeader();

        final int VERSION = 0;
        final int WELLINFO = 1;
        final int CRVINFO = 2;
        final int CRVDATA = 3;
        final int CRVOTHER = 4;
        final int CRVPARMS = 5;

        int delimiter = SPACE;

        int section = -1;
        boolean foundMd = false;

        try
        {
			while(true)
			{

				String line = bufRdr.readLine();
                if(line == null) return null;
                line = StsStringUtils.detabString(line);
                line.trim();
                if(line.length() < 2) continue;
                if(line.indexOf("#") == 0) continue;             // Comment

                if(isSection(line, VER_HDR))
                    section = VERSION;
                else if(isSection(line, WELL_HDR))
                    section = WELLINFO;
                else if(isSection(line, CRV_HDR))
                    section = CRVINFO;
                else if(isSection(line, OTHER_HDR))
                    section = CRVOTHER;
                else if(isSection(line, PARAM_HDR))
                    section = CRVPARMS;
                else if(isSection(line, DATA_HDR))
                {
                    section = CRVDATA;
                    return fileHeader;
                }
                else
                {
                    StringTokenizer stok = new StringTokenizer(line,":");
                    String keyWord = null;
                    switch(section)
                    {
                        case VERSION:
                            keyWord = stok.nextToken().trim();
                            if(keyWord.indexOf(FMT_VERSION) == 0)
                            {
                                fileHeader.formatVersion = keyWord.substring(keyWord.indexOf(".")+1, keyWord.length()).trim();
                            }
                            else if(keyWord.indexOf(LINEWRAP) == 0)
                            {
                                if((keyWord.indexOf("NO") != 0) ||
                                                (keyWord.indexOf("No") != 0) ||
                                                               (keyWord.indexOf("no") != 0))
                                    fileHeader.lineWrap = false;
                                else
                                    fileHeader.lineWrap = true;
                            }
                            break;
                        case WELLINFO:
                            keyWord = stok.nextToken().trim();
                            if(keyWord.indexOf(WELL_NAME) == 0)
                            {
                                keyWord = keyWord.substring(keyWord.indexOf("."),keyWord.length());
                                if(keyWord.length() > 1) fileHeader.wellName = keyWord.substring(keyWord.indexOf(" "), keyWord.length()).trim();
                            }
                            else if(keyWord.indexOf(API) == 0)
                            {
                                keyWord = keyWord.substring(keyWord.indexOf("."),keyWord.length());
                                if(keyWord.length() > 1) fileHeader.uwiNumber = keyWord.substring(keyWord.indexOf(" "), keyWord.length()).trim();
                            }
                            else if(keyWord.indexOf(COMPANY) == 0)
                            {
                                keyWord = keyWord.substring(keyWord.indexOf("."),keyWord.length());
                                if(keyWord.length() > 1) fileHeader.company = keyWord.substring(keyWord.indexOf(" "), keyWord.length()).trim();
                            }
                            else if(keyWord.indexOf(DATE) == 0)
                            {
                                keyWord = keyWord.substring(keyWord.indexOf("."),keyWord.length());
                                if(keyWord.length() > 1) fileHeader.spudDate = keyWord.substring(keyWord.indexOf(" "), keyWord.length()).trim();
                            }
                            else if(keyWord.indexOf(WELL_LABEL) == 0)
                            {
                                keyWord = keyWord.substring(keyWord.indexOf("."),keyWord.length());
                                if(keyWord.length() > 1) fileHeader.wellLabel = keyWord.substring(keyWord.indexOf(" "), keyWord.length()).trim();
                            }
                            else if(keyWord.indexOf(WELL_NUMBER) == 0)
                            {
                                keyWord = keyWord.substring(keyWord.indexOf("."),keyWord.length());
                                if(keyWord.length() > 1) fileHeader.wellNumber = keyWord.substring(keyWord.indexOf(" "), keyWord.length()).trim();
                            }
                            else if(keyWord.indexOf(OPERATOR) == 0)
                            {
                                keyWord = keyWord.substring(keyWord.indexOf("."),keyWord.length());
                                if(keyWord.length() > 1) fileHeader.operator = keyWord.substring(keyWord.indexOf(" "), keyWord.length()).trim();
                            }
                            else if(keyWord.indexOf(X_ORIGIN) == 0)
                            {
                                keyWord = keyWord.substring(keyWord.indexOf("."),keyWord.length());
                                if(keyWord.length() > 1) fileHeader.xOrigin = new Double(keyWord.substring(keyWord.indexOf(" "), keyWord.length())).doubleValue();
                            }
                            else if(keyWord.indexOf(Y_ORIGIN) == 0)
                            {
                                keyWord = keyWord.substring(keyWord.indexOf("."),keyWord.length());
                                if(keyWord.length() > 1) fileHeader.yOrigin = new Double(keyWord.substring(keyWord.indexOf(" "),keyWord.length())).doubleValue();
                            }
                            else if(keyWord.indexOf(NULL_VALUE) == 0)
                            {
                                keyWord = keyWord.substring(keyWord.indexOf("."),keyWord.length());
                                if(keyWord.length() > 1) fileHeader.nullValue = new Float(keyWord.substring(keyWord.indexOf(" "),keyWord.length())).floatValue();
                            }
                            else if(keyWord.indexOf(Z_START) == 0)
                            {
                                keyWord = keyWord.substring(keyWord.indexOf("."),keyWord.length());
                                if(keyWord.length() > 1) fileHeader.startZ = new Float(keyWord.substring(keyWord.indexOf(" "),keyWord.length())).floatValue();
                            }
                            else if(keyWord.indexOf(Z_STOP) == 0)
                            {
                                keyWord = keyWord.substring(keyWord.indexOf("."),keyWord.length());
                                if(keyWord.length() > 1) fileHeader.stopZ = new Float(keyWord.substring(keyWord.indexOf(" "),keyWord.length())).floatValue();
                            }
                            else if(keyWord.indexOf(Z_STEP) == 0)
                            {
                                keyWord = keyWord.substring(keyWord.indexOf("."),keyWord.length());
                                if(keyWord.length() > 1) fileHeader.stepZ = new Float(keyWord.substring(keyWord.indexOf(" "),keyWord.length())).floatValue();
                            }
                            else if(keyWord.indexOf(WID) == 0)
                            {
                                keyWord = keyWord.substring(keyWord.indexOf("."),keyWord.length());
                                if(keyWord.length() > 1) fileHeader.wellIdNumber = keyWord.substring(keyWord.indexOf(" "),keyWord.length()).trim();
                            }
                            else if(keyWord.indexOf(ELEV) == 0)
                            {
                                keyWord = keyWord.substring(keyWord.indexOf("."),keyWord.length());
                                if(keyWord.length() > 1) fileHeader.datumElevation = new Float(keyWord.substring(keyWord.indexOf(" "),keyWord.length())).floatValue();
                            }
                            break;
                        case CRVINFO:
                            String[] curveNames = new String[0];
                            while(!lineHasKeyword(line, DATA_HDR) && !lineHasKeyword(line, OTHER_HDR) && !lineHasKeyword(line, PARAM_HDR))
                            {
                                if((line.substring(0) != "#") && (line.length() > 0) && (line.indexOf(".") > 0))
                                {
                                    line = line.substring(0,line.indexOf(".")).trim();
                                    curveNames = (String[])StsMath.arrayAddElement(curveNames, line);
                                    if(isMeasuredDepth(line))
                                        foundMd = true;
                                }
                                else
                                    break;
                                bufRdr.mark(80);
                                line = bufRdr.readLine();
                                if(line == null) break;
                                line.trim();
                            }
                            bufRdr.reset();
                            fileHeader.curveNames = curveNames;
						// Jaguar-Pemex hack for now.  For td curve, we don't want mdepth, so don't switch depth to mdepth.
						/*
                            if(!foundMd)
                            {
                                StsMessageFiles.infoMessage("Unable to find MDEPTH in file, assuming DEPTH is measured depth.");
                                setDepthToMDepth(curveNames);
                            }
	                    */
                            break;
                        case CRVOTHER:
                            break;
                        case CRVPARMS:
                            break;
                        default:
                            System.out.println("Unrecognized Line in LAS File.");
                            break;
                    }
                }
			}
        }
        catch(Exception e)
        {
            StsException.outputException("StsLasKeywordIO.readFileHeader() failed.",e, StsException.WARNING);
			return null;
        }
    }

	static private boolean isSection(String line, String sectionString)
	{
		if(line.length() < sectionString.length())
           sectionString = sectionString.substring(0,line.length());
		return line.substring(0, sectionString.length()).toUpperCase().equals(sectionString);
	}

    private static void setDepthToMDepth(String[] names)
    {
        for(int i=0; i<names.length; i++)
        {
            if(findType(names[i]) == StsLogVector.DEPTH)
            {
                names[i] = "MDEPTH";
                break;
            }
        }
    }

    private static byte findType(String name)
    {
        String clippedType = null;
        int nTypes = StsLogVector.types.length;

        for(int n = 0; n < nTypes; n++)
        {
            if(name.length() < StsLogVector.types[n].length())
                clippedType = StsLogVector.types[n].substring(0,name.length());
            else
                clippedType = StsLogVector.types[n];

            if (name.equals(clippedType))
            {
                return (byte)n;
            }
        }
        return (byte) -1;
    }

    private static boolean isMeasuredDepth(String name)
    {
        String clippedType = null;
        int nTypes = StsLogVector.types.length;

        for(int n = 0; n < nTypes; n++)
        {
            if(name.length() < StsLogVector.types[n].length())
                clippedType = StsLogVector.types[n].substring(0,name.length());
            else
                clippedType = StsLogVector.types[n];

            if (name.equals(clippedType))
            {
                if(n == StsLogVector.MDEPTH)
                    return true;
                else
                    return false;
            }
        }
        return false;
    }

	// used for las files to capture header info
	static public class FileHeader
	{
		String formatVersion = "unknown";
        boolean lineWrap = false;
        String wellName = "unknown";
        String wellLabel = "unknown";
        String wellNumber = "unknown";
        String operator = "operator";
        String spudDate = "unknown";
        float startZ = 0.0f;
        float stopZ = 0.0f;
        float stepZ = 0.0f;
        String company = null;
        String wellIdNumber = null;
		double xOrigin = 0.0;
		double yOrigin = 0.0;
        float datumElevation = 0.0f;
        String uwiNumber = "00000000000000";
        String field = "unknown";
		float nullValue = StsParameters.nullValue;
		String[] curveNames;
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

    static private String[] readCurveNames(BufferedReader bufRdr)
    {
		String line;
        String[] curveNames = new String[0];

        try
        {
			if(!findKeyword(bufRdr, CRV_HDR)) return new String[0];
			bufRdr.mark(1024);

			line = bufRdr.readLine().trim();
			while(!lineHasKeyword(line, DATA_HDR))
			{
                if(line.substring(0) != "#")
                {
                    curveNames = (String[])StsMath.arrayAddElement(curveNames, line);
                    bufRdr.mark(1024);
                }
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
    static private boolean readCurveValues(BufferedReader bufRdr, FileHeader fileHdr, String binaryDataDir, String[] curveNames,
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
            StsException.outputException("StsWellKeywordIO.readCurveValues() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    static private boolean readAsciiCurveValues(BufferedReader bufRdr, FileHeader fileHdr, String[] curveNames,
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

		    String line;
			int nLines = 0;
            while ((line = bufRdr.readLine()) != null)
            {
                line = line.trim();
                if(line.equals("")) continue;  // blank line
                line = StsStringUtils.detabString(line);
                if(line.substring(0).equals("#")) continue; // comment line

                StringTokenizer stok = new StringTokenizer(line);
				int nTokens = stok.countTokens();

                for(int col = 0; col < nTokens; col++)
                {
                    double value = Double.parseDouble(stok.nextToken());
					int nVector = colAssignments[col];
					if(nVector == -1) continue;

                    if((float)value == fileHdr.nullValue)
                        value = logCurveNull;

                    if(nVector == depthVectorIdx)
                        curveVectors[nVector].append(value + datumShift);
                    else if(nVector == owtVectorIdx)
                        curveVectors[nVector].append(value * 2.0f);
                    else
                        curveVectors[nVector].append(value);
                }
				if(++nLines%1000 == 0)
                    StsMessageFiles.logMessage("File: " + filename + ": " + nLines + " values read.");
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
            curves =  StsLasKeywordIO.readDeviation(path, path, filename, StsParameters.nullValue, StsWellImport.WELL_DEV_PREFIX);
        }
        catch (Exception e)
        {
            StsException.outputException("StsWellKeywordIO.main: ",
                                e, StsException.WARNING);
        }
  	}
}
