
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
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;

import java.awt.*;
import java.io.*;
import java.util.*;

abstract public class StsWellIO extends StsKeywordIO
{
    public StsModel model;
    boolean deleteBinaries = false;

    StsProject project;
    String binaryDataDirectory;

    public String wellname;
    protected StsWell well;

    double xOrigin = 0.0;
    double yOrigin = 0.0;
    float nullValue = StsParameters.nullValue;
    float datumShift = 0.0f;
    String[] curveNames;

    String formatVersion = "unknown";
    boolean lineWrap = false;
    String wellName = "unknown";
    String wellLabel = "unknown";
    String wellNumber = "unknown";
    String operator = "operator";
    String date = "unknown";
    float startZ = 0.0f;
    float stopZ = 0.0f;
    float stepZ = 0.0f;
    String company = null;
    String wellId = null;
    float elevation = 0.0f;
    String apiNumber = "00000000000000";
    String field = "unknown";

    public float vScalar = 1.0f;
    public float hScalar = 1.0f;
    public byte verticalUnits = StsParameters.DIST_NONE;
    public byte horizontalUnits = StsParameters.DIST_NONE;

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
    static public final String[] Z_KEYWORDS = { DEPTH, MDEPTH, SSDEPTH, TIME };

    abstract public long getWellFileDate();
    abstract public StsLogCurve readTdCurve(float logCurveNull, byte vUnits);
    abstract public StsLogVector[] readDeviationVectors(float logCurveNull, byte vUnits, float datumShift);
    abstract public StsLogCurve[] readLogCurves(float logCurveNull, byte vUnits);
    abstract public void addWellMarkers(StsWell well);

    public StsWellIO(StsModel model, boolean deleteBinaries)
	{
        this.model = model;
        this.deleteBinaries = deleteBinaries;

        this.project = model.getProject();
        binaryDataDirectory = project.getBinaryDirString();
    }

    public StsWell createWell(float logNull, float datumShift, StsProgressPanel progressPanel)
    {
        StsLogVector[] xyzmtVectors;
        StsLogCurve[] logCurves;
        StsLogCurve tdCurve;

        try
        {
            nullValue = logNull;
            this.datumShift = datumShift;

            tdCurve = readTdCurve(logNull, verticalUnits);

            xyzmtVectors = readDeviationVectors(logNull, verticalUnits, datumShift);
            if (xyzmtVectors == null)
            {
            	progressPanel.appendLine("Failed to read value vectors (" + wellname + "), check file header.");
            	progressPanel.setDescriptionAndLevel("Failed to load " + wellname + ".", StsProgressBar.ERROR);
                return null;
            }

            // Build the Well
            well = buildWell(xyzmtVectors, tdCurve, progressPanel);

            if (well == null)
            {
            	progressPanel.appendLine("Failed to build well from supplied information (" + wellname + "), check file contents.");
            	progressPanel.setDescriptionAndLevel("Failed to build" + wellname + ".", StsProgressBar.ERROR);
                return null;
            }
            // Add the TD Curve
            if (tdCurve != null)
                well.addLogCurve(tdCurve);

            // Set units to the ones from the binaries, whether just created or old.
            byte vUnits = StsLogVector.getVectorOfType(xyzmtVectors, StsLogVector.DEPTH).getUnits();
            byte hUnits = StsLogVector.getVectorOfType(xyzmtVectors, StsLogVector.X).getUnits();

            logCurves = readLogCurves(logNull, vUnits);
            if (logCurves == null || logCurves.length == 0)
                return well;
            well.addLogCurves(logCurves);

            if (progressPanel != null)
            {
                progressPanel.appendLine("Processing well: " + wellname + "...");
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsWellWizard failed.", e, StsException.WARNING);
            StsMessageFiles.logMessage("Unable to read well deviation files.");
            progressPanel.setDescriptionAndLevel("Unable to read well deviation files.\n", StsProgressBar.ERROR);

            return null;
        }
        return well;
    }

  /** build the well entity with its well line displayed in the 3d window */
    private StsWell buildWell(StsLogVector[] xyzmtVectors, StsLogCurve tdCurve, StsProgressPanel progressPanel)
    {
        StsWell well = new StsWell(wellname, true);

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
            addHeaderInfo(well);
            return well;
        }
        catch (Exception e)
        {
            StsException.outputException("StsWellWizard.buildWell() failed.", e, StsException.WARNING);
            progressPanel.setDescriptionAndLevel("StsWellWizard.buildWell() failed.\n", StsProgressBar.ERROR);
            return null;
        }
    }

    /** override this if well data has additional info to add. */
    protected void addHeaderInfo(StsWell well)
    {
    }
	/** check that binary file dates are all newer than ascii file date */
	 protected boolean binaryFileDatesOK(String[] curveNames, String filePrefix)
	 {
		 File file;

		 try
		 {
             long sourceFileDate = getWellFileDate();
			 boolean binaryDatesOK = true;
			 int nNames = curveNames.length;
			 for(int n = 0; n < nNames; n++)
			 {
				 String curveName = curveNames[n];
				 if(curveName.equals("") || curveName.equalsIgnoreCase("ignore")) continue;
				 String binaryFilename = filePrefix + ".bin." + wellname + "." + curveName + "." + version;
				 file = new File(binaryDataDirectory + File.separator + binaryFilename);
				 long binaryFileDate = file.lastModified();
				 if(binaryFileDate < sourceFileDate) binaryDatesOK = false;
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
    protected boolean checkWriteBinaryFiles(StsLogVector[] vectors)
    {
        int nVectors = vectors.length;
        boolean writeOK = true;
        for(int n = 0; n < nVectors; n++)
            if(!vectors[n].checkWriteBinaryFile(binaryDataDirectory)) writeOK = false;

        return writeOK;
    }

     /** check that binary files exist */
      public boolean binaryFileExist(String name, String filePrefix, byte type)
      {
          File file;

          try
          {
              String binaryFilename = filePrefix + ".bin." + name + "." + StsLogVector.getStringFromType(type) + "." + version;
              file = new File(binaryDataDirectory + File.separator + binaryFilename);
              return file.exists();
          }
          catch (Exception e)
          {
              StsException.outputException("StsWellKeywordIO.logVectorDatesOK() failed.",
                  e, StsException.WARNING);
              return false;
          }
      }

	  protected void deleteBinaryFiles(String[] curveNames, String filePrefix)
      {
		  File file;

		  try
		  {
			  int nNames = curveNames.length;
			  for(int n = 0; n < nNames; n++)
			  {
				  String curveName = curveNames[n];
				  if(curveName.equals("") || curveName.equalsIgnoreCase("ignore")) continue;
				  String binaryFilename = filePrefix + ".bin." + wellname + "." + curveName + "." + version;
				  file = new File(binaryDataDirectory + File.separator + binaryFilename);
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
    protected StsLogVector[] constructLogVectors(String[] curveNames, String filePrefix)
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

    /* read the well name */
    protected boolean readFileHeader(BufferedReader bufRdr)
    {
        try
        {
			while(true)
			{
				String line = bufRdr.readLine().trim();
                line = StsStringUtils.detabString(line);
				if(line.endsWith(WELL_NAME))
					wellname =  new String(bufRdr.readLine().trim());
				else if(line.indexOf(ORIGIN) >= 0)  // is origin keyword there?
				{
					boolean xyOrder = true;
					if (line.indexOf(YX) >= 0) xyOrder = false;  // determine x-y order
					line = bufRdr.readLine().trim();  // get the next line

					// tokenize the x-y values and convert to a point object
					StringTokenizer stok = new StringTokenizer(line);
					xOrigin = Double.valueOf(stok.nextToken()).doubleValue();
					yOrigin = Double.valueOf(stok.nextToken()).doubleValue();

					if(!xyOrder)
					{
						double temp = xOrigin;
						xOrigin = yOrigin;
						yOrigin = temp;
					}
				}
				else if(line.endsWith(NULL_VALUE))
				{
					line = bufRdr.readLine().trim();  // get the next line
					StringTokenizer stok = new StringTokenizer(line);
					nullValue = Float.valueOf(stok.nextToken()).floatValue();
				}
				else if(line.endsWith(CURVE))
				{
					String[] curveNames = new String[0];
					line = bufRdr.readLine().trim();
					while(!lineHasKeyword(line, VALUE))
					{
						curveNames = (String[])StsMath.arrayAddElement(curveNames, line);
						line = bufRdr.readLine().trim();
					}
					return true;
				}
				else if(line.endsWith(VALUE))
					return false;
				else
					return false;
			}
        }
        catch(Exception e)
        {
            StsException.outputException("StsWellKeywordIO.readFileHeader() failed.",
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
}
