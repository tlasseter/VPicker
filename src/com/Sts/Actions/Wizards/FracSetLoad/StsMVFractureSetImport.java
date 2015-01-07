//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.FracSetMVLoad;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;

import java.text.*;

public class StsMVFractureSetImport
{
    static StsModel model;
    static String currentDirectory = null;
    static boolean success = false;

    static private boolean reloadAscii = false;

    static byte vUnits = StsParameters.DIST_NONE;
    static byte hUnits = StsParameters.DIST_NONE;
 
    static public final boolean debug = false;
    static StsTimer timer = new StsTimer();

    //static double modelXOrigin, modelYOrigin;

    static public void initialize(StsModel model_)
    {
        model = model_;
        StsProject project = model.getProject();
        currentDirectory = project.getRootDirString();
        hUnits = project.getXyUnits();
        vUnits = project.getDepthUnits();
    }

    static public String getCurrentDirectory()
    {
        return currentDirectory;
    }

    static public void setCurrentDirectory(String dirPath)
    {
        currentDirectory = dirPath;
    }

    /** Create frac sets */
    static public StsMVFractureSet createFractureSet(StsModel model, StsProgressPanel progressPanel, StsMVFracSetFile selectedFile)
    {
    	initialize(model);
        StsMVFractureSet fracSet = null;

        try
        {
            if ((model.getObjectWithName(StsMVFractureSet.class, selectedFile.name) == null) || reloadAscii)
            {
                StsMVFractureSet temp = (StsMVFractureSet) model.getObjectWithName(StsMVFractureSet.class, selectedFile.name);
                if ( (temp != null) && reloadAscii)
                {
                    temp.delete();
                }

                 fracSet = createFractureSetFromAsciiFile(model.getProject().getRotatedBoundingBox(), currentDirectory,
                             selectedFile, progressPanel, model.getProject().getTimeDateFormat());
                	
            }
            else
            {
                if (progressPanel != null)
                {
                    progressPanel.appendLine("FracSet: " + selectedFile.name + " already loaded...\n");
                }
            }
            return fracSet;

        }
        catch (Exception e)
        {
            StsException.outputException("StsSensorImport.createSensors() failed.", e, StsException.WARNING);
            progressPanel.setDescriptionAndLevel("StsSensorImport.createWells() failed.\n", StsProgressBar.WARNING);
            return null;
        }
    }

    /** read in multiple well deviations from a list of Ascii files */
    static public StsMVFractureSet createFractureSetFromAsciiFile(StsRotatedBoundingBox rBox, String dataDir,
       StsMVFracSetFile fracSetFile, StsProgressPanel progressPanel, SimpleDateFormat dateFormat)
    {
        StsMVFractureSet fracSet = StsMVFractureSet.constructor();

        String fracSetName = fracSetFile.name;
        fracSet.setName(fracSetName);

        if (dataDir == null)
        {
            return null;
        }

        try
        {
            //boolean deleteBinaries = reloadAscii;
            fracSetFile.analyzeFile();
            double easting = StsParameters.nullValue;
            double northing = StsParameters.nullValue;
            float depth = StsParameters.nullValue;
            float dip = StsParameters.nullValue;
            float azimuth = StsParameters.nullValue;
            float strike = StsParameters.nullValue;
            float aperture = StsParameters.nullValue;
            float aspectRatio = StsParameters.nullValue;
            float length = StsParameters.nullValue;

            String[] attributeNames = null;
            float[] attributeVals = null;
            if(fracSetFile.numAttributes > 0)
            {
                attributeNames = new String[fracSetFile.numAttributes];
                attributeVals = new float[fracSetFile.numAttributes];
            }
            while (fracSetFile.getPropertyValues(model, model.getProject().getDateOrder()))
            {
            	int numProps = fracSetFile.validCurves.length;
                int cnt = 0;
            	for (int i=0; i<numProps; i++)
            	{           		
            		String unit = fracSetFile.currentUnits[i];  
            		if ( unit.equalsIgnoreCase("Metre") || unit.equalsIgnoreCase("Meter")) unit = "Meters";
            		float hScalar = model.getProject().getXyScalar(unit);
            		float vScalar = model.getProject().getDepthScalar(unit);
            		boolean isMetric = (unit.equalsIgnoreCase("Meter") || unit.equalsIgnoreCase("Metre"));
            		boolean isEnglish = unit.equalsIgnoreCase("Feet");
            		if ( fracSetFile.validCurves[i].equalsIgnoreCase("East") ||
            				fracSetFile.validCurves[i].equalsIgnoreCase("X"))
            		{
            			easting = (float)fracSetFile.currentValues[i];
            			easting *= hScalar;
//            			if (hUnits == StsParameters.DIST_FEET && isMetric )
//            				easting *= StsParameters.DIST_FEET_SCALE;
//            			if (hUnits == StsParameters.DIST_METER && isEnglish )
//            				easting /= StsParameters.DIST_FEET_SCALE;
            		}
            		else if ( fracSetFile.validCurves[i].equalsIgnoreCase("North") ||
            				fracSetFile.validCurves[i].equalsIgnoreCase("Y"))
            		{
            			northing = (float)fracSetFile.currentValues[i];
            			northing *= hScalar;
//            			if (hUnits == StsParameters.DIST_FEET && isMetric )
//            				northing *= StsParameters.DIST_FEET_SCALE;
//            			if (hUnits == StsParameters.DIST_METER && isEnglish )
//            				northing /= StsParameters.DIST_FEET_SCALE;
            		}
            		else if ( fracSetFile.validCurves[i].equalsIgnoreCase("Z") ||
            				fracSetFile.validCurves[i].equalsIgnoreCase("depth"))
            		{
                        depth = (float)fracSetFile.currentValues[i];
                        boolean isElevation = !fracSetFile.validCurves[i].equalsIgnoreCase("depth");
                        if(isElevation) depth *= -1.0;
                        depth *= vScalar;
//                        if (vUnits == StsParameters.DIST_FEET && isMetric )
//            				depth *= StsParameters.DIST_FEET_SCALE;
//            			if (vUnits == StsParameters.DIST_METER && isEnglish )
//            				depth /= StsParameters.DIST_FEET_SCALE;            			
                    }
            		else if ( fracSetFile.validCurves[i].equalsIgnoreCase("dip"))
            		{
            			dip = (float)fracSetFile.currentValues[i];
            		}
            		else if ( fracSetFile.validCurves[i].equalsIgnoreCase("azimuth"))
            		{
            			azimuth = (float)fracSetFile.currentValues[i];
            		}
            		else if ( fracSetFile.validCurves[i].equalsIgnoreCase("strike"))
            		{
            			strike = (float)fracSetFile.currentValues[i];
            		}
            		else if ( fracSetFile.validCurves[i].equalsIgnoreCase("aperture"))
            		{
            			aperture = (float)fracSetFile.currentValues[i];
            		}
            		else if ( fracSetFile.validCurves[i].equalsIgnoreCase("aspect_ratio"))
            		{
            			aspectRatio = (float)fracSetFile.currentValues[i];
            		}
            		else if (fracSetFile.validCurves[i].equalsIgnoreCase("length"))
            		{
            			length = (float)fracSetFile.currentValues[i];
            			length *= vScalar;
//            			if (vUnits == StsParameters.DIST_FEET && isMetric )
//            				length *= StsParameters.DIST_FEET_SCALE;
//            			if (vUnits == StsParameters.DIST_METER && isEnglish )
//            				length /= StsParameters.DIST_FEET_SCALE;
            		}
                    else
                    {
                        attributeNames[cnt] = fracSetFile.validCurves[i];
                        attributeVals[cnt] = (float)fracSetFile.currentValues[i];
                        cnt++;
                        // Other attributes.
                    }
            	}
                if(attributeNames != null)
                {
                    attributeNames = (String[])StsMath.trimArray(attributeNames, cnt);
                    attributeVals = (float[])StsMath.trimArray(attributeVals, cnt);
                }

                StsProject project = model.getProject();
                project.checkSetOrigin(easting, northing);
                float[] xy = project.getRotatedRelativeXYFromUnrotatedAbsoluteXY(easting, northing);
                float[] center = new float[] {xy[0], xy[1], depth};
                azimuth -= project.getAngle();
                StsMVFracture fracture = new StsMVFracture(fracSet, center, azimuth, dip, length, aspectRatio,
                                            aperture, fracSetFile.currentTime, attributeVals, attributeNames);
                fracSet.addFracture(fracture);
            }

            if (progressPanel != null)
            {
                progressPanel.appendLine("Processing ASCII Frac Set: " + fracSetName + "...");
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsSensorImport failed.", e, StsException.WARNING);
            StsMessageFiles.logMessage("Unable to read Ascii sensor files.");
            progressPanel.appendLine("Unable to read Ascii sensor files.");

            return null;
        }
        return fracSet;
    }

 }