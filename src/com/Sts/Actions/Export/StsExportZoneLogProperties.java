
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Export;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.RFUI.*;
import com.Sts.Types.*;
import com.Sts.UI.*;

import java.awt.*;
import java.io.*;

public class StsExportZoneLogProperties extends StsAction
{
    static final protected String NL = System.getProperty("line.separator");
    static final protected int CONTINUOUS_TYPE = 0;
    static final protected int CATEGORICAL_TYPE = 1;

    StsZoneClass zoneClass;
    protected int nZones;
    protected StsZone[] zones;

    protected StsGridDefinition gridDefinition;
//    protected String logName;
    protected int propertyType = CONTINUOUS_TYPE;
    protected BufferedWriter out = null;

 	public StsExportZoneLogProperties(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public boolean start()
    {
        statusArea.setTitle("Export Continuous Properties:");
        zoneClass = (StsZoneClass)model.getCreateStsClass(StsZone.class);
        if(!zoneClass.checkBuildWellZones())
        {
            statusArea.textOnly();
            return false;
        }
        String logName = userQuery();
        if (logName == null)
        {
            statusArea.textOnly();
            return false;
        }

        File f = buildExportFile(logName);
        if( f != null )
        {
            logMessage("Export file " + f.getPath() + " created successfully.");
			actionManager.endCurrentAction();
            return true;
        }

        return false;
    }

	public File export(String logName)
    {
        if (retrieveModelObjects())
        {
            return buildExportFile(logName);
//            if( f != null ) return calculatePropertyVolume(f, logName);
        }
        return null;
    }


    public boolean calculatePropertyVolume(File scatterDataFile, String logName)
    {
        try
        {
        /*
            CFunction createDGrid = new CFunction("S2SStatLib.dll", "createDGrid");
            CPointer dg = createDGrid.callCPointer(new Object[] {scatterDataFile.getPath()});
            if( dg != null )
            {
                CFunction createInverseDistance = new CFunction("S2SStatLib.dll", "createInverseDistance");
                Float power = new Float(2.0f);
                CPointer id = createInverseDistance.callCPointer(new Object[] { dg, power });


                String path = ".";
                try { path = model.getProject().getSts_3DModelsFullDirString(); }
                catch (Exception e) { }
                String filename = new String("Mod3DT1_" +
                                            model.getGridDefinition().getNCols() + "_" +
                                            model.getGridDefinition().getNRows() + "_" +
                                            model.getNLayers() + "_" +
                                            logName);
                File propertyVolumeFile = new File(path, filename);
                CFunction createPropertyVolume = new CFunction("S2SStatLib.dll", "createPropertyVolume");
                createPropertyVolume.callVoid(new Object[] { dg, id, propertyVolumeFile.getPath() });

                CFunction deleteInverseDistance = new CFunction("S2SStatLib.dll", "deleteInverseDistance");
                deleteInverseDistance.callVoid(new Object[] { id });
                CFunction deleteDGrid = new CFunction("S2SStatLib.dll", "deleteDGrid");
                deleteDGrid.callVoid(new Object[] { dg });
            }
        */
        }
        catch(Exception e) { e.printStackTrace(); }
        return true;
    }


    public boolean end()
    {
        statusArea.textOnly();
        return true;
    }

    protected boolean retrieveModelObjects()
    {
        nZones = zoneClass.getSize();
        if(nZones == 0)
        {
            StsMessageFiles.logMessage("No zones found.");
 			actionManager.endCurrentAction();
            return false;
        }

        // check to see if we have well zone sets for all zones
        for (int i=0; i<nZones; i++)
        {
            StsZone zone = (StsZone)zoneClass.getElement(i);
            if (zone == null)
            {
                StsMessageFiles.logMessage("Null zone detected! Terminating export ...");
                return false;
            }
            if (zone.getWellZoneSet() == null)
            {
                StsMessageFiles.logMessage("Zone: " + zone.getName() + " has no corresponding well zone set. Terminating export...");
                return false;
            }
        }

        // build array of zones for convenience
        zones = new StsZone[nZones];
        for (int i=0; i<nZones; i++) zones[i] = (StsZone)zoneClass.getElement(i);

        // Get the modelGrid geometry
        gridDefinition = model.getGridDefinition();
        if(gridDefinition == null)
        {
            logMessage("Must define horizons grids first. Terminating export ...");
            return false;
        }

        return true;
    }

    protected String userQuery()
    {
        return askForLogName("Select a continuous log curve: ", "No continuous log curve selected.");
    }

    protected String askForLogName(String title, String notFoundText)
    {
        String logName = null;
        try
        {
            // ask for a log curve name
            StsSelectStsObjects selector = StsSelectStsObjects.constructor(model,
                    StsLogCurveType.class, model.getName(), title, true);
			if(selector == null) return null;
            logName = selector.selectName();
            if (logName==null) logMessage(notFoundText);
        }
        catch(Exception e)
        {
            logMessage("Unable to select a log curve.");
        }
        return logName;
    }

    protected void writeHeader(String logName) throws IOException
    {
        out.write("# property name and type (0=continuous, 1=categorical)" + NL);
        out.write(logName + "  " + propertyType + NL);
        out.write("# Max#Wells Max#Zones" + NL);
        int nWells = model.getNObjects(StsWell.class);
//        StsWellClass wellClass = (StsWellClass)model.getCreateStsClass(StsWell.class);
//        int nWells = wellClass.getNWells(StsParameters.WELL);
        out.write(nWells + "  " + nZones + NL);
    }


    protected File buildExportFile(String logName)
    {
        StsCursor cursor = null;
        String filename = null;
        File f = null;
        try
        {
            // set the wait cursor
            cursor = new StsCursor(model.win3d, Cursor.WAIT_CURSOR);

            // analyze subzones
            int totalSubZones = 0;
            int maxNSubZones = 0;
            for (int i=0; i<nZones; i++)
            {
                int nSubZones = zones[i].getNLayers();
                totalSubZones += nSubZones;
                if (nSubZones>maxNSubZones) maxNSubZones = nSubZones;
            }
            if (totalSubZones==0)
            {
                logMessage("No sub-zones found!  "
                        + "Terminating export ...");
                cursor.restoreCursor();
                return null;
            }

            // open ascii export file
            filename = "DB." + model.getName() + ".P." + logName + ".scat";
            f = new File(model.getProject().getDataFullDirString(), filename);
            if (f.exists() && !f.canWrite())
            {
                logMessage("Unable to write to existing file "
                        + filename + ".  Terminating export ...");
                cursor.restoreCursor();
                return null;
            }
            try
            {
                out = new BufferedWriter(new FileWriter(f));
            }
            catch (Exception e)
            {
                logMessage("Unable to open file "
                        + filename + ".  Terminating export ...");
                cursor.restoreCursor();
                return null;
            }

            // create export file header
            writeHeader(logName);

            // get xy grid extent info
            float xMin = gridDefinition.getXMin();
            float yMin = gridDefinition.getYMin();
            float xSize = gridDefinition.getXSize();
            float ySize = gridDefinition.getYSize();
            float xMax = xMin + xSize;
            float yMax = yMin + ySize;
            float xInc = gridDefinition.getXInc();
            float yInc = gridDefinition.getYInc();
            final double absoluteXOrigin = model.getXOrigin();
            final double absoluteYOrigin = model.getYOrigin();
            xMin += absoluteXOrigin;  // change from relative to absolute X
            yMin += absoluteYOrigin;  // change from relative to absolute Y
            float angle = 0.0f;
            out.write("# Total#Layers XOrigin Y0rigin XWidth YHeight dX dY Angle" + NL);
            out.write(totalSubZones + " " + xMin + " " + yMin + " " + xSize + " " +
                      ySize + " " + xInc + " " + yInc + " " + angle + NL);

            // output zone names
            out.write("# Zone# ZoneName" + NL);
            for (int i=0; i<nZones; i++)
            {
                out.write((i+1) + "    " + zones[i].getName() + "   " + zones[i].getNSubZones() + NL);
            }

            // calculate average thickness for subZones
            out.write("# Zone# Layer# AvgLayerThickness" + NL);
            for (int i=0; i<nZones; i++)  // zones are reordered top->bottom
            {
                int nSubZones = zones[i].getNLayers();
                StsObjectRefList wellZoneList = zones[i].getWellZoneSet().getWellZones();
                if (wellZoneList == null) continue;
                int nWellZones = wellZoneList.getSize();

                // go thru well zones & compute subZone thicknesses
                int nValues = 0;
                float[] avgThicknesses = new float[nSubZones];
                for (int j=0; j<nWellZones; j++)
                {
                    StsWellZone wellZone = (StsWellZone)wellZoneList.getElement(j);
                    float[] subZoneThicknesses = wellZone.getLayerThicknesses();
                    if (subZoneThicknesses==null) continue;

                    StsLogCurve logCurve = wellZone.getWell().getLastLogCurveOfType(logName);
                    if (logCurve==null) continue;

                    for (int k=0; k<nSubZones; k++)
                    {
                        avgThicknesses[k] += subZoneThicknesses[k];
                    }
                    nValues++;
                }
                if (nValues > 0)
                {
                    for (int k=0; k<nSubZones; k++)
                    {
                        avgThicknesses[k] /= nValues;
                        out.write((i+1) + "    " + (k+1) + "    "
                                + avgThicknesses[k] + NL);
                    }
                }
            }

            // write out subzone averaged values
            out.write("# Well# GlobalLayer# Zone# Layer# X Y Value dZ" + NL);
            int nLayersAbove = 0;
            int nSubZones = 0;
            for (int i=0; i<nZones; i++)
            {
                if (zones[i].getWellZoneSet() == null) continue;
                StsObjectRefList wellZoneList = zones[i].getWellZoneSet().getWellZones();

                nLayersAbove += nSubZones;
                nSubZones = zones[i].getNLayers();
                int nWellZones = wellZoneList.getSize();
                for (int j=0; j<nWellZones; j++)
                {
                    StsWellZone wellZone = (StsWellZone)wellZoneList.getElement(j);
                    StsWell well = wellZone.getWell();
                    if (well == null) continue;

                    StsLogCurve logCurve = well.getLastLogCurveOfType(logName);
                    if (logCurve==null)
                    {
                        logMessage("Unable to access log " +
                                logName + " in well " + well.getName() +
                                ".  Continuing...");
                        continue;
                    }
                    float[] subZoneLogAverages = getLayerValues(logCurve, wellZone);
                    if (subZoneLogAverages==null)
                    {
                        logMessage("Unable to compute layer "
                                + "average for log " + logName + " in well "
                                + well.getName() + ".  Continuing...");
                        continue;
                    }

                    float[] subZoneThicknesses = wellZone.getLayerThicknesses();
                    if (subZoneThicknesses==null)
                    {
                        logMessage("Unable to compute layer "
                                + " thicknesses in well "
                                + well.getName() + ".  Continuing...");
                        continue;
                    }
                    float[] subZoneCenterZs = wellZone.getLayerCenterZs();
                    if (subZoneCenterZs==null)
                    {
                        logMessage("Unable to compute layer "
                                + " center Zs in well "
                                + well.getName() + ".  Continuing...");
                        continue;
                    }

                    // write subzone properties
                    for (int k=0; k<nSubZones; k++)
                    {
                        StsPoint p = well.getPointAtZ(subZoneCenterZs[k], true);
                        // iw, iwz, i, k, x, y, val, dzs
                        out.write(well.getIndex()+1 + "    " + (nLayersAbove+k+1)
                                + "    " + (i+1) + "    " + (k+1)
                                + "    " + (p.getX()+absoluteXOrigin)
                                + "    " + (p.getY()+absoluteYOrigin)
                                + "    " + subZoneLogAverages[k]
                                + "    " + subZoneThicknesses[k] + NL);
                    }
                }
            }
        }
        catch (Exception e)
        {
            closeFile(out);
            System.out.println("Exception in StsExportZoneLogProperties.start()" + NL + e);
            logMessage("Export Zone Log Properties failed!");
            if (cursor!=null) cursor.restoreCursor();
            return null;
        }

        cursor.restoreCursor();
        if (closeFile(out))
        {
            return f;
        }
        else return null;
    }

    private boolean closeFile(BufferedWriter out)
    {
        if (out!=null)
        {
            try
            {
                out.flush();
                out.close();
            }
            catch (IOException e)
            {
                logMessage("Error!  Unable to Close export file!");
                return false;
            }
        }
        return true;
    }

    protected float[] getLayerValues(StsLogCurve logCurve, StsWellZone wellZone)
    {
        float[] values = wellZone.getLayerAverages(logCurve);
        if (values == null) return null;

        // fix nulls
    /*
        for (int i=0; i<values.length; i++)
        {
            if (values[i] == StsParameters.nullValue)
            {
                values[i] = StsParameters.HDFnullValue;
            }
        }
    */
        return values;
    }

    protected float[] getSubZoneValues(StsLogCurve logCurve, StsWellZone wellZone)
    {
        float[] values = wellZone.getSubZoneAverages(logCurve);
        if (values == null) return null;

        // fix nulls
    /*
        for (int i=0; i<values.length; i++)
        {
            if (values[i] == StsParameters.nullValue)
            {
                values[i] = StsParameters.HDFnullValue;
            }
        }
    */
        return values;
    }

}
