package com.Sts.Actions.Wizards.OSWell;

import com.Sts.Actions.Import.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jun 20, 2008
 * Time: 10:24:41 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsOSWellIO extends StsWellIO
{
    StsOSWell osWell;

    public StsOSWellIO(StsModel model, StsOSWell osWell, boolean deleteBinaries)
    {
        super(model, deleteBinaries);
        this.osWell = osWell;
        this.wellname = osWell.getName();
    }

    public long getWellFileDate()
    {
        return 0L;
    }

    public StsLogCurve readTdCurve(float logCurveNull, byte vUnits)
    {
    	// TD curves from OpenSpirit do not contain null values
    	StsLogCurve tdCurve = osWell.getTdCurve();
        if(tdCurve != null)
        {
            StsWellKeywordIO.checkWriteBinaryFiles(new StsLogVector[] {tdCurve.getDepthVector()}, model.getProject().getBinaryFullDirString());
            StsWellKeywordIO.checkWriteBinaryFiles(new StsLogVector[] {tdCurve.getMDepthVector()}, model.getProject().getBinaryFullDirString());
            StsWellKeywordIO.checkWriteBinaryFiles(new StsLogVector[] {tdCurve.getTimeVector()}, model.getProject().getBinaryFullDirString());
            StsWellKeywordIO.checkWriteBinaryFiles(new StsLogVector[] {tdCurve.getValueVector()}, model.getProject().getBinaryFullDirString());
        }
        return tdCurve;
    }

    public StsWell createWell(float logNull, float datumShift, StsProgressPanel progressPanel)
    {
        StsWell well = super.createWell(logNull, datumShift, progressPanel);
        if(well != null)
        {
            well.setName(osWell.getName());
            well.setWellLabel(osWell.getName());
            well.setApi(osWell.getApi());
            well.setUwi(osWell.getUwi());
            well.setOperator(osWell.getOperator());
            well.setField(osWell.getField());
            well.setDate(osWell.getDate());
            well.setKbElev(osWell.getKbElev());
            well.setElev(osWell.getElev());
            well.setElevDatum(osWell.getElevDatum());
            well.setSpudDate(osWell.getSpudDate());
            well.setCompletionDate(osWell.getCompletionDate());
            well.setXOrigin(osWell.getxOrigin());
            well.setYOrigin(osWell.getyOrigin());
        }
        return well;
    }

    public StsLogVector[] readDeviationVectors(float logCurveNull, byte vUnits, float datumShift)
    {
    	StsLogVector[] vectors = null;
    	int numVec = 5;
    	if (osWell.getDevTWT() == null)
    		numVec = 4;

    	String[] names = new String[numVec];
        names[0] = StsLogVector.getStringFromType(StsLogVector.X);
        names[1] = StsLogVector.getStringFromType(StsLogVector.Y);
        names[2] = StsLogVector.getStringFromType(StsLogVector.MDEPTH);
        names[3] = StsLogVector.getStringFromType(StsLogVector.DEPTH);
        if (osWell.getDevTWT() != null)
            names[4] = StsLogVector.getStringFromType(StsLogVector.TWT);

        vectors = StsWellKeywordIO.constructLogVectors(names, wellname, StsLogVector.WELL_DEV_PREFIX);

        vectors[0].setNullValue(logCurveNull);
        vectors[0].setValues(osWell.getDevXOffset());
        vectors[1].setNullValue(logCurveNull);
        vectors[1].setValues(osWell.getDevYOffset());
        vectors[2].setNullValue(logCurveNull);
        vectors[2].setValues(osWell.getDevMD());
        vectors[3].setNullValue(logCurveNull);
        vectors[3].setValues(osWell.getDevTVD());
        if (osWell.getDevTWT() != null)
        {
            vectors[4].setNullValue(logCurveNull);
            vectors[4].setValues(osWell.getDevTWT());
        }
        StsWellKeywordIO.checkWriteBinaryFiles(vectors, model.getProject().getBinaryFullDirString());
        return vectors;
    }

    public StsLogCurve[] readLogCurves(float logCurveNull, byte vUnits)
    {
    	// null values are already set in each log
        StsLogCurve[] curves = osWell.getLogs();
        if(curves != null)
        {
            StsLogVector[] vectors = new StsLogVector[curves.length];
            StsWellKeywordIO.checkWriteBinaryFiles(new StsLogVector[] {curves[0].getDepthVector()}, model.getProject().getBinaryFullDirString());
            StsWellKeywordIO.checkWriteBinaryFiles(new StsLogVector[] {curves[0].getMDepthVector()}, model.getProject().getBinaryFullDirString());
            for(int i=0; i<curves.length; i++)
            {
                vectors[i] = curves[i].getValueVector();
            }
            StsWellKeywordIO.checkWriteBinaryFiles(vectors, model.getProject().getBinaryFullDirString());
        }
        return curves;
    }

    public void addWellMarkers(StsWell well)
    {
        StsOSWellMarker[] osWellMarkers = osWell.getMarkers();
        if(osWellMarkers == null) // No markers found
            return;
        int nMarkers = osWellMarkers.length;
        if(nMarkers == 0) return;
        StsWellMarker[] wellMarkers = new StsWellMarker[nMarkers];

        for(int n = 0; n < nMarkers; n++)
        {
            StsPoint loc = well.getPointAtMDepth(osWellMarkers[n].getMdepth(), true);
            if(loc == null)
            {
                StsMessageFiles.errorMessage("Unable to locate marker (" + osWellMarkers[n].getName() + ") at md=" + osWellMarkers[n].getMdepth() + " for well " + well.getName());
                continue;
            }
            wellMarkers[n] = StsWellMarker.constructor(osWellMarkers[n].getName(), well, StsMarker.GENERAL, loc);
        }
    }
}
