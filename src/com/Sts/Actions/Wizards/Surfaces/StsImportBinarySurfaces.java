
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Surfaces;

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.net.*;

public class StsImportBinarySurfaces extends StsImportSurfaces
{
    /** read a generic binary surface file, create the map point array and
        build a new surface; return true if created, false otherwise */
    static protected StsSurface createSurface(StsModel model, StsAbstractFile file, String name, StsProgressPanel panel)
    {
        StsSurface newSurface;
        URLConnection urlConnection;

        if(file == null) return null;

      	StsMessageFiles.logMessage("Loading grid from: " + file.getFilename() + " ...");

		StsColor color = model.getCurrentSpectrumColor("Basic");
        newSurface = constructSurface(model, file, name, color, panel);
        if(newSurface == null) return null;

        StsMessageFiles.logMessage("Adding surface to model ...");
//        newSurface.persistPoints();
        model.incrementSpectrumColor("Basic");

        // Translate Surface to Model Units
        newSurface.toProjectUnits();

        if(!model.getProject().addToProject(newSurface))
        {
            newSurface.delete();
            new StsMessage(model.win3d, StsMessage.WARNING, "Failed to load surface: " + name);
            return null;
        }

        if(debug) System.out.println("Loaded binary grid from file: " + file.getFilename());
        StsMessageFiles.logMessage("Loaded grid: " + name);
        return newSurface;
    }

    static private StsSurface constructSurface(StsModel model, StsAbstractFile file, String name, StsColor stsColor, StsProgressPanel panel)
    {
        try
        {
            StsMessageFiles.logMessage("Reading grid coordinates from binary file...");
            return readBinaryFile(model, file, name, stsColor, panel);
        }
        catch(Exception e)
        {
            StsException.outputException("StsImportBinarySurface.constructSurface() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

	static private StsSurface readBinaryFile(StsModel model, StsAbstractFile file, String name, StsColor stsColor, StsProgressPanel progressPanel)
	{
		StsBinaryFile binaryFile = null;
        float xInc, yInc;
        int nCols, nRows, nPnts;
        float rowNumMin = StsParameters.nullValue;
        float rowNumMax = StsParameters.nullValue;
        float rowNumInc = StsParameters.nullValue;
        float colNumMin = StsParameters.nullValue;
        float colNumMax = StsParameters.nullValue;
        float colNumInc = StsParameters.nullValue;
        float angle;
        double xOrigin, yOrigin;
        float zMin, zMax;
        float nullZValue;
        boolean hasNulls;
        byte zDomain;
        float[][] pointsZ;

		try
		{
            binaryFile = new StsBinaryFile(file);
            binaryFile.openRead();
            byte[] units = binaryFile.getByteValues();
            int[] ints = binaryFile.getIntegerValues();
            nRows = ints[0]; nCols = ints[1];
//            System.out.println("rows: " + nRows + " ncols: " + nCols);
            double[] doubles = binaryFile.getDoubleValues();
            xOrigin = doubles[0];
            yOrigin = doubles[1];
//            System.out.println("xOrigin: " + xOrigin + " yOrigin: " + yOrigin);

            float[] floats;
            floats = binaryFile.getFloatValues();
            xInc = floats[0];
            yInc = floats[1];
            angle = floats[2];
            zMin = floats[3];
            zMax = floats[4];
            nullZValue = floats[5];
            if(floats.length > 6)
            {
                rowNumMin = floats[6];
                rowNumMax = floats[7];
                rowNumInc = floats[8];
                colNumMin = floats[9];
                colNumMax = floats[10];
                colNumInc = floats[11];
            }
            boolean[] booleans = binaryFile.getBooleanValues();
            hasNulls = booleans[0];
            byte[] bytes = binaryFile.getByteValues();
            zDomain = bytes[0];

//            System.out.println("StsGrid.readBinaryGrid() header initialized");

            pointsZ = new float[nRows][];
            for(int row = 0; row < nRows; row++)
                pointsZ[row] = binaryFile.getFloatValues();

//            xMax = xMin + (nCols-1)*xInc;
//            yMax = yMin + (nRows-1)*yInc;

            StsSurface surface = StsSurface.constructSurface(name, stsColor, StsSurface.IMPORTED,
				nCols, nRows, xOrigin, yOrigin, xInc, yInc, 0.0f, 0.0f, angle, pointsZ, hasNulls,
                StsParameters.nullValue, zDomain, units[0], units[1], progressPanel);
            if(surface == null) return null;
            if(rowNumMin != StsParameters.nullValue)
            {
                surface.rowNumMin = rowNumMin;
                surface.rowNumMax = rowNumMax;
                surface.rowNumInc = rowNumInc;
                surface.colNumMin = colNumMin;
                surface.colNumMax = colNumMax;
                surface.colNumInc = colNumInc;
            }
            surface.setZMin(zMin);
            surface.setZMax(zMax);
//            surface.setUnits(units);

//            grid.writeBinaryFile(getDirectory(), createFilename(name));
            return surface;
		}
		catch(Exception e)
		{
			StsException.outputException("StsGrid.readBinaryFile() failed.",
					e, StsException.WARNING);
			return null;
		}
	}
}
