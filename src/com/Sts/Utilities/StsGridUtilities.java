
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Utilities;



public class StsGridUtilities
{
    static public final float nullValue = StsParameters.nullValue;
/*
   	static public float[][][] computeGridNormals(float[][] zValues, StsBoundingBox boundingBox,
                                               float[][][] normals)
    {
		int row, col;
        float z, zRite, zUp;
        float[] normal = null;
        int nRows, nCols;
        float xInc, yInc;
        float[] difRite, difUp;

        try
        {
            nRows = boundingBox.getNRows();
            nCols = boundingBox.getNCols();
            xInc = boundingBox.xInc;
            yInc = boundingBox.yInc;

            if(normals == null) normals = new float[nRows][nCols][];

            difRite = new float[3];
            difUp = new float[3];

            difUp[0] = 0.0f; difUp[1] = yInc;
            difRite[0] = xInc; difRite[1] = 0.0f;

            for(row = 0; row < nRows-1; row++)
            {
                zRite = zValues[row][0];
                for(col = 0; col < nCols-1; col++)
                {
                    z = zRite;
                    zRite = zValues[row][col+1];
                    if(z == nullValue || zRite == nullValue) continue;
                    zUp = zValues[row+1][col];
                    if(zUp == nullValue) continue;
                    difUp[2] = zUp - z;
                    difRite[2] = zRite - z;
                    normal = StsGridPoint.leftCrossProduct(difUp, difRite);
                    normals[row][col] = normal;
                }
                // Copy normal for last column from second to last column
                normals[row][nCols-1] = normal;
            }

            // Copy normals for top row from second to top row
            for(col = 0; col < nCols; col++)
                normals[nRows-1][col] = normals[nRows-2][col];

            return normals;
        }
        catch(Exception e)
        {
            StsException.outputException("StsGridUtilities.computeGridNormals() failed.",
                e, StsException.WARNING);
            return null;
        }
	}
*/
}
