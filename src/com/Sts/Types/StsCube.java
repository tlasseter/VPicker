

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

public class StsCube
{
	public float xMin = -50.0f;
	public float xMax = 50.0f;
    public float dX = 10.0f;
	public float yMin = -50.0f;
	public float yMax = 50.0f;
    public float dY = 10.0f;
	public float zMin = 0.0f;
	public float zMax = 100.0f;
    public float dZ = 10.0f;
    public float angle = 0.0f;

	public double originX;
	public double originY;

    private float xSize, ySize, zSize;
    private float fnjx, fniy;

    private StsMatrix4f cubeMatrix = null;
    private StsMatrix4f cubeInverseMatrix = null;

    public StsCube()
    {
        xSize = xMax - xMin;
        ySize = yMax - yMin;

        fnjx = xSize/dX;
        fniy = ySize/dY;

        createTransformMatrices();
    }

    private void createTransformMatrices()
    {
        // make matrix to transform from world to index cube coordinates

        cubeMatrix = new StsMatrix4f();
        cubeMatrix.setIdentity();
        cubeMatrix.translate(-xMin, -yMin, -zMin);
        if(angle != 0.0) cubeMatrix.rotZ(-angle);
        cubeMatrix.scale(1.0f/dX, 1.0f/dY, 1.0f/dZ);

        // make inverse-matrix to transform from index cube to world coordinates

        cubeInverseMatrix = new StsMatrix4f();
        cubeInverseMatrix.setIdentity();
        cubeInverseMatrix.scale(dX, dY, dZ);
        if(angle != 0.0) cubeInverseMatrix.rotZ(angle);
        cubeInverseMatrix.translate(xMin, yMin, zMin);
    }

    // Given an index cube from 0 to xSize and 0 to ySize, clip line from start to end

    public static boolean clipLineToCube(StsPoint2D start, StsPoint2D end, float xSize, float ySize)
    {
	    float max, min;

	    /* clip ends between x = 0.0 and x = xSize 				*/

	    if(start.x < end.x)
		    clipXLineToCube(xSize, start, end);
	    else
		    clipXLineToCube(xSize, end, start);

	    /* check if clipped line is below y = 0 or above y = ySize	*/

	    max = Math.max(start.y, end.y);
	    min = Math.min(start.y, end.y);

	    if(max <= 0.0 || min >= ySize)
		    return false;

	    /* clip ends between y = 0.0 and y = ySize 					*/

	    if(start.y < end.y)
		    clipYLineToCube(ySize, start, end);
    	else
		    clipYLineToCube(ySize, end, start);

	    /* check if clipped line is below x = 0 or above x = xSize	*/

	    max = Math.max(start.x, end.x);
	    min = Math.min(start.x, end.x);

	    if(max <= 0.0 || min >= xSize)
		    return false;
	    else
		    return true;
    }

    // Line is increasing in X-index from start to end and extends beyond
    // an X-direction boundary: clip it.

    static void clipXLineToCube(float xMax, StsPoint2D start, StsPoint2D end)
    {
	    float f;

	    if(start.x < 0.0 && end.x > 0.0)
	    {
		    f = -start.x/(end.x - start.x);
		    start.y = start.y + f*(end.y - start.y);
		    start.x = 0.0f;
	    }

	    if(end.x > xMax && start.x < xMax)
	    {
		    f = (xMax - start.x)/(end.x - start.x);
		    end.y = start.y + f*(end.y - start.y);
		    end.x = xMax;
	    }
    }

    static void clipXLineToCube(float xMin, float xMax, float[] start, float[] end)
    {
	    float f;

	    if(start[0] < xMin && end[0] > xMin)
	    {
		    f = (xMin - start[0])/(end[0] - start[0]);
		    start[0] = xMin;
		    start[1] = start[1] + f*(end[1] - start[1]);
		    start[2] = start[2] + f*(end[2] - start[2]);
	    }

	    if(end[0] > xMax && start[0] < xMax)
	    {
		    f = (xMax - start[0])/(end[0] - start[0]);
		    end[0] = xMax;
		    end[1] = start[1] + f*(end[1] - start[1]);
		    end[2] = start[2] + f*(end[2] - start[2]);
	    }
    }

    // Line is increasing in Y-index from start to end and extends beyond
    // a Y-direction boundary: clip it.

    static void clipYLineToCube(float yMax, StsPoint2D start, StsPoint2D end)
    {
	    float f;

	    if(start.y < 0.0f && end.y > 0.0f)
	    {
		    f = -start.y/(end.y - start.y);
		    start.x = start.x + f*(end.x - start.x);
		    start.y = 0.0f;
	    }

	    if(end.y > yMax && start.y < yMax)
	    {
		    f = (yMax - start.y)/(end.y - start.y);
		    end.x = start.x + f*(end.x - start.x);
		    end.y = yMax;
	    }
    }

    static void clipYLineToCube(float yMin, float yMax, float[] start, float[] end)
    {
	    float f;

	    if(start[1] < yMin && end[1] > yMin)
	    {
		    f = (yMin - start[1])/(end[1] - start[1]);
		    start[1] = yMin;
		    start[0] = start[0] + f*(end[0] - start[0]);
		    start[2] = start[2] + f*(end[2] - start[2]);
	    }

	    if(end[1] > yMax && start[1] < yMax)
	    {
		    f = (yMax - start[1])/(end[1] - start[1]);
		    end[1] = yMax;
		    end[0] = start[0] + f*(end[0] - start[0]);
		    end[2] = start[2] + f*(end[2] - start[2]);
	    }
    }

    static void clipZLineToCube(float zMin, float zMax, float[] start, float[] end)
    {
	    float f;

	    if(start[2] < zMin && end[2] > zMin)
	    {
		    f = (zMin - start[2])/(end[2] - start[2]);
		    start[2] = zMin;
		    start[0] = start[0] + f*(end[0] - start[0]);
		    start[1] = start[1] + f*(end[1] - start[1]);
	    }

	    if(end[2] > zMax && start[2] < zMax)
	    {
		    f = (zMax - start[2])/(end[2] - start[2]);
		    end[2] = zMax;
		    end[0] = start[0] + f*(end[0] - start[0]);
		    end[1] = start[1] + f*(end[1] - start[1]);
	    }
    }

    public float getFnjx()
    {
        return fnjx;
    }

    public float getFniy()
    {
        return fniy;
    }

    public static boolean clipLineToCube(float[] start, float[] end, float xMin, float xMax,
                                         float yMin, float yMax, float zMin, float zMax)
    {
	    float max, min;

	    /* clip ends between x = xMin and x = xMax 				*/

	    max = Math.max(start[0], end[0]);
	    min = Math.min(start[0], end[0]);

	    if(max <= xMin || min >= xMax) return false;

	    if(start[0] < end[0])
		    clipXLineToCube(xMin, xMax, start, end);
	    else
		    clipXLineToCube(xMin, xMax, end, start);

	    /* check if clipped line is below yMin or above yMax	*/

	    max = Math.max(start[1], end[1]);
	    min = Math.min(start[1], end[1]);

	    if(max <= yMin || min >= yMax) return false;

	    /* clip ends between y = yMin and y = yMax 					*/

	    if(start[1] < end[1])
		    clipYLineToCube(yMin, yMax, start, end);
    	else
		    clipYLineToCube(yMin, yMax, end, start);

	    /* clip ends between z = zMin and z = zMax  */

	    max = Math.max(start[2], end[2]);
	    min = Math.min(start[2], end[2]);
	    if(max <= zMin || min >= zMax) return false;

	    if(start[2] < end[2])
		    clipZLineToCube(zMin, zMax, start, end);
    	else
		    clipZLineToCube(zMin, zMax, end, start);

        return true;
    }
}
