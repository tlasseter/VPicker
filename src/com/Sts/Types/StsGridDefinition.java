
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;


import com.Sts.DB.DBCommand.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

public class StsGridDefinition extends StsRotatedGridBoundingBox implements Cloneable, StsXYGridable
{
    public double rowXInc, rowYInc, colXInc, colYInc;

    static public final double doubleNullValue = StsParameters.doubleNullValue;

    static final float nullValue = StsParameters.nullValue;

    public StsGridDefinition()
    {
    }

    public StsGridDefinition(int nRows, int nCols, double xOrigin, double yOrigin,
                             float xMin, float yMin, float xInc, float yInc, StsModel model)
    {
        this.nRows = nRows;
        this.nCols = nCols;
        this.xOrigin = xOrigin;
        this.yOrigin = yOrigin;
        this.xMin = xMin;
        this.yMin = yMin;
        this.xInc = xInc;
        this.yInc = yInc;
 //       initializeGrid();
        StsChangeCmd cmd = new StsChangeCmd(currentModel.getProject(), this, "gridDefinition", false);
        model.addTransactionCmd("StsGridDefinition", cmd);
    }

    public StsGridDefinition(StsXYSurfaceGridable grid, StsModel model)
     {
         xOrigin = grid.getXOrigin();
         yOrigin = grid.getYOrigin();
         xInc = grid.getXInc();
         yInc = grid.getYInc();
         zInc = grid.getZInc();
         nRows = grid.getNRows();
         nCols = grid.getNCols();
         angle = grid.getAngle();
		 xMin = grid.getXMin();
		 yMin = grid.getYMin();
		 xMax = xMin + grid.getXSize();
		 yMax = yMin + grid.getYSize();
         setAngle();
         StsChangeCmd cmd = new StsChangeCmd(currentModel.getProject(), this, "gridDefinition", false);
         model.addTransactionCmd("StsGridDefinition", cmd);
    }

    public boolean initialize(StsModel model)
    {
        return true;
    }

    public float[][] getPointsZ()
    {
        StsException.systemError("StsGridDefinition.getPointsZ() called: not implementable. ");
        return null;
    }

    public float[][] getAdjPointsZ()
    {
        StsException.systemError("StsGridDefinition.getAdjPointsZ() called: not implementable. ");
        return null;
    }

    public float interpolateBilinearZ(StsPoint point, boolean computeIfNull, boolean setPoint)
    {
        StsException.systemError("StsGridDefinition.interpolateBilinearZ(StsPoint...) called: not implementable. ");
        return nullValue;
    }

    public float interpolateBilinearZ(StsGridPoint gridPoint, boolean computeIfNull, boolean setPoint)
    {
        StsException.systemError("StsGridDefinition.interpolateBilinearZ(StsGridPoint) called: not implementable. ");
        return nullValue;
    }

    public float getComputePointZ(int row, int col)
    {
        StsException.systemError("StsGridDefinition.getComputePointZ() called: not implementable. ");
        return nullValue;
    }

    public StsPoint getPoint(int row, int col)
    {
        StsException.systemError("StsGridDefinition.getStsPoint() called: not implementable. ");
        return null;
    }

    public float[] getXYZ(int row, int col)
    {
        StsException.systemError("StsGridDefinition.getXYZ() called: not implementable. ");
        return null;
    }

    public void checkConstructGridNormals()
    {
        StsException.systemError("StsGridDefinition.checkConstructGridNormals() called: not implementable. ");
        return;
    }

    public float[] getNormal(int row, int col)
    {
        StsException.systemError("StsGridDefinition.getNormal() called: not implementable. ");
        return null;
    }

    public float[] getNormal(float rowF, float colF)
    {
        StsException.systemError("StsGridDefinition.getNormal() called: not implementable. ");
        return null;
    }

    public boolean clipLineToGrid(StsPoint point0, StsPoint point1)
    {
        StsGridPoint start = new StsGridPoint(point0, this);
        StsGridPoint end = new StsGridPoint(point1, this);
        return StsGridDefinition.clipLineToGrid(start, end, this, false);
    }

    static public boolean clipLineToGrid(StsGridPoint start, StsGridPoint end, StsXYGridable grid, boolean outsideGridOK)
    {
        if(!clipColLine(start, end, grid, outsideGridOK)) return false;
        if(!clipRowLine(start, end, grid, outsideGridOK)) return false;

        return true;
    }

    // Line is increasing in X-index from start to end and extends beyond
    // an X-direction boundary: clip it.

    static private boolean clipColLine(StsGridPoint start, StsGridPoint end, StsXYGridable grid, boolean outsideGridOK)
    {
	    float f;

	    float max = Math.max(start.colF, end.colF);
	    float min = Math.min(start.colF, end.colF);


        float colMin = grid.getColMin();
        float colMax = grid.getColMax();

        if(outsideGridOK)
        {
            float addRange = 0.1f*(colMax - colMin);
            colMin -= addRange;
            colMax += addRange;
        }

        if(max < colMin || min > colMax) return false;
        if(min >= colMin && max <= colMax) return true;

	    if(start.colF < colMin)
	    {
		    f = (colMin - start.colF)/(end.colF - start.colF);
            start.interpolate(start, end, f);
	    }
        else if(end.colF < colMin)
	    {
		    f = (colMin - start.colF)/(end.colF - start.colF);
            end.interpolate(start, end, f);
	    }

	    if(start.colF > colMax)
	    {
		    f = (colMax - start.colF)/(end.colF - start.colF);
		    start.interpolate(start, end, f);
	    }
	    else if(end.colF > colMax)
	    {
		    f = (colMax - start.colF)/(end.colF - start.colF);
		    end.interpolate(start, end, f);
	    }

        return true;
    }

    static private boolean clipRowLine(StsGridPoint start, StsGridPoint end, StsXYGridable grid, boolean outsideGridOK)
    {
	    float f;

	    float max = Math.max(start.rowF, end.rowF);
	    float min = Math.min(start.rowF, end.rowF);

        float rowMin = grid.getRowMin();
        float rowMax = grid.getRowMax();

        if(outsideGridOK)
        {
            float addRange = 0.1f*(rowMax - rowMin);
            rowMin -= addRange;
            rowMax += addRange;
        }
        if(max < rowMin || min > rowMax) return false;
        if(min >= rowMin && max <= rowMax) return true;

	    if(start.rowF < rowMin)
	    {
		    f = (rowMin - start.rowF)/(end.rowF - start.rowF);
            start.interpolate(start, end, f);
	    }
        else if(end.rowF < rowMin)
	    {
		    f = (rowMin - start.rowF)/(end.rowF - start.rowF);
            end.interpolate(start, end, f);
	    }

	    if(start.rowF > rowMax)
	    {
		    f = (rowMax - start.rowF)/(end.rowF - start.rowF);
		    start.interpolate(start, end, f);
	    }
	    else if(end.rowF > rowMax)
	    {
		    f = (rowMax - start.rowF)/(end.rowF - start.rowF);
		    end.interpolate(start, end, f);
	    }

        return true;
    }

    static private boolean clipLineToZLimits(StsGridPoint start, StsGridPoint end, StsXYSurfaceGridable grid)
    {
	    float f;

        float zMin = grid.getZMin();
        float zMax = grid.getZMax();
		if(zMax - zMin < 1.0f)
		zMax += 0.5f;
	    zMin -= 0.5f;

        float startZ, endZ;
    /*
        if(isDepth)
        {
            startZ = start.pxyz[3];
            endZ = end.pxyz[3];
        }
        else
        {
    */
            startZ = start.getZorT();
            endZ = end.getZorT();
 //       }

	    float max = Math.max(startZ, endZ);
	    float min = Math.min(startZ, endZ);

        if(max < zMin || min > zMax) return false;
        if(min >= zMin && max <= zMax) return true;

	    if(startZ < zMin)
	    {
		    f = (zMin - startZ)/(endZ - startZ);
            start.interpolate(start, end, f);
            startZ = start.getZorT();
        }
        else if(endZ < zMin)
	    {
		    f = (zMin - startZ)/(endZ - startZ);
            end.interpolate(start, end, f);
            endZ = end.getZorT();
        }
	    if(startZ > zMax)
	    {
		    f = (zMax - startZ)/(endZ - startZ);
		    start.interpolate(start, end, f);
	    }
        else if(endZ > zMax)
        {
		    f = (zMax - startZ)/(endZ - startZ);
		    end.interpolate(start, end, f);
	    }
    /*
        if(isDepth)
         {
             startZ = start.pxyz[3];
             endZ = end.pxyz[3];
         }
         else
         {
    */
             startZ = start.getZorT();
             endZ = end.getZorT();
//        }
        return true;
    }

    static public StsGridPoint getSurfacePosition(double[] p0, double[] p1, StsXYSurfaceGridable grid, boolean outsideGridOK)
    {
        boolean debug = false;

        StsGridPoint point0 = new StsGridPoint(p0, grid);
        StsGridPoint point1 = new StsGridPoint(p1, grid);

        if(!clipLineToZLimits(point0, point1, grid)) return null;
        if(!clipLineToGrid(point0, point1, grid, outsideGridOK)) return null;

        int difRows = Math.abs(point0.row - point1.row);
        int difCols = Math.abs(point0.col - point1.col);

        // calculate a increment to step along the near-far line
        float df = 0.99999f / StsMath.max3(difRows, difCols, 1);

        // step along the near-far line looking for an intersection
        float lastError = nullValue;
        float lastF = 0.0f;

        StsGridPoint point = new StsGridPoint(grid);
//        StsGridPoint lastPoint = new StsGridPoint();

        float f = 0.0f;
        while (f < 1.0f)
        {
            getBoundedGridPoint(grid, point0, point1, f, point);
//            point.interpolate(point0, point1, f);
            if(debug) System.out.println("interp Pt:  x = " + point.getX() +
                ", y = " + point.getY() + ", z = " + point.getZorT());
            float lineZ = point.getZorT();
            float gridZ = grid.interpolateBilinearZ(point, false, false);

            if(debug) System.out.println("interpolated z = " + gridZ);
            if(gridZ != nullValue)
            {
                float error = gridZ - lineZ;
                if (error == 0.0f) return point;

                if (lastError != nullValue)
                {
                    if (lastError<0.0f && error>0.0f || lastError>0.0f && error<0.0f)
                    {
                        // interpolate between lastF (f - df) and f
                        float ff = error/(error - lastError);
                        f = f - ff*df;
                        getBoundedGridPoint(grid, point0, point1, f, point);
//                        point.interpolate(point0, point1, f);

                        float pickedZ = grid.interpolateBilinearZ(point, false, false);

                        if(debug) System.out.println("interp Pt:  x = " + point.getX() +
                            ", y = " + point.getY() + ", z = " + pickedZ);

//                        StsGridSectionPoint gridPoint = new StsGridSectionPoint(point.pxyz, false);
                        return point;
                    }
                }
                lastError = error;
                lastF = f;
            }
            f += df;
        }

        return null;
    }

    static private void getBoundedGridPoint(StsXYSurfaceGridable grid, StsGridPoint point0, StsGridPoint point1,
                                     float f, StsGridPoint point)
    {
        point.interpolate(point0, point1, f);

        int rowMin = grid.getRowMin();
        int rowMax = grid.getRowMax();
        if(point.rowF < rowMin)
            point.setRow(rowMin);
        else if(point.rowF > rowMax)
            point.setRow(rowMax);

        int colMin = grid.getColMin();
        int colMax = grid.getColMax();
        if(point.colF < colMin)
            point.setCol(colMin);
        else if(point.colF > colMax)
            point.setCol(colMax);
    }

    // interface compatability
    public String getLabel()
    {
        return new String("grid Definition");
    }
}