
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Utilities.Interpolation;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

import java.util.*;

public class StsBlockGridInterpolationLstSqPlaneFit extends StsBlockGridInterpolation
{
/* ================================================================================================*/
/* Least Squares Gridding.                                  	                                    */
/* z = a + b*dx + c*dy                                                                             */
/* Which means use determinants to solve for constant term a                                        */
/* We need to solve a system of equations defined as follows                                        */
/*     wz  =  | w    wx    wy |  a                                                                  */
/*     wxz =  | wx   wxx  wxy |  b                                                                  */
/*     wyz =  | wy   wxy  wyy |  c                                                                  */
/* Using Cramer's rule:                                                                             */
/*                                                                                                  */
/* determinant = w*(wxx*wyy - wxy*wxy) + wx*(wy*wxy - wx*wyy) + wy*(wx*wxy - wxx*wy)                */
/* z0 = a = (wz*(wxx*wyy - wxy*wxy) + wxz*(wxy*wy - wx*wyy) + wyz*(wx*wxy - wxx*wy))/determinant    */
/* dz/dx = b = w*(wxz*wyy - wyz*wxy) + wz*(wxy*wy - wx*wyy) + wy*(wx*wyz - wxz*wy))/determinant     */
/* dz/dy = c = w*(wxx*wyz - wxz*wxy) + wx*(wxz*wy - wx*wyz) + wz*(wx*wxy - wxx*wy))/determinant     */
/*                                                                                                  */
/* =================================================================================================*/



    static StsBlockGridInterpolationLstSqPlaneFit interpolator = null;

    float w, wx, wy, wz, wxx, wxy, wyy, wxz, wyz;

    float determinant;

    static final float minDeterminant = 1.0e-5f;

    static public StsBlockGridInterpolationLstSqPlaneFit getInstance(StsModelSurface surface, StsBlockGrid blockGrid)
     {
         if(interpolator == null) interpolator = new StsBlockGridInterpolationLstSqPlaneFit();
         interpolator.initializeGrid(surface, blockGrid);
         return interpolator;
     }

     public void initializeGrid(StsModelSurface surface, StsBlockGrid blockGrid)
     {
         if(this.surface == surface && this.blockGrid == blockGrid) return;
         if(surface.faultDistances == null && StsBlock.hasModelAuxiliarySections())
             surface.computeFaultArcLengthWeights();
         if(surface.dZdX == null) surface.computeGradients();
         super.initializeGrid(surface, blockGrid);
     }

    static public void constructor(StsModelSurface surface, StsBlockGrid blockGrid)
    {
        if(surface.faultDistances == null && StsBlock.hasModelAuxiliarySections())
            surface.computeFaultArcLengthWeights();

        if(!(interpolator instanceof StsBlockGridInterpolationLstSqPlaneFit)) interpolator = new StsBlockGridInterpolationLstSqPlaneFit();
        interpolator.initializeGrid(surface, blockGrid);
    }

    public float interpolate(int iCenter, int jCenter, boolean useGradient, boolean isInside, boolean debugInterpolate)
    {
        if(interpolator == null)
        {
            StsException.outputException(new StsException(StsException.FATAL, "StsInterpolateLstSqPlaneFit.interpolate() failed.",
                                            "interpolator has not been constructed."));
            return nullValue;
        }
        interpolator.initialize(iCenter, jCenter, useGradient, isInside, debugInterpolate);
        return interpolator.interpolatePoint(iCenter, jCenter);
    }
    /* ================================================================================================*/
    /* Least Squares Gridding.                                  	                                    */
    /* z = a + b*dx + c*dy                                                                             */
    /* Which means use determinants to solve for constant term a                                        */
    /* We need to solve a system of equations defined as follows                                        */
    /*     wz  =  | w    wx    wy |  a                                                                  */
    /*     wxz =  | wx   wxx  wxy |  b                                                                  */
    /*     wyz =  | wy   wxy  wyy |  c                                                                  */
    /* Using Cramer's rule:                                                                             */
    /*                                                                                                  */
    /* determinant = w*(wxx*wyy - wxy*wxy) + wx*(wy*wxy - wx*wyy) + wy*(wx*wxy - wxx*wy)                */
    /* z0 = a = (wz*(wxx*wyy - wxy*wxy) + wxz*(wxy*wy - wx*wyy) + wyz*(wx*wxy - wxx*wy))/determinant    */
    /* dz/dx = b = w*(wxz*wyy - wyz*wxy) + wz*(wxy*wy - wx*wyy) + wy*(wx*wyz - wxz*wy))/determinant     */
    /* dz/dy = c = w*(wxx*wyz - wxz*wxy) + wx*(wxz*wy - wx*wyz) + wz*(wx*wxy - wxx*wy))/determinant     */
    /*                                                                                                  */
    /* =================================================================================================*/
    
    float getZ(ArrayList sortedPoints)
    {
        StsBlockGridInterpolation.Point point;
        float z, wt, dx, dy;
        int nPointsUsed;

        try
        {
            int nTotalPoints = sortedPoints.size();
            nPointsUsed = Math.min(nTotalPoints, currentMinNPoints);
            point = (StsBlockGridInterpolation.Point)sortedPoints.get(nPointsUsed-1);
            float weightMin = point.wt;
            for(int n = nPointsUsed+1; n < nTotalPoints; n++)
            {
                point = (StsBlockGridInterpolation.Point)sortedPoints.get(n);
                if(point.wt < weightMin) break;
                nPointsUsed++;
            }

            // remove the points we are not going to use
            for(int n = nTotalPoints-1; n >= nPointsUsed; n--)
                sortedPoints.remove(n);

            if(!useGradient)
            {
                w = 0.0f;
                wz = 0.0f;

                for(int n = 0; n < nPointsUsed; n++)
                {
                    point = (StsBlockGridInterpolation.Point)sortedPoints.get(n);
                    z = point.z;
                    wt = point.wt;
                    w += wt;
                    wz += wt*z;
                }
                z = wz/w;
                if(debugInterpolate) StsMessageFiles.infoMessage("      not using gradient. z: " + z);
            }
            else
            {
                w = 0.0f;
                wx = 0.0f;
                wy = 0.0f;
                wz = 0.0f;
                wxx = 0.0f;
                wxy = 0.0f;
                wyy = 0.0f;
                wxz = 0.0f;
                wyz = 0.0f;

                for(int n = 0; n < nPointsUsed; n++)
                {
                    point = (StsBlockGridInterpolation.Point)sortedPoints.get(n);
                    z = point.z;
                    wt = point.wt - weightMin;
                    dx = point.dx;
                    dy = point.dy;

                    w += wt;
                    wz += wt*z;
                    wx += wt*dx;
                    wy += wt*dy;
                    wxx += wt*dx*dx;
                    wxy += wt*dx*dy;
                    wyy += wt*dy*dy;
                    wxz += wt*dx*z;
                    wyz += wt*dy*z;
                }
                if(computeDeterminantOK())
                {
                    z = (wz*(wxx*wyy - wxy*wxy) + wxz*(wxy*wy - wx*wyy) + wyz*(wx*wxy - wxx*wy))/determinant;
                    if(debugInterpolate) StsMessageFiles.infoMessage("      determinant OK; using gradient. z: " + z);
                }
                else
                {
                    z = wz/w;
                    if(debugInterpolate) StsMessageFiles.infoMessage("      determinant too small, not using gradient. z: " + z);
                }
                if(debugInterpolate) StsMessageFiles.infoMessage("      determinant: " + determinant);
            }

            if(debugInterpolate)
            {
                for(int n = 0; n < nPointsUsed; n++)
                {
                    point = (StsBlockGridInterpolation.Point)sortedPoints.get(n);
                    StsMessageFiles.infoMessage("     " + n + " interpolatePoint row: " + point.row +
                                  " col: " + point.col + " z: " + point.z + " wt: " + point.wt);
                }
            }
            return z;
        }
        catch(Exception e)
        {
            StsException.outputException("StsINterpolationLstSqPlaneFit.getZ() failed.",
                e, StsException.WARNING);
            return nullValue;
        }
    }

    boolean computeDeterminantOK()
    {
        if(nPoints == 0)
        {
            determinant =  0.0f;
            return false;
        }
        else
        {
            determinant = w*(wxx*wyy - wxy*wxy) + wx*(wy*wxy - wx*wyy) + wy*(wx*wxy - wxx*wy);
            return Math.abs(determinant) >= minDeterminant;
        }
    }
}
