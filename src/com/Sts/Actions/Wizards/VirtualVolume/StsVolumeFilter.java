package com.Sts.Actions.Wizards.VirtualVolume;

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Aug 13, 2007
 * Time: 9:24:43 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsVolumeFilter implements StsVolumeFilterFace
{
    private String name;
    private int xSize = 3;
    private int ySize = 3;
    private int zSize = 3;

    transient StsJPanel filterPanel=null;

    public static final String BARTLETT = StsFilterVVolumePanel.BARTLETT;
    public static final String GAUSSIAN = StsFilterVVolumePanel.GAUSSIAN;
    public static final String LAPLACIAN = StsFilterVVolumePanel.LAPLACIAN;
    public static final String MEAN = StsFilterVVolumePanel.MEAN;
    public static final String MEDIAN = StsFilterVVolumePanel.MEDIAN;
    public static final String MAX = StsFilterVVolumePanel.MAX;
    public static final String MIN = StsFilterVVolumePanel.MIN;
    public static final String VARIANCE = StsFilterVVolumePanel.VARIANCE;
    
    public static final String KUWAHARA = StsEPFVVolumePanel.KUWAHARA;

    public static final int XDIR = StsParameters.XDIR;
    public static final int YDIR = StsParameters.YDIR;
    public static final int ZDIR = StsParameters.ZDIR;

    private static final byte CONVOLUTION = 1;
    private static final byte RANK = 2;
    private static final byte MULTI_WINDOW = 3;

    public StsVolumeFilter()
    {
    }

    public StsVolumeFilter(String name)
    {
        this.name = name;
        initialize();
    }

    static public StsVolumeFilter constructor(String name)
    {
        return new StsVolumeFilter(name);
    }

    public void initialize()
    {
        //Add filter panels for filter specific parameters
    }

    public void setName(String name) {this.name = name;}
    public String getName() {return name;}
    public void setXSize(int size) {xSize = size;}
    public void setYSize(int size) {ySize = size;}
    public void setZSize(int size) {zSize = size;}
    public int getXSize() { return xSize; }
    public int getYSize() { return ySize; }
    public int getZSize() { return zSize; }

    public StsJPanel getFilterPanel() { return filterPanel; }

    public byte[] processBytePlaneData(int dir, int nCenterPlane, StsSeismicVolume inputVolume, StsVirtualVolume outputVolume)
    {
        byte[] planeByteData;
        float[] outPlaneFloatData;
        float[][] planeFloatData;

        try
        {
            if (inputVolume == null || outputVolume == null) return null;
            int nPlanes = getNPlanes(dir);
            int nSidePlanes = nPlanes/2;
            int nPlaneMin = Math.max(0, nCenterPlane - nSidePlanes);
            int nDirSize = outputVolume.getBoxDimension(dir);
            int nPlaneMax = Math.min(nCenterPlane + nSidePlanes, nDirSize-1);
            nPlanes = nPlaneMax - nPlaneMin + 1;
            int nPlaneRows = outputVolume.getNCursorRows(dir);
            int nPlaneCols = outputVolume.getNCursorCols(dir);
            int nPlanePoints = nPlaneRows*nPlaneCols;

            planeFloatData = new float[nPlanes][nPlanePoints];
            for (int plane=0; plane<nPlanes; plane++)
            {
                planeByteData = inputVolume.readBytePlaneData(dir, nPlaneMin + plane);
                if (planeByteData == null) return null;

                nPlanePoints = planeByteData.length;
                for(int i=0; i<nPlanePoints; i++)
                    planeFloatData[plane][i] = inputVolume.getScaledValue(planeByteData[i]);
            }
            nCenterPlane -= nPlaneMin;
            outPlaneFloatData =  processPlaneData(dir, nPlaneRows, nPlaneCols, nPlaneMin, nPlaneMax,
                                                  planeFloatData, nPlanePoints, outputVolume);
            return outputVolume.computeUnscaledByteValues(outPlaneFloatData);
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "processBytePlaneData", "StsSeismicVolume.readPlaneData() failed.", e);
            return null;
        }
    }

    public float[] processFloatPlaneData(int dir, int nCenterPlane, StsSeismicVolume inputVolume, StsVirtualVolume outputVolume)
     {
         float[][] planeFloatData;

         try
         {
             if (inputVolume == null || outputVolume == null) return null;
             int nPlanes = getNPlanes(dir);
             int nSidePlanes = nPlanes/2;
             int nPlaneMin = Math.max(0, nCenterPlane - nSidePlanes);
             int nDirSize = outputVolume.getBoxDimension(dir);
             int nPlaneMax = Math.min(nCenterPlane + nSidePlanes, nDirSize-1);
             nPlanes = nPlaneMax - nPlaneMin + 1;
             int[] nRowCols = outputVolume.getCursorDisplayNRowCols(dir);
             int nPlaneRows = nRowCols[0];
             int nPlaneCols = nRowCols[1];
             int nPlanePoints = nPlaneRows*nPlaneCols;

             planeFloatData = new float[nPlanes][nPlanePoints];
             for (int n = 0, nPlane = nPlaneMin; nPlane < nPlanes; nPlane++, n++)
                 planeFloatData[n] = inputVolume.readRowPlaneFloatData(nPlane);

             nCenterPlane -= nPlaneMin;
             return processPlaneData(dir, nPlaneRows, nPlaneCols, nPlaneMin, nPlaneMax, planeFloatData, nPlanePoints, outputVolume);
         }
         catch (Exception e)
         {
             StsException.outputWarningException(this, "processBytePlaneData", "StsSeismicVolume.readPlaneData() failed.", e);
             return null;
         }
     }


    public float[] processPlaneData(int dir, int nPlaneRows, int nPlaneCols, int nPlaneMin, int nPlaneMax,
                                   float[][] planeData, int nPlanePoints, StsVirtualVolume outputVolume)
    {
        float[] outFloatValues = new float[nPlanePoints];
        byte filterType = getFilterType(name);
        if( filterType==CONVOLUTION)
            outFloatValues = StsConvolve.convolveFloat3D(dir, nPlaneRows, nPlaneCols, nPlaneMin, nPlaneMax, planeData,
                                                         name, xSize, ySize, zSize);
        else if( filterType==RANK)
            outFloatValues = StsRankFilters.rankFloat3D(dir, nPlaneRows, nPlaneCols, nPlaneMin, nPlaneMax, planeData,
                                                        name, xSize, ySize, zSize);
        
        else if( filterType==MULTI_WINDOW)
            outFloatValues = StsMultiWindowFilters.multiWinFloat3D(dir, nPlaneRows, nPlaneCols, nPlaneMin, nPlaneMax, planeData,
                                                        name, xSize, ySize, zSize);
        return outFloatValues;
    }

    private int getNPlanes(int dir)
    {
        switch(dir)
        {
            case XDIR:
                return ySize;
            case YDIR:
                return xSize;
            case ZDIR:
                return zSize;
            default:
                return 0;
        }
    }

    private byte getFilterType(String name)
    {
        if(name.equals(BARTLETT))
            return CONVOLUTION;

        else if(name.equals(GAUSSIAN))
            return CONVOLUTION;

        else if(name.equals(LAPLACIAN))
            return CONVOLUTION;

        else if(name.equals(MEDIAN))
            return RANK;

        else if(name.equals(MEAN))
            return RANK;

        else if(name.equals(MIN))
            return RANK;

        else if(name.equals(MAX))
            return RANK;

        else if(name.equals(VARIANCE))
            return RANK;
        
        else if(name.equals(KUWAHARA))
            return MULTI_WINDOW;

        return 0;
    }
}