package com.Sts.SeismicAttributes;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

import com.Sts.Actions.Wizards.VolumeFilter.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;


public class StsVolumeFilterConstructor extends StsSeismicVolumeConstructor
{
    StsSeismicVolume inputVolume;
    int nSamples;
    //Moving window dimensions
    int winColSize = 3;
    int winRowSize = 3;
    int winSliceSize = 1;
 
    StsTimer timer = null;
    boolean runTimer = false;

    static public final String KUWAHARA = "kuwahara";
    static public final String LUM = "lowerUpperMiddle";

    private StsVolumeFilterConstructor(StsModel model, StsSeismicVolume volume, boolean isDataFloat, StsProgressPanel panel,
    									int winX, int winY, int winZ, byte analysis)
    {
        this.model = model;
        this.panel = panel;
        winColSize = winX;
        winRowSize = winY;
        winSliceSize = winZ;
        StsMultiWindowFilters.winMethod = analysis;
        String analysisType = "MEAN";
        if (analysis == StsVolumeFilterWizard.ANALYSIS_MEDIAN)
        	analysisType = "MEDIAN";

		nInputRows = volume.nRows;
		nInputCols = volume.nCols;

        volumeName = KUWAHARA + "_" + analysisType + winColSize + winRowSize + winSliceSize;
        float dataMin = volume.dataMin;
        float dataMax = volume.dataMax;
        outputVolume = StsSeismicVolume.initializeAttributeVolume(model, volume, dataMin, dataMax, true, true, volume.stemname, volumeName, "rw");

		if(panel != null) panel.initialize(nInputRows);

        nSamples = volume.nSlices;
        inputVolume = volume;
        createOutputVolume();
    }

    static public StsVolumeFilterConstructor constructor(StsModel model, StsSeismicVolume data, boolean isDataFloat, StsProgressPanel panel,
			int winX, int winY, int winZ, byte analysis)
    {
        try
        {
        	StsVolumeFilterConstructor filtConstructor = new StsVolumeFilterConstructor(model, data, isDataFloat, panel, winX, winY, winZ, analysis);
            return filtConstructor;
        }
        catch(Exception e)
        {
            StsMessage.printMessage("StsVolumeFilterConstructor.constructor() failed.");
            return null;
        }
    }

    public boolean doProcessInputBlock(int nBlock, StsMappedBuffer[] inputBuffers)
    {
        if (runTimer)
        {
            timer.start();
        }
        processFloatPlaneData(nBlock);
		if (isCanceled()) return false;
		if(panel != null )  panel.setValue(nBlock);
		if (runTimer)
        {
            timer.stopPrint("process " + nInputBlockTraces + " traces for block " + nBlock + ":");

        }
        return true;
    }
    
    public void processFloatPlaneData(int row)
    {
        float[][][] floatData;

        try
        {
            if (inputVolume == null || outputVolume == null) return;
			int stRow = Math.max(0, row-winRowSize+1);
   			int enRow = Math.min(row+winRowSize-1, nInputRows-1);
            
   			for (int col=0; col<nInputCols; col++)
    		{
   				int stCol = Math.max(0, col-winColSize+1);
   				int enCol = Math.min(col+winColSize-1, nInputCols-1);
  				for (int slice=0; slice<nSamples; slice++)
   				{
  					int stSlice = Math.max(0, slice-winSliceSize+1);
  					int enSlice = Math.min(slice+winSliceSize-1, nSamples-1);
					floatData = new float[enSlice-stSlice+1][enRow-stRow+1][enCol-stCol+1];
					for (int s=stSlice; s<=enSlice; s++)
					{
	  					for (int r=stRow; r<=enRow; r++)
						{
				            for (int c=col; c<=enCol; c++)
				            {
				            	floatData[s-stSlice][r-stRow][c-stCol]= inputVolume.getRowFloatBlockValue(r, c, s);
				            }
	        			}
					}
					float val = StsMultiWindowFilters.applyKuwaharaFilter(floatData, row, stRow, enRow, winRowSize, col, stCol, enCol, winColSize,
							stSlice, enSlice, winSliceSize);
					outputFloatBuffer.put(val);
					int intVal = (int)val;
					outputVolume.accumulateHistogram((byte) intVal);
            	} 				
            }          	   
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "processBytePlaneData", "StsVolumeFilterConstructor.processFloatPlaneData() failed.", e);
            return;
        }
    }

}