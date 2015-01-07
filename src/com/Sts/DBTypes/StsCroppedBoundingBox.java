package com.Sts.DBTypes;

import com.Sts.Types.*;
import com.Sts.Utilities.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Nov 30, 2007
 * Time: 5:02:20 PM
 * To change this template use File | Settings | File Templates.
 */

/** croppedBoundingBox is initialized to original boundingBox and then reduced in size based on a cropped box passed in.
 *  Specifically, the row, col, and slice range variables in this class are adjusted from the originalBoundingBox to the croppedBoundingBox.
 *  On completion, nRows, nCols, and nCroppedSlices correspond with the croppedVolume
 */
public class StsCroppedBoundingBox extends StsRotatedGridBoundingBox
{
    public StsRotatedGridBoundingBox originalBoundingBox;
    // additional cropped parameters required in addition to rotatedGridBoundingBox parameters in superClass StsRotatedGridBoundingBox
    public int rowMin, rowMax, rowInc, colMin, colMax, colInc, sliceMin, sliceMax, sliceInc;
    public  boolean nullClip = false;
    public boolean isCropped = false;
    // number of samples in a row plane of the output cube
    public int nSamplesPerRow = 0;
    // number of samples in a col plane of the output cube
    public int nSamplesPerCol = 0;
    // number of samples in a slice plane of the output cube
    public int nSamplesPerSlice = 0;

    // row number of the first and last rows in the current block being processed
    public int blockRowMin = 0;
    public int blockRowMax = 0;

    // number of rows in this block
    public int nBlockRows = 0;

    // convenience values for the block
    public int nBlockColSamples;
    public int nBlockSliceSamples;
    public int nBlockSamples;

    public int outputRowMin = 0;
    public int outputRowMax = -1;

    public StsCroppedBoundingBox()
    {
    }

    public StsCroppedBoundingBox(boolean persistent)
    {
        super(persistent);
    }

    public StsCroppedBoundingBox(StsRotatedGridBoundingBox boundingBox, StsRotatedGridBoundingBox croppedBox, boolean persistent)
    {
        super(persistent);
        this.originalBoundingBox = boundingBox;
        initializeToBoundingBox(croppedBox);
        setCroppedBoxRange();
    }
    
    public void setCroppedBoxRange()
    {
        isCropped = false;

        // rowMin, rowMax, rowInc will be the range and step thru the original data
        // if rowIncRatio is < 1.0, then we are subdividing original data
        // same is true for col range
        rowMin = originalBoundingBox.getRowFromRowNum(rowNumMin);
        rowMax = originalBoundingBox.getRowFromRowNum(rowNumMax);
//        rowMin = StsMath.ceiling((rowNumMin - originalBoundingBox.rowNumMin) / originalBoundingBox.rowNumInc);
//        rowMax = StsMath.floor((rowNumMax - originalBoundingBox.rowNumMin) / originalBoundingBox.rowNumInc);
        float rowIncRatio = getIntegerRatio(rowNumInc, originalBoundingBox.rowNumInc);
        rowInc = Math.max(1, Math.round(rowIncRatio));

        // we are assuming that the largest rowNum is adjusted so that the difference is divisible by rowInc
        if(rowNumMax > rowNumMin)
            while(((rowMax - rowMin) % rowInc) != 0)
            {
                rowMax--;
                rowNumMax -= rowNumInc;
            }
        else
            while(((rowMax - rowMin) % rowInc) != 0)
            {
                rowMin++;
                rowNumMin += rowNumInc;
            }

        nRows = 1 + Math.round((rowMax - rowMin) / rowInc);

        colMin = StsMath.ceiling((colNumMin - originalBoundingBox.colNumMin) / originalBoundingBox.colNumInc);
        colMax = StsMath.floor((colNumMax - originalBoundingBox.colNumMin) / originalBoundingBox.colNumInc);
        float colIncRatio = getIntegerRatio(colNumInc, originalBoundingBox.colNumInc);
        colInc = Math.max(1, Math.round(colIncRatio));

        while(((colMax - colMin) % colInc) != 0)
        {
            colMax--;
            colNumMax -= colNumInc;
        }
        nCols = 1 + Math.round((colMax - colMin) / colInc);

        // Process and Clean the Crop Slices
        // hack for now: if boundingBox.zInc is 0.0, then volume originally has a 0.0 zInc
        // even though we've set the croppedBox zInc to a nonzero on the editRange panel.
        // So if boundingBox.zInc is 0.0, use cropped zRange for slice range calculation.

        if(zInc > 0.0f)
        {
            sliceMin = Math.round((zMin - originalBoundingBox.zMin) / originalBoundingBox.zInc);
            sliceMax = Math.round((zMax - originalBoundingBox.zMin) / originalBoundingBox.zInc);

            sliceInc = Math.round(zInc / originalBoundingBox.zInc);
        }
        else
        {
            sliceMin = 0;
            sliceMax = Math.round((zMax - zMin) / zInc);
            sliceInc = 1;
        }
        nSlices = (sliceMax - sliceMin)/sliceInc + 1;

        if(rowMin > 0)
        {
            isCropped = true;
            yMin = originalBoundingBox.getYCoor(rowMin);
        }
        if(rowMax < originalBoundingBox.nRows - 1)
        {
            isCropped = true;
            yMax = originalBoundingBox.getYCoor(rowMax);
        }
        if(rowIncRatio != 1.0f)
        {
            isCropped = true;
            yInc = rowIncRatio*originalBoundingBox.yInc;
        }
        if(colMin > 0)
        {
            isCropped = true;
            xMin = originalBoundingBox.getXCoor(colMin);
        }
        if(colMax < originalBoundingBox.nCols - 1)
        {
            isCropped = true;
            xMax = originalBoundingBox.getXCoor(colMax);
        }
        if(colInc != 1)
        {
            isCropped = true;
            xInc = colIncRatio*originalBoundingBox.xInc;
        }
        else if(sliceMin > 0 || sliceMax < originalBoundingBox.nSlices - 1 || sliceInc != 1)
            isCropped = true;

        nSamplesPerRow = nCols * nSlices;
        nSamplesPerCol = nRows * nSlices;
        nSamplesPerSlice = nRows * nCols;
    }

    private float getIntegerRatio(float f1, float f2)
    {
        float ratio = f1/f2;
        if(ratio == 1.0f)
            return ratio;
        else if(ratio > 1.0)
            return Math.round(ratio);
        else
        {
            ratio = f2/f1;
            ratio = Math.round(ratio);
            return 1.0f/ratio;
        }
    }

    public void adjustCropBlockRows(int uncropBlockRowMin, int uncropBlockRowMax, int nInputFirstRow)
    {
        if(isCropped)
        {
            blockRowMin = StsMath.intervalRoundUp(uncropBlockRowMin, nInputFirstRow, rowInc);
            blockRowMax = StsMath.intervalRoundDown(uncropBlockRowMax, nInputFirstRow, rowInc);
        }
        else
        {
            blockRowMin = uncropBlockRowMin;
            blockRowMax = uncropBlockRowMax;
        }

        nBlockRows = Math.max(0, (blockRowMax - blockRowMin) / rowInc + 1);
        nBlockSamples = nBlockRows * nSamplesPerRow;
        nBlockColSamples = nBlockRows * nSlices;
        nBlockSliceSamples = nBlockRows * nCols;

        outputRowMin = outputRowMax + 1;
        outputRowMax = outputRowMin + nBlockRows - 1;
    }

    public void setNullClip(boolean value)
    {
        nullClip = value;
    }

    public boolean getNullClip()
    {
        return nullClip;
    }
}
