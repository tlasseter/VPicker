package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.PreStack.*;
import com.Sts.Utilities.Interpolation.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.nio.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Feb 2, 2009
 * Time: 11:28:38 AM
 * To change this template use File | Settings | File Templates.
 */

/** sources and receivers can be shown in two ways: location, where they are displayed as points of a particular color; and
 *  as a texture showing interpolated elevation.  When shown as a location, the instance is called as an StsDisplayable (instance
 *  of a class which can be displayed).  As a texture, it is displayed on the z-slice StsCursorSection and is displayed generally
 *  if no other texture is being displayed there.
 */
public abstract class StsSeismicSourcesReceivers extends StsDisplayPreStackAttribute implements StsInstance3dDisplayable
{
    StsPreStackLine line = null;
    protected String xString;
    protected String yString;
    protected StsColor normalColor;
    double[] x;
    double[] y;
    int[] colorIndices;
    StsColor[] colors;
    int nGather = -1;
    int nFirstGatherTrace;
    int nLastGatherTrace;
    int nGatherTraces;

    String currentCulture = StsPreStackLineSet.DISPLAY_CULTURE_NONE;

    public StsSeismicSourcesReceivers(StsPreStackLineSet lineSet, StsModel model)
    {
        super(lineSet, model);
    }

    protected boolean computeBytes(String elevStg, String xStg, String yStg)
    {
        int nRows = lineSet.nRows;
        int nCols = lineSet.nCols;
        StsPreStackLine[] lines = lineSet.lines;
        double max = -StsParameters.largeDouble;
        double min = StsParameters.largeDouble;
        float[][] elevation = new float[nRows][nCols];
        double[] elev, x, y;
        float[] xy ;
        int rCol, rRow;

        try
        {
            StsSpiralRadialInterpolation interpolation = lineSet.getSpiralRadialInterpolation(1, 10);
            interpolation.initialize();

            for (int row = 0; row < nRows; row++)
            {
                for (int col = 0; col < nCols; col++)
                {
                    elevation[row][col] = StsParameters.nullValue;
                }
            }

            for (int row = 0; row < nRows; row++)
            {
                StsPreStackLine line = lines[row];
                if (line == null)
                {
                    continue;
                }

                elev = line.getAttributeArray(elevStg);
                x = line.getAttributeArray(xStg);
                y = line.getAttributeArray(yStg);
                /*
                   public float[] getRotatedRelativeXYFromUnrotatedAbsoluteXY(double x, double y)
                   {
                 return this.getRotatedRelativeXYfromUnrotatedRelativeXY((float)(x - xOrigin), (float)(y - yOrigin));
                  }
                 */
                for (int k = 0; k < elev.length; k++)
                {
                    xy = model.getProject().getRotatedRelativeXYFromUnrotatedAbsoluteXY(x[k], y[k]);
                    rCol = (int) lineSet.getColCoor(xy[0]);
                    rRow = (int) lineSet.getRowCoor(xy[1]);
                    if ( (rCol >= 0) && (rCol < nCols) &&
                        (rRow >= 0) && (rRow < nRows))
                    {
                        elevation[rRow][rCol] = (float) elev[k];
                        if (max < elev[k])
                        {
                            max = elev[k];
                        }
                        if (min > elev[k])
                        {
                            min = elev[k];
                        }
                    }
                }

            }
            if (max == 0)
            {
                return false;
            }
            int nDataPoints = 0;
            for (int row = 0; row < nRows; row++)
            {
                for (int col = 0; col < nCols; col++)
                {
                    if (elevation[row][col] != StsParameters.nullValue)
                    {
                        interpolation.addDataPoint(new Float(elevation[row][col]), row, col);
                        nDataPoints++;
                    }
                }
            }
            System.out.println("nDataPoints: " + nDataPoints + " out of " + (nRows + 1) * (nCols + 1));

            interpolation.run();
            byte[] byteElevation = new byte[nRows * nCols];
            float scale = 254.0f / (float) max;
            int i = 0;
            for (int row = 0; row < nRows; row++)
            {
                for (int col = 0; col < nCols; col++, i++)
                {
                    StsSpiralRadialInterpolation.Weights dataWeights = interpolation.getWeights(row, col);
                    if (dataWeights == null || dataWeights.nWeights == 0)
                    {
                        byteElevation[i] = StsParameters.nullByte;
                    }
                    else
                    {
                        int nWeights = dataWeights.nWeights;
                        double[] weights = dataWeights.weights;
                        Object[] dataObjects = dataWeights.dataObjects;
                        double value = 0.0;
                        double sumWeight = 0.0;
                        for (int n = 0; n < nWeights; n++)
                        {
                            float z = ( (Float) dataObjects[n]).floatValue();
                            if (z != StsParameters.largeFloat)
                            {
                                value += z * weights[n];
                                sumWeight += weights[n];
                            }
                        }
                        if (sumWeight != 0.0)
                        {
                            byteElevation[i] = StsMath.unsignedIntToUnsignedByte254( (int) (scale * value / sumWeight));
                        }
                        else
                        {
                            byteElevation[i] = StsParameters.nullByte;
                        }
                    }
                }
            }
            if(colorscale == null)
            {
                getInitializeColorscale(StsSpectrumClass.SPECTRUM_RAINBOW, 0.0f, (float) max);
            }
            
            byteBuffer = ByteBuffer.wrap(byteElevation);
            return true;
        }
        catch (Exception e)
        {
            StsMessageFiles.errorMessage("Failed to build elevation map for " + elevStg);
            return false;
        }
    }

    public void display(StsGLPanel3d glPanel3d)
    {
        StsPreStackLineSetClass lineSetClass = (StsPreStackLineSetClass)lineSet.getStsClass();
        int size = lineSetClass.getDisplayPointSize();
        StsPreStackLine currentLine = lineSet.currentLine;
        if (currentLine == null)
        {
            return;
        }
        if (currentLine != line)
        {
            line = currentLine;
            currentCulture = StsPreStackLineSet.DISPLAY_CULTURE_NONE;
            try
            {
                x = currentLine.getAttributeArray(xString);
                y = currentLine.getAttributeArray(yString);
            }
            catch (Exception e)
            {
                StsMessageFiles.errorMessage("Failed to find receiver locations for line " + currentLine.getName());
                return;
            }
        }
        StsSuperGather superGather = lineSet.getSuperGather(glPanel3d.window);
        nFirstGatherTrace = superGather.centerGather.nFirstGatherLineTrace;
        nLastGatherTrace = superGather.centerGather.nLastGatherLineTrace;
        nGatherTraces = superGather.centerGather.nGatherTraces;
//            System.out.println("Plotting sources: First trace=" + nFirstGatherTrace + " Last Trace=" + nLastGatherTrace);
//        checkComputeCultureColorIndices();

        GL gl = glPanel3d.getGL();
        float[] point = new float[3];
        float[] xy;

        model.getGlPanel3d().setViewShift(gl, 2.0f);
        point[2] = model.getCursor3d().getCurrentDirCoordinate(StsCursor3d.ZDIR);
        if ( (nFirstGatherTrace == 0) && (nLastGatherTrace == 0))
        {
            return;
        }
        for (int i = 0, n = nFirstGatherTrace; n <= nLastGatherTrace; i++, n++)
        {
            xy = model.getProject().getRotatedRelativeXYFromUnrotatedAbsoluteXY(x[n], y[n]);
            point[0] = xy[0];
            point[1] = xy[1];
            StsGLDraw.drawPoint(point, StsColor.BLACK, model.getGlPanel3d(), size + 2);
            if (currentCulture != StsPreStackLineSet.DISPLAY_CULTURE_NONE)
            {
                StsGLDraw.drawPoint(point, colors[colorIndices[i]], model.getGlPanel3d(), size, 2.0);
            }
            else
            {
                StsGLDraw.drawPoint(point, normalColor, model.getGlPanel3d(), size, 2.0);
            }
        }
        model.getGlPanel3d().resetViewShift(gl);
    }
/*
    private boolean checkComputeCultureColorIndices()
    {
        if (currentCulture == displayCulture)
        {
            return true;
        }
        currentCulture = displayCulture;
        if (currentCulture == StsPreStackLineSet.DISPLAY_CULTURE_NONE)
        {
            colors = null;
            colorIndices = null;
            return true;
        }
        if (colors == null)
        {
            colors = currentModel.getSpectrum(StsSpectrumClass.SPECTRUM_RWB).getStsColors();
        }
        try
        {
            double[] range = new double[2];
            range[0] = StsParameters.largeDouble;
            range[1] = -StsParameters.largeDouble;
            double[] attrVal = currentLine.getAttributeArray(currentCulture);
            for (int i = 0; i < attrVal.length; i++)
            {
                if (range[0] > attrVal[i])
                {
                    range[0] = attrVal[i];
                }
                if (range[1] < attrVal[i])
                {
                    range[1] = attrVal[i];
                }
            }
            colorIndices = initializeAttributeColors(attrVal, range);
            return true;
        }
        catch (Exception e)
        {
            new StsMessage(currentModel.win3d, StsMessage.WARNING, "Failed to find attribute " + currentCulture);
            displayCulture = StsPreStackLineSet.DISPLAY_CULTURE_NONE;
            currentCulture = StsPreStackLineSet.DISPLAY_CULTURE_NONE;
            colors = null;
            colorIndices = null;
            return false;
        }
    }

    private int[] initializeAttributeColors(double[] attrVal, double[] range)
    {
        double min = range[0];
        double max = range[1];
        int[] indices = new int[maxNTracesPerGather];
        int nColors = colors.length;
        float binIncrement = (float) (max - min) / nColors;
        for (int i = 0, n = nFirstGatherTrace; n <= nLastGatherTrace; i++, n++)
        {
            indices[i] = (int) ( (float) (attrVal[n] - min) / binIncrement);
            indices[i] = StsMath.minMax(indices[i], 0, nColors - 1);
        }
        return indices;
    }
*/
}
