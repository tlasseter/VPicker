package com.Sts.DBTypes;

import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;
import java.awt.event.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: May 27, 2009
 * Time: 8:10:32 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsSurfaceCurvatureAttribute extends StsSurfaceAttribute
{
	private static final long serialVersionUID = 1L;
	public byte curveType = CURVPos;
    public int filterSize = 5;

    transient StsColorscalePanel colorscalePanel = null;
    transient private StsProgressPanel progressPanel = null;

    static final float largeFloat = StsParameters.largeFloat;
    static final float nullValue = StsParameters.nullValue;

    static private boolean debug = false;
    static private boolean runTimer = false;
    static private StsTimer timer;

    static final byte FILTER_ON_CHI_SQ = 1;
    static final byte FILTER_ON_STD_DEV = 2;
    static final byte filterType = FILTER_ON_STD_DEV;

    // Multiply # pts in SVD to get ChiSqr limit
    static final private double chiSqrMultiplyer = 2;
    static final private double stdDevFactor = 2.0;

    static public final float badCurvature =  StsQuadraticCurvature.badCurvature;
    static public final float curvatureTest = StsQuadraticCurvature.curvatureTest;

    //Curvature Attribute Types
    static public final byte CURVDip 	= StsQuadraticCurvature.CURVDip;
	static public final byte CURVStrike = StsQuadraticCurvature.CURVStrike;
	static public final byte CURVMean 	= StsQuadraticCurvature.CURVMean;
	static public final byte CURVGauss  = StsQuadraticCurvature.CURVGauss;
	static public final byte CURVMax	= StsQuadraticCurvature.CURVMax;
	static public final byte CURVMin  	= StsQuadraticCurvature.CURVMin;
	static public final byte CURVPos  	= StsQuadraticCurvature.CURVPos;
	static public final byte CURVNeg  	= StsQuadraticCurvature.CURVNeg;

    public static final String CURVDipString 	= StsQuadraticCurvature.CURVDipString;
    public static final String CURVStrikeString = StsQuadraticCurvature.CURVStrikeString;
    public static final String CURVMeanString 	= StsQuadraticCurvature.CURVMeanString;
    public static final String CURVGaussString 	= StsQuadraticCurvature.CURVGaussString;
    public static final String CURVPosString 	= StsQuadraticCurvature.CURVPosString;
    public static final String CURVNegString 	= StsQuadraticCurvature.CURVNegString;
    public static final String CURVMinString 	= StsQuadraticCurvature.CURVMinString;
    public static final String CURVMaxString 	= StsQuadraticCurvature.CURVMaxString;

    public static final String CURVDipName	 	= StsQuadraticCurvature.CURVDipName;
    public static final String CURVStrikeName 	= StsQuadraticCurvature.CURVStrikeName;
    public static final String CURVMeanName 	= StsQuadraticCurvature.CURVMeanName;
    public static final String CURVGaussName 	= StsQuadraticCurvature.CURVGaussName;
    public static final String CURVPosName 		= StsQuadraticCurvature.CURVPosName;
    public static final String CURVNegName 		= StsQuadraticCurvature.CURVNegName;
    public static final String CURVMinName 		= StsQuadraticCurvature.CURVMinName;
    public static final String CURVMaxName 		= StsQuadraticCurvature.CURVMaxName;
    public static String[] CURV_ATTRIBUTE_STRINGS = {CURVDipString, CURVStrikeString, CURVMeanString, CURVGaussString, 
        CURVMaxString, CURVMinString, CURVPosString, CURVNegString };
    public static String[] CURV_ATTRIBUTE_NAMES = {CURVDipName, CURVStrikeName, CURVMeanName, CURVGaussName, 
        CURVMaxName, CURVMinName, CURVPosName, CURVNegName };
    public StsSurfaceCurvatureAttribute()
    {
    }
    
    public StsSurfaceCurvatureAttribute(StsSurface surface, byte curveType, int filterSize)
    {
        super(false);
        this.surface = surface;
        this.curveType = curveType;
        this.filterSize = filterSize;
        name = surface.getName() + "-" + CURV_ATTRIBUTE_STRINGS[curveType] + filterSize;
    }

    public String getAttributeString() { return CURV_ATTRIBUTE_STRINGS[curveType]; }

    public String toString() { return getName(); }

    public boolean createAttribute(StsProgressPanel progressPanel)
    {
        this.progressPanel = progressPanel;
        return createAttribute();
    }

    public boolean createAttribute()
    {
        int nRows = surface.nRows;
        int nCols = surface.nCols;
        byte[][] pointsNull = surface.pointsNull;
        float[][] pointsZ = surface.pointsZ;

        float[][] values = new float[nRows][nCols];
        byte[] textureData = new byte[nRows * nCols];
        if (progressPanel != null)
            progressPanel.initialize(nRows);

        dataMin = largeFloat;
        dataMax = -largeFloat;

        int winRows = filterSize;
        int winCols = filterSize;
        int halfWinRows = winRows/2;
        int halfWinCols = winCols/2;
        int num = 0;
        double sum = 0, sumSqr = 0;
        double mean;
        // Determine quadratic coefficients for this neighborhood

        if(runTimer)
        {
            timer = new StsTimer();
            timer.start();
        }

        float[][] fitPoints = new float[winRows*winCols][3];
        for (int row = 0; row < nRows; row++)
        {
            for (int col = 0; col < nCols; col++)
            {
            	int nFitPoints = 0;  // number of equations
                if (pointsNull[row][col] == StsSurface.SURF_PNT)
                {
                    float zc = pointsZ[row][col];
                    int col0 = (col-halfWinCols > 0) ? (col-halfWinCols) : 0;
                    int col1 = (col+halfWinCols < nCols) ? (col+halfWinCols) : nCols-1;
                    int row0 = (row-halfWinRows > 0) ? (row-halfWinRows) : 0;
                    int row1 = (row+halfWinRows < nRows) ? (row+halfWinRows) : nRows-1;
                    for (int winRow = row0; winRow <= row1; winRow++)
                    {
                        for (int winCol = col0; winCol <= col1; winCol++)
                        {
                            if (pointsNull[winRow][winCol] == StsSurface.SURF_PNT)
                            {
		                    	fitPoints[nFitPoints][0] = (winCol-col) * surface.xInc;
		                    	fitPoints[nFitPoints][1] = (winRow-row) * surface.yInc;
		                    	fitPoints[nFitPoints][2] = pointsZ[winRow][winCol]-zc;
		                    	nFitPoints++;
                            }
                        }
                    }
                    if (nFitPoints < StsQuadraticCurvature.nCoefs)
                    {
                        values[row][col] = nullValue;
                        continue;
                    }

                    if(!StsQuadraticCurvature.computeSVD(fitPoints, nFitPoints))continue;
                    
                    float val;
                    try
                    {
                        val = StsQuadraticCurvature.getCurvatureComponent(curveType);
                    }
                    catch(Exception e)
                    {
                        StsException.systemError(this, "computeCurvature", "getCurvatureComponent failed.");
                        continue;
                    }
                    
                    //if (filterType == FILTER_ON_CHI_SQ)
                    //Always do chiSqr test
                    {
                    	double chiSqrTest = chiSqrMultiplyer * nFitPoints;
                        double chiSqr = StsQuadraticCurvature.computeChiSquared();
                        if(chiSqr > chiSqrTest)
                        {
                            // if(StsPatchVolume.debug) StsException.systemDebug(this, "computeCurvature", "ChiSqr = " + chiSqr + " at row, col " + row + " " + col);
                            //continue;
                        	if (val > 0) val = badCurvature;
                        	if (val < 0) val = -badCurvature;
                        }
                    } 
                    values[row][col] = val;
                    if(Math.abs(val) > curvatureTest) continue;
                    // ChiSqr filtered vals not used for dataMin / dataMax & Statistics
                    num++;
                    sum += val;
                    dataMax = Math.max(dataMax, val);
                    dataMin = Math.min(dataMin, val);

                }
            	else
            	{
            		values[row][col] = nullValue;
            	}
            }
            if (progressPanel == null) continue;
            if (progressPanel.isCanceled())
            {
                progressPanel.setDescriptionAndLevel("Cancelled by user.", StsProgressBar.ERROR);
                textureData = null;
                values = null;
                return false;
            }
            progressPanel.setValue(row+1);
        }

        if (runTimer)
            timer.stopPrint("Time to compute curvature for a " + winRows + " x " + winRows + "window.");

        //compute variance
        mean = sum / num;
        num = 0;
       //Set filter at 2.5 * stdDev for now
        if (filterType == FILTER_ON_STD_DEV)
        {
            for (int row = 0; row < nRows; row++)
            {
                for (int col = 0; col < nCols; col++)
                {
                    if (values[row][col] != nullValue && Math.abs(values[row][col]) < curvatureTest)
                    {
                        sumSqr += ((double)values[row][col]-mean)*((double)values[row][col]-mean);
                        num++;
                    }
                }
            }
            double variance = (sumSqr) / num;
            double stdDev = Math.sqrt(variance) * stdDevFactor;
            float maxFilt = (float)(mean + stdDev);
            float minFilt = (float)(mean - stdDev);
            dataMin = Math.max(dataMin, minFilt);
            dataMax = Math.min(dataMax, maxFilt);
            float scale = StsMath.floatToUnsignedByteScale(dataMin, dataMax);
            float offset = StsMath.floatToUnsignedByteScaleOffset(scale, dataMin);

            int i = 0;
            for (int row = 0; row < nRows; row++)
            {
                for (int col = 0; col < nCols; col++)
                {
                    float value = values[row][col];
                    if (value != nullValue)
                    {
                        if (value > maxFilt) value = maxFilt;
                        else if (value < minFilt) value = minFilt;
                        textureData[i++] = StsMath.floatToUnsignedByte254WithScale(value, scale, offset);
                    }
                    else
                    {
                        textureData[i++] = -1;
                    }
                }
            }
        }
        else
        {
            float scale = StsMath.floatToUnsignedByteScale(dataMin, dataMax);
            float offset = StsMath.floatToUnsignedByteScaleOffset(scale, dataMin);
            int i = 0;
            for (int row = 0; row < nRows; row++)
            {
                for (int col = 0; col < nCols; col++)
                {
                    float value = values[row][col] ;
                    if (value != nullValue)
                        textureData[i++] = StsMath.floatToUnsignedByte254WithScale(value, scale, offset);
                    else
                        textureData[i++] = -1;
                }
            }
        }
        if (progressPanel != null)
            progressPanel.finished();
        initialize(dataMin, dataMax, values, textureData);
        surfaceTexture.setTextureData(textureData);
        progressPanel.appendLine("Attribute " + name + " completed successfully.");
        return true;
    }


    public StsColorscalePanel getColorscalePanel()
    {
        if (colorscalePanel != null) return colorscalePanel;

        StsColorscale colorscale = getColorscale();
        colorscalePanel = new StsColorscalePanel(colorscale, StsParameters.HORIZONTAL, false, false);
        colorscalePanel.setMinimumSize(new Dimension(200, 350));
        colorscalePanel.setPreferredSize(new Dimension(200, 350));
        colorscalePanel.addActionListener(new ActionListener()
        {
            public void actionPerformed(ActionEvent e)
            {
                curvColorscalePanelChanged();
            }
        });
        return colorscalePanel;
    }

    protected void curvColorscalePanelChanged()
    {
        if (surface == null)
        {
            return;
        }
        surface.setTextureChanged();
        surface.colorListChanged();
        StsSerialize.currentModel.win3dDisplay();
    }

    public float[] getHistogram()
    {
        if(surfaceTexture == null) return null;
        byte[] textureData = surfaceTexture.getTextureData();
        return StsToolkit.buildHistogram(textureData, dataMin, dataMax);
    }

    protected void attributeSelected()
    {
        surface.curvatureAttribute = this;
    }
}
