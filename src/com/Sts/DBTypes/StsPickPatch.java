package com.Sts.DBTypes;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version beta 1.0
 */

import com.Sts.Actions.Wizards.Horpick.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

import java.awt.*;

public class StsPickPatch extends StsMainObject
{
    /** seed point */
    protected StsGridPoint seedPoint = null;
    /** horizon this patch is part of */
    protected StsHorpick horpick;
    /** type of pick to make: see pickTypes */
    public byte pickType;
    /** window size for correlation around pick point */
    public int windowSize = 21;
    /** max allowable dif in sample spacing between picks on adjacent traces */
    public float maxPickDif = 1.0f;
    /** color of seedPoint and patch */
    public StsColor stsColor;
    /** number used as name */
//    public int numberName;
    /** range of correlations around seed point */
    public float correlTestMin = 0.0f, correlTestMax = 1.0f;
    /** avg distance between picks */
    public float avgTestPickDif;
    /** min wavelength for test area */
    public float minWaveLength;
    /** max wavelength for test area */
    public float maxWaveLength = 0;
    /** max distance between picks found on test area */
    public float maxTestPickDif;
    /** number of rows and cols for surface */
    transient int nRows, nCols;
    /** number of traces allocated for this patch */
    transient int nPatchTraces = 0;

    /** color of seed points and resulting patch */
//    static protected int spectrumColorIndex = -1;

    public StsPickPatch()
    {
    }

    private StsPickPatch(StsHorpick horpick, byte pickType, int windowSize, float maxPickDif)
    {
        this.horpick = horpick;
        this.pickType = pickType;
        this.windowSize = windowSize;
        this.maxPickDif = maxPickDif;
		stsColor = currentModel.getSpectrum("Basic").getColor(getIndex());
//        numberName = ++spectrumColorIndex;
//        stsColor = currentModel.getSpectrum("Basic").getColor(++spectrumColorIndex);
    }

    static public StsPickPatch constructor(StsHorpick horpick, byte pickType, StsHorpickWizard wizard)
    {
        try
        {
            int windowSize = wizard.getWindowSize();
            float maxPickDif = wizard.getMaxPickDif();
            return new StsPickPatch(horpick, pickType, windowSize, maxPickDif);
        }
        catch(Exception e)
        {
            StsException.outputException("StsPickPatch.constructor() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public boolean initialize(StsModel model)
    {
//        spectrumColorIndex = Math.max(spectrumColorIndex, numberName);
        return true;
    }
/*
	static public void initializeNumber()
	{
		spectrumColorIndex = -1;
	}
*/
    public void createGridPoint(StsCursorPoint cursorPoint, StsSurface surface)
    {
        seedPoint = new StsGridPoint(cursorPoint.point, surface);
        seedPoint.setNullType(StsGridPoint.SURF_PNT);
    }

    public StsGridPoint getSeedPoint() { return seedPoint; }
	public void setGridPoint(StsGridPoint gridPoint) { this.seedPoint = gridPoint; }
/*
    public void setStsColor(Color color)
    {
        Color currentColor = getColor();
        if(currentColor != null && currentColor.equals(color)) return;
        stsColor.setBeachballColors(color);
        StsChangeCmd cmd = new StsChangeCmd(this, stsColor, "stsColor", false);
        currentModel.getActionManager().getCreateTransactionAddCmd("colorChange", cmd, currentModel);
        currentModel.win3dDisplayAll();
    }
*/
    public Color getColor()
    {
        if(stsColor == null) return null;
        return stsColor.getColor();
    }

    public String getName() { return Integer.toString(getIndex()); }

    public StsColor getStsColor()
    {
        if(stsColor == null) return null;
        return stsColor;
    }
/*
    public void initializeGrid(StsSurface surface)
    {
        nRows = surface.getNRows();
        nCols = surface.getNCols();
        if(pointsZ == null) pointsZ = new float[nRows][nCols];
        if(pointsNull == null)
        {
            pointsNull = new byte[nRows][nCols];
            setPointsNull(StsSurface.NULL_GAP_NOT_FILLED);
        }
    }

    public void deleteGrid()
    {
        pointsZ = null;
        pointsNull = null;
        corCoefs = null;
    }

    public float[][] getCorCoefs(int nRows, int nCols)
    {
        if(corCoefs == null) corCoefs = new float[nRows][nCols];
        return corCoefs;
    }

    public void setPointsNull(byte nullType)
    {
		if(pointsNull == null) return;
        for(int row = 0; row < nRows; row++)
            for(int col = 0; col < nCols; col++)
                pointsNull[row][col] = nullType;
    }
*/
    public void displaySeed(StsGLPanel3d glPanel3d)
    {
        if(seedPoint == null || seedPoint.getPoint() == null) return;
        StsGLDraw.drawPoint(seedPoint.getXYZorT(), StsColor.BLACK, glPanel3d, 8, 1.0);
        StsGLDraw.drawPoint(seedPoint.getXYZorT(), stsColor, glPanel3d, 4, 2.0);
    }

    public void setPickZ(float pickZ)
    {
        if(seedPoint == null) return;
        seedPoint.setZ(pickZ);
        seedPoint.setNullType(StsSurface.SURF_PNT);
        horpick.getSurface().setGridPoint(seedPoint);
    }

    public float getPickZ()
    {
        if(seedPoint == null) return StsParameters.nullValue;
        return seedPoint.getZorT();
    }

    public void setPickType(byte pickType) { this.pickType = pickType; }
    public byte getPickType() { return pickType; }
/*
    public void setPickTypeName(String pickTypeName)
    {
        for(byte n = 0; n < 4; n++)
        {
            if(pickTypeName.equals(pickTypeNames[n]))
            {
                pickType = n;
                break;
            }
        }
    }
    public String getPickTypeName() { return pickTypeNames[pickType]; }

    public void setStopCriteriaName(String stopCriteriaName)
    {
        for(byte n = 0; n < 4; n++)
        {
            if(stopCriteriaName.equals(stopCriteriaNames[n]))
            {
                stopCriteria = n;
                break;
            }
        }
    }
    public String getStopCriteriaName() { return stopCriteriaNames[stopCriteria]; }
*/
    public void setWindowSize(int windowSize) { this.windowSize = windowSize; }
    public int getWindowSize() { return windowSize; }

    public float getMaxPickDif() { return maxPickDif; }
    public void setMaxPickDif(float maxPickDif) { this.maxPickDif = maxPickDif; }
//    public byte getStopCriteria() { return stopCriteria; }
//    public void setStopCriteria(byte stopCriteria) { this.stopCriteria = stopCriteria; }

    public StsHorpick getHorpick() { return horpick; }

//    public byte[][] getPointsNull() { return pointsNull; }
//    public float[][] getPointsZ() { return pointsZ; }
//    public void setPointsZ(float[][] pointsZ) { this.pointsZ = pointsZ; }

    public int getNPatchTraces() { return nPatchTraces; }
    public void setNPatchTraces(int nPatchTraces) { this.nPatchTraces = nPatchTraces; }
//    public int getNumberName() { return numberName; }
//    public void setNumberName(int numberName) { this.numberName = numberName; }
}
