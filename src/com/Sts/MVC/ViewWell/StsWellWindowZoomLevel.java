package com.Sts.MVC.ViewWell;

/**
 * Title:        Sts Well Viewer Description:  Well Model-Viewer Copyright:    Copyright (c) 2001
 * Company:      4D Systems LLC
 * @author T.J.Lasseter
 * @version 1.0
 */
public class StsWellWindowZoomLevel
{
    static public StsWellWindowZoomLevel[] zoomLevels;
    static public int nZoomLevels;
    static int pixelsPerTick = 10;
    //static double[] unitsPerTickArray = new double[]
    //{
     //   1000.0, 500.0, 250.0, 100.0, 50.0, 25.0, 10.0, 5.0, 2.0, 1.0, 0.5
    //};
	static double[] unitsPerTickArray = new double[]
	{
		500.0, 400., 300., 200., 100., 50., 30., 15, 7. , 3.
	};
    static int ticksPerLabel = 5;

    static
    {
        nZoomLevels = unitsPerTickArray.length;
        zoomLevels = new StsWellWindowZoomLevel[nZoomLevels];
        for (int n = 0; n < nZoomLevels; n++)
            zoomLevels[n] = new StsWellWindowZoomLevel(unitsPerTickArray[n], pixelsPerTick, ticksPerLabel);
    }

    public double unitsPerTick;
    public double pixelsPerUnit;
    public double unitsPerPixel;
    public double unitsPerLabel;

    public StsWellWindowZoomLevel(double unitsPerTick, int pixelsPerTick, int ticksPerLabel)
    {
        this.unitsPerTick = unitsPerTick;
        this.pixelsPerTick = pixelsPerTick;
        this.ticksPerLabel = ticksPerLabel;
        pixelsPerUnit = pixelsPerTick / unitsPerTick;
        unitsPerPixel = unitsPerTick / pixelsPerTick;
        unitsPerLabel = unitsPerTick * ticksPerLabel;
    }
}
