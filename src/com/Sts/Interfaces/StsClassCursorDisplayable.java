package com.Sts.Interfaces;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;

public interface StsClassCursorDisplayable
{
    public void drawOnCursor3d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, StsPoint[] planePoints, boolean isDragging);
    public void drawOnCursor2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate,
                               boolean axesFlipped, boolean xAxisReversed, boolean axisReversed);
    public void drawOn3dCurtain(StsGLPanel3d glPanel3d, StsSeismicCurtain seismicCurtain);
}