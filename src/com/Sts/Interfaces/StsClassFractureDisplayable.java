package com.Sts.Interfaces;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version beta 1.0
 */

import com.Sts.Types.StsRotatedGridBoundingBox;

import java.util.ArrayList;

/** An interface representing a class whose instances can be displayed on the 3d cursors */
public interface StsClassFractureDisplayable
{
    public ArrayList<StsFractureDisplayable> getDisplayableFractures();
	StsRotatedGridBoundingBox getFractureBoundingBox();
}
