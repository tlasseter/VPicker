package com.Sts.Interfaces;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version beta 1.0
 */

import com.Sts.DBTypes.StsObject;
import com.Sts.MVC.StsModel;
import com.Sts.MVC.View3d.StsCursor3d;
import com.Sts.MVC.View3d.StsCursor3dTexture;
import com.Sts.Types.StsBoundingBox;
import com.Sts.Types.StsGolderFracture;
import com.Sts.Types.StsPoint;

/** An interface representing a class whose instances can be displayed on the 3d cursors */
public interface StsFractureDisplayable
{
    public StsBoundingBox getBoundingBox();
	public boolean getIsVisible();
	public StsPoint[] getSectionIntersectionAtZ(float z);
	public boolean intersects(StsGolderFracture fracture);
}
