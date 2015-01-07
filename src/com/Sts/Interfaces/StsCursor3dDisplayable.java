package com.Sts.Interfaces;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version beta 1.0
 */

import com.Sts.DBTypes.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;

public interface StsCursor3dDisplayable
{
    public boolean setDirCoordinate(float dirCoordinate, int nPlane);
    public void clearTextureDisplay();
    public void cropChanged();
    public boolean isDisplayableObject(StsObject object);
    public boolean isDisplayingObject(StsObject object);
    public void display(StsGLPanel3d glPanel3d, boolean is3d, byte[] subVolumePlane);
    public String logReadout(StsPoint point);
    public void setObject(StsObject object);
}
