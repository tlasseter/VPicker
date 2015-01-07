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

/** An interface representing instances of this class that can be displayed on the 3d cursors */
public interface StsClassCursor3dDisplayable
{
    public StsObject getCurrentObject();
    public StsCursor3dDisplayable constructDisplayableSection(StsCursor3dTexture cursor3dTexture);
    public boolean drawLast();
}
