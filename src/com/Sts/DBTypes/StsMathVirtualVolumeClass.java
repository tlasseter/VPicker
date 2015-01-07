package com.Sts.DBTypes;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.UI.Beans.StsBooleanFieldBean;
import com.Sts.UI.Beans.StsFieldBean;

public class StsMathVirtualVolumeClass extends StsVirtualVolumeClass implements StsSerializable, StsClassSurfaceDisplayable, StsClassCursor3dTextureDisplayable //, StsClassDisplayable
{
    public StsMathVirtualVolumeClass()
    {
        userName = "Virtual Volumes from Math Operations on others";
    }

    public void setIsVisibleOnCursor(boolean isVisibleOnCursor)
    {
        if(this.isVisibleOnCursor == isVisibleOnCursor) return;
        this.isVisibleOnCursor = isVisibleOnCursor;
 //       setDisplayField("isVisibleOnCursor", isVisibleOnCursor);
//        currentModel.win3dDisplayAll();
    }

    public boolean getIsVisibleOnCursor()
	{
		return isVisibleOnCursor;
	}
    
    public void initializeDisplayFields()
	{
		displayFields = new StsFieldBean[]
							 {
							 new StsBooleanFieldBean(this, "isVisibleOnCursor", "On 3D Cursors"),
                             new StsBooleanFieldBean(this, "isPixelMode", "Pixel Display Mode"),
							 new StsBooleanFieldBean(this, "contourColors", "Contoured Seismic Colors")
        };
	}


    public void setContourColors(boolean contour)
    {
        if(this.contourColors == contour)
           return;
        this.contourColors = contour;
        currentModel.win3dDisplayAll();
    }
    
}
