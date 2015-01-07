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
import com.Sts.UI.Beans.*;

public class StsSingleSurfaceSubVolumeClass extends StsSubVolumeClass implements StsSerializable
{
    public StsSingleSurfaceSubVolumeClass()
    {
        userName = "SubVolume from Single Surface";
    }

	public void initializeDisplayFields()
   {
//	   initColors(StsSingleSurfaceSubVolume.displayFields);

	   displayFields = new StsFieldBean[]
	   {
		   new StsBooleanFieldBean(this, "isVisible", "Visible"),
		   new StsBooleanFieldBean(this, "isApplied", "Applied")
	   };
   }
}
