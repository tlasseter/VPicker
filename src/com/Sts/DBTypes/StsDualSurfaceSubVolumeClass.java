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

public class StsDualSurfaceSubVolumeClass extends StsSubVolumeClass implements StsSerializable
{
    public StsDualSurfaceSubVolumeClass()
    {
        userName = "SubVolumes between 2 surfaces";                                               
    }

	public void initializeDisplayFields()
	{
		displayFields = new StsFieldBean[]
		{
		    new StsBooleanFieldBean(this, "isVisible", "Visible"),
		    new StsBooleanFieldBean(this, "isApplied", "Applied")
		};
	}
}
