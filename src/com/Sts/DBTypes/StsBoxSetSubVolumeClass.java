package com.Sts.DBTypes;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;

import java.util.*;

public class StsBoxSetSubVolumeClass extends StsSubVolumeClass implements StsSerializable, StsClassDisplayable
{
	public StsBoxSetSubVolumeClass()
	{
        userName = "Set of 3D Boxes SubVolume";
	}

	public void initializeDisplayFields()
	{
//		initColors(StsBoxSetSubVolume.displayFields);

		displayFields = new StsFieldBean[]
			{new StsBooleanFieldBean(this, "isVisible", "Visible"), new StsBooleanFieldBean(this, "isApplied", "Applied")
		};
	}

	public void displayClass(StsGLPanel3d glPanel3d)
	{
		if (!isVisible)
		{
			return;
		}

		Iterator iter = getObjectIterator();
		while (iter.hasNext())
		{
			StsBoxSetSubVolume boxSetSubVolume = (StsBoxSetSubVolume) iter.next();
			boolean isCurrentObject = (boxSetSubVolume == currentObject);
			boxSetSubVolume.display(glPanel3d, isCurrentObject);
		}
	}
}
