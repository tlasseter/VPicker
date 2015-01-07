package com.Sts.UI;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

public class StsAGCPreStackProperties extends StsAGCProperties implements StsSerializable
{
	boolean applyAGCPoststack = false;
	boolean hasAppliedAGCPoststack = false;
	int postStackWindowWidth = 1500; //in milliseconds

	transient StsBooleanFieldBean applyAGCPoststackBean;
	transient StsIntFieldBean postStackWindowWidthBean;

	static private final String title = "AGC Properties";

	public StsAGCPreStackProperties()
	{
	}

	public StsAGCPreStackProperties(String fieldName)
	{
		super(title, fieldName);
	}

	public StsAGCPreStackProperties(StsObject parentObject, StsAGCPreStackProperties defaultProperties, String fieldName)
	{
            super(parentObject, title, fieldName);
            initializeDefaultProperties(defaultProperties);
    }
    public void initializeDefaultProperties(Object defaultProperties)
    {
        if(defaultProperties != null)
            StsToolkit.copyObjectNonTransientFields(defaultProperties, this);
    }
	public void initializeBeans()
	{
		int zMax = (int) (StsPreStackLineSetClass.currentProjectPreStackLineSet.zMax - 1);
		propertyBeans = new StsFieldBean[]
		{
            applyAGCBean = new StsBooleanFieldBean(this, "applyAGC", "Apply AGC Prestack"),
            windowWidthBean = new StsIntFieldBean(this, "windowWidth", 1, zMax, "Prestack AGC Window Width (ms)", true),
            applyAGCPoststackBean = new StsBooleanFieldBean(this, "applyAGCPoststack", "Apply AGC Poststack"),
            postStackWindowWidthBean = new StsIntFieldBean(this, "postStackWindowWidth", 1, zMax, "Poststack AGC Window Width (ms)", true)
		};
	}

	public boolean getApplyAGCPoststack() { return applyAGCPoststack; }
	public void setApplyAGCPoststack(boolean apply)
	{ 
	    applyAGCPoststack = apply; 
	}
	public void setHasAppliedAGCPoststack(boolean applied)
	{ 
	    hasAppliedAGCPoststack = applied; 
	}
	public int getPostStackWindowWidth() { return postStackWindowWidth; }
	public void setPostStackWindowWidth(int postStackWindowWidth) 
	{
	    this.postStackWindowWidth = postStackWindowWidth;
	}
	
	public boolean getApplyAGC(byte dataType)
	{
		if (dataType == StsFilterProperties.PRESTACK)
		{
			return applyAGC;
		} else
		{
			return applyAGCPoststack;
		}
	}
	
	public int getWindowWidth(byte dataType)
	{
		if (dataType == StsFilterProperties.PRESTACK)
		{
			return windowWidth;
		} else
		{
			return postStackWindowWidth;
		}
	}
}
