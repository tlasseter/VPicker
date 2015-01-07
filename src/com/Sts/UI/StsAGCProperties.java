package com.Sts.UI;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

public class StsAGCProperties extends StsPanelProperties implements StsSerializable
{
	int windowWidth = 500; //in milliseconds
	boolean applyAGC = false;

	transient StsGroupBox groupBox = null;
	transient StsIntFieldBean windowWidthBean;
	transient StsBooleanFieldBean applyAGCBean;

	transient boolean hasAppliedAGC = false;
	transient public boolean recompute = true;
//	transient boolean rescaleRequired = true;

	static private final String title = "AGC Properties";

	public StsAGCProperties()
	{
	}

	public StsAGCProperties(String fieldName)
	{
		this(title, fieldName);
	}

	public StsAGCProperties(String title, String fieldName)
	{
		super(title, fieldName);
	}


    public StsAGCProperties(StsObject parentObject, String title, String fieldName)
	{
        super(parentObject, title, fieldName);
    }

    public StsAGCProperties(StsObject parentObject, StsAGCProperties defaultProperties, String fieldName)
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
		propertyBeans = new StsFieldBean[] { 
				windowWidthBean = new StsIntFieldBean(this, "windowWidth", 1, zMax, "AGC Window Width (ms)", true),
				applyAGCBean    = new StsBooleanFieldBean(this, "applyAGC", "Apply AGC") 
				};
	}

	public int getWindowWidth(byte dataType) { return windowWidth; }
	/** this is just for fieldBean use - use getWindowWidth(dataType) instead!! */
	public int getWindowWidth() { return windowWidth; }
	public boolean getApplyAGC(byte dataType) { return applyAGC; }
	/** this is just for fieldBean use - use getApplyAGC(dataType) instead!! */
	public boolean getApplyAGC() { return applyAGC; }

	public void setWindowWidth(int windowWidth)
	{ 
	    this.windowWidth = windowWidth;
	}
	
	public void setApplyAGC(boolean apply) 
	{ 
	    applyAGC = apply;
	}

	public void setHasAppliedAGC(boolean applied)
	{ 
	    hasAppliedAGC = applied; 
	}
}
