package com.Sts.DBTypes;

import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;

import java.util.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public class StsSeismicLineSetClass extends StsSeismicClass implements StsSerializable, StsClassSurfaceDisplayable, StsClassTextureDisplayable
{
	public StsSeismicLineSetClass()
	{
		super();
        userName = "Sets of Seismic Lines";
	}

	public void initializeDisplayFields()
	{
		displayFields = new StsFieldBean[]
        {
            new StsBooleanFieldBean(this, "isPixelMode", "Pixel Display Mode"),
            new StsBooleanFieldBean(this, "displayWiggles", "Wiggle Traces"),
            new StsBooleanFieldBean(this, "contourColors", "Contoured Prestack Colors"),
            new StsIntFieldBean(this, "wiggleToPixelRatio", 1, 100, "Wiggle to Pixel Ratio:")
		};
	}

	public boolean setCurrentObject(StsObject object)
	{
		StsSeismic oldCurrentVolume = (StsSeismic)currentObject;
		if (!super.setCurrentObject(object)) return false;
		if (oldCurrentVolume != null) oldCurrentVolume.deleteTransients();
		return true;
	}

    public void displayClass(StsGLPanel3d glPanel3d)
    {
        Iterator iter = getObjectIterator();
        while(iter.hasNext())
        {
            StsSeismicLineSet lineSet = (StsSeismicLineSet)iter.next();
            lineSet.display(glPanel3d);
        }
    }

   public void projectRotationAngleChanged()
   {
       forEach("setLinesXYs");
   }

    public void cropChanged()
    {        
    }

    public boolean textureChanged(StsObject object)
    {
        return ((StsSeismicLineSet)object).textureChanged();
    }
}
