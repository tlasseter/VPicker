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

public class StsWellSetSubVolumeClass extends StsSubVolumeClass implements StsSerializable, StsClassDisplayable
{
    public StsWellSetSubVolumeClass()
    {
        userName = "SubVolumes Around a Well";                                                                              
    }

	public void initializeDisplayFields()
   {
//	   initColors(StsWellSetSubVolume.displayFields);

	   displayFields = new StsFieldBean[]
	   {
		   new StsBooleanFieldBean(this, "isVisible", "Visible"),
		   new StsBooleanFieldBean(this, "isApplied", "Applied")
	   };
   }

    public void displayClass(StsGLPanel3d glPanel3d)
     {
         if(!isVisible) return;

         Iterator iter = getObjectIterator();
         while(iter.hasNext())
         {
             StsWellSetSubVolume wellSetSubVolume = (StsWellSetSubVolume)iter.next();
             boolean isCurrentObject = (wellSetSubVolume == currentObject);
             wellSetSubVolume.display(glPanel3d, isCurrentObject);
         }
     }

}
