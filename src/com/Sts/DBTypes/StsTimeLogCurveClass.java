package com.Sts.DBTypes;

import com.Sts.DB.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Feb 13, 2009
 * Time: 2:20:26 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsTimeLogCurveClass extends StsClass implements StsSerializable // , StsClassTimeDisplayable
{
    protected boolean isPixelMode = false; // Blended Pixels or Nearest
	protected boolean contourColors = true; // shader
    boolean enableTime = true;

    public StsTimeLogCurveClass()
    {
        userName = "Time Dependent Logs";
    }

   public void displayTimeClass(StsGLPanel3d glPanel3d, long time)
   {
       Iterator iter = getVisibleObjectIterator();
       while(iter.hasNext())
       {
           StsTimeLogCurve dataset = (StsTimeLogCurve)iter.next();
           if((enableTime && dataset.isAlive(time)) || (!enableTime))
              dataset.display(glPanel3d);
       }
   }

	public void setContourColors(boolean contourColors)
	{
		if (this.contourColors == contourColors) return;
		this.contourColors = contourColors;
		currentModel.win3dDisplayAll();
	}

	public boolean getContourColors()
	{
		return contourColors;
	}

	public void setIsPixelMode(boolean isPixelMode)
	{
		if (this.isPixelMode == isPixelMode)return;
		this.isPixelMode = isPixelMode;
		currentModel.win3dDisplayAll();
	}

	public boolean getIsPixelMode()
	{
		return isPixelMode;
	}

    public void initializeDisplayFields()
	{
        StsLine.initColors();

        displayFields = new StsFieldBean[]
		{
			new StsBooleanFieldBean(this, "isPixelMode", "Pixel Display Mode"),
			new StsBooleanFieldBean(this, "contourColors", "Contoured Prestack Colors")
		};
    }
}

