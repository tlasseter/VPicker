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
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;

import java.util.*;

public class StsMVFractureSetClass extends StsClass implements StsSerializable, StsClassTimeDisplayable, StsClassCursorDisplayable, StsClassTimeSeriesDisplayable
{
    /**
	 * 
	 */
	private static final long serialVersionUID = -2677014498211804642L;

    boolean enableTime = false;
    private String defaultSpectrumName = "RedWhiteBlue";
    private StsColor defaultFractureSetColor = new StsColor(StsColor.BLUE);

    public StsMVFractureSetClass()
    {
        userName = "Sets of Midland Valley Fractures";
    }

   public void initializeFields()
   {
       displayFields = new StsFieldBean[]
       {
           new StsBooleanFieldBean(this, "enableTime", "Enable Time"),
       };

       defaultFields = new StsFieldBean[]
       {
		   new StsColorComboBoxFieldBean(this, "defaultFractureSetColor", "Color:", StsColor.colors32),
		   new StsComboBoxFieldBean(this, "defaultSpectrumName", "Spectrum:", StsSpectrumClass.cannedSpectrums)
       };
   }

    public void displayTimeClass(StsGLPanel3d glPanel3d, long time)
    {
        Iterator<Object> iter = getVisibleObjectIterator();
        while (iter.hasNext())
        {
        	StsMVFractureSet fracSet = (StsMVFractureSet) iter.next();
            //if ((enableTime && sensor.isAlive(time)) || (!enableTime))
        	fracSet.display(glPanel3d);
        }
    }

 
    public void drawOnCursor3d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, StsPoint[] planePoints, boolean isDragging) { }

    public void drawOnCursor2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, boolean axesFlipped, boolean xAxisReversed, boolean yAxisReversed)
    {
//        long time = currentModel.getProject().getProjectTime();
//        Iterator iter = getVisibleObjectIterator();
//        while (iter.hasNext())
//        {
//            StsSensor sensor = (StsSensor) iter.next();
//            if ((enableTime && sensor.isAlive(time)) || (!enableTime))
//                sensor.display2d(glPanel3d, dirNo, dirCoordinate, axesFlipped, xAxisReversed, yAxisReversed);
//        }
    }

 
    public void drawOn3dCurtain(StsGLPanel3d glPanel3d, StsSeismicCurtain seismicCurtain) { }
    // fails if object is null
    // should simply use superclass default method:
    // update Properties should be handled when the object is added to model
    // TJL June 3 08

    /*
       public boolean setCurrentObject(StsObject object)
       {
            if (currentObject == object)
                return false;
            currentObject = object;
            ((StsSensor)object).updateProperties();
            return super.setCurrentObject(object);
       }
   }
   */
//    public void setDisplayNames(boolean displayNames)
//    {
//        if(this.displayNames == displayNames) return;
//        this.displayNames = displayNames;
////      setDisplayField("displayNames", displayNames);
//        currentModel.win3dDisplayAll();
//    }
//
//    public boolean getDisplayNames() {	return displayNames; }
//
    public void setEnableTime(boolean enable)
    {
        if (this.enableTime == enable) return;
        this.enableTime = enable;
//      setDisplayField("enableTime", enableTime);
        currentModel.win3dDisplayAll();
    }

    public boolean getEnableTime() { return enableTime; }

	public StsColor getNextColor()
	{
		int i = list.getSize();
	    return StsColor.colors16[i%16];
	}

	public void setDefaultFractureSetColor(StsColor defaultFractureSetColor) {
		this.defaultFractureSetColor = defaultFractureSetColor;
	}

	public StsColor getDefaultFractureSetColor() {
		return defaultFractureSetColor;
	}
    public String getDefaultSpectrumName()  { return defaultSpectrumName; }
   public void setDefaultSpectrumName(String value)
   {
       if(this.defaultSpectrumName.equals(value)) return;
       this.defaultSpectrumName = value;
   }
 }