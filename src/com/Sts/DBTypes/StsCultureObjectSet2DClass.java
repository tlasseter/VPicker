package com.Sts.DBTypes;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2002</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */

import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;

import java.util.*;

public class StsCultureObjectSet2DClass extends StsClass implements StsSerializable, StsClassDisplayable, StsClassCursorDisplayable
{
    private String defaultSpectrumName = StsSpectrumClass.SPECTRUM_RWB;

    public StsCultureObjectSet2DClass()
    {
        userName = "2D Sets of Culture";
    }

    public void initializeFields()
    {
        displayFields = new StsFieldBean[]
        {
            new StsComboBoxFieldBean(this, "defaultSpectrumName", "Spectrum:", StsSpectrumClass.cannedSpectrums)
        };

        defaultFields = null;
    }
   
    public void displayClass(StsGLPanel3d glPanel3d)
    {
        Iterator iter = getVisibleObjectIterator();
        while(iter.hasNext())
        {
        	StsCultureObjectSet2D set = (StsCultureObjectSet2D)iter.next();
            set.display(glPanel3d);
        }
    }
    
    public void drawOnCursor3d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, StsPoint[] planePoints, boolean isDragging) { }
    public void drawOn3dCurtain(StsGLPanel3d glPanel3d, StsSeismicCurtain seismicCurtain) { }
    public void drawOnCursor2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, boolean axesFlipped, boolean xAxisReversed, boolean yAxisReversed)
    {
        Iterator iter = getVisibleObjectIterator();
        while(iter.hasNext())
        {
            StsCultureObjectSet2D culture = (StsCultureObjectSet2D)iter.next();
            if(dirNo == StsCursor3d.ZDIR)
                culture.display2d(glPanel3d);
        }
   }
   
   public String getDefaultSpectrumName()  { return defaultSpectrumName; }
   public void setDefaultSpectrumName(String value)
   {
       if(this.defaultSpectrumName.equals(value)) return;
       this.defaultSpectrumName = value;
//       setDisplayField("defaultSpectrumName", defaultSpectrumName);
   }

}
