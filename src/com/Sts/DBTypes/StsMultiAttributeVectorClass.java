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

public class StsMultiAttributeVectorClass extends StsClass implements StsSerializable, StsClassDisplayable, StsClassCursorDisplayable
{
   private String defaultSpectrumName = StsSpectrumClass.SPECTRUM_RWB;

   public StsMultiAttributeVectorClass()
   {
       userName = "Vector Sets from multiple Volumes";
   }

   public void initializeFields()
   {
       displayFields = new StsFieldBean[]
           {
           new StsComboBoxFieldBean(this, "defaultSpectrumName", "Spectrum:", StsSpectrumClass.cannedSpectrums)
       };

       defaultFields = null;
   }

   public void drawOnCursor3d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, StsPoint[] planePoints, boolean isDragging) { }
   public void drawOn3dCurtain(StsGLPanel3d glPanel3d, StsSeismicCurtain seismicCurtain) { }
   public void drawOnCursor2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, boolean axesFlipped, boolean xAxisReversed, boolean yAxisReversed)
   {
       Iterator iter = getVisibleObjectIterator();
       while(iter.hasNext())
       {
           StsMultiAttributeVector vector = (StsMultiAttributeVector)iter.next();
           if(dirNo == StsCursor3d.ZDIR)
               vector.display2d(glPanel3d);
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
