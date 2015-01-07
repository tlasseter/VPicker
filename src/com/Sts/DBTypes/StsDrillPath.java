package com.Sts.DBTypes;


import com.Sts.Interfaces.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.Utilities.*;

public class StsDrillPath extends StsLine implements StsTreeObjectI
 {
     final static int DEVIATED = 0;
     final static int VERTICAL = 1;
     final static int SIDETRACK = 2;

 //    protected StsObjectRefList changeRate = null;
     protected float rateOfChange = 1.0f;

     static StsDateFieldBean bornField = new StsDateFieldBean(StsPlatform.class, "bornDate", "Born Date:");
     static StsDateFieldBean deathField = new StsDateFieldBean(StsPlatform.class, "deathDate", "Death Date:");
     // display fields
     static public StsFieldBean[] pathDisplayFields = null;

     public StsDrillPath()
     {
 //        super();
 //        setName(new String("DrillPath-" + index));
 //       changeRate = StsObjectRefList.constructor(2, 2, "changeRate", this);
     }

     public StsDrillPath(boolean persistent, String name)
     {
         super();
         setName(name);
     }

     static public StsDrillPath buildVerticalPath(StsGridPoint gridPoint) {
         try {
             return StsDrillPath.buildVerticalPath(gridPoint);
         }
         catch (Exception e) {
             StsException.systemError(
                 "StsDrillPath.buildVerticalPath(gridPoint) failed.");
             return null;
         }
     }
 /*
     static public StsDrillPath buildPath(StsGridSectionPoint[] gridPoints)
     {
         try
         {
             StsDrillPath line = new StsDrillPath();
             line.construct(gridPoints);
             return line;
         }
         catch (Exception e)
         {
             StsException.systemError(
                 "StsDrillPath.buildPath(gridPoints) failed.");
             return null;
         }
     }
 */
     public StsColor getStsColor()
     {
         return super.getStsColor();
     }

     public StsSurfaceVertex addLineVertex(StsPoint point, float rate, boolean computePath)
     {
         StsSurfaceVertex vertex = super.addLineVertex(point, computePath, false);

         // Add the rate of change to the reference list - SAJ

         return vertex;
     }

     public String getLabel()
     {
         String name = getName();
         if(name != null)
             return name;
         else
             return new String("DrillPath-" + getIndex() + " ");
     }

     public boolean delete()
     {
         if(!super.delete()) return false;
         return true;
     }

     /** returns true if changed */
     public boolean setIsVisibleNoDisplay(boolean isVisible)
     {
         if(this.isVisible == isVisible) return false;
         this.isVisible = isVisible;
         return true;
     }

     public boolean computePoints()
     {
         try
         {
             StsPoint[] projectPoints = super.computeRotatedCoorVertexPoints();
             StsPoint[] slopes = StsBezier.computeXYZLineSlopes(projectPoints);
             return computeRotatedPointsFromVertexPoints(projectPoints, slopes);
         }
         catch( Exception e )
         {
             StsException.outputException("Exception in StsLine.computeXYZPoints()", e, StsException.WARNING);
             return false;
         }
     }

     public boolean projectRotationAngleChanged()
     {
         return computePoints();
     }

     public StsFieldBean[] getDisplayFields()
     {
         if(pathDisplayFields == null)
         {
             pathDisplayFields = new StsFieldBean[]
             {
                 new StsBooleanFieldBean(StsDrillPath.class, "isVisible", "Enable"),
         bornField,
         deathField,
                 new StsColorComboBoxFieldBean(StsDrillPath.class, "stsColor", "Color", currentModel.getSpectrum("Basic").getStsColors()),
                 new StsFloatFieldBean(StsDrillPath.class, "topZ", false, "Min Depth"),
                 new StsFloatFieldBean(StsDrillPath.class, "botZ", false, "Max Depth"),
                 new StsDoubleFieldBean(StsDrillPath.class, "xOrigin", false, "X Origin"),
                 new StsDoubleFieldBean(StsDrillPath.class, "yOrigin", false, "Y Origin")
             };
         }
         return pathDisplayFields;
     }
     public StsFieldBean[] getPropertyFields() { return null; }
     public Object[] getChildren() { return new Object[0]; }
     public boolean anyDependencies() { return false; }

     public void setRate(float rate) { rateOfChange = rate; }
     public float getRate() { return rateOfChange; }

     static public StsFieldBean[] getStaticDisplayFields() { return pathDisplayFields; }

     public void display(StsGLPanel3d glPanel)
     {
         if (glPanel == null)
         {
             return;
         }
         super.display(glPanel);
     }

     public StsObjectPanel getObjectPanel()
     {
         if(objectPanel == null) objectPanel = StsObjectPanel.constructor(this, true);
         return objectPanel;
     }

     public void treeObjectSelected()
     {
         currentModel.getCreateStsClass("com.Sts.DBTypes.StsDrillPath").selected(this);
     }
     public void setBornDate(String born)
     {
         if(!StsDateFieldBean.validateDateInput(born))
         {
             bornField.setValue(StsDateFieldBean.convertToString(bornDate));
             return;
         }
         super.setBornDate(born);
     }
     public void setDeathDate(String death)
     {
         if(!StsDateFieldBean.validateDateInput(death))
         {
             deathField.setValue(StsDateFieldBean.convertToString(deathDate));
             return;
         }
         super.setDeathDate(death);
     }
 }