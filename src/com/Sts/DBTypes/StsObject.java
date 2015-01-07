//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.DB.*;
import com.Sts.Interfaces.StsDialogFace;
import com.Sts.Interfaces.StsTreeObjectI;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.StsGLPanel3d;
import com.Sts.UI.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;
import com.Sts.Types.StsMouse;

import javax.media.opengl.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;

/** Base class for all Sts Objects
  * @author TJLasseter
  * Index number of this geometric object uniquely identifies it.
  * The current lastIndex is maintained by the corresponding StsClass which manages the
  * instances of a particular type of StsObject (StsWell instances are managed by StsWellClass).
  * So this lastIndex is incremented and assigned to a new instance of that type.
  * If object is transient, index stays -1.
  */
abstract public class StsObject extends StsSerialize implements Cloneable, ActionListener, StsSerializable
{
    private transient int index = -1;

    final public void setIndex(int index) { this.index = index; }
    final public int getIndex() { return this.index; }

    public StsObject()
    {
        try
        {
            if (currentModel != null) currentModel.add(this);
        }
        catch (Exception e)
        {
            StsException.outputException(e, StsException.FATAL);
        }
    }

    public void addToModel()
    {
		if(currentModel == null) return;
		currentModel.add(this);
        refreshObjectPanel();
    }

    /** Override this in concrete subclass to respond when currentObject is set by corresponding StsClass. */
    public void setToCurrentObject()
    {

    }

	public void addCopyToModel()
	{
		StsObject object = (StsObject)StsToolkit.copyAllObjectFields(this, false);
	}

    public void refreshObjectPanel()
    {
        currentModel.refreshObjectPanel(this);
    }

    public static int getObjectIndex(StsObject object)
    {
        if (object == null)
        {
            return -99;
        }
        else
        {
            return object.getIndex();
        }
    }

    public static String getObjectIndexString(StsObject object)
    {
        if (object == null)
        {
            return new String("null");
        }
        int index = object.getIndex();
        if (index >= 0)
        {
            return new String("" + index);
        }
        else
        {
            return new String("transient");
        }
    }

    public StsObject(boolean persistent)
    {
        if (persistent)
        {
            try
            {
                if (currentModel != null)
                {
                    currentModel.add(this);
                }
            }
            catch (Exception e)
            {
                StsException.outputException(e, StsException.FATAL);
            }
        }
    }

    public boolean isPersistent()
    {
        return getIndex() >= 0;
    }

    public boolean indexOK()
    {
        return getIndex() >= 0;
    }

    public StsClass getStsClass()
    {
        try
        {
            return currentModel.getStsClass(this.getClass());
        }
        catch (Exception e)
        {
            return null;
        }
    }

    public StsClass getCreateStsClass()
    {
        try
        {
            return currentModel.getCreateStsClass(this.getClass());
        }
        catch (Exception e)
        {
            return null;
        }
    }
/*
    public StsObject getStsClassObjectWithIndex(Class c, int index)
    {
        StsClass instanceList = currentModel.getCreateStsClass(c);
        if (instanceList == null)
        {
            return null;
        }
        if (index < 0 || index >= instanceList.getSize())
        {
            return null;
        }
        return instanceList.getElementWithIndex(index);
    }
*/
    /*
     static public StsClass getCreateStsClass(Class c)
     {
      try	{ return currentModel.getCreateStsClass(c); }
      catch(Exception e) { return null; }
     }
     */
    public StsObject getStsClassObjectWithName(Class listClass, String name)
    {
        StsClass list = currentModel.getCreateStsClass(listClass);
        return getStsClassObjectWithName(list, name);
    }

    public StsObject getStsClassObjectWithName(String name)
    {
        StsClass list = currentModel.getCreateStsClass(this.getClass());
        return getStsClassObjectWithName(list, name);
    }

    public StsObject getStsClassObjectWithName(StsClass list, String name)
    {
        int nObjects = list.getSize();
        for (int n = 0; n < nObjects; n++)
        {
            StsObject stsObject = (StsObject) list.getElement(n);
			String objectName = stsObject.getName();
			if(objectName == null) continue; // shouldn't happen, but safeguards against bugs
            if (objectName.equals(name))
                return stsObject;
        }
        return null;
    }

    /*
     public StsObject getStsClassObjectWithName(StsClass list, String name)
     {
      int nObjects = list.getSize();
      if(nObjects == 0) return null;
      Object object = list.getElement(0);
      if(object instanceof StsMainObject)
      {
       for(int n = 0; n < nObjects; n++)
       {
        StsMainObject mainObject = (StsMainObject)list.getElement(n);
        if(mainObject.getName().equals(name)) return mainObject;
       }
      }
      else if(object instanceof StsObject)
      {
       for(int n = 0; n < nObjects; n++)
       {
        StsObject stsObject = (StsObject)list.getElement(n);
        if(stsObject.getName().equals(name)) return stsObject;
       }
      }
      return null;
     }
     */
    /** Methods to be implemented by subclasses as needed */

    public void display(StsGLPanel glPanel)
    {
        System.err.println("System error: display(StsGLPanel) not implemented in called class " + getClassname());
    }

    public boolean initialize(StsModel model)
    {
        System.err.println("initialize(StsModel model) not implemented in class " + getClassname());
        return false;
    }

    public void pick(GL gl, StsGLPanel glPanel)
    {
        System.err.println("System error: pick(GL) not implemented in called class: " + getClassname());
    }

    public boolean setHighlight(boolean state)
    {
        System.err.println("System error: setHighlight(boolean) not implemented in called class.");
        return false;
    }

    /** clone this StsObject and add it to the appropriate instance list */
    public Object clone()
    {
        try
        {
            Object newObject = super.clone();
			StsObject newStsObject = (StsObject)newObject;
			newStsObject.setIndex(-1);
            if(currentModel != null) currentModel.add(newStsObject);
            return newObject;
        }
        catch (Exception e)
        {
            System.out.println("Exception in StsObject()\n" + e);
            return null;
        }
    }

    public boolean delete()
    {
        if (getIndex() == -1)
        {
            return false;
        }
        return currentModel.delete(this);
    }

    public float compare(StsObject object)
    {
        return (float) StsParameters.UNDEFINED;
    }

    public String getName()
    {
        return new String(getClass().toString() + "-" + getIndex());
    }

    static public String getObjectName(StsObject object) { return object.getName(); }

    public String getLabel()
    {
        return new String(getClass().toString() + "-" + getIndex());
    }

    public String toString()
    {
        return new String(getClass().toString() + "[" + getIndex()+ "]");
    }

    public String toDebugString()
    {
        return new String(getClass().toString() + "[" + getIndex()+ "]");
    }

    // generally redundant as getName returns className-index.  TJL 2/7/07
    // suggest just using toString()
    public String getClassAndNameString()
	{
		return getClass().toString() + "[" + getIndex()+ "]: " + getName();
    }

	public String getClassname()
	{
		return getClass().getName();
    }

    /** methods to override */

    public void setStsColor(StsColor color)
    {
        System.err.println("System error: setStsColor(StsColor) not implemented in called class " + getClassname());
    }

    public StsColor getStsColor()
    {
        System.err.println("System error: getStsColor() not implemented in called class " + getClassname());
        return null;        
    }

    //TODO remove all calls to getColor and setBeachballColors so we are consistently working with StsColor(s)
    public void setColor(Color color)
    {
        System.err.println("System error: setBeachballColors(Color) not implemented in called class " + getClassname());
    }

    public boolean getIsVisible()
    {
        return true;
    } // StsMainObject overrides this

    /*
     public void setBeachballColors(Color color)
     {
      StsSpectrum spectrum = null;
      try { spectrum = currentModel.getSpectrum("Basic"); }
      catch(Exception e) { }
      if( spectrum == null ) return;
      StsColor[] colors = spectrum.getStsColors();
      for( int i=0; i<colors.length; i++ )
      {
       if( colors[i].getColor().equals(color) )
       {
        setBeachballColors(colors[i]);
        break;
       }
      }
     }
     */

    public float getOrderingValue()
    {
        return StsParameters.nullValue;
    }

    public boolean isInList(StsObject[] objects)
    {
        if (objects == null)
        {
            return false;
        }
        int nObjects = objects.length;
        for (int n = 0; n < nObjects; n++)
        {
            if (this == objects[n])
            {
                return true;
            }
        }
        return false;
    }

    public void logMessage(String msg)
    {
        StsMessageFiles.logMessage(msg);
    }

    public void infoMessage(String msg)
    {
        StsMessageFiles.infoMessage(msg);
    }

    public void errorMessage(String msg)
    {
        StsMessageFiles.errorMessage(msg);
    }

    public boolean launch()
    {
        String className = getClassname();
        StsMessageFiles.infoMessage("Objects of class " + className + " cannot be launched.");
        return false;
    }
    public boolean goTo()
    {
        String className = getClassname();
        StsMessageFiles.infoMessage("Cannot change location to objects of class " + className + ".");
        return false;
    }    
    public boolean canExport() { return true; }

    public boolean export()
	{
		 String className = StsToolkit.getSimpleClassname(this);
		 try
		 {
			 String directory = currentModel.getProject().getDataFullDirString();
			 String filename = className + ".obj." + getName();
			 StsDBFileObjectTrader.exportStsObject(directory + File.separator + filename, this, null);
			 StsMessageFiles.infoMessage("Successfully exported file:" + directory + File.separator + filename);
			 return true;
		 }
		 catch(Exception e)
		 {
			 StsException.outputException("StsObject export failed for " + getName(), e, StsException.WARNING);
			 return false;
		 }
	}

    public boolean canLaunch() { return false; }

    public void close()
    {
    }

    public void instanceChanged(String reason)
    {
        currentModel.instanceChange(this, reason);
    }

	public void actionPerformed(ActionEvent e)
	{
		System.err.println("System error: actionPerformed(ActionEvent) not implemented in called class " + getClassname());
	}

	public void objectPropertiesChanged(Object object) { }

    public StsColorscale getColorscaleWithName(String name)
    {
        System.err.println("System error: getColorscaleWithName() not implemented in called class " + getClassname());
        return null;
    }

    // Can only be used in viewer since the beans are static. The object panel is not shown when in viewer only mode so not a problem.
    public void popupPropertyPanel()
    {
        if(!StsTreeObjectI.class.isAssignableFrom(getClass())) return;
		try
        {
            StsFieldBeanPanel box = buildPropertyPanel((StsTreeObjectI)this);
            box.gbc.fill = box.gbc.HORIZONTAL;
            new StsOkDialog(currentModel.win3d, box, getName() + " Properties", false);
        }
        catch(Exception ex)
        {
            StsException.outputException("Failed to construct property panel for " + this.getName(), ex, StsException.WARNING);
        }
    }

    public void popupClassListPanel()
    {
        getStsClass().popupClassListPanel();    
    }

    // Can only be used in viewer since the beans are static. The object panel is not shown when in viewer only mode so not a problem.
    public void popupClassPropertyPanel()
    {
        if(!StsTreeObjectI.class.isAssignableFrom(getClass())) return;
		try
        {
            StsFieldBeanPanel box = buildPropertyPanel((StsTreeObjectI)getStsClass());
            box.gbc.fill = box.gbc.HORIZONTAL;
            new StsOkDialog(currentModel.win3d, box, getName() + " Properties", false);
        }
        catch(Exception ex)
        {
            StsException.outputException("Failed to construct class property panel for " + this.getName(), ex, StsException.WARNING);
        }
    }

    public StsFieldBeanPanel buildPropertyPanel(StsTreeObjectI treeObject)
    {
        StsFieldBeanPanel propertiesBox = new StsFieldBeanPanel();
        StsFieldBean[] displayFields = treeObject.getDisplayFields();
        if(displayFields == null) return null;
		for( int i=0; i<displayFields.length; i++ )
        {
		    propertiesBox.addBeanPanel(displayFields[i]);
		}
        propertiesBox.setPanelObject(this);
        return propertiesBox;
    }

}