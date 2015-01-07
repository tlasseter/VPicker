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

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;

public class StsCultureObjectSet2D extends StsMainObject implements StsTreeObjectI
{
    StsCultureObject2D[] objects;
    boolean mapToSurface = true;
    StsObjectRefList cultureDisplayables;
    public StsCultureDisplayable cultureDisplayableObject = null;
    protected byte zDomainSupported = StsProject.TD_DEPTH;
    protected byte zDomainOriginal = StsParameters.TD_NONE;
    
    static public StsFieldBean[] displayFields = null;

    static public StsComboBoxFieldBean cultureDisplayableList;
    static public StsFieldBean[] propertyFields = null;

    static StsObjectPanel objectPanel = null;

    public StsCultureObjectSet2D()
    {
    }

    public StsCultureObjectSet2D(boolean persist)
    {
        super(persist);
    }

    public StsCultureObjectSet2D(String name)
    {
        super(false);
        if(name == null)
            setName("Culture-" + getIndex());
        else
            setName(name);
        setStsColor(StsColor.WHITE);
        cultureDisplayableObject = currentModel.getProject();
    }

    public boolean initialize(StsModel model)
    {
        cultureDisplayableObject = model.getProject();
        for(int n = 0; n < objects.length; n++)
        {
            objects[n].initialize();
        }
        return true;
    }


    public void addObject(StsCultureObject2D cultureObject)
    {
        if(objects == null)
            objects = new StsCultureObject2D[] { cultureObject };
        else
            objects = (StsCultureObject2D[])StsMath.arrayAddElement(objects, cultureObject);
    }

    public StsCultureObject2D[] getCultureObjects() { return objects; }

    public void addCultureDisplayable(StsObject cultureDisplayableObject)
    {
        if(cultureDisplayables == null) cultureDisplayables = StsObjectRefList.constructor(2, 2, "cultureDisplayables", this, true);
        cultureDisplayables.add(cultureDisplayableObject);
    }

    public void display(StsGLPanel3d glPanel3d)
    {
        if(objects == null)
            return;
        if(!currentModel.getProject().canDisplayZDomain(zDomainSupported))
            return;       
        if(!getIsVisible())
            return;
        GL gl = glPanel3d.getGL();

        for(int n = 0; n < objects.length; n++)
        {
            if(cultureDisplayableObject != null)
            {
                objects[n].draw(cultureDisplayableObject, mapToSurface, zDomainOriginal, glPanel3d);
            }
        }
    }

    public int getNumberOfObjects()
    {
    	if(objects == null) return 0;
    	return objects.length;
    }
    
    static public StsCultureObjectSet2DClass getCultureObject2DClass()
    {
        return (StsCultureObjectSet2DClass) currentModel.getCreateStsClass(StsCultureObjectSet2D.class);
    }

    public void display2d(StsGLPanel glPanel)
    {
        if(objects == null)
            return; 
        if(!currentModel.getProject().canDisplayZDomain(zDomainSupported))
            return;        
        if(!getIsVisible())
            return;

        for(int n = 0; n < objects.length; n++)
        {
            if(cultureDisplayableObject != null)
            {
                objects[n].draw2d(cultureDisplayableObject, mapToSurface, zDomainOriginal, glPanel);
            }
        }
    }

    public double getMaxDepth()
    {
        double maxDepth = currentModel.getProject().getDepthMax();
        if(objects == null) return maxDepth;       
         for(int n = 0; n < objects.length; n++)
        {
            if(maxDepth < objects[n].getMaxDepth(zDomainOriginal))
                maxDepth = objects[n].getMaxDepth(zDomainOriginal);
        }
        return maxDepth;
    }

    public double getMaxTime()
    {
        double maxTime = currentModel.getProject().getTimeMax();
        if(objects == null) return maxTime;        
         for(int n = 0; n < objects.length; n++)
        {
            if(maxTime < objects[n].getMaxTime(zDomainOriginal))
                maxTime = objects[n].getMaxTime(zDomainOriginal);
        }
        return maxTime;
    }

    public double getMinDepth()
    {
        double minDepth = currentModel.getProject().getDepthMin();
        if(objects == null) return minDepth;
         for(int n = 0; n < objects.length; n++)
        {
            if(minDepth > objects[n].getMinDepth(zDomainOriginal))
                minDepth = objects[n].getMinDepth(zDomainOriginal);
        }
        return minDepth;
    }

    public double getXMax()
    {
    	double xMax = -StsParameters.largeDouble;
        for(int n = 0; n < objects.length; n++)   	
        {
            if(xMax < objects[n].getXMax())
                xMax = objects[n].getXMax();
        }
        return xMax;
    }
    
    public double getYMax()
    {
    	double yMax = -StsParameters.largeDouble;
        for(int n = 0; n < objects.length; n++)   	
        {
            if(yMax < objects[n].getYMax())
                yMax = objects[n].getYMax();
        } 
        return yMax;
    }
    public double getMinTime()
    {
        double minTime = currentModel.getProject().getTimeMax();
        if(objects == null) return minTime;
         for(int n = 0; n < objects.length; n++)
        {
            if(minTime < objects[n].getMinTime(zDomainOriginal))
                minTime = objects[n].getMinTime(zDomainOriginal);
        }
        return minTime;
    }
    public void treeObjectSelected()
    {

    }

    public boolean anyDependencies()
    {
        return false;
    }

    public StsObjectPanel getObjectPanel()
    {
        if (objectPanel == null)
        {
            objectPanel = StsObjectPanel.constructor(this, true);
        }
        initializePanel();
        return objectPanel;
    }

    private void initializePanel()
    {
        StsObject[] modelSurfaces = (StsObject[])currentModel.getObjectList(StsModelSurface.class);
        StsObject[] surfaces = (StsObject[])currentModel.getObjectList(StsSurface.class);
        int nModelSurfaces = modelSurfaces.length;
        int nSurfaces = surfaces.length;
        Object[] listItems = new Object[nModelSurfaces + nSurfaces + 1];
        listItems[0] = currentModel.getProject();
        System.arraycopy(surfaces, 0, listItems, 1, nSurfaces);
        System.arraycopy(modelSurfaces, 0, listItems, nSurfaces+1, nModelSurfaces);
        cultureDisplayableList.setListItems(listItems);
    }

    public Object[] getChildren()
    {
        return new Object[0];
    }

    public StsFieldBean[] getDisplayFields()
    {
        if(displayFields != null) return displayFields;

        displayFields = new StsFieldBean[]
        {
            new StsStringFieldBean(StsCultureObjectSet2D.class, "name", "Name"),
            new StsBooleanFieldBean(StsCultureObjectSet2D.class, "isVisible", "Enable"),
            new StsBooleanFieldBean(StsCultureObjectSet2D.class, "mapToSurface", "Map To Surface"),

//            new StsIntFieldBean(StsCultureObjectSet2D.class, "width", true, "Line Width:", true)
        };
        return displayFields;
    }

    public boolean addToProject()
    {
        this.zDomainSupported = StsProject.TD_TIME_DEPTH;
        this.zDomainOriginal = currentModel.getProject().getZDomain();
        currentModel.getProject().addToProject(this);
        return true;
    }
    public byte getZDomainSupported() { return zDomainSupported; }
    public byte getZDomainOriginal() { return zDomainOriginal; }
    public boolean getMapToSurface() { return mapToSurface; }
    public void setMapToSurface(boolean map)
    {
        mapToSurface = map;
        dbFieldChanged("mapToSurface" ,mapToSurface);
        currentModel.win3dDisplay();
    }
    public StsFieldBean[] getPropertyFields()
    {
        if(propertyFields != null) return propertyFields;
        cultureDisplayableList = new StsComboBoxFieldBean(StsCultureObjectSet2D.class, "cultureDisplayable", "Culture Surface");
        propertyFields = new StsFieldBean[]
        {
            cultureDisplayableList
        };
        return propertyFields;
    }

    public void setCultureDisplayable(StsCultureDisplayable cultureDisplayable)
    {
        cultureDisplayableObject = cultureDisplayable;
        dbFieldChanged("cultureDisplayableObject", cultureDisplayableObject);
        currentModel.win3dDisplayAll();
    }
    public StsCultureDisplayable getCultureDisplayable() { return cultureDisplayableObject; }

//    public StsCultureDisplayable getCultureDisplayable() { return (StsCultureDisplayable)cultureDisplayables.getLast(); }
}
