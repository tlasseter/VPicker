package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Nov 9, 2007
 * Time: 5:45:18 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsFaultStickSet extends StsMainObject implements StsSelectable, StsTreeObjectI
{
    public StsObjectRefList faultSticks = null;
    static public StsFieldBean[] displayFields = null;
    static protected StsObjectPanel objectPanel = null;
    static StsFaultStickSetClass faultStickSetClass = null;
    
    protected byte zDomainSupported = StsProject.TD_DEPTH;
    protected byte zDomainOriginal = StsParameters.TD_NONE;

    static public final String lmkGrp = "faultSticks";
    static public final String rmsGrp = "faultSticks-rms";
    static public final String[] fileGroups = new String[] { lmkGrp, rmsGrp };
    public StsFaultStickSet()
	{
	}

    static public StsFaultStickSet constructor(String name)
    {
        try
        {
            StsFaultStickSet stickSet = new StsFaultStickSet();
            stickSet.initialize(name);
            return stickSet;
        }
        catch(Exception e)
        {
            return null;
        }
    }

    private void initialize(String name)
    {
        faultSticks = StsObjectRefList.constructor(4, 4, "faultSticks", this);
        if(name == null)
            name = "FaultStickSet-" + getIndex();
        setName(name);
    }

    public boolean initialize(StsModel model)
    {
        return true;
    }

    public int getNumberOfSticks()
    {
    	return faultSticks.getSize();
    }
 
    public byte getZDomainSupported() { return zDomainSupported; }
    public byte getZDomainOriginal() { return zDomainOriginal; }
    
    public boolean delete()
    {
        if(!super.delete()) return false;
        if(faultSticks != null)  faultSticks.deleteAll();
        return true;
    }

    public void addLine(StsFaultLine line)
    {
    	faultSticks.add(line);
    }
    
    public void setIsVisible(boolean value)
    {
        for(int i=0; i<faultSticks.getSize(); i++)
        	((StsFaultLine)faultSticks.getElement(i)).setIsVisible(value);
            	
        currentModel.viewObjectRepaint(this, this);
    }
    
    public boolean getIsVisible()
    {
    	if(faultSticks.getSize() < 1) return true;
        return faultSticks.getElement(0).getIsVisible();
    }
    
    public void setStsColor(StsColor color)
    {
        for(int i=0; i<faultSticks.getSize(); i++)
        {
        	faultSticks.getElement(i).setStsColor(color);
        }
        currentModel.viewObjectRepaint(this, this);
    }
    
    public void setDrawLabels(boolean state)
    {
    	if(faultSticks.getSize() < 1) return;    	
        for(int i=0; i<faultSticks.getSize(); i++)
        {
        	StsFaultLine stick = (StsFaultLine)faultSticks.getElement(i);
        	stick.setDrawLabels(state);
        }
        currentModel.viewObjectRepaint(this, this);
    }

    public boolean getDrawLabels()
    {
    	if(faultSticks.getSize() < 1) return false;    	
    	StsFaultLine stick = (StsFaultLine)faultSticks.getElement(0);
    	return stick.getDrawLabels();
    }
    
    public StsColor getStsColor() 
    { 
    	if(faultSticks.getSize() < 1) return StsColor.RED;
        return	faultSticks.getElement(0).getStsColor();
    }
	public StsFieldBean[] getDisplayFields()
	{
        if(displayFields != null) return displayFields;
        displayFields = new StsFieldBean[]
        {
            new StsBooleanFieldBean(StsFaultStickSet.class, "isVisible", "Enable"),
            new StsBooleanFieldBean(StsFaultStickSet.class, "drawLabels", "Draw Labels"),
            new StsColorComboBoxFieldBean(StsFaultStickSet.class, "stsColor", "Color", currentModel.getSpectrum("Basic").getStsColors())
        };
        return displayFields;
    }

    public StsFieldBean[] getPropertyFields()
	{
        return null;
    }

	public Object[] getChildren()
	{
        return new Object[0];
	}

	public boolean anyDependencies()
	{
		return false;
	}

	public StsObjectPanel getObjectPanel()
	{
        if (objectPanel != null) return objectPanel;
        objectPanel = StsObjectPanel.constructor(this, true);
        return objectPanel;
    }

	public void treeObjectSelected()
	{
		currentModel.getCreateStsClass("com.Sts.DBTypes.StsFractureSet").selected(this);
	}

	public String getLabel()
	{
		return new String("fracture set: " + getName());
	}

    public StsObjectRefList getFaultStickList()
    {
        return faultSticks;
    }

    public int getFaultStickIndex(StsFaultLine stick)
    {
        return faultSticks.getIndex(stick);
    }
}
