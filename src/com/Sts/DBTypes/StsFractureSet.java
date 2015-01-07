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
public class StsFractureSet extends StsMainObject implements StsSelectable, StsTreeObjectI
{
    protected StsObjectRefList fractures = null;
    protected StsObjectRefList zCursorSurfaces = null;    // StsCursorZSurfaces used to define this fracture set.
    protected boolean drawStimulatedVolume = true;
    private float stimulatedRadius = 500.0f;
    private transient float rowGridSize;
    private transient float colGridSize;

    static public StsFieldBean[] displayFields = null;
    static public StsFieldBean[] propertyFields = null;
    static protected StsObjectPanel objectPanel = null;
    static int newNameIndex = 0;
    static int nextColorIndex = 0;
    static StsFractureSetClass fractureSetClass = null;
    
    public StsFractureSet()
	{
	}

	public StsFractureSet(boolean persistent)
	{
		super(persistent);
    }

    public void initialize()
    {
        checkCreateFractureSetClass();
        initializeVolumeStimulated();
        fractures = StsObjectRefList.constructor(4, 4, "fractures", this);
        zCursorSurfaces = StsObjectRefList.constructor(4, 4, "zCursorSurfaces", this);
        setName(fractureSetClass.getNextName());
        StsColor color = new StsColor(fractureSetClass.getNextColor());
        setStsColor(color);
    }

	static public StsFractureSet constructor()
	{
        try
        {
            StsFractureSet fractureSet = new StsFractureSet(false);
            fractureSet.initialize();
            return fractureSet;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(StsFractureSet.class, "constructor", e);
            return null;
        }
    }

    public boolean initialize(StsModel model)
    {
        checkCreateFractureSetClass();
        reinitializeVolumeStimulated();
        return true;
    }

    public void initializeVolumeStimulated()
    {
        setStimulatedRadius(fractureSetClass.getDefaultStimulatedRadius());
        reinitializeVolumeStimulated();
    }
    public void reinitializeVolumeStimulated()
    {
        rowGridSize = fractureSetClass.getRowGridSize();
        colGridSize = fractureSetClass.getColGridSize();
    }

    public void checkCreateFractureSetClass()
    {
        if(fractureSetClass == null)
            fractureSetClass = (StsFractureSetClass)currentModel.getCreateStsClass(StsFractureSet.class);
    }
    
    public void setDrawStimulated(boolean draw)
    {
        for(int i=0; i<fractures.getSize(); i++)
        	((StsFracture)fractures.getElement(i)).setDrawStimulated(draw);

        currentModel.viewObjectRepaint(this, this);
    }
    public boolean getDrawStimulated()
    {
    	if(fractures.getSize() < 1) return true;
        return	((StsFracture)fractures.getElement(0)).getDrawStimulated();
    }
    
    public void setAreaScale(int scale) 
    { 
        for(int i=0; i<fractures.getSize(); i++)
        	((StsFracture)fractures.getElement(i)).setAreaScale(scale);

        currentModel.viewObjectRepaint(this, this);
    }
    public int getAreaScale()
    {
    	if(fractures.getSize() < 1) return 1;
        return	((StsFracture)fractures.getElement(0)).getAreaScale();
    }
    
    public boolean delete()
    {
        if(!super.delete()) return false;
        if(fractures != null)  fractures.deleteAll();
    	if(zCursorSurfaces != null) zCursorSurfaces.deleteAll();
        // if no fracture sets remain, delete the cursorZSurfaces used for construction
        if(!currentModel.getStsClass(StsFractureSet.class).hasObjects())
        {
        	StsClass cursorZSurfaceClass = currentModel.getStsClass(StsCursorZSurface.class);
        	if(cursorZSurfaceClass.hasObjects()) cursorZSurfaceClass.deleteAll();
        }
        return true;
    }

    public void addFracture(StsFracture fracture)
    {
    	fractures.add(fracture);
        currentModel.checkAddToCursor3d(fracture);
    }
    
    public void addCursorZSurface(StsCursorZSurface zCursor)
    {
    	zCursorSurfaces.add(zCursor);
    } 
    
    public StsCursorZSurface[] getCursorZSurfaces()
    { 
    	if(zCursorSurfaces.getSize() == 0)
    		return new StsCursorZSurface[0];
    	else
    		return (StsCursorZSurface[]) zCursorSurfaces.getCastList(StsCursorZSurface.class);
    }
	
    public void setCursorZSurfaces(StsCursorZSurface[] zSurfaces)
    {
    	zCursorSurfaces.add(zSurfaces);
    }

    public void setDrawLines(boolean value)
    {
        for(int i=0; i<fractures.getSize(); i++)
        	((StsFracture)fractures.getElement(i)).setDrawLines(value);

        currentModel.viewObjectRepaint(this, this);
    }
       
    public boolean getDrawLines()
    {
    	if(fractures.getSize() < 1) return true;
        return	((StsFracture)fractures.getElement(0)).getDrawLines();
    }
    
    public void setDrawEdges(boolean value)
    {
        for(int i=0; i<fractures.getSize(); i++)
        	((StsFracture)fractures.getElement(i)).setDrawEdges(value);
        
        currentModel.viewObjectRepaint(this, this);
    }   
    
    public boolean getDrawEdges()
    {
    	if(fractures.getSize() < 1) return true;
        return	((StsFracture)fractures.getElement(0)).getDrawEdges();
    }
    
    public void setIsVisible(boolean value)
    {
        for(int i=0; i<fractures.getSize(); i++)
        	((StsFracture)fractures.getElement(i)).setIsVisible(value);
            	
        currentModel.viewObjectChangedAndRepaint(this, this);
    }
    
    public boolean getIsVisible()
    {
    	if(fractures.getSize() < 1) return true;
        return fractures.getElement(0).getIsVisible();
    }
    
    public void setStsColor(StsColor color)
    {
        for(int i=0; i<fractures.getSize(); i++)
        {
        	((StsFracture)fractures.getElement(i)).setStsColor(color);
        }
        currentModel.viewObjectRepaint(this, this);
    }

    public StsColor getStsColor() 
    { 
    	if(fractures.getSize() < 1) return StsColor.RED;
        return	((StsFracture)fractures.getElement(0)).getStsColor();
    }

    public byte getZDomain()
    {
        StsCursorZSurface zCursorSurface = (StsCursorZSurface)zCursorSurfaces.getFirst();
        return zCursorSurface.zDomain;
    }

    public boolean export()
    {
        byte zDomain = getZDomain();
        return export(StsParameters.TD_ALL_STRINGS[zDomain]);
    }

    public boolean export(String timeOrDepth)
    {
        return StsFractureSetExportDialog.exportFractureSet(currentModel, currentModel.win3d, "FractureSet Export Utility", true, this, timeOrDepth);
    }

	public StsFieldBean[] getDisplayFields()
	{
        if(displayFields == null)
        {
            displayFields = new StsFieldBean[]
            {
                new StsBooleanFieldBean(StsFractureSet.class, "isVisible", "Enable Fractures"),
                new StsFloatFieldBean(this, "stimulatedRadius", 0.0f, 10000.0f, "Stimulated radius"),
                new StsBooleanFieldBean(StsFractureSet.class, "drawEdges", "Edges"),
                new StsBooleanFieldBean(StsFractureSet.class, "drawLines", "Lines"),
                new StsColorComboBoxFieldBean(StsFractureSet.class, "stsColor", "Color", currentModel.getSpectrum("Basic").getStsColors())
            };
        }
        return displayFields;
    }

    public boolean getStimVisible() { return drawStimulatedVolume; }
    public void setStimVisible(boolean val)
    {
        drawStimulatedVolume = val;
    }
    public StsFieldBean[] getPropertyFields()
	{
        return propertyFields;
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

    public StsObjectRefList getFractureList() { return fractures; }

    public float getStimulatedRadius()
    {
        return stimulatedRadius;
    }

    public void setStimulatedRadius(float stimulatedRadius)
    {
        this.stimulatedRadius = stimulatedRadius;
        dbFieldChanged("stimulatedRadius", stimulatedRadius);
        currentModel.viewObjectChangedAndRepaint(this, this);
    }

    public float getRowGridSize()
    {
        return rowGridSize;
    }

    public float getColGridSize()
    {
        return colGridSize;
    }
   
}
