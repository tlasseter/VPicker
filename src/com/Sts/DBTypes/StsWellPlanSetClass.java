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
import com.Sts.Utilities.*;

import java.util.*;

public class StsWellPlanSetClass extends StsClass implements StsSerializable, StsClassTimeDisplayable, StsClassCursorDisplayable
{
	//static final Class[] subClassClasses = new Class[] { StsWellPlanSet.class };
	//static StsWellPlanSetClass[] subClasses = null;
    boolean enableTime = true;

	public StsWellPlanSetClass()
	{
        userName = "Well Plan Sets";                                                                              
	}

    // for mainDebug purposes only
    public void add(StsObject obj) throws StsException
    {
        super.add(obj);
    }

    // for mainDebug purposes only
    public Object[] getChildren()
    {
        return (Object[])list.copyArrayList();
    }


    public void initializeFields()
    {
        //initializeSubClasses();
    }
/*
    private void initializeSubClasses()
    {
		int nSubClasses = subClassClasses.length;
		subClasses = new StsWellPlanSetClass[nSubClasses];
		int nActualInstances = 0;
		for(int n = 0; n < nSubClasses; n++)
		{
			StsWellPlanSetClass subClassInstance = (StsWellPlanSetClass) currentModel.getStsClass(subClassClasses[n]);
			if(subClassInstance != null) subClasses[nActualInstances++] = subClassInstance;
		}
		subClasses = (StsWellPlanSetClass[])StsMath.trimArray(subClasses, nActualInstances);
    }
*/
    public void initializeDisplayFields()
	{
//        initColors(StsWellPlanSet.displayFields);

		displayFields = new StsFieldBean[]
        {
			new StsBooleanFieldBean(this, "isVisible", "Visible"),
            new StsBooleanFieldBean(this, "enableTime", "Enable Time")
        };
	}
    public void setEnableTime(boolean enable)
    {
        if(this.enableTime == enable) return;
        this.enableTime = enable;
//        setDisplayField("enableTime", enableTime);
        currentModel.win3dDisplayAll();
    }
    public boolean getEnableTime() {	return enableTime; }
    public void selected(StsWellPlanSet planSet)
	{
		super.selected(planSet);
		setCurrentObject(planSet);
	}

	public boolean setCurrentObject(StsObject object)
	{
		boolean changed = super.setCurrentObject(object);
		if (changed && object != null) ( (StsWellPlanSet) object).treeObjectSelected();
		return changed;
	}

	public boolean setCurrentWellPlanSetName(String name)
	{
		StsWellPlanSet newPlanSet = (StsWellPlanSet)getObjectWithName(name);
		return setCurrentObject(newPlanSet);
	}

	public void close()
	{
		list.forEach("close");
	}

	/** called to set this superClass or one of its subTypes with
	 *  new isVisible value.  If superClass, subTypes are called and set.
	 */
	public void setIsVisible(boolean isVisible)
	{
		if (this.isVisible == isVisible) return;
		this.isVisible = isVisible;
		boolean changed = false;
		ArrayList<StsClass> subClasses = getSubClasses();
		if(subClasses != null)
		{
			for(StsClass subClass : subClasses)
				if(((StsWellPlanSetClass)subClass).setIsVisibleNoDisplay(isVisible)) changed = true;
		}
		else
		{
			// apply to instances of this subClass
			if(setIsVisibleNoDisplay(isVisible)) changed = true;
			// check superClass to see if its visibility is changed
			StsWellPlanSetClass wellPlanSetClass = getWellPlanSetClass();
			wellPlanSetClass.checkSetSuperClassIsVisible(isVisible);
		}
		if(changed) currentModel.win3dDisplay();
	}

	/** Called only on the superclass StsDrillPathClass.
	 *  A subClass visibility has changed to isVisible;
	 *  if subClass is true, then this superClass must be true
	 *  to be consistent.  Otherwise if all subTypes are false,
	 *  then set superClass to false..
	 */
	public void checkSetSuperClassIsVisible(boolean isVisible)
	{
		this.isVisible = isVisible;
		/*
		if(isVisible)
			this.isVisible = true;
		else
		{
			for (int n = 0; n < subClasses.length; n++)
				if (subClasses[n].getIsVisible()) return;
			this.isVisible = false;
		}
		*/
	}

	/** Called to set subClass instances visible values. */
	private boolean setIsVisibleNoDisplay(boolean isVisible)
	{
		this.isVisible = isVisible;
		boolean changed = false;
		StsList list = getList();
		for (int n = 0; n < list.getSize(); n++)
		{
			StsWellPlanSet wellPlanSet = (StsWellPlanSet) list.getElement(n);
			if(wellPlanSet.setIsVisibleNoDisplay(isVisible)) changed = true;
		}
		return changed;
	}

	/** An instance has been toggled to state isVisible.
	 * If this is a StsSubVolumeClass subClass instance,
	 * check and set toggle and then move up to superClass
	 * and do the set and check there.
	 */
	public void checkSetClassIsVisible(boolean isVisible)
	{
		if (isVisible) // if object below isVisible, then this must be visible as well
		{
			if(this.isVisible) return;
			this.isVisible = true;
		}
		else
		{
			if (getClass() != StsWellPlanSetClass.class)
			{
				StsList list = getList();
				for (int n = 0; n < list.getSize(); n++)
				{
					StsWellPlanSet wellPlanSet = (StsWellPlanSet) list.getElement(n);
					if (wellPlanSet.getIsVisible())return;
				}
				this.isVisible = false;
			}
		}
		// isVisible flag has been changed for this class: check superClass
		StsWellPlanSetClass wellPlanSetClass = getWellPlanSetClass();
		wellPlanSetClass.checkSetSuperClassIsVisible(isVisible);
	}
    public void displayTimeClass(StsGLPanel3d glPanel3d, long time)
    {
        Iterator iter = getVisibleObjectIterator();
        while(iter.hasNext())
        {
            StsWellPlanSet plan = (StsWellPlanSet)iter.next();
            if((enableTime && plan.isAlive(time)) || (!enableTime))
                plan.display(glPanel3d);
        }
    }
    public void displayClass(StsGLPanel3d glPanel3d)
    {
        Iterator iter = getVisibleObjectIterator();
        while(iter.hasNext())
        {
            StsWellPlanSet plan = (StsWellPlanSet)iter.next();
            plan.display(glPanel3d);
        }
    }

    public void drawOnCursor3d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, StsPoint[] planePoints, boolean isDragging)
    {
    }

    public void drawOn3dCurtain(StsGLPanel3d glPanel3d, StsSeismicCurtain seismicCurtain)
    {
    }

    public void displayClass2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, boolean axesFlipped, boolean xAxisReversed, boolean yAxisReversed)
    {
        Iterator iter = getVisibleObjectIterator();
        while(iter.hasNext())
        {
            StsWellPlanSet planSet = (StsWellPlanSet)iter.next();
			planSet.display2d(glPanel3d, dirNo, dirCoordinate, axesFlipped, xAxisReversed, yAxisReversed);
        }
    }

    public void drawOnCursor2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, boolean axesFlipped, boolean xAxisReversed, boolean yAxisReversed)
    {
        Iterator iter = getVisibleObjectIterator();
        while(iter.hasNext())
        {
            StsWellPlanSet planSet = (StsWellPlanSet)iter.next();
            planSet.display2d(glPanel3d, dirNo, dirCoordinate, axesFlipped, xAxisReversed, yAxisReversed);
        }
    }

	static public StsWellPlanSetClass getWellPlanSetClass()
	{
		return (StsWellPlanSetClass)currentModel.getCreateStsClass(StsWellPlanSet.class);
	}
}
