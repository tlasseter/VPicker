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

public class StsDrillPathClass extends StsClass implements StsSerializable, StsClassDisplayable
{

	boolean isApplied = false; // indicates subVolume is applied to cursor display

	// static final Class[] subClassClasses = new Class[] { StsDrillPath.class };
	// static StsDrillPathClass[] subClasses = null;

	public StsDrillPathClass()
	{
        userName = "Drill Paths";
	}

    public void initializeFields()
    {
        //initializeSubClasses();
    }
/*
    private void initializeSubClasses()
    {
		int nSubClasses = subClassClasses.length;
		subClasses = new StsDrillPathClass[nSubClasses];
		int nActualInstances = 0;
		for(int n = 0; n < nSubClasses; n++)
		{
			StsDrillPathClass subClassInstance = (StsDrillPathClass) currentModel.getStsClass(subClassClasses[n]);
			if(subClassInstance != null) subClasses[nActualInstances++] = subClassInstance;
		}
		subClasses = (StsDrillPathClass[])StsMath.trimArray(subClasses, nActualInstances);
    }
*/
    public void initializeDisplayFields()
	{
 //       initColors(StsDrillPath.pathDisplayFields);

		displayFields = new StsFieldBean[]
		{
			new StsBooleanFieldBean(this, "isVisible", "Visible")
		};
	}

	public void selected(StsDrillPath drillPath)
	{
		super.selected(drillPath);
		setCurrentObject(drillPath);
	}

	public StsDrillPath getCurrentDrillPath()
	{
		return (StsDrillPath) currentObject;
	}

	public boolean setCurrentObject(StsObject object)
	{
		boolean changed = super.setCurrentObject(object);
		if (changed && object != null) ( (StsDrillPath) object).treeObjectSelected();
		return changed;
	}

	public boolean setCurrentSubVolumeName(String name)
	{
		StsDrillPath newDrillPath = (StsDrillPath) getObjectWithName(name);
		return setCurrentObject(newDrillPath);
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
	/*
		if(getClass() == StsDrillPathClass.class) // apply to subTypes
		{
			for(int n = 0; n < subClasses.length; n++)
				if(subClasses[n].setIsVisibleNoDisplay(isVisible)) changed = true;
		}
		else
	*/
		{
			// apply to instances of this subClass
			if(setIsVisibleNoDisplay(isVisible)) changed = true;
			// check superClass to see if its visibility is changed
			StsDrillPathClass drillPathClass = getDrillPathClass();
			drillPathClass.checkSetSuperClassIsVisible(isVisible);
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
			StsDrillPath drillPath = (StsDrillPath) list.getElement(n);
			if(drillPath.setIsVisibleNoDisplay(isVisible)) changed = true;
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
			if (getClass() != StsDrillPathClass.class)
			{
				StsList list = getList();
				for (int n = 0; n < list.getSize(); n++)
				{
					StsDrillPath drillPath = (StsDrillPath) list.getElement(n);
					if (drillPath.getIsVisible())return;
				}
				this.isVisible = false;
			}
		}
		// isVisible flag has been changed for this class: check superClass
		StsDrillPathClass drillPathClass = getDrillPathClass();
		drillPathClass.checkSetSuperClassIsVisible(isVisible);
	}

    public void displayClass(StsGLPanel3d glPanel3d)
    {
        Iterator iter = getVisibleObjectIterator();
        while(iter.hasNext())
        {
            StsDrillPath path = (StsDrillPath)iter.next();
            path.display(glPanel3d);
        }
    }

	static public StsDrillPathClass getDrillPathClass()
	{
		return (StsDrillPathClass)currentModel.getCreateStsClass(StsDrillPath.class);
	}

    public StsDrillPath[] getDrillPaths()
    {
/*
        Object drillPathList;

        StsDrillPath[] drillPaths = new StsDrillPath[0];
        drillPathList = currentModel.getCastObjectList(StsDualSurfaceSubVolume.class);
        drillPaths = (StsDrillPath[])StsMath.arrayAddArray(subVolumes, subVolumeList);
        drillPathList = currentModel.getCastObjectList(StsSingleSurfaceSubVolume.class);
        drillPaths = (StsDrillPath[])StsMath.arrayAddArray(subVolumes, subVolumeList);
        drillPathList = currentModel.getCastObjectList(StsBoxSetSubVolume.class);
        drillPaths = (StsDrillPath[])StsMath.arrayAddArray(subVolumes, subVolumeList);
        return drillPaths;
*/
	return null;
    }

}
