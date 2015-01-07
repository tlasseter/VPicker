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
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import java.util.ArrayList;

public class StsSubVolumeClass extends StsClass implements StsSerializable
{

	boolean isApplied = false; // indicates subVolume is applied to cursor display
    boolean isUnion = true; // indicates whether union or intersection

	//static final Class[] subClassClasses = new Class[] { StsBoxSetSubVolume.class, StsDualSurfaceSubVolume.class, StsWellSetSubVolume.class };
	//static StsSubVolumeClass[] subClasses = null;

	public StsSubVolumeClass()
	{
	}

    public void initializeFields()
    {
        //initializeSubClasses();
    }
/*
    private void initializeSubClasses()
    {
		int nSubClasses = subClassClasses.length;
		subClasses = new StsSubVolumeClass[nSubClasses];
		int nActualInstances = 0;
		for(int n = 0; n < nSubClasses; n++)
		{
			StsSubVolumeClass subClassInstance = (StsSubVolumeClass) currentModel.getStsClass(subClassClasses[n]);
			if(subClassInstance != null) subClasses[nActualInstances++] = subClassInstance;
		}
		subClasses = (StsSubVolumeClass[])StsMath.trimArray(subClasses, nActualInstances);
    }
*/
    public void initializeDisplayFields()
	{
		displayFields = new StsFieldBean[]
            {
			// new StsBooleanFieldBean(this, "isVisible", "Visible:"),
			new StsBooleanFieldBean(this, "isApplied", "Applied"),
            new StsBooleanFieldBean(this, "isUnion", "Union")

		};
	}

	public void selected(StsSubVolume subVolume)
	{
		super.selected(subVolume);
		setCurrentObject(subVolume);
	}

	public StsSubVolume getCurrentSubVolume()
	{
		return (StsSubVolume) currentObject;
	}

    /** if this is not the currentObject, set it (changed = true).
     *  if changed, and object is != null set this as the treeObject selected.
     *  return true if changed.
     * @param object
     * @return true if currentObject is changed
     */
    public boolean setCurrentObject(StsObject object)
	{
		if(!super.setCurrentObject(object)) return false;
//        if(object != null) ( (StsSubVolume) object).treeObjectSelected();
		return true;
	}

	public boolean setCurrentSubVolumeName(String name)
	{
		StsSubVolume newSubVolume = (StsSubVolume) getObjectWithName(name);
		return setCurrentObject(newSubVolume);
	}

	public void close()
	{
		list.forEach("close");
	}

	/** called to set this superClass or one of its subTypes with
	 *  new isApplied value.  If superClass, subTypes are called and set.
	 */
	public void setIsApplied(boolean isApplied)
	{
        if(this.isApplied == isApplied) return;
        this.isApplied = isApplied;
 //      setDisplayField("isApplied", isApplied);

		boolean changed = false;
		ArrayList<StsClass> subClasses = getSubClasses();
		if(subClasses != null)
		{
			for(StsClass subClass : subClasses)
				changed = changed | ((StsSubVolumeClass)subClass).setIsAppliedNoChange(isApplied);
		}
		else
		{
			// apply to instances of this subClass
			if(setIsAppliedNoChange(isApplied)) changed = true;
			// check superClass to see if its visibility is changed
			StsSubVolumeClass subVolumeClass = getSubVolumeClass();
			subVolumeClass.checkSetSuperClassIsApplied(isApplied);
		}
		if(changed) currentModel.subVolumeChanged();
        currentModel.win3dDisplay();
	}

    public void setIsUnion(boolean isUnion)
    {
        if(this.isUnion == isUnion) return;
        this.isUnion = isUnion;
 //      setDisplayField("isUnion", isUnion);
        currentModel.subVolumeChanged();
        currentModel.win3dDisplay();
    }
	/** Called only on the superclass StsSubVolumeClass.
	 *  A subClass isApplied has changed:
	 *  if subClass is true, then this superClass must be true
	 *  to be consistent.  Otherwise if all subTypes are false,
	 *  then set superClass to false..
	 */
	public void checkSetSuperClassIsApplied(boolean isApplied)
	{
		if(subClasses != null)
		{
			for (StsClass subClass : subClasses)
				if (((StsSubVolumeClass) subClass).getIsApplied()) return;
			this.isApplied = false;
		}
		else if(isApplied)
			this.isApplied = true;
	}

	/** Called to set subClass instances applied value.
	 *  If no instances are changed because they don't exist
	 *  then change is false is returned.
	 */
	private boolean setIsAppliedNoChange(boolean isApplied)
	{
		this.isApplied = isApplied;
		boolean changed = false;
		StsList list = getList();
		for(int n = 0; n < list.getSize(); n++)
		{
			StsSubVolume subVolume = (StsSubVolume)list.getElement(n);
//			if(subVolume.setIsAppliedNoChange(isApplied)) changed = true;
            subVolume.setIsApplied(isApplied);
            changed = true;
		}
		return changed;
	}

    /** called to set this superClass or one of its subTypes with
	 *  new isVisible value.  If superClass, subTypes are called and set.
	 */
	public void setIsVisible(boolean isVisible)
	{
		if (this.isVisible == isVisible) return;
		this.isVisible = isVisible;
 //      setDisplayField("isVisible", isVisible);
		boolean changed = false;
		if(getClass() == StsSubVolumeClass.class) // apply to subTypes
		{
			if(subClasses != null)
			{
				for(StsClass subClass : subClasses)
					if(((StsSubVolumeClass)subClass).setIsVisibleNoDisplay(isVisible)) changed = true;
			}
		}
		else
		{
			// apply to instances of this subClass
			if(setIsVisibleNoDisplay(isVisible)) changed = true;
			// check superClass to see if its visibility is changed
			StsSubVolumeClass subVolumeClass = getSubVolumeClass();
			subVolumeClass.checkSetSuperClassIsVisible(isVisible);
		}
		if(changed) currentModel.win3dDisplay();
	}

	/** Called only on the superclass StsSubVolumeClass.
	 *  A subClass visibility has changed to isVisible;
	 *  if subClass is true, then this superClass must be true
	 *  to be consistent.  Otherwise if all subTypes are false,
	 *  then set superClass to false..
	 */
	public void checkSetSuperClassIsVisible(boolean isVisible)
	{
		if(isVisible)
			this.isVisible = true;
		else if(subClasses != null)
		{
			for (StsClass subClass : subClasses)
				if (subClass.getIsVisible()) return;
			this.isVisible = false;
		}
	}

	/** Called to set subClass instances visible values. */
	private boolean setIsVisibleNoDisplay(boolean isVisible)
	{
		this.isVisible = isVisible;
		boolean changed = false;
		StsList list = getList();
		for (int n = 0; n < list.getSize(); n++)
		{
			StsSubVolume subVolume = (StsSubVolume) list.getElement(n);
			if(subVolume.setIsVisibleNoDisplay(isVisible)) changed = true;
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
			if (getClass() != StsSubVolumeClass.class)
			{
				StsList list = getList();
				for (int n = 0; n < list.getSize(); n++)
				{
					StsSubVolume subVolume = (StsSubVolume) list.getElement(n);
					if (subVolume.getIsVisible())return;
				}
				this.isVisible = false;
			}
		}
		// isVisible flag has been changed for this class: check superClass
		StsSubVolumeClass subVolumeClass = getSubVolumeClass();
		subVolumeClass.checkSetSuperClassIsVisible(isVisible);
	}

	static public StsSubVolumeClass getSubVolumeClass()
	{
		return (StsSubVolumeClass)currentModel.getCreateStsClass(StsSubVolume.class);
	}

	/** An instance has been toggled to state isApplied.
	 * If this is a StsSubVolumeClass subClass instance,
	 * check and set toggle and then move up to superClass
	 * and do the set and check there.
	 */
	public void checkSetClassIsApplied(boolean isApplied)
	{
		if (isApplied) // if object below isVisible, then this must be visible as well
		{
			if(this.isApplied) return;
			this.isApplied = true;
		}
		else
		{
			if (getClass() != StsSubVolumeClass.class)
			{
				StsList list = getList();
				for (int n = 0; n < list.getSize(); n++)
				{
					StsSubVolume subVolume = (StsSubVolume) list.getElement(n);
					if (subVolume.getIsApplied())return;
				}
				this.isApplied = false;
			}
		}
		// isVisible flag has been changed for this class: check superClass
		StsSubVolumeClass subVolumeClass = getSubVolumeClass();
		subVolumeClass.checkSetSuperClassIsApplied(isApplied);
	}
    public boolean getIsUnion() { return isUnion; }
	public boolean getIsApplied() { return isApplied; }

    public byte[] getSubVolumePlane(int dir, float dirCoordinate, StsRotatedGridBoundingBox rotatedBoundingBox, byte zDomainData)
    {
 //       isVisible = currentModel.getBooleanProperty(getClass().getName() + ".isVisible");
 //       if(!isVisible) return null;
        StsSubVolume[] subVolumes = getSubVolumes();
        if(subVolumes.length == 0) return null;

        boolean noneApplied = true;
        for(int n = 0; n < subVolumes.length; n++)
            if(subVolumes[n].getIsApplied()) noneApplied = false;
        if(noneApplied) return null;

//        StsGridBoundingBox boundingBox = new StsGridBoundingBox();

        int nRows = rotatedBoundingBox.getNCursorRows(dir);
        int nCols = rotatedBoundingBox.getNCursorCols(dir);
        byte[] subVolumePlane = new byte[nRows*nCols];

        if(isUnion)
        {
            for(int n = 0; n < subVolumes.length; n++)
                subVolumes[n].addUnion(subVolumePlane, dir, dirCoordinate, rotatedBoundingBox, zDomainData);
        }
        else
        {
            byte[][] subVolumePlanes = new byte[subVolumes.length][];
            for(int n = 0; n < subVolumes.length; n++)
            {
                subVolumePlanes[n] = new byte[nRows*nCols];
                subVolumes[n].addUnion(subVolumePlanes[n], dir, dirCoordinate, rotatedBoundingBox, zDomainData);
            }
            for(int i=0; i<nRows*nCols; i++)
            {
                subVolumePlane[i] = 1;
                for(int j = 0; j < subVolumePlanes.length; j++)
                {
                    if(subVolumePlanes[j][i] == 0)
                        subVolumePlane[i] = 0;
                }
            }
        }
        return subVolumePlane;
    }

    public StsSubVolume[] getSubVolumes()
    {
        Object subVolumeList;

        StsSubVolume[] subVolumes = new StsSubVolume[0];
        subVolumeList = currentModel.getCastObjectList(StsDualSurfaceSubVolume.class);
        subVolumes = (StsSubVolume[])StsMath.arrayAddArray(subVolumes, subVolumeList);
        subVolumeList = currentModel.getCastObjectList(StsBoxSetSubVolume.class);
        subVolumes = (StsSubVolume[])StsMath.arrayAddArray(subVolumes, subVolumeList);
        return subVolumes;
    }

    public StsGridBoundingBox getBoundingBox()
    {
        if(!isApplied) return null;
        StsGridBoundingBox boundingBox = new StsGridBoundingBox();
        StsSubVolume[] subVolumes = getSubVolumes();
        for(int n = 0; n < subVolumes.length; n++)
        {
            StsGridBoundingBox subVolumeBoundingBox = subVolumes[n].getGridBoundingBox();
            if(subVolumeBoundingBox != null) boundingBox.addBoundingBox(subVolumeBoundingBox);
        }
        return boundingBox;
    }
}
