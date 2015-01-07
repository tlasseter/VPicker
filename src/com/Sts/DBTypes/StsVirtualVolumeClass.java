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
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Utilities.*;

public class StsVirtualVolumeClass extends StsSeismicVolumeClass implements StsSerializable, StsClassSurfaceDisplayable, StsClassCursor3dTextureDisplayable //, StsClassDisplayable
{
    static public StsVirtualVolume currentVirtualVolumeObject = null;
/*
    static final Class[] subClassClasses = new Class[] { StsMathVirtualVolume.class, StsCrossplotVirtualVolume.class,
        StsRGBAVirtualVolume.class, StsFilterVirtualVolume.class, StsBlendedVirtualVolume.class, StsSensorVirtualVolume.class};
	static StsVirtualVolumeClass[] subClasses = null;
*/
    public StsVirtualVolumeClass()
    {
        userName = "Virtual Volumes";
    }

    public void initializeFields()
    {
        //initializeSubClasses();
    }
/*
    private void initializeSubClasses()
    {
        int nSubClasses = subClassClasses.length;
        subClasses = new StsVirtualVolumeClass[nSubClasses];
        int nActualInstances = 0;
        for(int n = 0; n < nSubClasses; n++)
        {
            StsVirtualVolumeClass subClassInstance = (StsVirtualVolumeClass) currentModel.getStsClass(subClassClasses[n]);
            if(subClassInstance != null) subClasses[nActualInstances++] = subClassInstance;
        }
        subClasses = (StsVirtualVolumeClass[])StsMath.trimArray(subClasses, nActualInstances);
    }
*/
    public StsVirtualVolume[] getVirtualVolumes()
    {
        Object virtualVolumeList;
        StsVirtualVolume[] virtualVolumes = new StsVirtualVolume[0];

        if(subClasses == null) return new StsVirtualVolume[0];
        for(StsClass subClass : subClasses)
        {
            virtualVolumeList = subClass.getCastObjectList();
            virtualVolumes = (StsVirtualVolume[])StsMath.arrayAddArray(virtualVolumes, virtualVolumeList);
        }
        return virtualVolumes;
    }
    
    public StsVirtualVolume getVirtualVolumeWithName(String name)
    {
        StsVirtualVolume volume = null;
        if(subClasses == null) return volume;
        for(StsClass subClass : subClasses)
        {
            volume = (StsVirtualVolume) subClass.getObjectWithName(name);
            if(volume != null) break;
        }
        return volume;
    }
    public void setIsVisibleOnCursor(boolean isVisibleOnCursor)
    {
        if(subClasses == null) return;
        for(StsClass subClass : subClasses)
        {
            StsVirtualVolumeClass vvClass = (StsVirtualVolumeClass) subClass;
            vvClass.setIsVisibleOnCursor(isVisibleOnCursor);
        }
    }
	public boolean getContourColors() {	return contourColors; }
    public void setContourColors(boolean contour)
    {
        if(subClasses == null) return;
        for(StsClass subClass : subClasses)
        {
            StsVirtualVolumeClass vvClass = (StsVirtualVolumeClass) subClass;
            vvClass.setContourColors(contour);
        }
    }
    public void selected(StsVirtualVolume virtualVolume)
    {
        super.selected(virtualVolume);
        setCurrentObject(virtualVolume);
    }
    
 	public boolean getIsVisibleOnCursor()
	{
		return isVisibleOnCursor;
	}

    public boolean setCurrentObject(StsObject object)
    {
        currentVirtualVolumeObject = (StsVirtualVolume)object;
        return super.setCurrentObject(object);
    }

    public void toggleOn(StsWin3dBase win3d)
    {
        if(debug) System.out.println("toggleOn called.");
        win3d.getCursor3d().toggleOn(currentVirtualVolumeObject);
        win3d.repaint();
    }

    public void toggleOff(StsWin3dBase win3d)
    {
        if(debug) System.out.println("toggleOff called.");
        win3d.getCursor3d().toggleOff(currentVirtualVolumeObject);
        win3d.repaint();
    }

    public void toggleOn(StsWin3d win3d)
    {
        if(debug) System.out.println("toggleOn called.");
        win3d.getCursor3d().toggleOn(currentVirtualVolumeObject);
        win3d.repaint();
    }

    public void toggleOff(StsWin3d win3d)
    {
        if(debug) System.out.println("toggleOff called.");
        win3d.getCursor3d().toggleOff(currentVirtualVolumeObject);
        win3d.repaint();
    }

    public void toggleOn(StsWin3dFull win3d)
    {
        if(debug) System.out.println("toggleOn called.");
        win3d.getCursor3d().toggleOn(currentVirtualVolumeObject);
        win3d.repaint();
    }

    public void toggleOff(StsWin3dFull win3d)
    {
        if(debug) System.out.println("toggleOff called.");
        win3d.getCursor3d().toggleOff(currentVirtualVolumeObject);
        win3d.repaint();
    }

    public StsCursor3dTexture constructDisplayableSection(StsModel model, StsCursor3d cursor3d, int dir)
    {
        return new StsSeismicCursorSection(model, (StsVirtualVolume)currentObject, cursor3d, dir);
    }

	public void setDisplayOnSubVolumes(boolean displayOnSubVolumes)
	{
		if(this.displayOnSubVolumes == displayOnSubVolumes) return;
        if(subClasses == null) return;
		boolean changed = false;
        for(StsClass subClass : subClasses)
        {
            StsVirtualVolumeClass vvClass = (StsVirtualVolumeClass) subClass;
            changed = changed | vvClass.setSubclassDisplayOnSubVolumes(displayOnSubVolumes);
        }
		this.displayOnSubVolumes = displayOnSubVolumes;
		if(!changed) return;
		currentModel.subVolumeChanged();
		currentModel.win3dDisplayAll();
	}

	public boolean setSubclassDisplayOnSubVolumes(boolean displayOnSubVolumes)
	{
		if(this.displayOnSubVolumes == displayOnSubVolumes) return false;
		this.displayOnSubVolumes = displayOnSubVolumes;
		return true;
	}
}
