package com.Sts.Interfaces;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version beta 1.0
 */

import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

abstract public class StsCursorSectionDisplayable
{
    // These four members are initialized by initializeTransients call from viewPersistManager
    transient protected StsGLPanel3d glPanel3d;
    transient public StsModel model;
    transient public int dirNo;
    transient public float dirCoordinate = -StsParameters.largeFloat;
    transient public StsRotatedGridBoundingBox cursorBoundingBox;

    transient public int nTextureRows, nTextureCols;

    abstract public boolean setDirCoordinate(float dirCoordinate);
//    public void clearTextureDisplay();
    abstract public void cropChanged();
    abstract public boolean isDisplayableObject(Object object);
    abstract public boolean isDisplayingObject(Object object);
    abstract public void displayTexture(StsGLPanel3d glPanel3d, boolean is3d, byte[] subVolumePlane);
    abstract public void display(StsGLPanel3d glPanel3d, boolean is3d);
    abstract public String propertyReadout(StsPoint point);
    abstract public void setObject(Object object);
    abstract public Object getObject();
    abstract public Class getDisplayableClass();
    abstract public boolean canDisplayClass(Class c);

    public boolean initialize(StsCursorSection cursorSection, StsGLPanel3d glPanel3d)
    {
//        this.glPanel3d = cursorSection.cursor3d.glPanel3d;
        this.dirNo = cursorSection.dir;
        this.dirCoordinate = cursorSection.dirCoordinate;
        this.cursorBoundingBox = cursorSection.rotatedBoundingBox;
        this.glPanel3d = glPanel3d;
        return true;
    }   

    protected boolean isInRange(int row, int col)
    {
        if(row < 0 || row >= nTextureRows) return false;
        if(col < 0 || col >= nTextureCols) return false;
        return true;
    }

	protected boolean isInRange(float rowF, float colF)
	{
		if(rowF < 0 || rowF >= nTextureRows) return false;
		if(colF < 0 || colF >= nTextureCols) return false;
		return true;
	}
}
