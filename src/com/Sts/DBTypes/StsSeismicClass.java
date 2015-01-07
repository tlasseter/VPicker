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
import com.Sts.Types.*;

abstract public class StsSeismicClass extends StsClass implements StsSerializable, StsClassSurfaceDisplayable //, StsClassDisplayable
{
	transient protected StsRotatedGridBoundingBox volume;
	transient protected StsRotatedGridBoundingBox cropVolume;

	protected boolean isPixelMode = false; // Blended Pixels or Nearest
	protected boolean contourColors = true; // shader
	protected boolean displayWiggles = false; // Display wiggles if data density allows - 2D Only
	protected int wiggleToPixelRatio = 4;
	protected String seismicSpectrumName = StsSpectrumClass.SPECTRUM_RWB;
    private boolean fillPlaneNulls = false;

    public StsSeismicClass()
	{
	}

	public void selected(StsSeismicVolume seismicVolume)
	{
		super.selected(seismicVolume);
		setCurrentObject(seismicVolume);
	}
// removed because it is identical to method in superclass which it overrides.  TJL.  2/2/08
/*
    public boolean setCurrentObject(StsObject object)
    {
//        if(currentObject == object) return false;
        currentModel.viewObjectChangedAndRepaint(object);
        currentObject = object;
//        setViewObject(object);
        return true;
    }
*/    
	/**
	 * Set the pixel mode for VSP displays. When set to off, pixels will be interpolated between samples.
	 * When set to on, there will be hard edges between samples.
	 * @param isPixelMode boolean
*/
	public void setIsPixelMode(boolean isPixelMode)
	{
		if (this.isPixelMode == isPixelMode)
			return;
		this.isPixelMode = isPixelMode;
//		setDisplayField("isPixelMode", isPixelMode);
		currentModel.win3dDisplayAll();
	}

	public boolean getIsPixelMode()
	{
		return isPixelMode;
	}

	/**
	 * Turn wiggle display on and off. When on it will use the wiggle to pixel ratio to determine when
	 * wigles are isVisible. Currently only works for 2D displays
	 * @param displayWiggles boolean
*/
	public void setDisplayWiggles(boolean displayWiggles)
	{
		if (this.displayWiggles == displayWiggles)
			return;
		this.displayWiggles = displayWiggles;
//		setDisplayField("displayWiggles", displayWiggles);
		currentModel.win3dDisplayAll();
	}

	public boolean getDisplayWiggles()
	{
		return displayWiggles;
	}

	public void setContourColors(boolean contourColors)
	{
		if (this.contourColors == contourColors)
			return;
		this.contourColors = contourColors;
//		setDisplayField("contourColors", contourColors);
		currentModel.win3dDisplayAll();
	}

	public boolean getContourColors()
	{
		return contourColors;
	}


	/**
	 * If the display wiggle flag is set to on, wiggle traces will be plotted on top of the
	 * texture maps and when they are isVisible is based on the wiggle to pixel ratio. The wiggle
	 * to pixel ratio is the number of pixels between traces before wiggles will be plotted. This is
	 * no avoid ridiculous plots where the number of traces plotted exceeds number of pixles resulting in
	 * solid black plots.
	 * @param wiggleToPixelRatio int
*/
	public void setWiggleToPixelRatio(int wiggleToPixelRatio)
	{
		if (this.wiggleToPixelRatio == wiggleToPixelRatio)
			return;
		this.wiggleToPixelRatio = wiggleToPixelRatio;
//		setDisplayField("wiggleToPixelRatio", wiggleToPixelRatio);
		currentModel.win3dDisplayAll();
	}

	public int getWiggleToPixelRatio()
	{
		return wiggleToPixelRatio;
	}
/*
	public void setViewObject(StsObject object)
	{
		if (currentModel.win3d != null)
			currentModel.win3d.glPanel3d.setObject(object);
	}
*/
/*
	  public boolean setCurrentSeismicVolumeName(String name)
	  {
		 StsSeismicVolume newSeismicVolume = (StsSeismicVolume)getObjectWithName(name);
		 return setCurrentObject(newSeismicVolume);
	  }

	  public StsSeismicVolume getSeismicVolumeDisplayableOnSurface()
	  {
		 if(currentObject == null) return null;
		 if(!isVisibleOnSurface || !isVisible) return null;
		 return (StsSeismicVolume)currentObject;
	  }
*/

/*
	  public void displayOnCursor(StsCursor3d.StsCursorSection cursorSection, StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, boolean is3d)
	  {
		 if(currentObject == null || !((StsSeismicVolume)currentObject).getIsVisibleOnCursor()) return;
		 StsSeismicCursorSection seismicCursorSection = cursorSection.hasCursor3dDisplayable(StsSeismicCursorSection.class);
		 if(seismicCursorSection == null)
		 {
		 seismicCursorSection = new StsSeismicCursorSection((StsSeismicVolume)currentObject, glPanel3d, dirNo, dirCoordinate);
		 cursorSection.addDisplayable(seismicCursorSection);
		 }
		 seismicCursorSection.display(this, glPanel3d, is3d);
	  }
*/
/*
	  public void displayOnCursor2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate)
	  {
		 if(currentObject == null || !((StsSeismicVolume)currentObject).getIsVisibleOnCursor()) return;
		 GL gl = glPanel3d.getGL();
		 GLU glu = glPanel3d.getGLU();
//        currentSeismicVolume.displayOnCursor2d(dirNo, dirCoordinate, gl, glu);
	  }
*/

/*
	  public boolean getApplyCrop() { return applyCrop; }
	  public void setApplyCrop(boolean crop)
	  {
		 if(crop == applyCrop) return;
		 applyCrop = crop;
		 currentModel.getProject().getCropVolume().valueChanged();
		 currentModel.clearDataDisplays();
		 currentModel.win3dDisplayAll();
	  }
*/


	public void setSeismicSpectrumName(String seismicSpectrumName)
	{
		if (this.seismicSpectrumName == seismicSpectrumName)
			return;
		this.seismicSpectrumName = seismicSpectrumName;
//		setDisplayField("seismicSpectrumName", seismicSpectrumName);
		currentModel.win3dDisplayAll();
	}

	public String getSeismicSpectrumName()
	{
		return seismicSpectrumName;
	}

	public void close()
	{
		list.forEach("close");
	}

	public boolean drawLast()
	{
		return false;
	}

    public boolean getFillPlaneNulls()
    {
        return fillPlaneNulls;
    }

    public void setFillPlaneNulls(boolean fillPlaneNulls)
    {
        this.fillPlaneNulls = fillPlaneNulls;
        currentModel.win3dDisplayAll();
    }
}
