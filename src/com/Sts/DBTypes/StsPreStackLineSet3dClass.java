package com.Sts.DBTypes;

import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Nov 9, 2006
 * Time: 2:30:36 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsPreStackLineSet3dClass extends StsPreStackLineSetClass
        implements StsSerializable, StsClassCursor3dTextureDisplayable, StsClassCursorDisplayable
{
    protected boolean isVisibleOnCursor = true;
    protected boolean isVisibleOnSection = false;
    protected boolean isVisibleOnSurface = false;
    protected boolean displayOnSubVolumes = true;

	public StsPreStackLineSet3dClass()
	{
		super();
        userName = "3D PreStack Sets of Lines";
	}

    public void initializePropertyFields()
	{
		propertyFields = new StsFieldBean[]
		{
			new StsButtonFieldBean("Wiggle Display Properties", "Edit Wiggle Display Properties.", this, "displayWiggleProperties"),
			new StsButtonFieldBean("Semblance Display Properties", "Edit Semblance Display Properties.", this, "displaySemblanceProperties"),
            new StsButtonFieldBean("CVS/VVS Display Properties", "Edit CVS Display Properties.", this, "displayCVSProperties"),
            new StsButtonFieldBean("AGC Properties", "Edit Automatic Gain Control Properties.", this, "displayAGCProperties"),
            new StsButtonFieldBean("Filter Properties", "Edit Filter properties.", this, "displayFilterProperties"),
            new StsButtonFieldBean("Datum Properties", "Edit Datim properties.", this, "displayDatumProperties")
		};
    }

    public void initializeDisplayFields()
	{
		displayFields = new StsFieldBean[]
		{
			new StsBooleanFieldBean(this, "isVisibleOnCursor", "On 3D Cursors"),
			new StsBooleanFieldBean(this, "isVisibleOnSection", "On Sections"),
			new StsBooleanFieldBean(this, "isVisibleOnSurface", "On Surfaces"),
			new StsIntFieldBean(this, "displayPointSize", true, "Point Size:"),
//			new StsBooleanFieldBean(this, "displayGridLines", "Grid Lines"),
			new StsBooleanFieldBean(this, "displayAxis", "Plot Axis"),
			new StsBooleanFieldBean(this, "isPixelMode", "Pixel Display Mode"),
			new StsBooleanFieldBean(this, "displayWiggles", "Wiggle Traces"),
			new StsBooleanFieldBean(this, "contourColors", "Contoured Prestack Colors"),
			new StsIntFieldBean(this, "wiggleToPixelRatio", 1, 100, "Wiggle to Pixel Ratio:"),
			//tpiBean,
			//new StsFloatFieldBean(this, "inchesPerSecond", true, "Inches Per Second:"),
			new StsFloatFieldBean(this, "maxZToStack", true, "Stacked max Z:")
		};
    }

    public void initializeDefaultFields()
    {
        defaultFields =	new StsFieldBean[]
		{
			new StsComboBoxFieldBean(this, "seismicSpectrumName", "Seismic Spectrum:", StsSpectrumClass.cannedSpectrums),
			new StsComboBoxFieldBean(this, "semblanceSpectrumName", "Semblance Spectrum:", StsSpectrumClass.cannedSpectrums),
			new StsBooleanFieldBean(this, "showVelStat", "Show Velocity Status"),
            new StsFloatFieldBean(this, "stackNeighborRadius", true, "Radius of stacks:",true),
// Not Implemented Yet            new StsIntFieldBean(this, "gatherIntervalToStack", true, "Gather Interval:",true),
// Not Implemented Yet            new StsFloatFieldBean(this, "sampleIntervalToStack", true, "Sample Interval:",true),
            new StsFloatFieldBean(this, "maxZToStack", true, "Maximum Z:",true),
			new StsFloatFieldBean(this, "maxOffsetToStack", true, "Maximum Offset:",true),
// Not Implemented Yet			new StsIntFieldBean(this, "xGatherExtent", true, "Number X Gathers:",true),
// Not Implemented Yet			new StsIntFieldBean(this, "yGatherExtent", true, "Number Y Gathers:",true)
		};
//	    initializeWigglePropertiesPanel();
	}

    public boolean getIsVisibleOnCursor()
	{
		return isVisibleOnCursor;
	}

	public void setIsVisibleOnCursor(boolean isVisibleOnCursor)
	{
		if(this.isVisibleOnCursor == isVisibleOnCursor)return;
		this.isVisibleOnCursor = isVisibleOnCursor;
//		setDisplayField("isVisibleOnCursor", isVisibleOnCursor);
		currentModel.win3dDisplayAll();
	}

	public void toggleVisibleOnCursor()
	{
		isVisibleOnCursor = !isVisibleOnCursor;
		currentModel.win3dDisplayAll();
	}

	public boolean getIsVisibleOnSection()
	{
		return isVisibleOnSection;
	}

	public void setIsVisibleOnSection(boolean isVisible)
	{
		if(this.isVisibleOnSection == isVisible)return;
		this.isVisibleOnSection = isVisible;
//		setDisplayField("isVisibleOnSection", isVisibleOnSection);
		currentModel.win3dDisplayAll();
	}

	public boolean getIsVisibleOnSurface()
	{
		return isVisibleOnSurface;
	}

	public void setIsVisibleOnSurface(boolean isVisible)
	{
		if(this.isVisibleOnSurface == isVisible)return;
		this.isVisibleOnSurface = isVisible;
//		setDisplayField("isVisibleOnSurface", isVisibleOnSurface);
		currentModel.win3dDisplayAll();
	}

    public boolean getDisplayOnSubVolumes()
    {
        return displayOnSubVolumes;
    }

    public void setDisplayOnSubVolumes(boolean b)
    {
        if(this.displayOnSubVolumes == b)return;
        this.displayOnSubVolumes = b;
//        setDisplayField("displayOnSubVolumes", displayOnSubVolumes);
        currentModel.subVolumeChanged();
        currentModel.win3dDisplayAll();
    }

    public StsPreStackLineSet3d getSeismicVolumeDisplayableOnSection()
	{
		if(currentObject == null)return null;
		if(!getIsVisibleOnSection() || !getIsVisible())return null;
		return(StsPreStackLineSet3d)currentObject;
	}

	public StsCursor3dTexture constructDisplayableSection(StsModel model, StsCursor3d cursor3d, int dir)
	{
		return new StsPreStackSeismicCursorSection((StsPreStackLineSet3d)currentObject, model, cursor3d, dir);
	}

    public void readoutOnCursor(StsCursorPoint cursorPoint)
	{
	}

    public void drawOnCursor2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, boolean axesFlipped, boolean xAxisReversed, boolean yAxisReversed)
	{
		Iterator iter = getVisibleObjectIterator();
		while(iter.hasNext())
		{
			StsPreStackLineSet3d psVol = (StsPreStackLineSet3d)iter.next();
			psVol.drawOnCursor2d(glPanel3d, dirNo, axesFlipped);
		}
	}

	public boolean setCurrentObject(StsObject object)
	{
		StsPreStackLineSet oldCurrentVolume = (StsPreStackLineSet) currentObject;
		if (!super.setCurrentObject(object)) return false;
		currentObject = object;
		if (oldCurrentVolume != null) oldCurrentVolume.deleteTransients();
    /*
        StsClass prestackVolume2dClass = currentModel.getStsClass(StsPreStackLineSet2d.class);
		if(prestackVolume2dClass != null) prestackVolume2dClass.setCurrentObject(null);
    */
		return true;
	}

    public void drawOnCursor3d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, StsPoint[] planePoints, boolean isDragging)
	{
	}

	public void drawOn3dCurtain(StsGLPanel3d glPanel3d, StsSeismicCurtain seismicCurtain)
	{
	}
}