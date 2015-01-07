package com.Sts.DBTypes;

import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Nov 9, 2006
 * Time: 3:31:15 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsPreStackLineSet2dClass extends StsPreStackLineSetClass implements StsSerializable, StsClassTextureDisplayable
{
	public StsPreStackLineSet2dClass()
	{
		super();
        userName = "2D PreStack Sets of Lines";
	}

    public void initializePropertyFields()
	{
		propertyFields = new StsFieldBean[]
		{
			new StsButtonFieldBean("Wiggle Display Properties", "Edit wiggle display properties.", this, "displayWiggleProperties"),
			new StsButtonFieldBean("Semblance Display Properties", "Edit semblanceBytes display properties.", this, "displaySemblanceProperties"),
			new StsButtonFieldBean("AGC Properties", "Edit automatic gain control properties.", this, "displayAGCProperties")
		};
    }

    public void initializeDisplayFields()
	{
        displayFields = new StsFieldBean[]
		{
			new StsIntFieldBean(this, "displayPointSize", true, "Point Size:"),
			new StsBooleanFieldBean(this, "displayAxis", "Plot Axis"),
			new StsBooleanFieldBean(this, "isPixelMode", "Pixel Display Mode"),
			new StsBooleanFieldBean(this, "displayWiggles", "Wiggle Traces"),
			new StsBooleanFieldBean(this, "contourColors", "Contoured Prestack Colors"),
			new StsIntFieldBean(this, "wiggleToPixelRatio", 1, 100, "Wiggle to Pixel Ratio:"),
			new StsFloatFieldBean(this, "maxZToStack", true, "Stacked max Z:")
		};
    }

    public void initializeDefaultFields()
    {
	    defaultFields =	new StsFieldBean[]
		{
			new StsComboBoxFieldBean(this, "seismicSpectrumName", "Seismic Spectrum:", StsSpectrumClass.cannedSpectrums),
			new StsComboBoxFieldBean(this, "semblanceSpectrumName", "Semblance Spectrum:", StsSpectrumClass.cannedSpectrums),
            new StsFloatFieldBean(this, "stackNeighborRadius", true, "Radius of stacks:",true),
            new StsIntFieldBean(this, "gatherIntervalToStack", true, "Gather Interval:",true),
            new StsFloatFieldBean(this, "sampleIntervalToStack", true, "Sample Interval:",true),
            new StsFloatFieldBean(this, "maxZToStack", true, "Maximum Z:",true),
            new StsFloatFieldBean(this, "maxOffsetToStack", true, "Maximum Offset:",true)
		};
	}

    public void displayClass(StsGLPanel3d glPanel3d)
    {
        if (currentObject instanceof StsPreStackLineSet2d)
        {
            StsPreStackLineSet2d lineSet = (StsPreStackLineSet2d) currentObject;
            if (lineSet != null) lineSet.display(glPanel3d);
        }
    }
    
	public boolean setCurrentObject(StsObject object)
	{
		StsPreStackLineSet oldCurrentVolume = (StsPreStackLineSet) currentObject;
		if (!super.setCurrentObject(object)) return false;
		if (oldCurrentVolume != null) oldCurrentVolume.deleteTransients();
    /*
        StsClass prestackVolume3dClass = currentModel.getStsClass(StsPreStackLineSet3d.class);
		if(prestackVolume3dClass != null) prestackVolume3dClass.setCurrentObject(null);
	*/
		return true;
	}
    public void cropChanged() { }

    public boolean textureChanged(StsObject object)
    {
        if(!(object instanceof StsPreStackLineSet2d)) return false;
        return ((StsPreStackLineSet2d)object).textureChanged();
    }
}
