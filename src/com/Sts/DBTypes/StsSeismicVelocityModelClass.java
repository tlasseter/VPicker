package com.Sts.DBTypes;

import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

//public class StsSeismicVelocityModelClass extends StsClass implements DBSerializable
public class StsSeismicVelocityModelClass extends StsSeismicVolumeClass implements StsSerializable, StsClassSurfaceDisplayable //, StsClassCursor3dTextureDisplayable
{
	private String velocityDisplay = V_AVG_STRING;
	static public final String V_AVG_STRING = StsParameters.V_AVG_STRING;
	static public final String V_INSTANT_STRING = StsParameters.V_INSTANT_STRING;
	static String[] velocityDisplays = new String[] { V_AVG_STRING, V_INSTANT_STRING };

    public StsSeismicVelocityModelClass()
    {
    }

	public void initializeDisplayFields()
	{
		//StsComboBoxFieldBean velocityDisplayBean = new StsComboBoxFieldBean(this, "velocityDisplay", "Velocity Display", velocityDisplays);
		//displayFields = new StsFieldBean[] { velocityDisplayBean };
	}

	/** Set the current velocity display type; if changed, clear cached data from seismic volumes of type VELOCITY */
	public void setVelocityDisplay(String velocityDisplay)
	{

		if(this.velocityDisplay == velocityDisplay) return;
		this.velocityDisplay = velocityDisplay;
		StsObject[] seismicVolumes = currentModel.getObjectList(StsSeismicVolume.class);
		for(int n = 0; n < seismicVolumes.length; n++)
		{
			StsSeismicVolume volume = (StsSeismicVolume)seismicVolumes[n];
			if(volume.type == StsParameters.VELOCITY) volume.clearCache();
		}
	}

    /** there should be only one instance of velocityModel in the project:  return it */
    public StsSeismicVelocityModel getInstance()
    {
        return (StsSeismicVelocityModel)getLast();
    }

    public String getVelocityDisplay() { return velocityDisplay; }
 /*
    public boolean setCurrentObject(StsObject object)
    {
        if(currentObject == object) return false;
        currentObject = object;
        setViewObject(object);
        return true;
    }
*/
}
