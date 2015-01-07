package com.Sts.SeismicAttributes;

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

/**
 *
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002-2003</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version beta 1.0
 */

/** Given a set of inputVolumes required for constructing this seismic attribute, construct the attribute volume.
 *  Currently implementd for Hilbert attributes which are subclasses of this abstract class.
 */
public abstract class StsVirtualVolumeConstructor extends StsSeismicVolumeConstructor
{
    public StsVirtualVolume virtualVolume;
    public String vvName = "vvContructor";
    boolean inputIsFloat = false;

    abstract public  boolean doProcessInputBlock(int nBlock, StsMappedBuffer[] inputBuffers);

    public StsVirtualVolumeConstructor()
	{
	}

	public StsVirtualVolumeConstructor(StsModel model, StsSeismicVolume[] inputVolumes, StsVirtualVolume virtualVolume)
	{
        this.model = model;
        this.inputVolumes = inputVolumes;
        this.virtualVolume = virtualVolume;
		this.volumeName = virtualVolume.getName();
	}

    public void initialize(StsSeismicVolume[] inputVolumes, StsVirtualVolume virtualVolume, String vvName)
    {
        this.inputVolumes = inputVolumes;
        this.virtualVolume = virtualVolume;        
        this.volumeName = vvName;
   }

   public void initialize(String vvName)
   {
       this.inputVolumes = null;
       this.virtualVolume = null;       
       this.volumeName = vvName;
   }

    static public StsSeismicVolume createVirtualVolume(StsModel model, StsSeismicVolume[] inputVolumes, StsVirtualVolume volume,
		boolean isDataFloat, StsProgressPanel panel)
	{
		StsVirtualVolumeConstructor attributeVolume = null;
		StsSeismicVolume seismicVolume;
		if(panel != null) panel.appendLine("Creating " + volume.getName() + " volume.");
		seismicVolume = checkGetAttributeVolume(model, volume.getName());
		if (seismicVolume != null) 
			return seismicVolume;

		switch(volume.getType())
		{
			case StsVirtualVolume.SEISMIC_MATH:
				attributeVolume = StsVirtualVolumeMathConstructor.constructor(model, inputVolumes, volume, isDataFloat, panel);
				break;
			case StsVirtualVolume.SEISMIC_XPLOT_MATH:
				// Not required for xplot volume since the vv will have the same data range as the original. On/off is all
				break;
			case StsVirtualVolume.SEISMIC_BLEND:
				attributeVolume = StsVirtualVolumeBlendedConstructor.constructor(model, inputVolumes, volume, isDataFloat, panel);
				break;
			case StsVirtualVolume.SENSOR_VOLUME:
				// Not implemented yet.
				//attributeVolume = StsVirtualVolumeSensorConstructor.constructor(model, inputSensors, volume, isDataFloat, panel);
				break;				
			case StsVirtualVolume.SEISMIC_FILTER:
				// Not implemented yet.
				//attributeVolume = StsVirtualVolumeFilterConstructor.constructor(model, inputVolumes, volume, isDataFloat, panel);;
				break;
		}

		seismicVolume = attributeVolume.getVolume();
		if(!seismicVolume.initialize(model))
		{
			seismicVolume.delete();
			return null;
		}
		if(panel != null) panel.appendLine("Successfully created " + volume.getName() + " volume.");		
		return seismicVolume;
	}

	static public StsSeismicVolume checkGetAttributeVolume(StsModel model, String volumeName)
	{
		StsSeismicVolume attributeVolume;
		try
		{
			attributeVolume = (StsSeismicVolume) model.getObjectWithName(StsSeismicVolume.class, volumeName);
			if (attributeVolume == null) return null;
			boolean deleteVolume = StsYesNoDialog.questionValue(model.win3d, "Volume " + volumeName + " already loaded. Delete and recreate?");
			if (!deleteVolume) return attributeVolume;
			attributeVolume.delete();
			attributeVolume = null;
			return null;
		}
		catch (Exception e)
		{
			StsException.outputException("StsVirtualVolumeConstructor.checkVolumeExistence() failed.", e, StsException.WARNING);
			return null;
		}
	}

	static public StsSeismicVolume getExistingVolume(StsModel model, String volumeName)
	{
		try
		{
			return (StsSeismicVolume) model.getObjectWithName(StsSeismicVolume.class, volumeName);
		}
		catch (Exception e)
		{
			StsException.outputException("StsVirtualVolumeConstructor.checkVolumeExistence() failed.", e, StsException.WARNING);
			return null;
		}
	}

    public void createVolume(StsSeismicVolume[] inputVolumes, StsVirtualVolume vVolume)
	{
		this.inputVolumes = inputVolumes;
		this.virtualVolume = vVolume;
        createOutputVolume();
    }
}
