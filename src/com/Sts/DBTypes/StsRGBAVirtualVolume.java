package com.Sts.DBTypes;

/**
 * <p>Title: Workflow development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2001</p>
 * <p>Company: 4D Systems LLC</p>
 * @author unascribed
 * @version 1.0
 */

import com.Sts.Interfaces.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.Utilities.*;

public class StsRGBAVirtualVolume extends StsVirtualVolume implements StsTreeObjectI
{
    static StsObjectPanel virtualVolumeObjectPanel = null;
    static StsEditableColorscaleFieldBean vcolorscaleBean = new StsEditableColorscaleFieldBean(StsRGBAVirtualVolume.class,"colorscale");

    static public final StsFieldBean[] virtualDisplayFields =
        {
        new StsBooleanFieldBean(StsVirtualVolume.class, "isVisible", "Enable"),
        new StsBooleanFieldBean(StsVirtualVolume.class, "readoutEnabled", "Mouse Readout"),
        vcolorscaleBean
    };
    static public final StsFieldBean[] virtualPropertyFields = new StsFieldBean[]
        {
        new StsStringFieldBean(StsRGBAVirtualVolume.class, "name", true, "Name"),
        new StsStringFieldBean(StsRGBAVirtualVolume.class, "volumeOneName", false, "Seismic PostStack3d One"),
        new StsStringFieldBean(StsRGBAVirtualVolume.class, "volumeTwoName", false, "Seismic PostStack3d Two"),
        new StsStringFieldBean(StsRGBAVirtualVolume.class, "volumeThreeName", false, "Seismic PostStack3d Three"),

    };

    public StsRGBAVirtualVolume()
    {
    }

    public StsRGBAVirtualVolume(StsObject[] volumeList, String name)
    {
        StsToolkit.copyDifferentStsObjectAllFields(volumeList[0], this);
		clearNonRelevantMembers();

        setVolumes(volumeList);

        this.type = RGB_BLEND;

        setName(name);
        colorscale = null;        
        initializeColorscale(); 
        initialize(currentModel);
        isVisible = true;
        getRGBAVirtualVolumeClass().setIsVisibleOnCursor(true);
    }
    
	public void setDataHistogram()
	{
		if (dataHist != null && colorscaleBean != null)
        {
            vcolorscaleBean.setHistogram(dataHist);
            vcolorscaleBean.revalidate();
        }
    }
	
    public byte[] processPlaneData(byte[][] planeValues)
    {
        try
        {
            int nPlanes = planeValues.length;
            int nPlanePoints = planeValues[0].length;

            StsVolumeDisplayable[] planeVolumes = new StsVolumeDisplayable[nPlanes];
            for (int n = 0; n < nPlanes; n++)
            {
                planeVolumes[n] = (StsVolumeDisplayable) volumes.getElement(n);
            }

            // Compute value based on these RGB values
            // Help TOM - Need to return the value derived from values[0-3]
            return planeValues[0];
        }
        catch (Exception e)
        {
            StsException.outputException("StsSeismicVolume.processBytePlaneData() failed.",
                                         e, StsException.WARNING);
            return null;
        }
    }
    public StsFieldBean[] getDisplayFields()
    {
        return virtualDisplayFields;
    }

    public StsFieldBean[] getPropertyFields()
    {
        return virtualPropertyFields;
    }

    public StsObjectPanel getObjectPanel()
    {
        if (virtualVolumeObjectPanel == null)
        {
            virtualVolumeObjectPanel = StsObjectPanel.constructor(this, true);
        }
        return virtualVolumeObjectPanel;
    }
    
    static public StsRGBAVirtualVolumeClass getRGBAVirtualVolumeClass()
    {
        return (StsRGBAVirtualVolumeClass) currentModel.getCreateStsClass(StsRGBAVirtualVolume.class);
    }
}
