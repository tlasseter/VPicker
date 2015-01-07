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

public class StsFilterVirtualVolume extends StsVirtualVolume implements StsTreeObjectI
{
    protected StsVolumeFilterFace filter = null;

    static StsObjectPanel virtualVolumeObjectPanel = null;
    static StsEditableColorscaleFieldBean vcolorscaleBean = new StsEditableColorscaleFieldBean(StsFilterVirtualVolume.class,"colorscale");

    static public final StsFieldBean[] virtualDisplayFields =
    {
        new StsBooleanFieldBean(StsVirtualVolume.class, "isVisible", "Enable"),
        new StsBooleanFieldBean(StsVirtualVolume.class, "readoutEnabled", "Mouse Readout"),
            vcolorscaleBean
    };
    static public final StsFieldBean[] virtualPropertyFields = new StsFieldBean[]
        {
        new StsStringFieldBean(StsFilterVirtualVolume.class, "name", true, "Name"),
        new StsStringFieldBean(StsFilterVirtualVolume.class, "volumeOneName", false, "Seismic PostStack3d One"),

    };

    public StsFilterVirtualVolume()
    {
    }

    public StsFilterVirtualVolume(StsSeismicVolume volume, StsVolumeFilterFace filter, String name)
    {
        super(false);
        StsToolkit.copyDifferentStsObjectAllFields(volume, this);
        setName(name);
        clearNonRelevantMembers();
        setVolume(volume);
        setName(name);
        this.type = SEISMIC_FILTER;
        this.filter = filter;
        colorscale = null;
        initializeColorscale();
        initialize(currentModel);
        isVisible = true;
        getFilterVirtualVolumeClass().setIsVisibleOnCursor(true);
        isDataFloat = true;
        addToModel();
    }
    
	public void setDataHistogram()
	{
		if (dataHist != null && colorscaleBean != null)
        {
            vcolorscaleBean.setHistogram(dataHist);
            vcolorscaleBean.revalidate();
        }
    }
	
    public byte[] readBytePlaneData(int dir, float dirCoordinate)
    {
        if(volumes == null) return null;
        StsSeismicVolume volume = (StsSeismicVolume)volumes.getFirst();
        int nPlane = volume.getCursorPlaneIndex(dir, dirCoordinate);
        return filter.processBytePlaneData(dir, nPlane, volume, this);
    }

	public byte[] readRowPlaneByteData(int nPlane)
	{
        StsSeismicVolume volume = (StsSeismicVolume)volumes.getFirst();
        return filter.processBytePlaneData(YDIR, nPlane, volume, this);
	}

	public float[] readRowPlaneFloatData(int nPlane)
	{
        StsSeismicVolume volume = (StsSeismicVolume)volumes.getFirst();
        return filter.processFloatPlaneData(YDIR, nPlane, volume, this);
	}

    /** included to satisfy superclass abstract method requirement. */
    public byte[] processPlaneData(byte[][] planeData)
    {
        return null;
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
    
    static public StsFilterVirtualVolumeClass getFilterVirtualVolumeClass()
    {
        return (StsFilterVirtualVolumeClass) currentModel.getCreateStsClass(StsFilterVirtualVolume.class);
    }
}
