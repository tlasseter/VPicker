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

public class StsCrossplotVirtualVolume extends StsVirtualVolume implements StsTreeObjectI
{
    // these members are persistent, but not loaded from seis3d.txt.name file
    protected StsCrossplot crossplot = null;
    protected boolean crossplotInclusive = true;

    static StsObjectPanel virtualVolumeObjectPanel = null;
    static StsEditableColorscaleFieldBean vcolorscaleBean = new StsEditableColorscaleFieldBean(StsCrossplotVirtualVolume.class,"colorscale");

    static public final StsFieldBean[] virtualDisplayFields =
        {
        new StsBooleanFieldBean(StsVirtualVolume.class, "isVisible", "Enable"),
        new StsBooleanFieldBean(StsVirtualVolume.class, "readoutEnabled", "Mouse Readout"),
                vcolorscaleBean
    };
    static public final StsFieldBean[] virtualPropertyFields = new StsFieldBean[]
        {
        new StsStringFieldBean(StsCrossplotVirtualVolume.class, "name", true, "Name"),
        new StsStringFieldBean(StsCrossplotVirtualVolume.class, "volumeOneName", false, "Seismic PostStack3d One"),
        new StsStringFieldBean(StsCrossplotVirtualVolume.class, "crossplotName", false, "Crossplot"),
        new StsStringFieldBean(StsCrossplotVirtualVolume.class, "isInclusive", false, "Inclusive"),      
    };

    public StsCrossplotVirtualVolume()
    {
    }
    
    public StsCrossplotVirtualVolume(boolean persistent)
    {
        super(persistent);    	
    }
    
    public StsCrossplotVirtualVolume(StsSeismicVolume volume, StsCrossplot crossplot, String name, boolean isInclusive)
    {
        StsToolkit.copyDifferentStsObjectAllFields(volume, this);
		clearNonRelevantMembers();

        setVolumes(new StsObject[] {volume, crossplot});

        this.crossplot = crossplot;
        crossplotInclusive = isInclusive;
        this.type = SEISMIC_XPLOT_MATH;

        setName(name);
        colorscale = null;
        initializeColorscale();
        initialize(currentModel);
        isVisible = true;
        getCrossplotVirtualVolumeClass().setIsVisibleOnCursor(true);
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
            if(nPlanes == 2)
            {
                StsVolumeDisplayable volume0 = planeVolumes[0];
                StsVolumeDisplayable volume1 = planeVolumes[1];

                byte[] planePoints = new byte[nPlanePoints];
                if (crossplotInclusive)
                {
                    for (int n = 0; n < nPlanePoints; n++)
                    {
                        if (volume1.isByteValueNull(planeValues[1][n]))
                            planePoints[n] = -1;
                        else
                            planePoints[n] = planeValues[0][n];
                    }
                }
                else
                {
                    for (int n = 0; n < nPlanePoints; n++)
                    {
                        if (volume1.isByteValueNull(planeValues[1][n]))
                            planePoints[n] = planeValues[0][n];
                        else
                            planePoints[n] = -1;
                    }
                }
                StsSeismicVolume seismicVolume = (StsSeismicVolume)volume0;
                this.dataMin = seismicVolume.dataMin;
                this.dataMax = seismicVolume.dataMax;
                resetDataRange(dataMin, dataMax, false);
                return planePoints;
            }
            else
            {
                return planeValues[0];
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsSeismicVolume.processBytePlaneData() failed.",
                                         e, StsException.WARNING);
            return null;
        }
    }

    public String getCrossplotName()
    {
        if (crossplot == null)
        {
            return "None";
        }
        else
        {
            return ( (StsCrossplot) crossplot).getName();
        }
    }

    public String getIsInclusive()
    {
        String rtn = "Not Applicable";
        if (crossplot == null)
        {
            return rtn;
        }
        else
        {
            if (crossplotInclusive)
            {
                rtn = "TRUE";
            }
            else
            {
                rtn = "FALSE";
            }
            return rtn;
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
    
    static public StsCrossplotVirtualVolumeClass getCrossplotVirtualVolumeClass()
    {
        return (StsCrossplotVirtualVolumeClass) currentModel.getCreateStsClass(StsCrossplotVirtualVolume.class);
    }
}
