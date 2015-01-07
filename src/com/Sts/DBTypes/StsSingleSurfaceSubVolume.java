
/**
 * <p>Title: S2S Development</p>
 * <p>Description: SubVolume Class instantiated by the SubVolume toolbar.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author Stuat A. Jackson
 * @version 1.0
 */
package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.StsProject;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsSingleSurfaceSubVolume extends StsSubVolume implements StsTreeObjectI
{
    public StsSurface surface;
    public float topOffset = 0;
    public float botOffset = 0;

    static protected StsObjectPanel objectPanel = null;

    static public StsFieldBean[] displayFields = null;
    static public StsFieldBean[] propertyFields = null;

    static public final float nullValue = StsParameters.nullValue;

    /**
     * Default constructor
     */
    public StsSingleSurfaceSubVolume()
    {
        super();
		initializeVisibleFlags();
   }

    /**
     * constructor
     */
    static public StsSingleSurfaceSubVolume constructor(String name, StsSurface surface, float topOffset, float botOffset)
    {
        try
        {
            StsSingleSurfaceSubVolume subVolume = new StsSingleSurfaceSubVolume();
            if (!subVolume.initialize(name, surface, topOffset, botOffset)) return null;
            return subVolume;
        }
        catch(Exception e)
        {
            StsException.outputException("StsSingleSurfaceSubVolume.constructor(win3d) failed.", e, StsException.WARNING);
            return null;
        }
    }

    public boolean initialize(String name, StsSurface surface, float topOffset, float botOffset)
    {
        try
        {
            if (surface == null)
            {
                new StsMessage(currentModel.win3d, StsMessage.ERROR, "Surface required for this subvolume type.");
                return false;
            }
            setName(name);
            this.type = SINGLE_SURFACE;
            this.surface = surface;
            this.topOffset = topOffset;
            this.botOffset = botOffset;
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsSingleSurfaceSubVolume.classInitialize() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    public boolean anyDependencies()
    {
        return false;
    }
    public StsFieldBean[] getDisplayFields()
    {
        if(displayFields == null)
        {
            displayFields = new StsFieldBean[]
            {
                new StsStringFieldBean(StsSingleSurfaceSubVolume.class, "name", "Name"),
                // new StsBooleanFieldBean(StsSingleSurfaceSubVolume.class, "isVisible", "Visible"),
                new StsBooleanFieldBean(StsSingleSurfaceSubVolume.class, "isApplied", "Applied"),
        //		new StsBooleanFieldBean(StsBoxSetSubVolume.class, "isInclusive", "Inclusive"),
                new StsColorComboBoxFieldBean(StsSingleSurfaceSubVolume.class, "stsColor", "Color", currentModel.getSpectrum("Basic").getStsColors())
            };
        }
        return displayFields;
    }

    public StsFieldBean[] getPropertyFields()
    {
        if(propertyFields != null) return propertyFields;
        propertyFields = new StsFieldBean[]
        {
            new StsStringFieldBean(StsSingleSurfaceSubVolume.class, "surfaceName", false, "Surface Name"),
            new StsFloatFieldBean(StsSingleSurfaceSubVolume.class, "topOffset", false, "Top Offset"),
            new StsFloatFieldBean(StsSingleSurfaceSubVolume.class, "botOffset", false, "Bottom Offset")
        };
        return propertyFields;
    }

    public Object[] getChildren() { return new Object[0]; }

    public StsObjectPanel getObjectPanel()
    {
        if(objectPanel == null) objectPanel = StsObjectPanel.constructor(this, true);
        return objectPanel;
    }

    public String getSurfaceName() { return surface.getName(); }
    public float getBotOffset() { return botOffset; }
    public void setBotOffset(float botOffset) { this.botOffset = botOffset; }
    public float getTopOffset() { return topOffset; }
    public void setTopOffset(float topOffset) { this.topOffset = topOffset; }

	public StsGridBoundingBox getGridBoundingBox()
	 {
		 StsGridBoundingBox boundingBox = new StsGridBoundingBox(false);
		 StsRotatedGridBoundingBox cursor3dBoundingBox = currentModel.getProject().getRotatedBoundingBox();
		 boundingBox.initialize(cursor3dBoundingBox);
		 float zMin = surface.getZMin() + topOffset;
		 boundingBox.zMin = zMin;
		 boundingBox.sliceMin = Math.max(0, StsMath.floor(cursor3dBoundingBox.getSliceCoor(zMin)));
		 float zMax = surface.getZMax() + botOffset;
		 boundingBox.zMax = zMax;
		 boundingBox.sliceMax = Math.min(cursor3dBoundingBox.getSliceMax(), StsMath.ceiling(cursor3dBoundingBox.getSliceCoor(zMax)));
		 return boundingBox;
	 }

    public void addUnion(byte[] subVolumePlane, int dir, float dirCoordinate, StsRotatedGridBoundingBox rotatedBoundingBox, byte zDomainData)
    {
        float zTop, zBot;
        int row, col;
        int nTopSlice, nBotSlice;
        int nTop, nBot;
        int n;
        float rowY, colX, rowYInc, colXInc;

        if(!isApplied) return;

        int nRows = rotatedBoundingBox.getNRows();
        int nCols = rotatedBoundingBox.getNCols();
        int nSlices = rotatedBoundingBox.getNSlices();

        boolean subVolumeIsDepth = surface.getZDomainOriginal() == StsProject.TD_DEPTH;
        // if subVolume domain and currentDomain are the same no velocity conversion necessary
        boolean velocityConvert = (isDepth != subVolumeIsDepth);
        StsSeismicVelocityModel velocityModel = null;
        if(velocityConvert)
        {
            velocityModel = currentModel.getProject().getSeismicVelocityModel();
            if(velocityModel == null) return;
        }
        StsProject project = currentModel.getProject();
        try
        {
            switch(dir)
            {
                case StsCursor3d.XDIR:
                    col = surface.getNearestColCoor(dirCoordinate);
                    if(col == -1) return;
                    int nRowStart = 0;
                    for(row = 0; row < nRows; row++, nRowStart += nSlices)
                    {
                        zTop = getSurfaceZorT(surface, row, col, topOffset, zDomainData, velocityModel);
                        if(zTop == nullValue) continue;
                        zBot = getSurfaceZorT(surface, row, col, botOffset, zDomainData, velocityModel);
                        if(zBot == nullValue) continue;

                        if(zBot < zTop)
                        {
                            float zTemp = zBot;
                            zBot = zTop;
                            zTop = zTemp;
                        }
                        nTopSlice = rotatedBoundingBox.getNearestSliceCoor(zTop);
                        nBotSlice = rotatedBoundingBox.getNearestSliceCoor(zBot);
                        nTop = nRowStart + nTopSlice;
                        nBot = nRowStart + nBotSlice;
                        for(n = nTop; n <= nBot; n++)
                            subVolumePlane[n] = (byte)1;
                    }
                    break;
                case StsCursor3d.YDIR:
                    row = surface.getNearestRowCoor(dirCoordinate);
                    if(row == -1) return;
                    int nColStart = 0;
                    for(col = 0; col < nCols; col++, nColStart += nSlices)
                    {
                        zTop = getSurfaceZorT(surface, row, col, topOffset, zDomainData, velocityModel);
                        if(zTop == nullValue) continue;
                        zBot = getSurfaceZorT(surface, row, col, botOffset, zDomainData, velocityModel);
                        if(zBot == nullValue) continue;

                        if(zBot < zTop)
                        {
                            float zTemp = zBot;
                            zBot = zTop;
                            zTop = zTemp;
                        }
                        nTopSlice = rotatedBoundingBox.getNearestSliceCoor(zTop);
                        nBotSlice = rotatedBoundingBox.getNearestSliceCoor(zBot);
                        nTop = nColStart + nTopSlice;
                        nBot = nColStart + nBotSlice;
                        for(n = nTop; n <= nBot; n++)
                            subVolumePlane[n] = (byte)1;
                    }
                    break;
                case StsCursor3d.ZDIR:
                    int i = 0;
                    for(row = 0; row < nRows; row++)
                        for(col = 0; col < nCols; col++, i++)
                        {
                             zTop = getSurfaceZorT(surface, row, col, topOffset, zDomainData, velocityModel);
                            if(zTop == nullValue) continue;
                            zBot = getSurfaceZorT(surface, row, col, botOffset, zDomainData, velocityModel);
                            if(zBot == nullValue) continue;
                            if(zBot < zTop)
                            {
                                float zTemp = zBot;
                                zBot = zTop;
                                zTop = zTemp;
                            }
                            if(dirCoordinate >= zTop && dirCoordinate <= zBot)
                                subVolumePlane[i] = (byte)1;
                        }
            }
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "addUnion", e);
        }
    }

    public boolean delete()
    {
        boolean success = super.delete();
        currentModel.viewObjectChangedAndRepaint(this, this);
        return success;
    }
    
    public StsGridBoundingBox getBoundingBox()
    {
        StsGridBoundingBox boundingBox = new StsGridBoundingBox(false);
        StsRotatedGridBoundingBox cursor3dBoundingBox = currentModel.getProject().getRotatedBoundingBox();
        boundingBox.initialize(cursor3dBoundingBox);
        float zMin = surface.getZMin() + topOffset;
        boundingBox.zMin = zMin;
        boundingBox.sliceMin = Math.max(0, StsMath.floor(cursor3dBoundingBox.getSliceCoor(zMin)));
        float zMax = surface.getZMax() + botOffset;
        boundingBox.zMax = zMax;
        boundingBox.sliceMax = Math.min(cursor3dBoundingBox.getSliceMax(), StsMath.ceiling(cursor3dBoundingBox.getSliceCoor(zMax)));
        return boundingBox;
    }
}
