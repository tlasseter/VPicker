
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
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsDualSurfaceSubVolume extends StsSubVolume implements StsTreeObjectI
{
    protected StsModelSurface topSurface, botSurface;
    protected float topOffset = 0.0f;
    protected float botOffset = 0.0f;

    static StsObjectPanel objectPanel = null;

    static public StsFieldBean[] displayFields = null;
    static public StsFieldBean[] propertyFields = null;

    static public final float nullValue = StsParameters.nullValue;

    /**
     * Default constructor
     */
    public StsDualSurfaceSubVolume()
    {
        super();
		initializeVisibleFlags();
    }

    /**
     * constructor
     */
    static public StsDualSurfaceSubVolume constructor(String name, StsModelSurface[] surfaces, float topOffset, float botOffset,
                                                      byte offsetDomain, StsSeismicVelocityModel velocityModel)
    {
        try
        {
            StsDualSurfaceSubVolume subVolume = new StsDualSurfaceSubVolume();
            if (!subVolume.initialize(name, surfaces, topOffset, botOffset, offsetDomain, velocityModel)) return null;
            return subVolume;
        }
        catch(Exception e)
        {
            StsException.outputException("StsDualSurfaceSubVolume.constructor(win3d) failed.", e, StsException.WARNING);
            return null;
        }

    }

    static public StsDualSurfaceSubVolume constructor(String name, StsModelSurface surface, float topOffset, float botOffset,
                                                      byte offsetDomain, StsSeismicVelocityModel velocityModel)
    {
        try
        {
            StsDualSurfaceSubVolume subVolume = new StsDualSurfaceSubVolume();
            if (!subVolume.initialize(name, surface, topOffset, botOffset, offsetDomain, velocityModel)) return null;
            return subVolume;
        }
        catch(Exception e)
        {
            StsException.outputException("StsDualSurfaceSubVolume.constructor(win3d) failed.", e, StsException.WARNING);
            return null;
        }
    }

    private boolean initialize(String name, StsModelSurface[] surfaces, float topOffset, float botOffset, byte offsetDomain, StsSeismicVelocityModel velocityModel)
    {
        try
        {
            if (surfaces == null || surfaces.length < 2)
            {
                new StsMessage(currentModel.win3d, StsMessage.ERROR, "Minimum of 2 surfaces required for this subvolume type.");
                return false;
            }
            setName(name);
            this.type = DUAL_SURFACE;
            topSurface = surfaces[0];
            botSurface = surfaces[1];
        /*
            if(topOffset != 0.0f)
                topSurface = StsModelSurface.constructModelSurfaceFromSurfaceAndOffset(surfaces[0], topOffset, offsetDomain, velocityModel);
            else
                topSurface = surfaces[0];
            if(botOffset != 0.0f)
                botSurface = StsModelSurface.constructModelSurfaceFromSurfaceAndOffset(surfaces[1], botOffset, offsetDomain, velocityModel);
            else
                botSurface = surfaces[1];
        */

            this.topOffset = topOffset;
            this.botOffset = botOffset;
            this.zDomainOffset = offsetDomain;
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsDualSurfaceSubVolume.classInitialize() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    private boolean initialize(String name, StsModelSurface surface, float topOffset, float botOffset, byte offsetDomain, StsSeismicVelocityModel velocityModel)
    {
        try
        {
            if (surface == null)
            {
                new StsMessage(currentModel.win3d, StsMessage.ERROR, "Surface required for this subvolume type.");
                return false;
            }
            if (topOffset == 0.0f && botOffset == 0.0f)
            {
                new StsMessage(currentModel.win3d, StsMessage.ERROR, "Both offsets cannot be zero as this creates a zero thickness subVolume.");
                return false;
            }
            setName(name);
            this.type = DUAL_SURFACE;
            topSurface = surface;
            botSurface = surface;
            this.topOffset = topOffset;
            this.botOffset = botOffset;
            this.zDomainOffset = offsetDomain;
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsDualSurfaceSubVolume.classInitialize() failed.",
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
        if(displayFields != null) return displayFields;
        displayFields = new StsFieldBean[]
        {
             new StsStringFieldBean(StsDualSurfaceSubVolume.class, "name", "Name"),
             // new StsBooleanFieldBean(StsDualSurfaceSubVolume.class, "isVisible", "Visible"),
             new StsBooleanFieldBean(StsDualSurfaceSubVolume.class, "isApplied", "Applied")
//		new StsBooleanFieldBean(StsBoxSetSubVolume.class, "isInclusive", "Inclusive"),
//             new StsColorListFieldBean(StsBoxSetSubVolume.class, "color", "Color", currentModel.getSpectrum("Basic").getColors())
         };
        return displayFields;
    }

    public StsFieldBean[] getPropertyFields()
    {
        if(propertyFields != null) return propertyFields;
        propertyFields = new StsFieldBean[]
        {
            new StsStringFieldBean(StsDualSurfaceSubVolume.class, "topSurfaceName", false, "Top Surface Name"),
            new StsStringFieldBean(StsDualSurfaceSubVolume.class, "botSurfaceName", false, "Bottom Surface Name"),
            new StsFloatFieldBean(StsDualSurfaceSubVolume.class, "topOffset", false, "Top Offset"),
            new StsFloatFieldBean(StsDualSurfaceSubVolume.class, "botOffset", false, "Bottom Offset"),
            zDomainStringBean
        };
        return propertyFields;
    }

    public Object[] getChildren() { return new Object[0]; }

    public StsObjectPanel getObjectPanel()
    {
        if(objectPanel == null) objectPanel = StsObjectPanel.constructor(this, true);
        return objectPanel;
    }

    public float getBotOffset() { return botOffset; }
    public void setBotOffset(float botOffset) { this.botOffset = botOffset; }
//    public void setBotSurfaceName(String name) { botSurfaceName = name; }
    public String getBotSurfaceName() { return botSurface.getName(); }
    public float getTopOffset() { return topOffset; }
    public String getTopSurfaceName() { return topSurface.getName(); }
    public void setTopOffset(float topOffset) { this.topOffset = topOffset; }
//    public void setTopSurfaceName(String topSurfaceName) { this.topSurfaceName = topSurfaceName; }


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

        // if no velocityModel is available, all zDomains (data, surfaces, offset) must be the same
        boolean velocityConvert = zDomainData != zDomainOffset || zDomainData != topSurface.getZDomainOriginal() || zDomainData != botSurface.getZDomainOriginal();
        StsSeismicVelocityModel velocityModel = null;
        if(velocityConvert)
        {
            velocityModel = currentModel.getProject().getSeismicVelocityModel();
            if(velocityModel == null) return;
        }
        StsProject project = currentModel.getProject();
        boolean surfacesCongruent = areSurfacesCongruent(rotatedBoundingBox);
        try
        {
            if(surfacesCongruent)
            {
                switch(dir)
                {
                    case StsCursor3d.XDIR:
                        col = topSurface.getNearestColCoor(dirCoordinate);
                        if(col == -1) return;
                        int nRowStart = 0;
                        for(row = 0; row < nRows; row++, nRowStart += nSlices)
                        {
                            zTop = getSurfaceZorT(topSurface, row, col, topOffset, zDomainData, velocityModel);
                            if(zTop == nullValue) continue;
                            zBot = getSurfaceZorT(botSurface, row, col, botOffset, zDomainData, velocityModel);
                            if(zBot == nullValue) continue;
                         /*
                            if(zBot < zTop)
                            {
                                float zTemp = zBot;
                                zBot = zTop;
                                zTop = zTemp;
                            }
                        */
                            nTopSlice = rotatedBoundingBox.getCeilingBoundedSliceCoor(zTop);
                            if(nTopSlice == -1) continue;
                            nBotSlice = rotatedBoundingBox.getFloorBoundedSliceCoor(zBot);
                            if(nBotSlice == -1) continue;
                            nTop = nRowStart + nTopSlice;
                            nBot = nRowStart + nBotSlice;
                            for(n = nTop; n <= nBot; n++)
                                subVolumePlane[n] = (byte)1;
                        }
                        break;
                    case StsCursor3d.YDIR:
                        row = topSurface.getNearestRowCoor(dirCoordinate);
                        if(row == -1) return;
                        int nColStart = 0;
                        for(col = 0; col < nCols; col++, nColStart += nSlices)
                        {
                            zTop = getSurfaceZorT(topSurface, row, col, topOffset, zDomainData, velocityModel);
                            if(zTop == nullValue) continue;
                            zBot = getSurfaceZorT(botSurface, row, col, botOffset, zDomainData, velocityModel);
                            if(zBot == nullValue) continue;
                        /*
                            if(zBot < zTop)
                            {
                                float zTemp = zBot;
                                zBot = zTop;
                                zTop = zTemp;
                            }
                        */
                            nTopSlice = rotatedBoundingBox.getCeilingBoundedSliceCoor(zTop);
                            if(nTopSlice == -1) continue;
                            nBotSlice = rotatedBoundingBox.getFloorBoundedSliceCoor(zBot);
                            if(nBotSlice == -1) continue;
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
                                zTop = getSurfaceZorT(topSurface, row, col, topOffset, zDomainData, velocityModel);
                                if(zTop == nullValue) continue;
                                zBot = getSurfaceZorT(botSurface, row, col, botOffset, zDomainData, velocityModel);
                                if(zBot == nullValue) continue;
                           /*
                                if(zBot < zTop)
                                {
                                    float zTemp = zBot;
                                    zBot = zTop;
                                    zTop = zTemp;
                                }
                            */
                                dirCoordinate = checkConvertZCoordinate(row, col, dirCoordinate, zDomainData, velocityModel);
                                if(dirCoordinate >= zTop && dirCoordinate <= zBot)
                                    subVolumePlane[i] = (byte)1;
                            }
                }
            }
            else
            {
                switch(dir)
                {
                    case StsCursor3d.XDIR:
                        colX = dirCoordinate;
                        rowY = rotatedBoundingBox.getYMin();
                        rowYInc = rotatedBoundingBox.getYInc();
                        int nRowStart = 0;
                        for(row = 0; row < nRows; row++, rowY += rowYInc, nRowStart += nSlices)
                        {
                            zTop = getSurfaceZorT(topSurface, colX, rowY, topOffset, zDomainData, velocityModel);
                            if(zTop == nullValue) continue;
                            zBot = getSurfaceZorT(botSurface, colX, rowY, botOffset, zDomainData, velocityModel);
                            if(zBot == nullValue) continue;
                        /*
                            if(zBot < zTop)
                            {
                                float zTemp = zBot;
                                zBot = zTop;
                                zTop = zTemp;
                            }
                        */
                            nTopSlice = rotatedBoundingBox.getCeilingBoundedSliceCoor(zTop);
                            if(nTopSlice == -1) continue;
                            nBotSlice = rotatedBoundingBox.getFloorBoundedSliceCoor(zBot);
                            if(nBotSlice == -1) continue;
                            nTop = nRowStart + nTopSlice;
                            nBot = nRowStart + nBotSlice;
                            for(n = nTop; n <= nBot; n++)
                                subVolumePlane[n] = (byte)1;
                        }
                        break;
                    case StsCursor3d.YDIR:
                        colX = rotatedBoundingBox.getXMin();
                        rowY = dirCoordinate;
                        colXInc = rotatedBoundingBox.getXInc();
                        int nColStart = 0;
                        for(col = 0; col < nCols; col++, colX += colXInc, nColStart += nSlices)
                        {
                            zTop = getSurfaceZorT(topSurface, colX, rowY, topOffset, zDomainData, velocityModel);
                            if(zTop == nullValue) continue;
                            zBot = getSurfaceZorT(botSurface, colX, rowY, botOffset, zDomainData, velocityModel);
                        /*
                            if(zBot < zTop)
                            {
                                float zTemp = zBot;
                                zBot = zTop;
                                zTop = zTemp;
                            }
                        */
                            nTopSlice = rotatedBoundingBox.getCeilingBoundedSliceCoor(zTop);
                            if(nTopSlice == -1) continue;
                            nBotSlice = rotatedBoundingBox.getFloorBoundedSliceCoor(zBot);
                            if(nBotSlice == -1) continue;
                            nTop = nColStart + nTopSlice;
                            nBot = nColStart + nBotSlice;
                            for(n = nTop; n <= nBot; n++)
                                subVolumePlane[n] = (byte)1;
                        }
                        break;
                    case StsCursor3d.ZDIR:
                        rowY = rotatedBoundingBox.getYMin();
                        float colXStart = rotatedBoundingBox.getXMin();
                        rowYInc = rotatedBoundingBox.getYInc();
                        colXInc = rotatedBoundingBox.getXInc();
                        int i = 0;
                        for(row = 0; row < nRows; row++)
                        {
                            colX = colXStart;
                            for(col = 0; col < nCols; col++, i++, colX += colXInc)
                            {
                                zTop = getSurfaceZorT(topSurface, colX, rowY, topOffset, zDomainData, velocityModel);
                                if(zTop == nullValue) continue;
                                zBot = getSurfaceZorT(botSurface, colX, rowY, botOffset, zDomainData, velocityModel);
                                if(zBot == nullValue) continue;
                            /*
                                if(zBot < zTop)
                                {
                                    float zTemp = zBot;
                                    zBot = zTop;
                                    zTop = zTemp;
                                }
                            */
                                dirCoordinate = checkConvertZCoordinate(colX, rowY, dirCoordinate, zDomainData, velocityModel);
                                if(dirCoordinate >= zTop && dirCoordinate <= zBot)
                                    subVolumePlane[i] = (byte)1;
                            }
                            rowY += rowYInc;
                        }
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "addUnion", e);
        }
    }

    private float checkConvertZCoordinate(int row, int col, float dirCoordinate, byte zDomainData, StsSeismicVelocityModel velocityModel)
    {
        byte zDomainProject = currentModel.getProject().getZDomain();
        if(zDomainProject == zDomainData) return dirCoordinate;
        if(zDomainProject == TD_TIME)
            return (float)velocityModel.getZ(row, col, dirCoordinate);
        else
            return (float)velocityModel.getT(row, col, dirCoordinate, 0.0f);
    }

    private float checkConvertZCoordinate(float x, float y, float dirCoordinate, byte zDomainData, StsSeismicVelocityModel velocityModel)
    {
        byte zDomainProject = currentModel.getProject().getZDomain();
        if(zDomainProject == zDomainData) return dirCoordinate;
        if(zDomainProject == TD_TIME)
            return (float)velocityModel.getZ(x, y, dirCoordinate);
        else
            return (float)velocityModel.getT(x, y, dirCoordinate, 0.0f);
    }

    public boolean delete()
    {
        boolean success = super.delete();
        currentModel.viewObjectChangedAndRepaint(this, this);
        return success;
    }
    
    private boolean areSurfacesCongruent(StsRotatedGridBoundingBox boundingBox)
    {
        if(!boundingBox.sameAs(topSurface, false)) return false;
        if(!boundingBox.sameAs(botSurface, false)) return false;
        return true;
    }

    public StsGridBoundingBox getGridBoundingBox()
    {
        StsGridBoundingBox boundingBox = new StsGridBoundingBox(false);
        StsRotatedGridBoundingBox cursor3dBoundingBox = currentModel.getProject().getRotatedBoundingBox();
        boundingBox.initialize(cursor3dBoundingBox);
        float zMin = topSurface.getZMin();
        boundingBox.zMin = zMin;
        boundingBox.sliceMin = Math.max(0, StsMath.floor(cursor3dBoundingBox.getSliceCoor(zMin)));
        float zMax = botSurface.getZMax();
        boundingBox.zMax = zMax;
        boundingBox.sliceMax = Math.min(cursor3dBoundingBox.getSliceMax(), StsMath.ceiling(cursor3dBoundingBox.getSliceCoor(zMax)));
        return boundingBox;
    }
}
