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
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.util.*;

public class StsZoneClass extends StsClass implements StsSerializable, StsClassDisplayable, StsClassCursorDisplayable
{
    protected int currentLayer = 1;
    protected boolean displayZones = true;
    protected boolean displaySubZoneColors = false;
    protected boolean displayLayers = false;
    protected boolean isVisibleOnCursor = true;

    transient DefaultBoundedRangeModel rangeModel = new DefaultBoundedRangeModel();

    public static final byte MODEL = 1;
    public static final byte NON_MODEL = 0;

    public StsZoneClass()
    {
    }

    public void initializeDisplayFields()
    {
//        initColors(StsZone.displayFields);
        rangeModel.setRangeProperties(0, 0, 0, getNLayers(), false);

        displayFields = new StsFieldBean[]
        {
            new StsBooleanFieldBean(StsZoneClass.class, "displayZones", "Enable"),
            //new StsBooleanFieldBean(StsZoneClass.class, "displaySubZoneColors", "SubUnitColors"),
            new StsBooleanFieldBean(this, "isVisibleOnCursor", "On 3D Cursors"),
            new StsBooleanFieldBean(StsZoneClass.class, "displayLayers", "Layer"),
            new StsIntFieldBean(StsZoneClass.class, "currentLayer", "Current"),
            new StsSliderFieldBean(StsZoneClass.class, "rangeModel", "rangeChanged")
        };
    }

    public StsFieldBean[] getDisplayFields()
    {
        rangeModel.setRangeProperties(0, 0, 0, getNLayers(), false);
        return super.getDisplayFields();
    }

    public boolean getDisplayZones() { return displayZones; }
    public boolean getDisplaySubZoneColors() { return displaySubZoneColors; }
    public boolean getDisplayLayers() { return displayLayers; }

    public void setDisplayZones(boolean b) { displayZones = b; }
    public void setDisplayLayers(boolean b) { displayLayers = b; }

    public void displayClass(StsGLPanel3d glPanel3d)
    {
        if(currentModel == null) return;
        if(!hasObjects()) return;
        String displayMode = StsBuiltModelClass.staticGetDisplayModeString();
        int nZones = getSize();
        for(int n = 0; n < nZones; n++)
        {
            StsZone zone = (StsZone)getElement(n);
            zone.display(glPanel3d, displayMode);
        }
        displayZoneBlocks(glPanel3d);
    }
    
    private void displayZoneBlocks(StsGLPanel3d glPanel3d)
    {
        if(!displayLayers) return;
        int nZones = getSize();
        if(nZones == 0) return;

        // display stratigraphic slice corresponding to current layer
        // determine which zone / layer we are displaying
        StsZone currentZone = null;
        int currentSublayer = 0;
        int layersAbove = 0;
        for (int i=0; i<nZones; i++)
        {
            StsZone zone = (StsZone)getElement(i);
            int subLayer = currentLayer - layersAbove;
            if( subLayer >= 0 && subLayer <= zone.getNLayers() )
            {
                currentZone = zone;
                currentSublayer = subLayer;
                break;
            }
            layersAbove = layersAbove + zone.getNLayers();
        }

        if( currentZone != null )
            currentZone.displayLayer(glPanel3d, currentSublayer, true);
    }

    /* Draw any map edges on all sections */
    public void drawOnCursor3d(StsGLPanel3d glPanel3d, int dir, float dirCoordinate, StsPoint[] planePoints, boolean isDragging)
    {
        if(!isVisibleOnCursor) return;
        if(currentModel.getCurrentObject(StsBuiltModel.class) == null) return;
        int nZones = getSize();
        String displayMode = StsBuiltModelClass.staticGetDisplayModeString();
        for (int i = 0; i < nZones; i++)
        {
            StsZone zone = (StsZone)getElement(i);
            zone.drawOnCursor3d(glPanel3d, dir, dirCoordinate, planePoints[0], planePoints[3], isDragging, displayMode);
        }
    }

    /* Draw any map edges on all 2d sections */
    public void drawOnCursor2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate,
                               boolean axesFlipped, boolean xAxisReversed, boolean axisReversed)
    {
    }

    public void drawOn3dCurtain(StsGLPanel3d glPanel3d, StsSeismicCurtain seismicCurtain)
    {
    }

    public StsZone getZone(StsModelSurface top, StsModelSurface base)
    {
        int nZones = getSize();
        for (int i = 0; i < nZones; i++)
        {
            StsZone zone = (StsZone)getElement(i);
            if (zone != null && zone.getTopSurface() == top && zone.getBotSurface() == base)
                return zone;
        }

        return null; // no match
    }

    public StsZone getZoneFromWellZone(String wellZoneName) throws StsException
    {
        if (wellZoneName == null) return null;

        int nZones = getSize();
        for (int i = 0; i < nZones; i++)
        {
            StsZone zone = (StsZone)getElement(i);
            StsWellZoneSet wellZoneSet = zone.getWellZoneSet();
            if (wellZoneSet == null) continue;
            if (wellZoneName.equals(wellZoneSet.getName()))
                return zone;
        }
        return null; // no match
    }

    public BoundedRangeModel getLayerRange()
    {
        int nLayers = getNLayers();
        if (nLayers == 0) return null;
        DefaultBoundedRangeModel rangeModel = new DefaultBoundedRangeModel();
        try
        {
            rangeModel.setRangeProperties(currentLayer, 0, 0, nLayers, false);
        }
        catch (Exception e) { e.printStackTrace(); }
        return rangeModel;
    }

    public void rangeChanged()
    {
        if(rangeModel.getValueIsAdjusting()) return;
        currentLayer = rangeModel.getValue();
    }

    public void setLayerRange(BoundedRangeModel range) { currentLayer = range.getValue(); }

    public int getCurrentLayer() { return currentLayer; }
    public void setCurrentLayer(int layer) { currentLayer = layer; }

    public int getNLayers()
    {
        int nZones = getSize();
        int nLayers = 0;
        for (int n = 0; n < nZones; n++)
        {
            StsZone zone = (StsZone)getElement(n);
            nLayers += zone.getNLayers();
        }
        return nLayers;
    }

    public void setDisplaySubZoneColors(boolean b)
    {
        if (displaySubZoneColors == b) return;
        displaySubZoneColors = b;
        forEach("displaySubZoneColorsChanged", null);
    }

    public boolean checkBuildWellZones()
    {
        int nZones = getSize();
        if (nZones == 0)
        {
            StsMessageFiles.logMessage("No zones found.");
            return false;
        }

        Iterator iterator = currentModel.getObjectIterator(StsWell.class);
        int nWells = 0;
        while (iterator.hasNext())
        {
            nWells++;
            StsWell well = (StsWell)iterator.next();
            well.checkConstructWellZones();
        }
        if (nWells > 0) return true;

        StsMessageFiles.logMessage("No wells found.");
        return false;
    }
/*
    public StsPropertyVolumeOld getProperty()
    {
        int nZones = getSize();
        if(nZones == 0) return null;
        return ((StsZone)getElement(0)).getPropertyVolume();
    }
*/

    public StsZone getExistingZone(StsModelSurface top, StsModelSurface base)
    {
        try
        {
            int nZones = getSize();
            for (int i=0; i<nZones; i++)
            {
                StsZone zone = (StsZone)getElement(i);
                if (zone==null) continue;  // shouldn't happen
                if (top==zone.getTopSurface() &&  base==zone.getBotSurface()) return zone;  // already used as a base of zone
            }
            return null;
        }
        catch(Exception e)
        {
            StsException.outputException(e, StsException.WARNING);
            return null;
        }
    }

    public boolean getIsVisibleOnCursor()
	{
		return isVisibleOnCursor;
	}

	public void setIsVisibleOnCursor(boolean isVisibleOnCursor)
	{
		if(this.isVisibleOnCursor == isVisibleOnCursor)return;
		this.isVisibleOnCursor = isVisibleOnCursor;
		currentModel.win3dDisplayAll();
	}
}