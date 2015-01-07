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
import com.Sts.RFUI.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.util.*;

public class StsSurfaceClass extends StsClass implements StsSerializable, StsClassTimeDisplayable, StsClassCursorDisplayable, StsClassTextureDisplayable
{
    /** display all surfaces */
    protected boolean displaySurfaces = true;
    /** display surface with lighting */
    protected boolean lighting = true;
    /** null points displayed in null color */
    protected boolean nullsFilled = false;
    /** null points displayed in null color */
    protected boolean enableTime = true;
    /** current null fill color */
    protected String nullColorName = "gray";

    private boolean defaultDisplayGrid = false;
    private boolean defaultDisplayFill = true;
    private StsColor defaultSurfaceColor;

    static final String[] nullColorNames = new String[]{"surface", "gray", "black", "white"};
    static StsColor[] nullColors = new StsColor[]{null, StsColor.GRAY, StsColor.BLACK, StsColor.WHITE};

    public StsSurfaceClass()
    {
        userName = "Surfaces";
    }

    public StsSurface[] getSurfaces()
    {
        Object surfaceList;

        StsSurface[] surfaces = new StsSurface[0];
        surfaceList = currentModel.getCastObjectList(StsSurface.class);
        surfaces = (StsSurface[]) StsMath.arrayAddArray(surfaces, surfaceList);
        return surfaces;
    }

    public void initializeDisplayFields()
    {
        displayFields = new StsFieldBean[]
                {
                        new StsBooleanFieldBean(this, "displaySurfaces", "Display"),
                        new StsBooleanFieldBean(this, "enableTime", "Enable Time"),
                        new StsBooleanFieldBean(this, "lighting", "Lighting"),
                        new StsBooleanFieldBean(this, "nullsFilled", "Fill Nulls"),
                        new StsComboBoxFieldBean(this, "nullColorName", "Null Color:", nullColorNames)
                };
    }

    public void initializeDefaultFields()
    {
        defaultFields = new StsFieldBean[]
                {
                        new StsBooleanFieldBean(this, "defaultDisplayGrid", "Show Grid"),
                        new StsBooleanFieldBean(this, "defaultDisplayFill", "Fill Surface")
                };
    }

    public void displayTimeClass(StsGLPanel3d glPanel3d, long time)
    {
        Iterator iter = getObjectIterator();
        while (iter.hasNext())
        {
            StsSurface surface = (StsSurface) iter.next();
            if ((enableTime && surface.isAlive(time)) || (!enableTime))
                surface.display(glPanel3d);
        }
    }

    public boolean getDisplaySurfaces() { return displaySurfaces; }

    public void setDisplaySurfaces(boolean b)
    {
        if (this.displaySurfaces == b) return;
        displaySurfaces = b;
 //      setDisplayField("displaySurfaces", displaySurfaces);

        int nSurfaces = getSize();
        for (int n = 0; n < nSurfaces; n++)
        {
            StsSurface surface = (StsSurface) getElement(n);
            surface.setIsVisible(displaySurfaces);
        }
        currentModel.win3dDisplayAll();
    }

    public boolean getDisplaySeismic() { return currentModel.getBooleanProperty("seismicOnSurface"); }

    public void setDisplaySeismic(boolean b) { currentModel.getProperties().set("seismicOnSurface", b); }

    /* Draw any map edges on all sections */
    public void drawOnCursor3d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, StsPoint[] planePoints, boolean isDragging)
    {
        int nSurfaces = getSize();
        for (int i = 0; i < nSurfaces; i++)
        {
            StsSurface s = (StsSurface) getElement(i);
            if (s != null && s.getIsVisible())
                s.drawOnCursor3d(glPanel3d, planePoints);
        }
    }

    public void drawOn3dCurtain(StsGLPanel3d glPanel3d, StsSeismicCurtain seismicCurtain)
    {
        int nSurfaces = getSize();
        for (int i = 0; i < nSurfaces; i++)
        {
            StsSurface s = (StsSurface) getElement(i);
            if (s != null && s.getIsVisible())
                s.drawOn3dCurtain(glPanel3d, seismicCurtain);
        }
    }

    /* Draw any map edges on all 2d sections */
    public void drawOnCursor2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate,
                               boolean axesFlipped, boolean xAxisReversed, boolean axisReversed)
    {
        if (!displaySurfaces) return;

        int nSurfaces = getSize();
        for (int i = 0; i < nSurfaces; i++)
        {
            StsSurface s = (StsSurface) getElement(i);
            if (s != null && s.getIsVisible())
                s.drawOnCursor2d(glPanel3d, dirNo, dirCoordinate);
        }
    }

    public StsSurface selectSurface(boolean visibleFirst, Class c)
    {
        try
        {
            StsSurface[] surfaces = (StsSurface[]) getCastObjectList(c);
            if (surfaces.length == 0) return null;

            if (surfaces.length == 1) return surfaces[0];

            StsSelectStsObjects selector = StsSelectStsObjects.constructor(currentModel, surfaces,
                    name, "Select a surface:", true);
            if (selector == null) return null;
            selector.setVisibleOnly(visibleFirst);
            StsSurface surface = (StsSurface) selector.selectObject();
            if (!visibleFirst) return surface;
            if (surface == null && selector.getNVisibleObjects() == 0)
            {
                selector.setVisibleOnly(false);
                surface = (StsSurface) selector.selectObject();
            }
            return surface;
        }
        catch (Exception e)
        {
            StsException.outputException("StsModel.selectSurface() failed.",
                    e, StsException.WARNING);
            return null;
        }
    }

    public void unloadSurfaces()
    {
        int nSurfaces = getSize();
        for (int n = 0; n < nSurfaces; n++)
        {
            StsSurface surface = (StsSurface) getElement(n);
            surface.unload();
        }
        currentModel.win3dDisplayAll();
    }

    public void setLighting(boolean lighting)
    {
        if (this.lighting == lighting) return;
        this.lighting = lighting;
        list.forEach("setLighting", new Boolean(lighting));
 //      setDisplayField("lighting", lighting);
        currentModel.win3dDisplayAll();
    }

    public boolean getLighting() { return lighting; }

    public void setNullsFilled(boolean nullsFilled)
    {
        if (this.nullsFilled == nullsFilled) return;
        this.nullsFilled = nullsFilled;
        list.forEach("nullsFilledChanged");
 //      setDisplayField("nullsFilled", nullsFilled);
        currentModel.win3dDisplayAll();
    }

    public boolean getNullsFilled() { return nullsFilled; }

    public void setEnableTime(boolean enable)
    {
        if (this.enableTime == enable) return;
        this.enableTime = enable;
 //      setDisplayField("enableTime", enableTime);
        currentModel.win3dDisplayAll();
    }

    public boolean getEnableTime() { return enableTime; }

    public void deleteDisplayLists()
    {
        forEach("deleteDisplayLists");
    }

    public void setNullColorName(String nullColorName)
    {
        if (nullColorName.equals(name)) return;
        this.nullColorName = nullColorName;
        list.forEach("nullsFilledChanged");
//        setDisplayField("nullColorName", nullColorName);
        currentModel.win3dDisplayAll();
    }

    public String getNullColorName() { return nullColorName; }

    public StsColor getNullColor(StsColor surfaceColor)
    {
        for (int n = 1; n < nullColorNames.length; n++)
            if (nullColorName.equals(nullColorNames[n])) return nullColors[n];
        return surfaceColor;
    }

    public StsSurface getTopVisibleSurface()
    {
        StsSurface surface = null;
        try
        {
            // try using visible surfaces only
            int nSurfaces = getSize();
            if (nSurfaces == 0) return null;

            for (int n = 0; n < nSurfaces; n++)
            {
                surface = (StsSurface) getElement(n);
                if (surface.getIsVisible()) return surface;
            }
            if (surface == null)
            {
                surface = (StsSurface) getElement(0); // if no visible ones available, use first (top) surface
                surface.setIsVisible(true);
            }

            if (!surface.checkIsLoaded()) return null;
            return surface;
        }
        catch (Exception e)
        {
            StsException.outputException("StsModelSurfaceClass.getTopVisibleSurface() failed.",
                    e, StsException.WARNING);
            return null;
        }
    }

    /** Given an array of surfaces, compute a common grid */
    public StsGridDefinition checkComputeUnionGrid(StsSurface[] surfaces)
    {
        if (!anglesOk(surfaces)) return null;
        StsGridDefinition grid = checkConstructSameGrid(surfaces);
        if (grid != null) return grid;
        return computeUnionGrid(surfaces);
    }

    private boolean anglesOk(StsSurface[] surfaces)
    {
        // surfaces are not congruent.  Check if rotation angles are different
        int nSurfaces = surfaces.length;

        float angle = surfaces[0].getAngle();
        float minAngle = angle;
        float maxAngle = angle;
        boolean anglesSame = true;
        for (int n = 1; n < nSurfaces; n++)
        {
            float angle0 = surfaces[n].getAngle();
            minAngle = Math.min(angle, minAngle);
            maxAngle = Math.max(angle, maxAngle);
            anglesSame = StsMath.sameAsTol(angle, angle0, 0.05f);
            angle = angle0;
        }

//        boolean anglesSame = StsMath.sameAsTol(minAngle, maxAngle);
        if (!anglesSame)
        {
            angle = (minAngle + maxAngle) / 2;
            new StsMessage(currentModel.win3d, StsMessage.WARNING, "Surfaces have different rotation angles: not supported yet.\n" +
                    "Cannot create horizons.");
        }
        return anglesSame;
    }

    private StsGridDefinition checkConstructSameGrid(StsSurface[] surfaces)
    {
        if (!StsSurface.surfacesGridsAreSame(surfaces)) return null;
        StsSurface surface = surfaces[0];
        return new StsGridDefinition(surface, currentModel);
    }

    public StsGridDefinition computeUnionGrid(StsSurface[] surfaces)
    {
        StsProject project = currentModel.getProject();
        StsGridDefinition gridDefinition = new StsGridDefinition();
        gridDefinition.setOrigin(project.getXOrigin(), project.getYOrigin());
        gridDefinition.setAngle(surfaces[0].getAngle());
        int nSurfaces = surfaces.length;
        for (int n = 0; n < nSurfaces; n++)
            gridDefinition.addVolumeBoundingBox(surfaces[n]);
        return gridDefinition;
    }

    public void setDefaultDisplayGrid(boolean defaultDisplayGrid)
    {
        if (this.defaultDisplayGrid == defaultDisplayGrid) return;
        this.defaultDisplayGrid = defaultDisplayGrid;
 //      setDisplayField("defaultDisplayGrid", defaultDisplayGrid);
        currentModel.win3dDisplayAll();
    }

    public boolean getDefaultDisplayGrid() { return defaultDisplayGrid; }

    public void setDefaultDisplayFill(boolean defaultDisplayFill)
    {
        if (this.defaultDisplayFill == defaultDisplayFill) return;
        this.defaultDisplayFill = defaultDisplayFill;
 //      setDisplayField("defaultDisplayFill", defaultDisplayFill);
        currentModel.win3dDisplayAll();
    }

    public boolean getDefaultDisplayFill() { return defaultDisplayFill; }

    public void cropChanged()
    {
        list.forEach("cropChanged");
    }

    public boolean textureChanged(StsObject object)
    {
        return ((StsSurface)object).textureChanged();
    }
}