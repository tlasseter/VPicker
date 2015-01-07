package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jun 12, 2009
 * Time: 12:00:01 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsPatchVolumeClass extends StsClass implements StsClassCursorDisplayable, StsClassDisplayable
{
    protected boolean displayCurvature = false;
    protected boolean displayVoxels = false;
    protected boolean displaySurfs = false;
    protected boolean displayChildPatches = true;
    protected int edgeWidth = 2;

    public void initializeDisplayFields()
    {
        displayFields = new StsFieldBean[]
            {
                new StsBooleanFieldBean(this, "displayCurvature", "Display curvature colorscale on edges and patches"),
                new StsBooleanFieldBean(this, "displaySurfs", "Display as Surfaces"),
                new StsIntFieldBean(this, "edgeWidth", 1, 5, "Line width"),
                new StsBooleanFieldBean(this, "displayVoxels", "Display as Voxel cloud"),
                new StsBooleanFieldBean(this, "displayChildPatches", "Display all connected patches")
            };
    }

    public void initializeDefaultFields()
    {

    }

    public void drawOnCursor3d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, StsPoint[] planePoints, boolean isDragging)
    {
        if(!displaySurfs) return;
        if (currentObject == null) return;
        StsPatchVolume patchVolume = (StsPatchVolume) currentObject;
        patchVolume.drawOnCursor3d(glPanel3d, dirNo, dirCoordinate);
    }

    public void drawOn3dCurtain(StsGLPanel3d glPanel3d, StsSeismicCurtain seismicCurtain)
    {
    }

    /* Draw any map edges on all 2d sections */
    public void drawOnCursor2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, boolean axesFlipped, boolean xAxisReversed, boolean axisReversed)
    {
        StsPatchVolume patchVolume = (StsPatchVolume) currentObject;
        if (currentObject == null) return;
        patchVolume.drawOnCursor2d(glPanel3d, dirNo, dirCoordinate, axesFlipped, xAxisReversed, axisReversed);
    }

    public void setDisplayVoxels(boolean displayVoxels)
    {
        if (this.displayVoxels == displayVoxels) return;
        this.displayVoxels = displayVoxels;
        int nPatchVolumes = getSize();
        for (int i = 0; i < nPatchVolumes; i++)
        {
            StsPatchVolume volume = (StsPatchVolume) getElement(i);
            volume.setDisplayVoxels(displaySurfs);
        }
        currentModel.win3dDisplayAll();
    }

    public boolean getDisplayVoxels()
    {
        return displayVoxels;
    }

    public void setDisplaySurfs(boolean displaySurfs)
    {
        if (this.displaySurfs == displaySurfs) return;
        this.displaySurfs = displaySurfs;
        if(currentObject == null) return;
        StsSeismicVolume seismicVolume = ((StsPatchVolume)currentObject).seismicVolume;
        if(displaySurfs)
            currentModel.toggleOffCursor3dObject(seismicVolume, StsCursor3d.ZDIR);
        else
            currentModel.toggleOnCursor3dObject(seismicVolume, StsCursor3d.ZDIR);
        int nPatchVolumes = getSize();
        for (int i = 0; i < nPatchVolumes; i++)
        {
            StsPatchVolume volume = (StsPatchVolume) getElement(i);
            volume.setDisplaySurfs(displaySurfs);
        }
        currentModel.win3dDisplayAll();
    }

    public boolean getDisplaySurfs()
    {
        return displaySurfs;
    }

    public void setDisplayCurvature(boolean display)
    {
        displayCurvature = display;
    }
    
    public boolean getDisplayCurvature() { return displayCurvature; }

    public void displayClass(StsGLPanel3d glPanel3d)
    {
        if(currentObject != null) currentObject.display(glPanel3d);
    }

    public void setEdgeWidth(int width) { edgeWidth = width; }
    public int getEdgeWidth() { return edgeWidth; }

    public void setDisplayChildPatches(boolean display) { displayChildPatches = display; }
    public boolean getDisplayChildPatches() { return displayChildPatches; }
}
