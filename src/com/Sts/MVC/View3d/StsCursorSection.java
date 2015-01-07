package com.Sts.MVC.View3d;

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;
import java.io.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Feb 4, 2008
 * Time: 3:45:28 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsCursorSection implements Serializable, StsSurfaceGridable
{
    StsCursor3dTexture[] displayableSections = new StsCursor3dTexture[0];
    public int dir; /** cursor direction */
    public float dirCoordinate = StsParameters.largeFloat; /** cursor position in coordinates */
    boolean isSelected = false;

    transient public StsCursor3d cursor3d;
    transient public StsRotatedGridBoundingBox rotatedBoundingBox;
    transient StsPoint[] planePoints = null;
    transient float[] totalHorizontalRange, totalVerticalRange;
    transient float[] initTotalHorizontalRange, initTotalVerticalRange;
    transient String horizontalAxisLabel, verticalAxisLabel, titleLabel;

    transient byte[] subVolumePlane = null;

    /** static variable NONE = -1 - No cursors on */
     public static final int NONE = StsParameters.NO_MATCH;
     /** static variable XDIR = 0 - X oriented cursor section */
     public static final int XDIR = StsParameters.XDIR;
     /** static variable YDIR = 1 - Y oriented cursor section */
     public static final int YDIR = StsParameters.YDIR;
     /** static variable ZDIR = 2 - Z oriented cursor section */
     public static final int ZDIR = StsParameters.ZDIR;

    static final long serialVersionUID = 1l;
    static final boolean debug = false;

    public StsCursorSection()
    {
    }

    /** CursorSection constructor */
    public StsCursorSection(StsCursor3d cursor3d, int dirNo, float dirCoordinate)
    {
        this.cursor3d = cursor3d;
        this.dir = dirNo;
        this.dirCoordinate = dirCoordinate;
        rotatedBoundingBox = cursor3d.rotatedBoundingBox;
        if (rotatedBoundingBox == null) return;
        rangeChanged();
 /*
        if (rotatedBoundingBox == null) return;
        planePoints = new StsPoint[4];
        rangeChanged();
        initializeDisplayableSections();
        computeSubVolumePlane();
*/
    }

    public void initialize()
    {
//        computeSubVolumePlane();
    }

    public void initializeTransients(StsModel model, StsCursor3d cursor3d, int dir)
    {
        this.cursor3d = cursor3d;
        rotatedBoundingBox = cursor3d.rotatedBoundingBox;
        rangeChanged();
        initializeTextureCursorSections(model, cursor3d, dir);
//        checkAddDisplayableSections();
    }

    public void initializeCursorSection()
    {
//        initializeTextureCursorSections();
//        computeSubVolumePlane();
    }

    private void initializeTextureCursorSections(StsModel model, StsCursor3d cursor3d, int dir)
    {
        for(int n = displayableSections.length - 1; n >= 0 ; n--)
            if(!displayableSections[n].initialize(model, cursor3d, dir))
                displayableSections = (StsCursor3dTexture[])StsMath.arrayDeleteElement(displayableSections, displayableSections[n]);
    }

    protected void rangeChanged()
    {
        if(rotatedBoundingBox == null) return;

        if(planePoints == null) planePoints = new StsPoint[4];
        switch(dir)
        {
            case XDIR:
                dirCoordinate = StsMath.minMax(dirCoordinate, rotatedBoundingBox.xMin, rotatedBoundingBox.xMax);
                planePoints[0] = new StsPoint(dirCoordinate, rotatedBoundingBox.yMin, rotatedBoundingBox.zMin);
                planePoints[1] = new StsPoint(dirCoordinate, rotatedBoundingBox.yMin, rotatedBoundingBox.zMax);
                planePoints[2] = new StsPoint(dirCoordinate, rotatedBoundingBox.yMax, rotatedBoundingBox.zMax);
                planePoints[3] = new StsPoint(dirCoordinate, rotatedBoundingBox.yMax, rotatedBoundingBox.zMin);
                initTotalHorizontalRange = new float[]
                    {rotatedBoundingBox.yMin, rotatedBoundingBox.yMax};
                totalHorizontalRange = new float[]
                    {rotatedBoundingBox.yMin, rotatedBoundingBox.yMax};
                initTotalVerticalRange = new float[]
                    {rotatedBoundingBox.zMax, rotatedBoundingBox.zMin};
                totalVerticalRange = new float[]
                    {rotatedBoundingBox.zMax, rotatedBoundingBox.zMin};
                break;
            case YDIR:
                dirCoordinate = StsMath.minMax(dirCoordinate, rotatedBoundingBox.yMin, rotatedBoundingBox.yMax);
                planePoints[0] = new StsPoint(rotatedBoundingBox.xMin, dirCoordinate, rotatedBoundingBox.zMin);
                planePoints[1] = new StsPoint(rotatedBoundingBox.xMin, dirCoordinate, rotatedBoundingBox.zMax);
                planePoints[2] = new StsPoint(rotatedBoundingBox.xMax, dirCoordinate, rotatedBoundingBox.zMax);
                planePoints[3] = new StsPoint(rotatedBoundingBox.xMax, dirCoordinate, rotatedBoundingBox.zMin);
                initTotalHorizontalRange = new float[]
                    {rotatedBoundingBox.xMin, rotatedBoundingBox.xMax};
                totalHorizontalRange = new float[]
                    {rotatedBoundingBox.xMin, rotatedBoundingBox.xMax};
                initTotalVerticalRange = new float[]
                    {rotatedBoundingBox.zMax, rotatedBoundingBox.zMin};
                totalVerticalRange = new float[]
                    {rotatedBoundingBox.zMax, rotatedBoundingBox.zMin};
                break;
            case ZDIR:
                dirCoordinate = StsMath.minMax(dirCoordinate, rotatedBoundingBox.zMin, rotatedBoundingBox.zMax);
                planePoints[0] = new StsPoint(rotatedBoundingBox.xMin, rotatedBoundingBox.yMin, dirCoordinate);
                planePoints[1] = new StsPoint(rotatedBoundingBox.xMax, rotatedBoundingBox.yMin, dirCoordinate);
                planePoints[2] = new StsPoint(rotatedBoundingBox.xMax, rotatedBoundingBox.yMax, dirCoordinate);
                planePoints[3] = new StsPoint(rotatedBoundingBox.xMin, rotatedBoundingBox.yMax, dirCoordinate);
                initTotalHorizontalRange = new float[]
                    {rotatedBoundingBox.xMin, rotatedBoundingBox.xMax};
                totalHorizontalRange = new float[]
                    {rotatedBoundingBox.xMin, rotatedBoundingBox.xMax};
                initTotalVerticalRange = new float[]
                    {rotatedBoundingBox.yMin, rotatedBoundingBox.yMax};
                totalVerticalRange = new float[]
                    {rotatedBoundingBox.yMin, rotatedBoundingBox.yMax};
                break;
        }
        //    nPlane = rotatedBoundingBox.getCursorPlaneIndex(dirNo, dirCoordinate);
    }

    protected void zRangeChanged()
    {
         if(rotatedBoundingBox == null) return;

        if(planePoints == null) planePoints = new StsPoint[4];
        switch(dir)
        {
            case XDIR:
                planePoints[0] = new StsPoint(dirCoordinate, rotatedBoundingBox.yMin, rotatedBoundingBox.zMin);
                planePoints[1] = new StsPoint(dirCoordinate, rotatedBoundingBox.yMin, rotatedBoundingBox.zMax);
                planePoints[2] = new StsPoint(dirCoordinate, rotatedBoundingBox.yMax, rotatedBoundingBox.zMax);
                planePoints[3] = new StsPoint(dirCoordinate, rotatedBoundingBox.yMax, rotatedBoundingBox.zMin);
                initTotalHorizontalRange = new float[]
                    {rotatedBoundingBox.yMin, rotatedBoundingBox.yMax};
                totalHorizontalRange = new float[]
                    {rotatedBoundingBox.yMin, rotatedBoundingBox.yMax};
                initTotalVerticalRange = new float[]
                    {rotatedBoundingBox.zMax, rotatedBoundingBox.zMin};
                totalVerticalRange = new float[]
                    {rotatedBoundingBox.zMax, rotatedBoundingBox.zMin};
                break;
            case YDIR:
                planePoints[0] = new StsPoint(rotatedBoundingBox.xMin, dirCoordinate, rotatedBoundingBox.zMin);
                planePoints[1] = new StsPoint(rotatedBoundingBox.xMin, dirCoordinate, rotatedBoundingBox.zMax);
                planePoints[2] = new StsPoint(rotatedBoundingBox.xMax, dirCoordinate, rotatedBoundingBox.zMax);
                planePoints[3] = new StsPoint(rotatedBoundingBox.xMax, dirCoordinate, rotatedBoundingBox.zMin);
                initTotalHorizontalRange = new float[]
                    {rotatedBoundingBox.xMin, rotatedBoundingBox.xMax};
                totalHorizontalRange = new float[]
                    {rotatedBoundingBox.xMin, rotatedBoundingBox.xMax};
                initTotalVerticalRange = new float[]
                    {rotatedBoundingBox.zMax, rotatedBoundingBox.zMin};
                totalVerticalRange = new float[]
                    {rotatedBoundingBox.zMax, rotatedBoundingBox.zMin};
                break;
            case ZDIR:
                dirCoordinate = StsMath.minMax(dirCoordinate, rotatedBoundingBox.zMin, rotatedBoundingBox.zMax);
                planePoints[0] = new StsPoint(rotatedBoundingBox.xMin, rotatedBoundingBox.yMin, dirCoordinate);
                planePoints[1] = new StsPoint(rotatedBoundingBox.xMax, rotatedBoundingBox.yMin, dirCoordinate);
                planePoints[2] = new StsPoint(rotatedBoundingBox.xMax, rotatedBoundingBox.yMax, dirCoordinate);
                planePoints[3] = new StsPoint(rotatedBoundingBox.xMin, rotatedBoundingBox.yMax, dirCoordinate);
                initTotalHorizontalRange = new float[]
                    {rotatedBoundingBox.xMin, rotatedBoundingBox.xMax};
                totalHorizontalRange = new float[]
                    {rotatedBoundingBox.xMin, rotatedBoundingBox.xMax};
                initTotalVerticalRange = new float[]
                    {rotatedBoundingBox.yMin, rotatedBoundingBox.yMax};
                totalVerticalRange = new float[]
                    {rotatedBoundingBox.yMin, rotatedBoundingBox.yMax};
        }
    }

    public void setAxisLabel()
    {
        if(cursor3d == null) return;
        
        boolean isGridCoordinates = cursor3d.getIsGridCoordinates();
        switch(dir)
        {
            case XDIR:
                if(isGridCoordinates)
                    horizontalAxisLabel = "Inline";
                else
                    horizontalAxisLabel = "Inline Y";
                verticalAxisLabel = cursor3d.model.getProject().getZDomainString();
                titleLabel = "Crossline ";
                break;

            case YDIR:
                if(isGridCoordinates)
                    horizontalAxisLabel = "Crossline";
                else
                    horizontalAxisLabel = "Crossline Y";
                verticalAxisLabel = cursor3d.model.getProject().getZDomainString();
                titleLabel = "Inline ";
                break;
            case ZDIR:
                if(isGridCoordinates)
                    horizontalAxisLabel = "Crossline";
                else
                    horizontalAxisLabel = "Crossline Y";
                if(isGridCoordinates)
                    verticalAxisLabel = "Inline";
                else
                    verticalAxisLabel = "Inline Y";
                  titleLabel = "Slice ";
        }
    }

    public void resetZRange(float zMin, float zMax, float z)
    {
        switch(dir)
        {
            case XDIR:
                planePoints[0].setZ(zMin);
                planePoints[1].setZ(zMax);
                planePoints[2].setZ(zMax);
                planePoints[3].setZ(zMin);
                break;
            case YDIR:
                planePoints[0].setZ(zMin);
                planePoints[1].setZ(zMax);
                planePoints[2].setZ(zMax);
                planePoints[3].setZ(zMin);
                break;
            case ZDIR:
                planePoints[0].setZ(z);
                planePoints[1].setZ(z);
                planePoints[2].setZ(z);
                planePoints[3].setZ(z);
                dirCoordinate = z;
                break;
        }
    }
/*
    private void checkAddDisplayableSections()
    {
        StsArrayList displayableClasses = cursor3d.model.getCursor3dTextureDisplayableClasses();
        for(int n = 0; n < displayableClasses.size(); n++)
        {
            Object classObject = displayableClasses.get(n);
            StsClass displayableClass = (StsClass)classObject;
            if(!displayableClass.hasObjects())continue;
            StsClassCursor3dTextureDisplayable textureDisplayableClass = (StsClassCursor3dTextureDisplayable)classObject;
            checkAddDisplayableSection(textureDisplayableClass);
        }
    }
*/
    /**
     * Set the cursor to the provided coordinate
     * @param dirCoor cursor coordinate
     */
/*
    void viewChanged(StsWin3dBase window, float dirCoor)
    {
        if(dirCoordinate == dirCoor)return;
        if(debug) StsException.systemDebug(this, "setCursor", " dirNo: " + dir + " coordinate: " + dirCoor);
        dirCoordinate = dirCoor;

        if(planePoints != null)
        {
            for(int n = 0; n < 4; n++)
                planePoints[n].v[dir] = dirCoordinate;
        }

        StsView[] displayedViews = window.getDisplayedViews();
        for(StsView view : displayedViews)
        {
            if(view.glPanel3d.displayableSections == null || view.glPanel3d.displayableSections[dir] == null) continue;
            StsCursorSectionDisplayable[] displayableSections = view.glPanel3d.displayableSections[dir];
            for(StsCursorSectionDisplayable sectionDisplayable : displayableSections)
                sectionDisplayable.setDirCoordinate(dirCoordinate);
        }
        subVolumeChanged();
    }

        void subVolumeChanged()
     {
         if(computeSubVolumePlane()) cursor3d.clearTextureDisplays(dir);
//			computeSubVolumePlane();
//			clearTextureDisplays();
     }
*/
    public void setDirCoordinate(float dirCoor)
    {
        dirCoordinate = dirCoor;
        if(planePoints == null) return;
        for(int n = 0; n < 4; n++)
            planePoints[n].v[dir] = dirCoordinate;
        for(StsCursor3dTexture displaySection : displayableSections)
            displaySection.setDirCoordinate(dirCoor);
        // subVolumePlane = null;
    }
/*
    byte[] getSubVolumePlane()
    {
        if(subVolumePlane == null) computeSubVolumePlane();
        return subVolumePlane;
    }
*/
    void clearSubVolumePlane()
    {
        subVolumePlane = null;
        for(StsCursor3dTexture displaySection : displayableSections)
            displaySection.subVolumeChanged();
    }
 /*
    boolean computeSubVolumePlane()
    {
    	if(cursor3d == null) return false;
        StsSubVolumeClass subVolumeClass = (StsSubVolumeClass)cursor3d.model.getStsClass(StsSubVolume.class);
        if(subVolumeClass == null)return false;
        subVolumePlane = subVolumeClass.getSubVolumePlane(dir, dirCoordinate, rotatedBoundingBox);
        return subVolumePlane != null;
    }
*/
    public boolean hasSubVolumes()
    {
        if(cursor3d == null) return false;
        StsSubVolumeClass subVolumeClass = (StsSubVolumeClass)cursor3d.model.getStsClass(StsSubVolume.class);
        if(subVolumeClass == null)return false;
        return true;
    }

    /**
     * Get the index of the current cursor plane
     * @return plane index
     */
    /*
            int getPlaneIndex()
            {
                if(nPlane != -1) return nPlane;
//            if(seismicCursorSection == null) return -1;
                nPlane = rotatedBoundingBox.getCursorPlaneIndex(dirNo, dirCoordinate);
                return nPlane;
            }
     */
    /**
     * Get the cropped axis ranges
     * @return [horz][min],[horz][max],[vert][min],[vert][max]
     */
    float[][] getTotalAxisRanges()
    {
        rangeChanged(); // jbw make sure
        return new float[][]
            {totalHorizontalRange, totalVerticalRange};
    }

    /**
     * Reset the axis ranges to the extents of the data and get axis ranges
     * @return [horz][min],[horz][max],[vert][min],[vert][max]
     */

    float[][] resetAxisRanges()
    {
        totalHorizontalRange[0] = initTotalHorizontalRange[0];
        totalHorizontalRange[1] = initTotalHorizontalRange[1];
        totalVerticalRange[0] = initTotalVerticalRange[0];
        totalVerticalRange[1] = initTotalVerticalRange[1];
        return new float[][]
            {totalHorizontalRange, totalVerticalRange};
    }

    /** for each class displayable on cursor, get the corresponding displayableSection for this cursorSection and display it */
    void displayTexture(StsGLPanel3d glPanel3d, boolean is3d)
    {
        if(displayableSections == null) return;
        for(StsCursor3dTexture displayableSection : displayableSections)
            displayableSection.displayTexture(glPanel3d, is3d, this);
    }
/*
    public void setSubVolumeFlags(boolean displaySubVolume)
    {
        if(displaySubVolume && subVolumePlane == null)
            computeSubVolumePlane();
    }
*/
    /** for each class displayable on cursor, get the corresponding displayableSection for this cursorSection and display it */
    void display(StsGLPanel3d glPanel3d, boolean is3d)
    {
        if(displayableSections == null) return;
        for(StsCursor3dTexture displayableSection : displayableSections)
            displayableSection.display(glPanel3d, is3d);
    }

    StsCursor3dTexture checkAddDisplayableSection(StsModel model, StsCursor3d cursor3d, StsClassCursor3dTextureDisplayable displayableClass)
    {
        if(displayableClass.getCurrentObject() == null)return null;
        Class objectClass = displayableClass.getObjectClass();
        for(StsCursor3dTexture displayableSection : displayableSections)
            if(displayableSection.canDisplayClass(objectClass)) return displayableSection;

        StsCursor3dTexture displayableSection = displayableClass.constructDisplayableSection(model, cursor3d, dir);
        if(displayableSection == null) return null;
        addDisplayableSection(displayableSection);
        return displayableSection;
    }

    void addDisplayableSection(StsCursor3dTexture displayableSection)
    {
        displayableSections = (StsCursor3dTexture[])StsMath.arrayAddElement(displayableSections, displayableSection);
    }

    public StsCursor3dTexture getDisplayableSection(Class objectClass)
    {
        for(int n = 0; n < displayableSections.length; n++)
            if(displayableSections[n].canDisplayClass(objectClass)) return displayableSections[n];
        return null;
    }

    public StsCursor3dTexture[] getDisplayableSections()
    {
        return displayableSections;
    }

    public StsCursor3dTexture[] getVisibleDisplayableSections()
    {
        StsCursor3dTexture[] visibleSections = new StsCursor3dTexture[displayableSections.length];
        int nVisibleSections = 0;
        for(StsCursor3dTexture displayableSection : displayableSections)
            if(displayableSection.isVisible) visibleSections[nVisibleSections++] = displayableSection;
        return (StsCursor3dTexture[])StsMath.trimArray(visibleSections, nVisibleSections);
    }

    public void deleteCursor3dTexture(StsObject object)
    {
        for(int n = displayableSections.length - 1; n >= 0; n--)
        {
            if(displayableSections[n].isDisplayableObject(object))
            {
                displayableSections = (StsCursor3dTexture[])StsMath.arrayDeleteElement(displayableSections, displayableSections[n]);
                break;
            }
        }
    }

    public void deleteClassTextureDisplays(Class displayableClass)
	{
       for(int n = displayableSections.length - 1; n >= 0; n--)
        {
            if(displayableSections[n].canDisplayClass(displayableClass))
            {
                displayableSections = (StsCursor3dTexture[])StsMath.arrayDeleteElement(displayableSections, displayableSections[n]);
                break;
            }
        }
    }

    public void deleteAllTextures(GL gl)
    {
        for(StsCursor3dTexture displayableSection : displayableSections)
            displayableSection.deleteTexturesAndDisplayLists(gl);
    }

    public boolean setObject(StsObject object)
    {
        boolean set = false;
        for(StsCursor3dTexture displayableSection : displayableSections)
        {
            if(!displayableSection.isDisplayableObject(object)) continue;
            Object displayedObject = displayableSection.getObject();
            set = true;
            if(displayedObject != object)
            {
                displayableSection.setObject(object);
            }
        }
        return set;
    }
    
    public boolean toggleOn(StsObject object)
    {
        for(StsCursor3dTexture displayableSection : displayableSections)
        {
            if(!displayableSection.isDisplayableObject(object)) continue;
            Object displayedObject = displayableSection.getObject();
            if(displayedObject != object)
            {
                displayableSection.setObject(object);
                return true;
            }
        }
        return false;
    }

    public boolean toggleOff(StsObject object)
    {
        for(StsCursor3dTexture displayableSection : displayableSections)
        {
            if(!displayableSection.isDisplayableObject(object)) continue;
            Object displayedObject = displayableSection.getObject();
            if(displayedObject == object)
            {
                displayableSection.setObject(null);
                return true;
            }
        }
        return false;
    }

    /** Clear the texture for sections displayed on this cursor plane. */
    public void clearTextureDisplays()
    {
        if(debug)System.out.println("Clearing textureDisplays for dirNo " + dir);
        for(StsCursor3dTexture displayableSection : displayableSections)
            displayableSection.clearTextureDisplay();
    }

    public boolean clearTextureDisplays(Object object)
    {
        boolean cleared = false;
        if(debug)System.out.println("Clearing textureDisplays for dirNo " + dir);
        for(StsCursor3dTexture displayableSection : displayableSections)
        {
            if(displayableSection.isDisplayingObject(object))
            {
                cleared = true;
                displayableSection.clearTextureDisplay();
            }
        }
        return cleared;
    }

    public boolean clearTextureDisplays(Class objectClass)
    {
        boolean cleared = false;
        if(debug)System.out.println("Clearing textureDisplays for dirNo " + dir);
        for(StsCursor3dTexture displayableSection : displayableSections)
        {
            Class displayedClass = displayableSection.getDisplayableClass();
            if(displayedClass == objectClass)
            {
                cleared = true;
                displayableSection.clearTextureDisplay();
            }
        }
        return cleared;
    }

    public boolean viewObjectChanged(Object object)
    {
        for(StsCursor3dTexture displayableSection : displayableSections)
        {
            if(displayableSection.isDisplayableObject(object))
            {
                if(displayableSection.setObject(object))
                {
                    if(debug)StsException.systemDebug(this, "viewObjectChanged", " dirNo: " + dir);
                    return true;
                }
            }
        }
        return false;
    }

    public boolean isDisplayableObject(Object object)
    {
        boolean isDisplayable = false;
        for(StsCursor3dTexture displayableSection : displayableSections)
        {
            if(displayableSection.isDisplayableObject(object))
                isDisplayable = true;
        }
        return isDisplayable;
    }

    public boolean isDisplayingObject(Object object)
    {
        for(StsCursor3dTexture displayableSection : displayableSections)
        {
            if(displayableSection.isDisplayingObject(object))
                return true;
        }
        return false;
    }

    public boolean objectChanged(Object object)
    {
        for(StsCursor3dTexture displayableSection : displayableSections)
        {
            if(displayableSection.isDisplayingObject(object))
            {

                displayableSection.clearTextureDisplay();
                return true;
            }
            /*
                if(displayableSections[n].isDisplayableObject(object))
                {

                    ( (StsTextureSurfaceFace)displayableSections[n]).clearTextureTileSurface(glPanel3d);
                    if (!displayableSections[n].isDisplayingObject(object))displayableSections[n].setViewObject(object);
                    return true;
                }
             */
        }
        return false;
    }

    public void cropChanged()
    {
        for(StsCursor3dTexture displayableSection : displayableSections)
            displayableSection.cropChanged();
    }

    public void subVolumeChanged()
    {
        for(StsCursor3dTexture displayableSection : displayableSections)
            displayableSection.subVolumeChanged();
    }

    public String rowColNumReadout(StsPoint point)
    {
        String valueString;
        StsRotatedGridBoundingBox boundingBox = cursor3d.model.getProject().getRotatedBoundingBox();
        if(dir == XDIR)
        {
//                int row = Math.round(boundingBox.getRowCoor(point.v[1]));
//                int col = Math.round(boundingBox.getSliceCoor(point.v[2]));
//                valueString = getValueString(row, col);
            float yLabel = boundingBox.getNumFromCoor(YDIR, point.v[1]);
            float zLabel = boundingBox.getNumFromCoor(ZDIR, point.v[2]);
            valueString = "Line: " + yLabel + " Slice: " + zLabel;
        }
        else if(dir == YDIR)
        {
//                int row = Math.round(boundingBox.getColCoor(point.v[0]));
//                int col = Math.round(boundingBox.getSliceCoor(point.v[2]));
//                valueString = getValueString(row, col);
            float xLabel = boundingBox.getNumFromCoor(XDIR, point.v[0]);
            float zLabel = boundingBox.getNumFromCoor(ZDIR, point.v[2]);
            valueString = "Crossline: " + xLabel + " Slice: " + zLabel;
        }
        else
        {
//                int row = Math.round(boundingBox.getRowCoor(point.v[1]));
//                int col = Math.round(boundingBox.getColCoor(point.v[0]));
//                valueString = getValueString(row, col);
            float xLabel = boundingBox.getNumFromCoor(XDIR, point.v[0]);
            float yLabel = boundingBox.getNumFromCoor(YDIR, point.v[1]);
            valueString = "Line: " + yLabel + " Crossline: " + xLabel;
        }
        /*
            StsSeismicVolume[] svs = (StsSeismicVolume[])glPanel3d.model.getCastObjectList(StsSeismicVolume.class);
            for(int i=0; i< svs.length; i++)
            {
                if((svs[i] == line2d) || (!svs[i].getReadoutEnabled())) continue;
                valueString = valueString + " -- " + svs[i].getName() + ": " +
                    new Float(svs[i].getValueFromByte((byte)svs[i].getValue(point.v))).toString();
            }
         */
        return valueString;
    }

    /**
     * Output the coordinate point
     * @param point relative coordinate
     * @return the output string
     */

    public String propertyReadout(StsGLPanel3d glPanel3d, StsPoint point)
    {
        StringBuffer stringBuffer = null;
        for(StsCursor3dTexture displayableSection : displayableSections)
        {
            String valueString = displayableSection.propertyReadout(point);
            if(valueString != null)
            {
                if(stringBuffer == null)stringBuffer = new StringBuffer();
                stringBuffer.append(valueString + " ");
            }
        }
        if(stringBuffer == null)return null;
        else return stringBuffer.toString();
    }

    public boolean cursorIntersected(StsBoundingBox boundingBox)
    {
        switch(dir)
        {
            case XDIR:
                float xInc = 2 * rotatedBoundingBox.xInc;
                return dirCoordinate >= boundingBox.xMin - xInc && dirCoordinate <= boundingBox.xMax + xInc;
            case YDIR:
                float yInc = 2 * rotatedBoundingBox.yInc;
                return dirCoordinate >= boundingBox.yMin - yInc && dirCoordinate <= boundingBox.yMax + yInc;
            case ZDIR:

//                    float zInc = 2*rotatedBoundingBox.zInc;
//                    return dirCoordinate >= boundingBox.zMin - zInc && dirCoordinate <= boundingBox.zMax + zInc;
                return true;
            default:
                return false;
        }
    }
    
    static protected StsPolygon concavePolygon = null;
    static private float[][] testPolygonScaledPoints = new float[][]
     {
         {0.25f, 0.25f},
         {0.25f, 0.75f},
         {0.75f, 0.75f},
         {0.5f, 0.5f},
         {0.75f, 0.25f}
     };
 //    static private float[][] testPolygonScaledPoints = new float[][] {
 //        { 0.75f, 0.25f }, { 0.5f, 0.5f }, { 0.75f, 0.75f }, { 0.25f, 0.75f }, { 0.25f, 0.25f } };
     static private float[][] testPolygonPoints = new float[5][];

    public void displayTestConcavePolygon(GL gl, GLU glu)
    {

        if(concavePolygon == null)concavePolygon = new StsPolygon();
        int nPoints = testPolygonPoints.length;
        for(int n = 0; n < nPoints; n++)
            testPolygonPoints[n] = rotatedBoundingBox.getCursorPntNrmlFromScaledPoint(testPolygonScaledPoints[n], dir, dirCoordinate);
        concavePolygon.setPntNrmls(testPolygonPoints);
        float[] normal;
        if(dir == XDIR)
        {
            normal = new float[]
                { -1.0f, 0.0f, 0.0f};
            concavePolygon.setStsColor(StsColor.RED);
        }
        else if(dir == YDIR)
        {
            normal = new float[]
                {0.0f, -1.0f, 0.0f};
            concavePolygon.setStsColor(StsColor.GREEN);
        }
        else
        {
            normal = new float[]
                {0.0f, 0.0f, -1.0f};
            concavePolygon.setStsColor(StsColor.BLUE);
        }

        gl.glEnable(GL.GL_LIGHTING); // shouldn't have to do this..should be enabled already
        System.out.println("Start test polygon draw.");
//            concavePolygon.draw(gl, glu, false, null, true);
        System.out.println("End test polygon draw.");
        concavePolygon.draw(gl, glu, true, normal, true, false, null, false);
    }
    public String getLabel() { return toString(); }

    public String toString() { return "CursorSection "  + StsParameters.coorLabels[dir]; }
    
    public float getColCoor(float[] xyz)
    {
        return rotatedBoundingBox.getCursorNormalizedColCoor(dir, xyz);
    }

    public float getRowCoor(float[] xyz)
    {
        return rotatedBoundingBox.getCursorNormalizedRowCoor(dir, xyz);
    }

    public StsPoint getPoint(int row, int col)
    {
        StsException.systemError(this, "getPoint", "Not implemented.  Shouldn't be called.");
        return null;
    }

    public float[] getXYZorT(int row, int col)
    {
        StsException.systemError(this, "getXYZorT", "Not implemented.  Shouldn't be called.");
        return null;
    }

    public StsPoint getPoint(float row, float col)
    {
        StsException.systemError(this, "getPoint", "Not implemented.  Shouldn't be called.");
        return null;
    }

    public float[] getXYZorT(float row, float col)
    {
        StsException.systemError(this, "getXYZorT", "Not implemented.  Shouldn't be called.");
        return null;
    }

    public float[] getNormal(int row, int col)
    {
       StsException.systemError(this, "getNormal", "Not implemented.  Shouldn't be called.");
        return null;
    }

    public float[] getNormal(float row, float col)
    {
       StsException.systemError(this, "getNormal", "Not implemented.  Shouldn't be called.");
        return null;
    }

    public void checkConstructGridNormals()
    {
       StsException.systemError(this, "checkConstructGridNormals", "Not implemented.  Shouldn't be called.");
    }

    public int getNRows() { return 2; }
    public int getNCols() { return 2; }
    public int getRowMin() { return 0; }
    public int getRowMax() { return 1; }
    public int getColMin() { return 0; }
    public int getColMax() { return 0; }

    public StsGridBoundingBox getGridBoundingBox()
    {
        return new StsGridBoundingBox(0, 1, 0, 1);
    }
}
