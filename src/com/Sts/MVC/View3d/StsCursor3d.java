//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.MVC.View3d;

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;
import com.magician.fonts.*;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;
import java.awt.event.*;
import java.io.*;
import java.text.*;

// This is a controller.
// All draws should go back to model and then out to the views.

// NB: truly serializable, beware of incompatible serializations!

public class StsCursor3d implements Serializable
{
    private boolean display3dCursor;
    /** this is current direction being dragged or last dragged: drawn bold */
    private int currentDirection = NONE;
    /** this is direction which has been selected with X, Y, Z buttons on cursor3dPanel */
    protected int selectedDirection = NONE;
    private int dragDirection = NONE;
    private boolean drawCursorPlanes = true;
    /** all cursor sections - default is three, one in each direction */
    public StsCursorSection[] cursorSections;

    transient public StsRotatedGridBoundingBox rotatedBoundingBox;
    /** Classes other than texture classes which can be displayed on cursor3d */
    transient protected StsArrayList displayableClasses = null;

    transient public StsModel model = null;
    //     transient public StsGLPanel3d glPanel3d = null;
    transient public StsWin3dFull window;
    /** row of current Cursor position - row of current superGather is more accurate, but this is used for storing current position in database */
    private int row;
    /** col of current Cursor position - col of current superGather is more accurate, but this is used for storing current position in database */
    private int col;

    static StsTimer timer = new StsTimer();
    static boolean runTimer = false;
    static boolean debug = false;

    static final short DOTTED_LINE = (short)0xCCCC;
    static private boolean testConcavePolygon = false;

    static final String noCursorTitle = "No cursor section.";
    /** static variable NONE = -1 - No cursors on */
    public static final int NONE = StsParameters.NO_MATCH;
    /** static variable XDIR = 0 - X oriented cursor section */
    public static final int XDIR = StsParameters.XDIR;
    /** static variable YDIR = 1 - Y oriented cursor section */
    public static final int YDIR = StsParameters.YDIR;
    /** static variable ZDIR = 2 - Z oriented cursor section */
    public static final int ZDIR = StsParameters.ZDIR;
    /** Default orientation labels */
    public static final String[] coorLabels = new String[]{"X", "Y", "Z"};

    static StsColor[] lineColors = new StsColor[]{StsColor.RED, StsColor.GREEN, StsColor.BLUE}; // line colors of cursor planes

    static final long serialVersionUID = 1l;

    public StsCursor3d()
    {
    }

    public StsCursor3d(StsModel model, StsWin3dFull window)
    {
        this.model = model;
        this.window = window;
        display3dCursor = true;
    }

    /* classInitialize from a database restore */
    public void initializeTransients(StsModel model, StsWin3dFull window)
    {
        this.model = model;
        this.window = window;
        display3dCursor = true;
        rotatedBoundingBox = model.getProject().getRotatedBoundingBox();
        if(cursorSections == null) initialize();
        for(int n = 0; n < 3; n++)
            cursorSections[n].initializeTransients(model, this, n);
        initializeCursorSections();
        if (StsPreStackLineSetClass.currentProjectPreStackLineSet == null) return;
        StsPreStackLineSetClass.currentProjectPreStackLineSet.setSuperGatherToCursor(window, this);
    }

    public void initializeCursors()
    {
        if(model.win3d == null) return;
        for(int n = 0; n < 3; n++)
        {
            cursorSections[n].initializeCursorSection();
            model.win3d.cursor3dPanel.setSliderSelected(n, cursorSections[n].isSelected);
        }
    }

    /**
     * Reposition the specified cursor section position.
     *
     * @param dir   cursor direction
     * @param dirCoor cursor position in project coordinates
     */
    public void adjustCursor(int dir, float dirCoor)
    {
        if(cursorSections == null) return;
        if(dir == -1) return;
        cursorSections[dir].setDirCoordinate(dirCoor);
        StsWindowFamily family = window.getWindowFamily();
        family.adjustCursor(dir, dirCoor);
        if (rotatedBoundingBox == null ) return;
        if (dir == StsCursor3d.XDIR)
        {
            float colNum = rotatedBoundingBox.getNumFromCoor(dir, dirCoor);
            col = rotatedBoundingBox.getColFromColNum(colNum);
        }
        else if(dir == StsCursor3d.YDIR)
        {
            float rowNum = rotatedBoundingBox.getNumFromCoor(dir, dirCoor);
            row = rotatedBoundingBox.getRowFromRowNum(rowNum);
        }
        subVolumeChanged(dir);
    }

    public void subVolumeChanged()
    {
        if(cursorSections == null) return;
        StsSubVolumeClass subVolumeClass = (StsSubVolumeClass)model.getStsClass(StsSubVolume.class);
        boolean displaySubVolume = subVolumeClass.getIsApplied();

        if(cursorSections == null) return;
        for(int dir = 0; dir < 3; dir++)
        {
            boolean display = displaySubVolume && dir != dragDirection;
            subVolumeChanged(dir, display);
        }
    }

    public void subVolumeChanged(int dir)
    {
        if(cursorSections == null) return;
        StsSubVolumeClass subVolumeClass = (StsSubVolumeClass)model.getStsClass(StsSubVolume.class);
        boolean display = subVolumeClass.getIsApplied() && dir != dragDirection;
        subVolumeChanged(dir, display);
    }

    public void subVolumeChanged(int dir, boolean display)
    {
        cursorSections[dir].clearTextureDisplays();
        cursorSections[dir].clearSubVolumePlane();
        // if(display) cursorSections[dir].computeSubVolumePlane();
    }

    /** Reposition the x and y cursor section in one action. */
    public void adjustCursorXY(float xCoor, float yCoor)
    {
        adjustCursor(XDIR,  xCoor);
        adjustCursor(YDIR,  yCoor);
    }

    /**
     * Get the cropped horizontal range
     *
     * @param dirNo cursor section direction
     * @return minimum and maximum
     */
    public float[] getTotalHorizontalRange(int dirNo)
    {
        return cursorSections[dirNo].totalHorizontalRange;
    }

    /**
     * Get the cropped vertical range
     *
     * @param dirNo cursor section direction
     * @return minimum and maximum
     */
    public float[] getTotalVerticalRange(int dirNo)
    {
        return cursorSections[dirNo].totalVerticalRange;
    }

    /**
     * Get the Horizontal Axis label
     *
     * @param dirNo direction
     * @return horizontal label
     */
    public String getHorizontalAxisLabel(int dirNo)
    {
        return cursorSections[dirNo].horizontalAxisLabel;
    }

    /**
     * Get the Vertical Axis label
     *
     * @param dirNo direction
     * @return vertical label
     */
    public String getVerticalAxisLabel(int dirNo)
    {
        return cursorSections[dirNo].verticalAxisLabel;
    }

    /**
     * Construct the label to be displayed on any title
     *
     * @param dirNo direction
     * @param label additional labeling
     * @return
     */
    public String getTitleLabel(int dirNo, String label)
    {
        if(dirNo == NONE) return noCursorTitle + label;
        return cursorSections[dirNo].titleLabel + label;
    }

    /**
     * Set whether the user is working in grid coordinates or relative coordinates
     *
     * @param isGridCoordinates is true if in grid coordinates
     */
    public void setIsGridCoordinates(boolean isGridCoordinates)
    {
        window.setCursorDisplayXY(!isGridCoordinates);
    }

    public boolean getIsGridCoordinates()
    {
        return !window.getCursorDisplayXY();
    }

    /**
     * Turn the 3D cursor on and off. Change the state of the 3D cursor button in the default Toolbar
     *
     * @param display3dCursor true if turning on
     * @see StsToggleButton#setSelected(boolean)
     */
    public void setDisplay3dCursor(boolean display3dCursor)
    {
        if(this.display3dCursor == display3dCursor) return;
        this.display3dCursor = display3dCursor;
        //        model.win3dDisplayAll(glPanel3d.window);

        // Determine which window is current
        if(window == null) return;
        StsToggleButton toggleButton = (StsToggleButton)window.getToolbarComponentNamed("Default Toolbar", "cursor3d");
        toggleButton.setSelected(display3dCursor);
    }

    /**
     * Is the 3D cursor isVisible. Turned off and on from the Object Panel
     *
     * @return true if displayed
     */
    public boolean isDisplay3dCursor()
    {
        return display3dCursor;
    }

    /** Initialize the 3d cursors based on Project extents */
    public void copy(StsCursor3d copyCursor3d)
    {
        //		initializeSuperGather();
        if(copyCursor3d == null) return;
        currentDirection = copyCursor3d.currentDirection;
        if(currentDirection == NONE) currentDirection = 2;
        rotatedBoundingBox = model.getProject().getRotatedBoundingBox();
        cursorSections = new StsCursorSection[3];
        for(int n = 0; n < 3; n++)
        {
            float dirCoordinate = copyCursor3d.cursorSections[n].dirCoordinate;
            cursorSections[n] = new StsCursorSection(this, n, dirCoordinate);

        }
        initializeCursorSections();
        setSelectedDirection(currentDirection, true);
        rangeChanged();
    }

    /** Initialize the 3d cursors based on Project extents */
    public void initialize()
    {
        rotatedBoundingBox = model.getProject().getRotatedBoundingBox();

        if(cursorSections == null)
        {
            cursorSections = new StsCursorSection[3];
            for(int n = 0; n < 3; n++)
            {
                float dirCoordinate = getInitialDirCoordinate(n);
                cursorSections[n] = new StsCursorSection(this, n, dirCoordinate);
            }
            if(currentDirection == NONE) currentDirection = 2;
            setSelectedDirection(currentDirection, true);
        }
        else
            rangeChanged();

        initializeCursorSections();
    }

    public void initializeCursorSections()
    {
        for(int dir = 0; dir < 3; dir++)
            cursorSections[dir].initialize();
        window.checkAddDisplayableSections();
    }

    public void rangeChanged()
    {
        for(int n = 0; n < 3; n++)
        {
            cursorSections[n].rangeChanged();
            cursorSections[n].setAxisLabel();
        }
        try
        {
            if(window != null && window.cursor3dPanel != null)
                window.cursor3dPanel.setSliderValues();
        }
        catch(Exception e)
        {
            e.printStackTrace();
            System.out.println("Cursor could not be set");
        }
    }

    public void zRangeChanged()
    {
        for(int n = 0; n < 3; n++)
            cursorSections[n].zRangeChanged();
        cursorSections[ZDIR].setAxisLabel();

        try
        {
            if(window != null && window.cursor3dPanel != null)
                window.cursor3dPanel.setSliderValues();
        }
        catch(Exception e)
        {
            e.printStackTrace();
            System.out.println("Cursor could not be set");
        }
    }

    public void setAxisLabels()
    {
        for(int n = 0; n < 3; n++)
            cursorSections[n].setAxisLabel();
    }

    private float getInitialDirCoordinate(int dir)
    {
        switch(dir)
        {
            case XDIR:
                return rotatedBoundingBox.xMax;
            case YDIR:
                return rotatedBoundingBox.yMax;
            case ZDIR:
                return rotatedBoundingBox.zMax;
            default:
                return StsParameters.nullValue;
        }
    }

    public void resetInitialCursorPositions()
    {
        adjustCursor(XDIR, rotatedBoundingBox.xMin);
        adjustCursor(YDIR, rotatedBoundingBox.yMin);
        adjustCursor(ZDIR, rotatedBoundingBox.zMax);
        currentDirection = ZDIR;
    }

    public void resetZDomain(float zMin, float zMax, float z)
    {
        for(int n = 0; n < 3; n++)
            cursorSections[n].resetZRange(zMin, zMax, z);

    }

    /**
     * Reset depth labels based on user specified value. Label is a property of the project
     * since all data must be in same units.
     */
    public void resetDepthLabels(String domainString)
    {
        cursorSections[0].verticalAxisLabel = domainString;
        cursorSections[1].verticalAxisLabel = domainString;
        cursorSections[2].titleLabel = domainString + " Slice ";
    }

    /**
     * @param dirNo
     * @return
     */
    public StsPoint[] getPlanePoints(int dirNo)
    {
        return cursorSections[dirNo].planePoints;
    }

    /**
     * Get the cropped cursor extents
     *
     * @return minimum and maximum in both dimensions [a][min][a][max][b][min][b][max]
     */
    public float[][] getTotalAxisRanges()
    {
        if(currentDirection == NONE) return null;
        else return cursorSections[currentDirection].getTotalAxisRanges();
    }

    /**
     * Reset the current cursor section ranges to the full project extents
     *
     * @return minimum and maximum in both dimensions [a][min][a][max][b][min][b][max]
     */

    public float[][] resetAxisRanges()
    {
        if(currentDirection == NONE) return null;
        return cursorSections[currentDirection].resetAxisRanges();
    }

    /**
     * Get the current cursor section direction
     *
     * @return current direction (XDIR, YDIR, ZDIR)
     */
    public int getCurrentDirNo()
    {
        return currentDirection;
    }

    public int getSelectedDirNo()
    {
        return selectedDirection;
    }

    /**
     * Set the current cursor section direction
     *
     * @param point new curspr position
     */
    public void setCoordinates(StsPoint point)
    {
        window.cursor3dPanel.setSliderValues(point.getX(), point.getY(), point.getZ());
    }

    /**
     * Set the current cursor section direction
     *
     * @param newDirNo new direction (XDIR, YDIR or ZDIR)
     */
    public void setCurrentDirNo(int newDirNo)
    {
        if(currentDirection == newDirNo) return;
        currentDirection = newDirNo;
        if(window.cursor3dPanel == null) return;
        window.cursor3dPanel.setCurrentDirNo(newDirNo);
    }

    private void setSelectedDirection(int dir)
    {
        if(selectedDirection == dir) return;
        selectedDirection = dir;
        StsWindowFamily family = window.getWindowFamily();
        family.viewChanged();

    }

    /*
       public int getCurrentCursorSectionNRows()
       {
           if(currentDirNo == NONE) return 0;
           return rotatedBoundingBox.getCursorRowMax(currentDirNo) + 1;
       }
    */
    /**
     * Get the index of the current cursor section. The current cursor section is the last one moved.
     * @param dir cursor direction (XDIR, YDIR, ZDIR)
     * @return index
     */
    /*
       public int getCurrentPlaneIndex()
           return getPlaneIndex(currentDirNo);
       }
    */
    /**
     * Get the index of the specified cursor section.
     * @param dir cursor direction (XDIR, YDIR, ZDIR)
     * @return index
     */
    /*
       public int getPlaneIndex(int dirNo)
       {
           if(dirNo == NONE) return -1;
           return cursorSectionDisplayables[dirNo].getPlaneIndex();
       }
    */
    /**
     * Get the position of the current cursor section. The current cursor section is the last one moved.
     *
     * @return position in relative project coordinates
     */
    public float getCurrentDirCoordinate()
    {
        return getCurrentDirCoordinate(currentDirection);
    }


    /**
     * Get the position of the specified cursor
     *
     * @param dir cursor direction (XDIR, YDIR, ZDIR)
     * @return position in relative project coordinates
     */
    public float getCurrentDirCoordinate(int dir)
    {
        if(dir == NONE) return 0.0f;
        else return cursorSections[dir].dirCoordinate;
    }

    public int getCurrentGridCoordinate(int dir)
    {
        float dirCoor = getCurrentDirCoordinate(dir);
        StsRotatedGridBoundingBox boundingBox = model.getProject().getRotatedBoundingBox();
        return boundingBox.getCursorPlaneIndex(dir, dirCoor);
    }
/*
    public StsCursor3dTexture[] getDisplayableSections(StsGLPanel3d glPanel3d)
    {
        return getDisplayableSections(glPanel3d, currentDirection);
    }

    public StsCursor3dTexture getDisplayableSection(int dirNo, StsClassCursor3dTextureDisplayable displayableClass)
    {
       return cursorSections[dirNo].checkAddDisplayableSection(displayableClass);
    }

    public StsCursor3dTexture[] getDisplayableSections(StsGLPanel3d glPanel3d, int dir)
    {
        if(dir == NONE) return null;
        return glPanel3d.getDisplayableSections(dir);
    }

    public StsCursor3dTexture[] getVisibleDisplayableSections(StsGLPanel3d glPanel3d, int dir)
    {
        if(dir == NONE) return new StsCursor3dTexture[0];
        return glPanel3d.getVisibleDisplayableSections(dir);
    }
*/
    public StsCursor3dTexture getDisplayableSection(int dir, Class displayableStsClass)
    {
        if(dir == NONE) return null;
        return cursorSections[dir].getDisplayableSection(displayableStsClass);
    }

    public void deleteClassTextureDisplays(Class displayableClass)
	{
        for(StsCursorSection cursorSection : cursorSections)
            cursorSection.deleteClassTextureDisplays(displayableClass);
    }

    public void deleteAllTextures(GL gl)
    {
        for(StsCursorSection cursorSection : cursorSections)
            cursorSection.deleteAllTextures(gl);
    }

    public void setSliderState(int direction, boolean isDragging)
    {
        currentDirection = direction;

         if(!isDragging)
         {
            this.dragDirection = NONE;
             if(debug) System.out.println("SliderState set: not dragging");
         }
         else
         {
             this.dragDirection = direction;
             if(debug) System.out.println("SliderState set: dragging dir: " + direction);
         }
    }

    /**
     * <Not currently used>
     *
     * @return
     */
    public boolean getIsDragging(int dir)
    {
        return dragDirection == dir;
    }

    public boolean getIsDragging()
    {
        return currentDirection != NONE && dragDirection == currentDirection;
    }

    /**
     * Is cursor section display enabled?
     *
     * @param dir direction (XDIR, YDIR, ZDIR)
     * @return true if displayable
     */
    public boolean isDisplayEnable(int dir)
    {
        if(dir == -1) return false;
        return cursorSections[dir].isSelected;
    }

    /**
     * Set the cursor section to displayable. Will not be drawn if not true
     *
     * @param dir    cursor direction (XDIR, YDIR, ZDIR)
     * @param enable true to display
     */
    public void setSelectedDirection(int dir, boolean enable)
    {
        if(dir == -1) return;
        if(cursorSections != null)
            cursorSections[dir].isSelected = enable;
        if(enable)
        {
            setCurrentDirNo(dir);
            setSelectedDirection(dir);
        }
    }

    public int getSelectedDirection() { return selectedDirection; }

    /**
     * Set the draw planes state. If set to true, redraw will be forced on next cycle
     *
     * @param draw true if stale
     */
    public void setDrawCursorPlanes(boolean draw)
    {
        drawCursorPlanes = draw;
    }

    /**
     * Are the cursor planes stale, used to determine if redraw is required
     *
     * @return true if stale
     */
    public boolean getDrawCursorPlanes()
    {
        return drawCursorPlanes;
    }

    /**
     * Display the cursor sections and all related cursor displayable objects
     *
     * @param glPanel3d
     * @param view3d
     * @param gl        graphics context
     * @param glu
     */
    public void display3d(StsGLPanel3d glPanel3d, StsView3d view3d, GL gl, GLU glu)
    {
        if(!drawCursorPlanes) return;

        displayableClasses = model.getCursorDisplayableClasses();

        gl.glDisable(GL.GL_LIGHTING);

        for(int i = 0; i < 3; i++)
            if(cursorSections[i].isSelected) display3d(glPanel3d, view3d, gl, i, cursorSections[i].dirCoordinate);

        if(model.getProject().getShowIntersection()) drawIntersections(glPanel3d, gl);
        gl.glEnable(GL.GL_LIGHTING);
    }

    private void drawIntersections(StsGLPanel3d glPanel3d, GL gl)
    {
        gl.glDisable(GL.GL_LIGHTING);

        StsPoint[] xPlane = cursorSections[XDIR].planePoints;
        StsPoint[] yPlane = cursorSections[YDIR].planePoints;
        StsPoint[] zPlane = cursorSections[ZDIR].planePoints;

        // view shift towards the user so always in front of cursor section
        glPanel3d.setViewShift(gl, StsGraphicParameters.edgeShift);

        StsPoint[] points = new StsPoint[2];

        // XY planes intersection
        points[0] = new StsPoint(xPlane[0].v[0], yPlane[0].v[1], xPlane[0].v[2]);
        points[1] = new StsPoint(xPlane[1].v[0], yPlane[0].v[1], xPlane[1].v[2]);
        StsGLDraw.drawDottedLine(gl, StsColor.WHITE, 1, points);

        // XZ planes intersection
        points[0] = new StsPoint(xPlane[0].v[0], xPlane[0].v[1], zPlane[0].v[2]);
        points[1] = new StsPoint(xPlane[2].v[0], xPlane[2].v[1], zPlane[0].v[2]);
        StsGLDraw.drawDottedLine(gl, StsColor.WHITE, 1, points);

        // YZ planes intersection
        points[0] = new StsPoint(yPlane[0].v[0], yPlane[0].v[1], zPlane[0].v[2]);
        points[1] = new StsPoint(yPlane[2].v[0], yPlane[2].v[1], zPlane[0].v[2]);
        StsGLDraw.drawDottedLine(gl, StsColor.WHITE, 1, points);

        glPanel3d.resetViewShift(gl);
        gl.glEnable(GL.GL_LIGHTING);
    }

    private void display3d(StsGLPanel3d glPanel3d, StsView3d view3d, GL gl, int dirNo, float dirCoordinate)
    {
        int i;
        boolean thisDirDragging;

        GLBitmapFont horizontalFont = GLHelvetica12BitmapFont.getInstance(gl);
        DecimalFormat labelFormat = new DecimalFormat("###0");

        if(dirNo == NONE) return;

        StsPoint[] planePoints = cursorSections[dirNo].planePoints;

        if(dirNo == currentDirection)
        {
            gl.glLineWidth(3.0f);
            //           thisDirDragging = isDragging;
        }
        else
        {
            gl.glLineWidth(1.0f);
            //            thisDirDragging = false;
        }
        thisDirDragging = (dragDirection == currentDirection);

        switch(dirNo)
        {
            case XDIR:
                StsColor.RED.setGLColor(gl);
                break;
            case YDIR:
                StsColor.GREEN.setGLColor(gl);
                break;
            case ZDIR:
                StsColor.BLUE.setGLColor(gl);
                break;
        }

        gl.glDisable(GL.GL_LIGHTING);

        /* Draw bounding lines around section plane */
        gl.glBegin(GL.GL_LINE_STRIP);

        for(i = 0; i < 4; i++)
            gl.glVertex3f(planePoints[i].v[0], planePoints[i].v[1], planePoints[i].v[2]);
        gl.glVertex3f(planePoints[0].v[0], planePoints[0].v[1], planePoints[0].v[2]);

        gl.glEnd();

        // Draw Labels
        String cursorLabel;
        float vec[] = new float[3];
        float x = 0, xMin = rotatedBoundingBox.xMin, xMax = rotatedBoundingBox.xMax;
        float y = 0, yMin = rotatedBoundingBox.yMin, yMax = rotatedBoundingBox.yMax;
        float z = 0, zMin = rotatedBoundingBox.zMin, zMax = rotatedBoundingBox.zMax;
        if(model.getProject().getShowLabels())
        {
            // Draw Edge Labels
            float sliceCoor = getCurrentDirCoordinate();
            float xExt = (rotatedBoundingBox.xMax - rotatedBoundingBox.xMin) * 0.05f;
            float yExt = (rotatedBoundingBox.yMax - rotatedBoundingBox.yMin) * 0.05f;
            float zExt = (rotatedBoundingBox.zMax - rotatedBoundingBox.zMin) * 0.05f;

            DecimalFormat coorFormat = model.getProject().getLabelFormat();

            if(!getIsGridCoordinates()) // Actual Coordinates
                cursorLabel = new String(coorFormat.format(sliceCoor));
            else // Relative Coordinates
            {
                float cursorNum = rotatedBoundingBox.getNumFromCoor(currentDirection, sliceCoor);
                cursorLabel = new String(coorFormat.format(cursorNum));
            }
            switch(currentDirection)
            {
                case XDIR:
                    StsColor.RED.setGLColor(gl);
                    vec[0] = sliceCoor;
                    vec[1] = rotatedBoundingBox.yMax + yExt;
                    vec[2] = rotatedBoundingBox.zMin - zExt;
                    break;
                case YDIR:
                    StsColor.GREEN.setGLColor(gl);
                    vec[0] = rotatedBoundingBox.xMax + xExt;
                    vec[1] = sliceCoor;
                    vec[2] = rotatedBoundingBox.zMin - zExt;
                    break;
                case ZDIR:
                    StsColor.BLUE.setGLColor(gl);
                    vec[0] = rotatedBoundingBox.xMax + xExt;
                    vec[1] = rotatedBoundingBox.yMax + yExt;
                    vec[2] = sliceCoor;
                    break;
            }
            StsGLDraw.fontOutput(gl, vec, cursorLabel, horizontalFont);

            double coord = 0.0f;
            float angle = view3d.getViewAzimuth();
            float gridDX = (xMax - xMin) / 10.0f, gridDY = (yMax - yMin) / 10.0f, gridDZ = (zMax - zMin) / 10.0f;
            if(model.getProject().getLabelFrequency() != model.getProject().AUTO)
            {
                gridDX = model.getProject().getLabelFrequency();
                gridDY = model.getProject().getLabelFrequency();
                gridDZ = model.getProject().getLabelFrequency();
            }
            switch(dirNo)
            {
                case ZDIR:
                    StsColor.BLUE.setGLColor(gl);
                    // Draw the Z slice with grid and labels
                    vec[2] = getCurrentDirCoordinate(ZDIR);
                    for(x = xMin + gridDX; x < xMax; x += gridDX)
                    {
                        if(!getIsGridCoordinates())
                            coord = x;
                        else
                            coord = rotatedBoundingBox.getNumFromCoor(StsCursor3d.XDIR, x);

                        cursorLabel = labelFormat.format(coord);
                        vec[0] = x;
                        if((angle > 90) && (angle < 270))
                            vec[1] = yMax + yExt;
                        else
                            vec[1] = yMin - yExt;
                        StsGLDraw.fontOutput(gl, vec, cursorLabel, horizontalFont);
                    }
                    for(y = yMin + gridDY; y < yMax; y += gridDY)
                    {
                        if(!getIsGridCoordinates())
                            coord = y;
                        else
                            coord = rotatedBoundingBox.getNumFromCoor(StsCursor3d.YDIR, y);

                        cursorLabel = labelFormat.format(coord);
                        if((angle > 180) && (angle < 360))
                            vec[0] = xMin - xExt;
                        else
                            vec[0] = xMax + xExt;
                        vec[1] = y;
                        StsGLDraw.fontOutput(gl, vec, cursorLabel, horizontalFont);
                    }
                    break;
                case XDIR:
                    // draw the X Slice with labels
                    StsColor.RED.setGLColor(gl);
                    vec[0] = getCurrentDirCoordinate(XDIR);
                    for(z = zMin + gridDZ; z < zMax; z += gridDZ)
                    {
                        if(!getIsGridCoordinates())
                            coord = z;
                        else
                            coord = rotatedBoundingBox.getNumFromCoor(StsCursor3d.ZDIR, z);

                        cursorLabel = labelFormat.format(coord);
                        vec[2] = z;
                        if((angle > 90) && (angle < 270))
                            vec[1] = yMax + yExt;
                        else
                            vec[1] = yMin - yExt;
                        StsGLDraw.fontOutput(gl, vec, cursorLabel, horizontalFont);
                    }
                    for(y = yMin + gridDY; y < yMax; y += gridDY)
                    {
                        if(!getIsGridCoordinates())
                            coord = y;
                        else
                            coord = rotatedBoundingBox.getNumFromCoor(StsCursor3d.YDIR, y);

                        cursorLabel = labelFormat.format(coord);
                        vec[2] = zMax + zExt;
                        vec[1] = y;
                        StsGLDraw.fontOutput(gl, vec, cursorLabel, horizontalFont);
                    }
                    break;
                case YDIR:
                    StsColor.GREEN.setGLColor(gl);
                    // draw the Y Slice with labels
                    vec[1] = getCurrentDirCoordinate(YDIR);
                    for(z = zMin + gridDZ; z < zMax; z += gridDZ)
                    {
                        if(!getIsGridCoordinates())
                            coord = z;
                        else
                            coord = rotatedBoundingBox.getNumFromCoor(StsCursor3d.ZDIR, z);

                        cursorLabel = labelFormat.format(coord);
                        vec[2] = z;
                        if((angle > 90) && (angle < 270))
                            vec[0] = xMax + xExt;
                        else
                            vec[0] = xMin - xExt;
                        StsGLDraw.fontOutput(gl, vec, cursorLabel, horizontalFont);
                    }
                    for(x = xMin + gridDX; x < xMax; x += gridDX)
                    {
                        if(!getIsGridCoordinates())
                            coord = x;
                        else
                            coord = rotatedBoundingBox.getNumFromCoor(StsCursor3d.XDIR, x);

                        cursorLabel = labelFormat.format(coord);
                        vec[2] = zMax + zExt;
                        vec[0] = x;
                        StsGLDraw.fontOutput(gl, vec, cursorLabel, horizontalFont);
                    }
                    break;
                default:
                    break;
            }
        }

        if(model.getProject().getShow3dGrid())
        {
            StsColor.LIGHT_GRAY.setGLColor(gl);
            float gridDX = (xMax - xMin) / 10.0f, gridDY = (yMax - yMin) / 10.0f, gridDZ = (zMax - zMin) / 10.0f;
            if(model.getProject().getGridFrequency() != model.getProject().AUTO)
            {
                gridDX = model.getProject().getGridFrequency();
                gridDY = model.getProject().getGridFrequency();
            }
            if(model.getProject().getZGridFrequency() != model.getProject().AUTO)
                gridDZ = model.getProject().getZGridFrequency();

            gl.glLineWidth(StsGraphicParameters.gridLineWidth);
            gl.glLineStipple(1, DOTTED_LINE);
            gl.glEnable(GL.GL_LINE_STIPPLE);
            glPanel3d.setViewShift(gl, StsGraphicParameters.edgeShift);

            switch(dirNo)
            {
                case ZDIR:
                    StsColor.BLUE.setGLColor(gl);
                    vec[2] = getCurrentDirCoordinate(ZDIR);
                    gl.glBegin(GL.GL_LINES);
                {
                    /* Draw the grid lines in the x-direction */
                    for(x = xMin + gridDX; x < xMax; x += gridDX)
                    {
                        vec[0] = x;
                        vec[1] = yMin;
                        gl.glVertex3fv(vec, 0);
                        vec[0] = x;
                        vec[1] = yMax;
                        gl.glVertex3fv(vec, 0);
                    }

                    /* Draw the grid lines in the y-direction */
                    for(y = yMin + gridDY; y < yMax; y += gridDY)
                    {
                        vec[0] = xMin;
                        vec[1] = y;
                        gl.glVertex3fv(vec, 0);
                        vec[0] = xMax;
                        vec[1] = y;
                        gl.glVertex3fv(vec, 0);
                    }
                }
                gl.glEnd();
                break;
                case XDIR:
                    StsColor.RED.setGLColor(gl);
                    vec[0] = getCurrentDirCoordinate(XDIR);
                    gl.glBegin(GL.GL_LINES);
                {
                    for(y = yMin + gridDY; y < yMax; y += gridDY)
                    {
                        vec[2] = zMin;
                        vec[1] = y;
                        gl.glVertex3fv(vec, 0);
                        vec[2] = zMax;
                        vec[1] = y;
                        gl.glVertex3fv(vec, 0);
                    }
                    for(z = zMin + gridDZ; z < zMax; z += gridDZ)
                    {
                        vec[1] = yMin;
                        vec[2] = z;
                        gl.glVertex3fv(vec, 0);
                        vec[1] = yMax;
                        vec[2] = z;
                        gl.glVertex3fv(vec, 0);
                    }
                }
                break;
                case YDIR:
                    StsColor.GREEN.setGLColor(gl);
                    vec[1] = getCurrentDirCoordinate(YDIR);
                    gl.glBegin(GL.GL_LINES);
                {
                    for(x = xMin + gridDX; x < xMax; x += gridDX)
                    {
                        vec[2] = zMin;
                        vec[0] = x;
                        gl.glVertex3fv(vec, 0);
                        vec[2] = zMax;
                        vec[0] = x;
                        gl.glVertex3fv(vec, 0);
                    }
                    for(z = zMin + gridDZ; z < zMax; z += gridDZ)
                    {
                        vec[0] = xMin;
                        vec[2] = z;
                        gl.glVertex3fv(vec, 0);
                        vec[0] = xMax;
                        vec[2] = z;
                        gl.glVertex3fv(vec, 0);
                    }
                }
                break;
                default:
                    break;
            }
            gl.glEnd();
            glPanel3d.resetViewShift(gl);
            gl.glDisable(GL.GL_LINE_STIPPLE);
        }

        // Display the texture displayableSections on this 3d cursor plane
        cursorSections[dirNo].displayTexture(glPanel3d, true);
        // display nontexture displayableSections on this 3d cursor plane
        cursorSections[dirNo].display(glPanel3d, true);

        // Draw any map edges on section
        gl.glDisable(GL.GL_LIGHTING);
  //      if(dirNo != ZDIR)
        if (true)
        {
            int nDisplayableClasses = displayableClasses.size();
            for(int n = 0; n < nDisplayableClasses; n++)
            {
                StsClassCursorDisplayable displayableClass = (StsClassCursorDisplayable)displayableClasses.get(n);
                displayableClass.drawOnCursor3d(glPanel3d, dirNo, dirCoordinate, planePoints, thisDirDragging);
            }
        }
        gl.glEnable(GL.GL_LIGHTING);

        if(testConcavePolygon && dirNo == currentDirection)
        {
            GLU glu = glPanel3d.getGLU();
            if(dirNo == currentDirection) cursorSections[dirNo].displayTestConcavePolygon(gl, glu);
        }
        /*
           StsSurfaceClass surfaceClass = (StsSurfaceClass)model.getCreateStsClass(StsSurface.class);
           boolean displayBlockColors = false;
           surfaceClass.drawSurfaceCursorEdges(glPanel3d, planePoints[0], planePoints[3], thisDirDragging, displayBlockColors);
           StsModelSurfaceClass modelSurfaceClass = (StsModelSurfaceClass)model.getCreateStsClass(StsModelSurface.class);
           modelSurfaceClass.drawSurfaceCursorEdges(glPanel3d, planePoints[0], planePoints[3], thisDirDragging, displayBlockColors);
        */
    }

    /** Display the seismic and crossplot data on the 2d cursor view */
    public void display2d(StsGLPanel3d glPanel3d, boolean axesFlipped, boolean xAxisReversed, boolean yAxisReversed, float pixelsPerXunit, float pixelsPerYunit)
    {
        if(selectedDirection == NONE) return;
        StsCursorSection cursorSection = cursorSections[selectedDirection];
        // display texture displayableSections on this 2d cursor plane
        cursorSection.displayTexture(glPanel3d, false);
        // display nontexture displayableSections on this 2d cursor plane
        cursorSection.display(glPanel3d, false);

        float dirCoordinate = cursorSections[selectedDirection].dirCoordinate;
        StsPoint[] planePoints = cursorSections[selectedDirection].planePoints;
        displayableClasses = model.getCursorDisplayableClasses();
        int nDisplayableClasses = displayableClasses.size();
        for(int n = 0; n < nDisplayableClasses; n++)
        {
            StsClassCursorDisplayable displayableClass = (StsClassCursorDisplayable)displayableClasses.get(n);
            displayableClass.drawOnCursor2d(glPanel3d, selectedDirection, dirCoordinate, axesFlipped, xAxisReversed, yAxisReversed);
        }

        // Draw plane intersection highlight lines
        // @JKF - glPanel3d.window.isCursorTracked is never updated by the GUI so is always false.
        //        Using getShowIntersection() from 3d checkbox instead. 16May06
        drawIntersections2d(glPanel3d, axesFlipped);
        // if(glPanel3d.window.isCursorTracked) drawIntersections2d(glPanel3d.getGL(), axesFlipped);
    }

    private void drawIntersections2d(StsGLPanel3d glPanel3d, boolean axesFlipped)
    {
        GL gl = glPanel3d.getGL();

        if(model.getProject().getShowIntersection())
        {
            StsPoint[] plane;
            int[] coorIndexes;
            glPanel3d.setViewShift(gl, 2.0);

            switch(selectedDirection)
            {
                case XDIR:
                    coorIndexes = new int[]
                            {1, 2};
                    plane = cursorSections[YDIR].planePoints;
                    StsGLDraw.drawDottedLine2d(gl, StsColor.GREEN, StsColor.BLACK, 2, plane[0], plane[2], coorIndexes);
                    plane = cursorSections[ZDIR].planePoints;
                    StsGLDraw.drawDottedLine2d(gl, StsColor.BLUE, StsColor.BLACK, 2, plane[0], plane[2], coorIndexes);
                    break;
                case YDIR:
                    coorIndexes = new int[]
                            {0, 2};
                    plane = cursorSections[XDIR].planePoints;
                    StsGLDraw.drawDottedLine2d(gl, StsColor.RED, StsColor.BLACK, 2, plane[0], plane[2], coorIndexes);
                    plane = cursorSections[ZDIR].planePoints;
                    StsGLDraw.drawDottedLine2d(gl, StsColor.BLUE, StsColor.BLACK, 2, plane[0], plane[2], coorIndexes);
                    break;
                case ZDIR:
                    if(!axesFlipped)
                        coorIndexes = new int[]
                                {0, 1};
                    else
                        coorIndexes = new int[]
                                {1, 0};
                    plane = cursorSections[XDIR].planePoints;
                    StsGLDraw.drawDottedLine2d(gl, StsColor.RED, StsColor.BLACK, 2, plane[0], plane[2], coorIndexes);
                    plane = cursorSections[YDIR].planePoints;
                    StsGLDraw.drawDottedLine2d(gl, StsColor.GREEN, StsColor.BLACK, 2, plane[0], plane[2], coorIndexes);
                    break;
            }

            glPanel3d.resetViewShift(gl);
        }
    }

    /**
     * <<<< Need to complete documentation >>>>
     *
     * @param dirNo
     * @param mouse
     * @return
     */
    public StsPoint getPointInCursorPlane(StsGLPanel3d glPanel3d, int dirNo, StsMouse mouse)
    {
        if(dirNo == NONE) return null;
        StsPoint planePoint = cursorSections[dirNo].planePoints[0];
        return getPointInPlaneAtMouse(glPanel3d, dirNo, planePoint, mouse);
    }

    /**
     * <<<< Need to complete documentation >>>>
     *
     * @param planePoint
     * @param mouse
     * @return
     */
    public StsPoint getPointInPlaneAtMouse(StsGLPanel3d glPanel3d, int nDir, StsPoint planePoint, StsMouse mouse)
    {
        return getPointInPlaneAtMouse(glPanel3d, nDir, planePoint.v, mouse);
    }

    public StsPoint getPointInPlaneAtMouse(StsGLPanel3d glPanel3d, int nDir, float[] xyz, StsMouse mouse)
    {
        double[][] points = glPanel3d.getViewLineAtMouse(mouse);
        double[] pointNear = points[0];
        double[] pointFar = points[1];

        double dS = pointFar[nDir] - pointNear[nDir];
        if(Math.abs(dS) < StsParameters.roundOff) return null;
        double s = (double)xyz[nDir];
        double f = (s - pointNear[nDir]) / dS;
        StsPoint point = new StsPoint(4);
        point.interpolatePoints(3, pointNear, pointFar, f);
        point.setF((float)f);
        if(pointInsideCube(point) == true)
            return point;
        else
            return null;
    }

    public StsCursorPoint getNearestPointInCursorPlane(StsGLPanel3d glPanel3d, StsMouse mouse, StsRotatedGridBoundingBox boundingBox)
    {
        StsCursorPoint cursorPoint = getNearestPointInCursorPlane(glPanel3d, mouse);
        if(cursorPoint == null) return null;
        cursorPoint.rowNum = boundingBox.getRowCoor(cursorPoint.point.getY());
        cursorPoint.colNum = boundingBox.getColCoor(cursorPoint.point.getX());
        return cursorPoint;
    }

    /**
     * nearest point is determined in screen z coordinates which range from
     * 0 at the near plane to 1 at the far plane
     */
    public StsCursorPoint getNearestPointInCursorPlane(StsGLPanel3d glPanel3d, StsMouse mouse)
    {
        int nearestDir = -1;
        float nearestDistance = StsParameters.largeFloat;
        StsPoint nearestPoint = null;

        for(int n = 0; n < 3; n++)
        {
            StsCursorSection cursorSection = cursorSections[n];
            if(cursorSections[n].isSelected)
            {
                StsPoint point = getPointInCursorPlane(glPanel3d, n, mouse);
                if(point == null) continue;
                float f = point.v[3];
                if(f < nearestDistance)
                {
                    nearestPoint = point;
                    nearestDir = n;
                    nearestDistance = point.v[3];
                }
            }
        }
        if(nearestPoint == null) return null;
        return new StsCursorPoint(nearestDir, nearestPoint);
    }

    /**
     * Is the supplied point inside the cube?
     *
     * @param point point in relative coordinates
     * @return true if inside
     */
    public boolean pointInsideCube(StsPoint point)
    {
        return rotatedBoundingBox.isInsideXYZ(point);
    }

    /** Get the cursor minimum X */
    public float getXMin()
    {
        return rotatedBoundingBox.xMin;
    }

    /** Get the cursor maximum X */
    public float getXMax()
    {
        return rotatedBoundingBox.xMax;
    }

    /** Get the cursor minimum Y */
    public float getYMin()
    {
        return rotatedBoundingBox.yMin;
    }

    /** Get the cursor maximum Y */
    public float getYMax()
    {
        return rotatedBoundingBox.yMax;
    }

    /** Get the cursor minimum Z */
    public float getZMin()
    {
        return rotatedBoundingBox.zMin;
    }

    /** Get the cursor maximum Z */
    public float getZMax()
    {
        return rotatedBoundingBox.zMax;
    }

    public boolean cursorIntersected(StsBoundingBox boundingBox, int dirNo)
    {
        if(cursorSections == null || dirNo == -1) return false;
        return cursorSections[dirNo].cursorIntersected(boundingBox);
    }
/*
    public byte[] getSubVolumePlane()
    {
        if(currentDirection == NONE) return null;
        return cursorSections[currentDirection].getSubVolumePlane();
    }
*/
    public boolean hasSubVolumes()
    {
        if(currentDirection == NONE) return false;
        return cursorSections[currentDirection].hasSubVolumes();
    }

    public byte[] getSubVolumePlane(int dir, float dirCoordinate, byte zDomainData)
    {
        return getSubVolumePlane(dir, dirCoordinate, rotatedBoundingBox, zDomainData);
    }

    public byte[] getSubVolumePlane(int dir, float dirCoordinate, StsRotatedGridBoundingBox rotatedBoundingBox, byte zDomainData)
    {
        StsSubVolumeClass subVolumeClass = (StsSubVolumeClass)model.getStsClass(StsSubVolume.class);
        if(subVolumeClass == null) return null;
        if(!subVolumeClass.getIsVisible()) return null;
        return subVolumeClass.getSubVolumePlane(dir, dirCoordinate, rotatedBoundingBox, zDomainData);
    }

    /**
     * Get the cursor point that is closest to the user from a mouse position. A vector from the user
     * \may intersect multiple cursor sections, need to determine the closest.
     *
     * @param mouse mouse position StsMouse
     */
    public StsCursorPoint getCursorPoint(StsGLPanel3d glPanel3d, StsMouse mouse)
    {
        //    	int leftButtonState = mouse.getButtonState(StsMouse.LEFT);
        //       	if(leftButtonState != StsMouse.PRESSED && leftButtonState != StsMouse.DRAGGED) return null;
        return getNearestPointInCursorPlane(glPanel3d, mouse);
    }

    /** Output the position of the mouse in relative coordinates for the 3D cursor. */
    public void logReadout(StsGLPanel3d glPanel3d, StsCursorPoint cursorPoint)
    {
        logReadout(glPanel3d, "", cursorPoint);
    }

    public void logReadout(StsGLPanel3d glPanel3d, String message, StsCursorPoint cursorPoint)
    {
        if(cursorPoint == null) return;
        StsPoint point = cursorPoint.point;
        int dirNo = cursorPoint.dirNo;
        StsMessageFiles.infoMessage(message + buildCursorString(point, dirNo));
        String propertyString = cursorSections[dirNo].propertyReadout(glPanel3d, point);
        if(propertyString != null) StsMessageFiles.infoMessage("    " + propertyString);
    }

    public String buildCursorString(StsPoint point, int dirNo)
    {
        String valueString = cursorSections[dirNo].rowColNumReadout(point);
        double[] xy = model.getProject().getAbsoluteXYCoordinates(point);
        DecimalFormat coorFormat = model.getProject().getLabelFormat();
        return "X: " + coorFormat.format(xy[0]) + " Y: " + coorFormat.format(xy[1]) +
                " Z: " + point.v[2] + " " + valueString;
    }

    /** moves cursor position of two planes orthogonal to the point on which this cursorPoint has been picked */
    public void moveCursor3d(StsGLPanel3d glPanel3d, StsCursorPoint cursorPoint)
    {
        logReadout(glPanel3d, "", cursorPoint);
        StsPoint point = cursorPoint.point;
        int dirNo = cursorPoint.dirNo;
        StsWindowFamily windowFamily = window.getWindowFamily();
        for(int n = 0; n < 3; n++)
        {
            if(n == dirNo) continue;
            cursorSections[n].setDirCoordinate(point.v[n]);
            windowFamily.adjustCursorAndSlider(n, point.v[n]);
        }
        model.win3dDisplayAll();
    }

    /**
     * Get the cursor point at the mouse position in the 2D view.
     *
     * @param mouse mouse position StsMouse
     * @return cursor point in relative coordinates
     * @see StsCursorPoint(int, StsPoint)
     */
    public StsCursorPoint getCursorPoint2d(StsGLPanel3d glPanel3d, StsMouse mouse, boolean axesFlipped)
    {
        StsPoint point = glPanel3d.getPointInPlaneAtMouse(mouse);
        if(point == null) return null;
        return new StsCursorPoint(currentDirection, cursorSections[currentDirection].dirCoordinate, axesFlipped, point);
    }

    /**
     * Output the position of the mouse in project coordinates for the 2D cursor. The 2D cursor is
     * the one used in the 2D Cursor View.
     */
    public void logReadout2d(StsGLPanel3d glPanel3d, StsCursorPoint cursorPoint, boolean axesFlipped)
    {
        StsPoint point = cursorPoint.point;
        String valueString = cursorSections[currentDirection].rowColNumReadout(point);
        if(!getIsGridCoordinates())
            StsMessageFiles.infoMessage("X: " + point.v[0] + " Y: " + point.v[1] + " Z: " + point.v[2] + valueString);
        else
        {
            DecimalFormat coorFormat = model.getProject().getLabelFormat();
            StsMessageFiles.infoMessage(" X: " + coorFormat.format(rotatedBoundingBox.getNumFromCoor(XDIR, point.v[0]))
                    + " Y: " + coorFormat.format(rotatedBoundingBox.getNumFromCoor(YDIR, point.v[1]))
                    + " Z: " + coorFormat.format(rotatedBoundingBox.getNumFromCoor(ZDIR, point.v[2]))
                    + " " + valueString);
        }
        String propertyString = cursorSections[currentDirection].propertyReadout(glPanel3d, point);
        if(propertyString != null) StsMessageFiles.infoMessage("    " + propertyString);

    }

    StsCursor3dTexture checkAddDisplayableSection(StsClassCursor3dTextureDisplayable displayableClass, int dir)
    {
        if(displayableClass.getCurrentObject() == null)return null;
        return cursorSections[dir].checkAddDisplayableSection(model, this, displayableClass);
    }
    /**
     * Set the current cursor object for valid cursor object types (StsSeismicVolume and StsCrossplot)
     *
     * @param object the desired object
     * @return true if changed is successful
     */
    public boolean setObject(StsObject object)
    {
        if(!model.isStsClassCursor3dTextureDisplayable(object)) return false;
        boolean set = false;
        for(StsCursorSection cursorSection : cursorSections)
            set = set | cursorSection.setObject(object);
        return set;
    }
/*
    public void deleteCursor3dTextureSection(StsObject object)
    {
        if(!model.isClassCursor3dTextureDisplayable(object)) return;
        for(StsCursorSection cursorSection : cursorSections)
            cursorSection.deleteCursor3dTexture(object);
    }


    public void deleteCursor3dTextureSection(StsObject object, int dir)
    {
        if(!model.isClassCursor3dTextureDisplayable(object)) return;
        cursorSections[dir].deleteCursor3dTexture(object);
    }
*/
    /**
     * Clears both the background texture and the subImage data from the cursorSectionDisplayables.
     * Results in the both being rebuilt.
     */
    public void clearTextureDisplays()
    {
        for(StsCursorSection cursorSection : cursorSections)
            cursorSection.clearTextureDisplays();
    }

    public void clearSubVolumePlanes()
    {
        for(StsCursorSection cursorSection : cursorSections)
            cursorSection.clearSubVolumePlane();
    }

    public void clearTextureDisplays(int dir)
    {
        cursorSections[dir].clearTextureDisplays();
        cursorSections[dir].clearSubVolumePlane();
    }

    public boolean viewObjectChanged(Object source, Object object)
    {
        if(!(object instanceof StsObject)) return false;
        boolean changed = false;
        for(StsCursorSection cursorSection : cursorSections)
            changed = changed | cursorSection.viewObjectChanged(object);
        return changed;
    }

    public boolean clearTextureDisplays(Object object)
    {
        if(!(object instanceof StsObject)) return false;
        boolean cleared = false;
        if(!model.isStsClassCursor3dTextureDisplayable(object)) return false;
        for(StsCursorSection cursorSection : cursorSections)
            cleared = cleared | cursorSection.clearTextureDisplays(object);
        return cleared;
    }

    public boolean clearTextureClassDisplays(Class objectClass)
    {
        boolean cleared = false;
        for(StsCursorSection cursorSection : cursorSections)
            cleared = cleared | cursorSection.clearTextureDisplays(objectClass);
        return cleared;
    }

    public boolean isDisplayableObject(Object object)
    {
        if(!(object instanceof StsObject)) return false;
        return model.isStsClassCursor3dTextureDisplayable((StsObject)object);
    }

    public boolean isDisplayingObject(Object object)
    {
        for(StsCursorSection cursorSection : cursorSections)
            if(cursorSection.isDisplayingObject(object))
                return true;
        return false;
    }

    public boolean objectChanged(Object object)
    {
        boolean viewChanged = false;
        for(StsCursorSection cursorSection : cursorSections)
            viewChanged = viewChanged | cursorSection.objectChanged(object);
        return viewChanged;
    }

    public void deleteCursor3dTextures(StsObject object)
    {
        for(StsCursorSection cursorSection : cursorSections)
            cursorSection.deleteCursor3dTexture(object);
    }

    public boolean toggleOn(StsObject object)
    {
        boolean toggled = false;
        for(StsCursorSection cursorSection : cursorSections)
            toggled = toggled | cursorSection.toggleOn(object);
        return toggled;
    }

    public boolean toggleOff(StsObject object)
    {
        boolean toggled = false;
        for(StsCursorSection cursorSection : cursorSections)
            toggled = toggled | cursorSection.toggleOff(object);
        return toggled;
    }

    public boolean toggleOffCursor(StsObject object, int dirNo)
    {
        return cursorSections[dirNo].toggleOff(object);
    }

    public boolean toggleOnCursor(StsObject object, int dirNo)
    {
        return cursorSections[dirNo].toggleOn(object);
    }

    /** adjust the texture coordinates after a crop change. */
    public void cropChanged()
    {
        for(StsCursorSection cursorSection : cursorSections)
            cursorSection.cropChanged();
    }

    /**
     * Based on the input class type, get the current object
     *
     * @param c class type (StsSeismicVolume.class or StsCrossplot.class)
     * @return the current object of class type
     */
    public StsObject getCurrentObject(Class c)
    {
        StsClass stsClass = model.getCreateStsClass(c);
        if(stsClass == null) return null;
        return stsClass.getCurrentObject();
    }

    /*
       public boolean rowOrColChanged(int[] currentRowCol)
       {
           int currentRow = currentRowCol[0];
           int currentCol = currentRowCol[1];
           int col = getCurrentGridCoordinate(StsParameters.XDIR);
           int row = getCurrentGridCoordinate(StsParameters.YDIR);

           boolean changed = row != currentRow || col != currentCol;
           if (!changed) return false;
           System.out.println("position changed from row/col " + currentRow + "," + currentCol + " to " + row + "," + col);
           currentRowCol[0] = row;
           currentRowCol[1] = col;
           return true;
       }
    */
    static public void toggleConcavePolygonTest(ItemEvent e)
    {
        testConcavePolygon = (e.getStateChange() == ItemEvent.SELECTED);
    }

    public void resurrect(StsCursor3d o)
    {
        setDisplay3dCursor(o.display3dCursor);
        setCurrentDirNo(o.currentDirection);
        drawCursorPlanes = o.drawCursorPlanes;
    }

    /**
     * x = xyz[0]
     * y = xyz[1]
     * z = xyz[2]
     * @return
     */
    public float[] getCurrentCoordinates()
    {
        float[] xyz = new float[3];
        xyz[0] = cursorSections[XDIR].dirCoordinate;// - cursorSections[XDIR].initTotalHorizontalRange[0];
        xyz[1] = cursorSections[YDIR].dirCoordinate;// - cursorSections[YDIR].initTotalHorizontalRange[0];
        xyz[2] = cursorSections[ZDIR].dirCoordinate;
        return xyz;
    }

    public int getRow()
    {
        return row;
    }

    public int getCol()
    {
        return col;
    }

    public void setRow(int row)
    {
        this.row = row;
    }

    public void setCol(int col)
    {
        this.col = col;
    }

    public StsCursorSection getColCursorSection() { return cursorSections[XDIR]; }
    public StsCursorSection getRowCursorSection() { return cursorSections[YDIR]; }
    public StsCursorSection getSliceCursorSection() { return cursorSections[ZDIR]; }
}