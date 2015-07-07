package com.Sts.MVC.View3d;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: Class defining two-dimensional cursor view. The cursor is a reference to the
 * three planes that are displayed in the 3D view. An object based on this class wuold result in a 2D
 * view of any one of the three cursor planes.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Types.PreStack.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;
import com.Sts.Actions.*;

import javax.media.opengl.*;
import javax.swing.*;
import java.awt.event.*;

// TODO should really be subclassed into ViewCursorSeismic and viewCursorPreStack to avoid display confusion
public class StsViewCursor extends StsView2d implements StsSerializable
{
    /** compass angle displayed in legend at top of 2D view */
    public float displayAngle;
    /** legend display type at top of 2D view */
    public byte legendType = S2S_LEGEND;
    /** Retain current axis ranges for each of the 3 axis ranges */
    public float[][][] cursorAxisRanges = new float[3][][];
    /** Retain current axis ranges for each of the 3 axis ranges */
    public float[][][] cursorTotalAxisRanges = new float[3][][];
    /** Current direction being viewed (XDIR, YDIR, ZDIR) */
//    public int selectedDirNo = -1;
    /** coordinate of current direction being viewed */
    //    public float selectedDirCoordinate = 0.0f;

    /** Is labeling in Grid or relative coordinates */
//    public boolean isGridCoordinates = false;
    /** seismic volume being displayed if any */
//    private transient StsSeismicVolume seismicVolume = null;
    /** pre-stack seismic volume being displayed if any */
    //    private transient StsPreStackLineSet3d prestackLineSet3d = null;
    /** cursor3d of this window */
//    private transient StsCursor3d cursor3d;
    /** cursor3d of parent window which is controlling this view */
//    transient StsCursor3d parentCursor3d;
	private transient JPopupMenu tp = null;
    private transient StsMenuItem aspectBtn = new StsMenuItem();
    /** vertical plane at a constant x going y+ and z+ directions */
    static final int XDIR = StsCursor3d.XDIR;
    /** vertical plane at a constant y going x+ and z+ directions */
    static final int YDIR = StsCursor3d.YDIR;
    /** horizontal plane at a constant z going x+ and y+ directions */
    static final int ZDIR = StsCursor3d.ZDIR;

    static final byte CORE_LABS_LEGEND = 1;
    static final byte S2S_LEGEND = 2;

    static final String[] compassDirections = new String[]{"E", "NE", "N", "NW", "W", "SW", "S", "SE"};
    /**
     * reverse and flip XY axes so they are closest to looking in a north direction with east to the right.
     * See setViewOrientation() for details of reverse and flip depending on view quadrant (N, S, E, W)
     */
    static final boolean reorient = false;

    static public final String viewName2d = "2D View";
    static public final String shortViewName2d = "2D";

    static final long serialVersionUID = 1l;

    /** Default constructor */
    public StsViewCursor()
    {
        System.out.println("viewCursor constructor");
    }

    /**
     * StsViewCursor constructor
     *
     * @param glPanel3d the graphics context
     */
    public StsViewCursor(StsGLPanel3d glPanel3d)
    {
        super(glPanel3d);
//        cursor3d = glPanel3d.getCursor3d();
        initialize();
    }

    static public String getStaticViewName()
    {
        return viewName2d;
    }

    static public String getStaticShortViewName()
    {
        return shortViewName2d;
    }

    // classInitialize after a db custom restore
    public void initializeTransients(StsGLPanel3d glPanel3d)
    {
        super.initializeTransients(glPanel3d);
//		totalAxisRanges = StsMath.copyFloatArray(glPanel3d.cursor3d.getTotalAxisRanges());
        viewChanged();
        initializeCurrentCursor();
        // set mouse mode to zoom
//		StsMouseActionToolbar toolbar = glPanel3d.window.getMouseActionToolbar();
//		if(toolbar != null) toolbar.zoom();

    }

    public void initialize()
    {
        initializeCursor();
        aspectBtn.setMenuActionListener("1:1 Aspect Ratio", this, "aspectRatio", null);
        // set mouse mode to zoom
//        StsMouseActionToolbar toolbar = glPanel3d.window.getMouseActionToolbar();
//        if(toolbar != null) toolbar.zoom();
    }

    private void initializeCursor()
    {
        initializeAxisRanges();
        setLabels();
        getCursor3d().setAxisLabels();

        setViewOrientation();
        setLabels();
    }

    private void initializeCurrentCursor()
    {
        getCursor3d().setAxisLabels();
        setViewOrientation();
        setLabels();
    }

    public void initializeAxisRanges()
    {
        for(int n = 0; n < 3; n++)
        {
            cursorAxisRanges[n] = model.getProject().getRotatedBoundingBox().getCursorAxisRange(n);
            cursorTotalAxisRanges[n] = StsMath.copyFloatArray(cursorAxisRanges[n]);
        }
        StsCursor3d cursor3d = glPanel3d.getCursor3d();
        int selectedDirection = cursor3d.selectedDirection;
        axisRanges = cursorAxisRanges[selectedDirection];
        totalAxisRanges = cursorTotalAxisRanges[selectedDirection];
        glPanel3d.viewChanged = true;
    }

    public void computeProjectionMatrix()
    {
        if(axisRanges == null) return;
        gl.glMatrixMode(GL.GL_PROJECTION);
//        glPanel3d.debugPrintViewport("StsViewCursor.computeProjectionMatrix");
        gl.glLoadIdentity();
        glu.gluOrtho2D(axisRanges[0][0], axisRanges[0][1], axisRanges[1][0], axisRanges[1][1]);
        gl.glGetDoublev(GL.GL_PROJECTION_MATRIX, glPanel3d.projectionMatrix, 0);
        gl.glMatrixMode(GL.GL_MODELVIEW);
    }
    /**
     * Set the parameters of the current view based on another viewCursor object
     * @param otherView the other viewCursor object
     */
/*
    public void setParameters(StsViewCursor otherView)
    {
        setRangesFromOtherView(otherView);
        axisRanges = StsMath.copyFloatArray(otherView.axisRanges);
        totalAxisRanges = StsMath.copyFloatArray(otherView.totalAxisRanges);
        isGridCoordinates = otherView.isGridCoordinates;
        setViewOrientation();
    }
*/
    /**
     * The total axis ranges regardless of zoom
     * @return X and Y data ranges [2][2] (xMin,xMax,yMin,yMax)
     */
    //    public float[][] getTotalAxisRanges() { return cursor3d.getTotalAxisRanges(); }
    /**
     * Is the cursor viewable
     *
     * @return true if viewable
     */
    public boolean isViewable()
    {
        return getCursor3d().selectedDirection != StsParameters.NO_MATCH;
    }

    public void viewChanged()
    {
        if(glPanel3d == null) return;
        glPanel3d.viewChanged = true;
    }

    /**
     * Display the current cursor plane view
     *
     * @param component the drawable component
     * @see StsSeismicVolume#getSeismicVolumeClass()
     * @see StsSeismicVolumeClass#getCurrentSeismicVolume()
     * @see #drawHorizontalAxis
     * @see #drawVerticalAxis
     */
    public void display(GLAutoDrawable component)
    {
          if(glPanel.panelViewChanged)
          {
              initializeView();
              glPanel.panelViewChanged = false;
          }
        // viewPort may be inset from last display call (needed for picking)
        resetViewPort();
        // glPanel3d.viewPortChanged = false;
        /*
        if(glPanel3d.viewPortChanged)
        {
            glPanel3d.setViewPort();
            glPanel3d.viewPortChanged = false;
        }
        */
        int selectedDirection = getCursor3d().selectedDirection;

        // no current cursor: display black screen
        if(selectedDirection == StsParameters.NO_MATCH)
        {
            clearToBackground(GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);
            //gl.glDrawBuffer(GL.GL_BACK);
            //glPanel3d.applyClearColor();
            //gl.glClear(GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);
            gl.glDisable(GL.GL_LIGHTING);
            return;
        }

        if(glPanel3d.viewChanged)
        {
            setAxisRanges();
            setLabels();
            computeProjectionMatrix();
            computeModelViewMatrix();
            glPanel3d.viewChanged = false;
        }

        // Display axis ranges; different from axis ranges if coordinates are different
        boolean isGridCoordinates = getCursor3d().getIsGridCoordinates();

        float[][] axisDisplayRanges;
        if(isGridCoordinates)
            axisDisplayRanges = model.getProject().getRotatedBoundingBox().getGridCoordinateRanges(axisRanges, selectedDirection, axesFlipped);
        else
            axisDisplayRanges = resetDisplayRanges(axisRanges, totalAxisRanges, selectedDirection);
        try
        {
            boolean displaySeismicAxis = true;
            setInsets(displaySeismicAxis);
            if(displaySeismicAxis) insetViewPort();

            if(cursorButtonState != StsMouse.CLEARED)
            {
                // if this is window where cursor is being dragged and we have focus, draw foreground cursor.
                // If not the window where cursor is being dragged, but we are displaying cursor here,
                // draw the windows;
                if(isCursorWindow && glPanel3d.hasFocus() || !isCursorWindow)
                {
                    drawForeground();
                    if(cursorButtonState != StsMouse.CLEARED) return;
                }
            }
            clearToBackground(GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);
            /*
            gl.glDrawBuffer(GL.GL_BACK);
            glPanel3d.applyClearColor();
            gl.glClear(GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);
            */
            // 2d display: objects drawn last are displayed on top of others
            gl.glDisable(GL.GL_DEPTH_TEST);
//            gl.glDepthFunc(GL.GL_ALWAYS);

//            if(displayAxes) resetViewPort();

            if(displaySeismicAxis)
            {
                resetViewPort();
                gl.glMatrixMode(GL.GL_PROJECTION);
                gl.glPushMatrix();
                gl.glLoadIdentity();
                // StsException.systemDebug(this, "display", "gluOrtho2d width: " + glPanel3d.getWidth() + " height: " + glPanel3d.getHeight());
                glu.gluOrtho2D(0, glPanel3d.getWidth(), 0, glPanel3d.getHeight());
                gl.glMatrixMode(GL.GL_MODELVIEW);

                String titleLabel, label;

                StsRotatedGridBoundingBox boundingBox = model.getProject().getRotatedBoundingBox();
                labelFormat = model.getProject().getLabelFormat();

                drawVerticalAxis(axisLabels[1], axisDisplayRanges[1], model.getProject().getShow2dGrid());

                float displayedDirCoordinate = getCursor3d().getCurrentDirCoordinate(selectedDirection);
                if(isGridCoordinates)
                    label = getLabelFormat().format(boundingBox.getNearestNumFromCoor(selectedDirection, displayedDirCoordinate));
                else
                    label = getLabelFormat().format(displayedDirCoordinate);
                titleLabel = getCursor3d().getTitleLabel(selectedDirection, label);
                drawHorizontalAxis(titleLabel, axisLabels[0], axisDisplayRanges[0], model.getProject().getShow2dGrid());
                drawTopLegend(gl);

                gl.glMatrixMode(GL.GL_PROJECTION);
                gl.glPopMatrix();
                gl.glMatrixMode(GL.GL_MODELVIEW);

                insetViewPort(); // sets up viewport for picking
                // glPanel3d.viewPortChanged = true;
            }

            // draw non-texture objects on cursor view
            getCursor3d().display2d(glPanel3d, axesFlipped, xAxisReversed, yAxisReversed, pixelsPerXunit, pixelsPerYunit);

            gl.glEnable(GL.GL_DEPTH_TEST);
            gl.glEnable(GL.GL_LIGHTING);
        }
        catch(Exception e)
        {
            StsException.outputException("StsViewCursor.display() failed.",
                    e, StsException.FATAL);
        }
    }
/*
    public void clearToBackground()
    {
          glPanel3d.setClearColor(model.getProject().getBackgroundStsColor());
          gl.glDrawBuffer(GL.GL_BACK);
          glPanel3d.applyClearColor();
          gl.glClear(GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);
    }
*/
    protected void initializeView()
    {
        gl.glShadeModel(GL.GL_FLAT);
        gl.glEnable(GL.GL_ALPHA_TEST);
        gl.glAlphaFunc(GL.GL_GREATER, 0.1f);
        initializeFont(gl);
        initializeInsets();
        computePixelScaling();
        computeProjectionMatrix();
        computeModelViewMatrix();
    }

   /** Set the aspect ratio to 1:1. */
   public void aspectRatio()
   {
	   if(getCursor3d().selectedDirection != getCursor3d().ZDIR) return;

       float halfXRange = Math.abs(axisRanges[0][1] - axisRanges[0][0])/2.0f;
       float yPos = axisRanges[1][0] + Math.abs(axisRanges[1][1] - axisRanges[1][0])/2.0f;
       axisRanges[1][0] = yPos - halfXRange;
       axisRanges[1][1] = yPos + halfXRange;
	   computePixelScaling();
       viewChangedRepaint();
   }
    /**
     * Display the popup menu)
     *
     * @param mouse
     */
    public void showPopupMenu(StsMouse mouse)
    {
        // chec if well was picked
        JPopupMenu tp = new JPopupMenu("2D View Popup");
        glPanel3d.add(tp);

        tp.add(aspectBtn);
        tp.show(glPanel3d, mouse.getX(), mouse.getY());
        clearCurrentKeyPressed();
    }
	public void cancelPopupMenu(StsMouse mouse)
	{
		   if (tp != null) tp.setVisible(false);
    }
    private float[][] resetDisplayRanges(float[][] axisRanges, float[][] totalAxisRanges, int selectedDirection)
    {
        float[][] axisDisplayRanges = new float[2][2];
        float xMin = totalAxisRanges[0][0];
        axisDisplayRanges[0][0] = axisRanges[0][0] - xMin;
        axisDisplayRanges[0][1] = axisRanges[0][1] - xMin;
        if(selectedDirection != ZDIR)
        {
            axisDisplayRanges[1][0] = axisRanges[1][0];
            axisDisplayRanges[1][1] = axisRanges[1][1];
        }
        else
        {
            float yMin = totalAxisRanges[1][0];
            axisDisplayRanges[1][0] = axisRanges[1][0] - yMin;
            axisDisplayRanges[1][1] = axisRanges[1][1] - yMin;
        }
        return axisDisplayRanges;
    }
/*
    public void displayForeground()
    {
       	if(cursorButtonState == StsMouse.PRESSED)
        {
//            System.out.println("Pressed");
            previousCursorPoint = null;
            currentCursorPoint = cursorPoint;
            drawForeground();
        }
        else if(leftButtonState == StsMouse.DRAGGED)
        {
//            System.out.println("Dragged");
            currentCursorPoint = cursorPoint;
            drawForeground();
        }
        else // StsMouse.RELEASED
        {
//            System.out.println("Released");
            drawForeground();
            cursorButtonState = StsMouse.CLEARED;
        }
    }
*/

    /**
     * Draws in foreground XORed against current view.  Must be called
     * again and draw same objects to erase.
     */
    public void drawForeground()
    {
        gl.glDisable(GL.GL_BLEND);
        gl.glDrawBuffer(GL.GL_FRONT);
        gl.glEnable(GL.GL_COLOR_LOGIC_OP);
        gl.glLogicOp(GL.GL_XOR);
        gl.glDepthMask(false);
        gl.glDepthFunc(GL.GL_ALWAYS);
        drawForegroundCursor();
        gl.glFlush();
        if(cursorButtonState == StsMouse.RELEASED)
        {
            gl.glLogicOp(GL.GL_COPY);
            gl.glDisable(GL.GL_COLOR_LOGIC_OP);
            gl.glDepthMask(true);
            gl.glDepthFunc(GL.GL_LESS);
            gl.glDrawBuffer(GL.GL_BACK);
            cursorButtonState = StsMouse.CLEARED;
        }
    }

    /**
     * Draw the cursor in the front buffer.  If previously drawn,
     * draw previous one again to erase it; then drawn new one.
     */
    private void drawForegroundCursor()
    {
        if(previousXYZ != null)
            drawCursorPoint(previousXYZ, this.glPanel3d);
        if(currentXYZ != null)
            drawCursorPoint(currentXYZ, glPanel3d);
        previousXYZ = currentXYZ;
        currentXYZ = null;
    }

    /**
     * Draw vertical and horizontal lines thru the cursor point.
     * because the crossed area is drawn twice it will be clear.
     */
    private void drawCursorPoint(float[] xyz, StsGLPanel3d glPanel3d)
    {
        switch(getCursor3d().selectedDirection)
        {
            case XDIR:
                StsGLDraw.drawLine2d(xyz[1], xyz[2], StsColor.GRAY, glPanel3d, 3, 31, StsParameters.NORTH, 0.0, -15.0);
                StsGLDraw.drawLine2d(xyz[1], xyz[2], StsColor.GRAY, glPanel3d, 3, 31, StsParameters.EAST, -15.0, 0.0);
                break;
            case YDIR:
                StsGLDraw.drawLine2d(xyz[0], xyz[2], StsColor.GRAY, glPanel3d, 3, 31, StsParameters.NORTH, 0.0, -15.0);
                StsGLDraw.drawLine2d(xyz[0], xyz[2], StsColor.GRAY, glPanel3d, 3, 31, StsParameters.EAST, -15.0, 0.0);
                break;
            case ZDIR:
                if(!axesFlipped)
                {
                    StsGLDraw.drawLine2d(xyz[0], xyz[1], StsColor.GRAY, glPanel3d, 3, 31, StsParameters.NORTH, 0.0, -15.0);
                    StsGLDraw.drawLine2d(xyz[0], xyz[1], StsColor.GRAY, glPanel3d, 3, 31, StsParameters.EAST, -15.0, 0.0);
                }
                else
                {
                    StsGLDraw.drawLine2d(xyz[1], xyz[0], StsColor.GRAY, glPanel3d, 3, 31, StsParameters.NORTH, 0.0, -15.0);
                    StsGLDraw.drawLine2d(xyz[1], xyz[0], StsColor.GRAY, glPanel3d, 3, 31, StsParameters.EAST, -15.0, 0.0);
                }
        }
    }
/*
    public void drawForeground()
    {
//        System.out.println("Drawing foreground. XOR: " + isXOR);
        gl.glDrawBuffer(GL.GL_FRONT);
        gl.glEnable(GL.GL_COLOR_LOGIC_OP);
        gl.glLogicOp(GL.GL_XOR);
        gl.glDepthMask((byte)GL.GL_FALSE);
        drawForegroundCursor();
        if(cursorButtonState == StsMouse.RELEASED)
        {
            gl.glLogicOp(GL.GL_COPY);
            gl.glDisable(GL.GL_COLOR_LOGIC_OP);
            gl.glDepthMask((byte)GL.GL_TRUE);
            gl.glDrawBuffer(GL.GL_BACK);
            cursorButtonState = StsMouse.CLEARED;
        }
    }

    private void drawForegroundCursor()
    {
        if(previousXYZ != null)
            drawCursorPoint2d(previousXYZ, glPanel3d);
        if(currentXYZ != null)
            drawCursorPoint2d(currentXYZ, glPanel3d);
        previousXYZ = currentXYZ;
        currentXYZ = null;
    }

    private void drawCursorPoint2d(float[] xyz, StsGLPanel3d glPanel3d)
    {
        switch(currentDirNo)
        {
            case XDIR:
//                StsGLDraw.drawPoint2d(xyz[1], xyz[2], StsColor.BLUE, glPanel3d, 10, 0.0);
                StsGLDraw.drawPoint2d(xyz[1], xyz[2], StsColor.GRAY, glPanel3d, 10, 0.0);
                break;
            case YDIR:
//                StsGLDraw.drawPoint2d(xyz[0], xyz[2], StsColor.BLUE, glPanel3d, 10, 0.0);
                StsGLDraw.drawPoint2d(xyz[0], xyz[2], StsColor.GRAY, glPanel3d, 10, 0.0);
                break;
            case ZDIR:
//                StsGLDraw.drawPoint2d(xyz[0], xyz[1], StsColor.BLUE, glPanel3d, 10, 0.0);
                if(!axesFlipped)
                    StsGLDraw.drawPoint2d(xyz[0], xyz[1], StsColor.GRAY, glPanel3d, 10, 0.0);
                else
                    StsGLDraw.drawPoint2d(xyz[1], xyz[0], StsColor.GRAY, glPanel3d, 10, 0.0);
        }
    }
*/

    /**
     * For Corelabs legend type, display the line direction for a vertical section at the top;
     * and the view direction for a horizontal section.  Line direction is the
     * compass direction to the left and right.  View direction is the compass direction
     * in the direction from bottom to top of the view.  Compass direction is one
     * of the 8 compass points.
     * <p>For S2S legend type, display the view direction and angle for both horizontal
     * and vertical sections.
     */
    public void drawTopLegend(GL gl)
    {
        int x, y;

        if(legendType == CORE_LABS_LEGEND)
        {
            String[] legends = getClbView2dLegend();
            y = glPanel3d.getHeight() - insets.top / 2 + halfWidth;
            if(legends.length == 2)
            {
                x = insets.left;
                StsGLDraw.fontOutput(gl, x, y, legends[0], horizontalFont);
                x = glPanel3d.getWidth() - insets.right - StsGLDraw.getFontStringLength(horizontalFont, legends[1]);
                StsGLDraw.fontOutput(gl, x, y, legends[1], horizontalFont);
            }
            else if(legends.length == 1)
            {
                x = glPanel3d.getWidth() - insets.right - StsGLDraw.getFontStringLength(horizontalFont, legends[0]) - 150;
                StsGLDraw.fontOutput(gl, x, y, legends[0], horizontalFont);
            }
        }
        else // S2S_LEGEND
        {
            String legend = getStsView2dLegend();
            int axisMaxX = glPanel3d.getWidth() - insets.right;
            x = axisMaxX - StsGLDraw.getFontStringLength(horizontalFont, legend);
            y = glPanel3d.getHeight() - insets.top / 2 + halfWidth;
            StsGLDraw.fontOutput(gl, x, y, legend, horizontalFont);
        }
    }

    /** if parent (main window) has a different selectedDirection, then set the corresponding axis ranges */
    public void setAxisRanges()
    {
        int selectedDirection = getCursor3d().selectedDirection;
//        if(parentCursor3d.selectedDirection == selectedDirection) return;
        axisRanges = cursorAxisRanges[selectedDirection];
        totalAxisRanges = cursorTotalAxisRanges[selectedDirection];
		// keep axisRanges within totalAxisRanges
		for(int i = 0; i < 2; i++)
			for(int j = 0; j < 2; j++)
				axisRanges[i][j] = StsMath.limitBetweenInclusive(axisRanges[i][j], totalAxisRanges[i][0], totalAxisRanges[i][1]);
//        glPanel3d.viewChanged = true;
//        boolean isGridCoordinates = cursor3d.isGridCoordinates();
    }
/*
    private void setRangesFromOtherView(StsViewCursor otherView)
    {
        axisRanges = StsMath.copyFloatArray(otherView.axisRanges);
        totalAxisRanges = StsMath.copyFloatArray(otherView.totalAxisRanges);
//        setDisplayRanges(otherView.isGridCoordinates);
        setViewOrientation();
    }
*/

    private void setLabels()
    {
        StsCursor3d cursor3d = getCursor3d();

        int selectedDirection = cursor3d.selectedDirection;
        if(!axesFlipped)
        {
            axisLabels[0] = cursor3d.getHorizontalAxisLabel(selectedDirection);
            axisLabels[1] = cursor3d.getVerticalAxisLabel(selectedDirection);
        }
        else
        {
            axisLabels[1] = cursor3d.getHorizontalAxisLabel(selectedDirection);
            axisLabels[0] = cursor3d.getVerticalAxisLabel(selectedDirection);
        }
    }

    /**
     * Defines in which quadrant the +Inline direction is which determines the layout of the window.
     * possible orientations and angles: 0: East (315-45), 1: North (45-135),
     * 2: West (135-225), 3: South (225-315).  Angles are measured counterclockwise from
     * the +X axis.  The displayAngle is compass angle measured clockwise from North.
     */
    private void setViewOrientation()
    {
        float angle = model.getProject().getAngle();
        displayAngle = angle;
        if(reorient)
        {
            xAxisReversed = false; // range on project x axis is reversed
            yAxisReversed = false; // range on project y axis is reversed

            int orientation = (int)((angle + 45.0f) / 90.0f);
            orientation = orientation % 4;

            switch(orientation)
            {
                case 0: // East
                    axesFlipped = false;
                    xAxisReversed = false;
                    yAxisReversed = false;
                    break;
                case 1: // North
                    axesFlipped = true;
                    xAxisReversed = false;
                    yAxisReversed = true;
                    break;
                case 2: // West
                    axesFlipped = false;
                    xAxisReversed = true;
                    yAxisReversed = true;
                    break;
                case 3: // South
                    axesFlipped = true;
                    xAxisReversed = true;
                    yAxisReversed = false;
                    break;
            }

            switch(getCursor3d().selectedDirection)
            {
                case XDIR: // vertical section defined by x value in y direction

                    if(yAxisReversed)
                        flipAxisDirection(0);
                    else
                        displayAngle = (displayAngle + 180f) % 360f;
                    break;
                case YDIR: // vertical section defined by y value in x direction
                    if(xAxisReversed)
                    {
                        flipAxisDirection(0);
                        displayAngle = (displayAngle + 270f) % 360f;
                    }
                    else
                        displayAngle = (displayAngle + 90f) % 360f;

                    break;
                case ZDIR: // horizontal section
                    if(yAxisReversed) flipAxisDirection(1);
                    if(xAxisReversed) flipAxisDirection(0);
                    if(axesFlipped)
                    {
                        flipAxes(totalAxisRanges);
                        flipAxes(axisRanges);
                        //                    flipAxes(axisDisplayRanges);
                        flipLabels(axisLabels);
                    }
                    displayAngle = (displayAngle + 90f - orientation * 90f) % 360f;
            }
        }
        else
        {
            switch(getCursor3d().selectedDirection)
            {
                case XDIR: // vertical section defined by x value in y direction
                    displayAngle = (displayAngle + 180f) % 360f;
                    break;
                case YDIR: // vertical section defined by y value in x direction
                    displayAngle = (displayAngle + 90f) % 360f;

                    break;
                case ZDIR: // horizontal section
                    displayAngle = (displayAngle + 90f) % 360f;
            }
        }
    }

    private String getStsView2dLegend()
    {
        int compassAngle = Math.round((450f - displayAngle) % 360);
        String compassDirection = get8PointCompassDirection(displayAngle);
        return new String(compassDirection + " (" + compassAngle + ")");
    }

    private String[] getClbView2dLegend()
    {
        String[] directionInfo;

        if(getCursor3d().selectedDirection == ZDIR)
        {
            directionInfo = new String[1];
            directionInfo[0] = get8PointCompassDirection(displayAngle);
        }
        else
        {
            directionInfo = new String[2];
            float leftAngle = (displayAngle + 90f) % 360f;
            directionInfo[0] = get8PointCompassDirection(leftAngle);
            float riteAngle = (displayAngle + 270f) % 360f;
            directionInfo[1] = get8PointCompassDirection(riteAngle);
        }
        return directionInfo;
    }

    private String get8PointCompassDirection(float angle)
    {
        if(angle < 0f)
        {
            StsException.systemError("Angle passed to get9PointCompassDirection() is < 0");
            return compassDirections[0];
        }
        int octant = ((Math.round(angle + 22.5f)) % 360) / 45;
        return compassDirections[octant];
    }

    private void flipAxisDirection(int index)
    {
        flipAxisDirection(totalAxisRanges[index]);
        flipAxisDirection(axisRanges[index]);
//        flipAxisDirection(axisDisplayRanges[index]);
    }

    private void flipAxisDirection(float[] axisRanges)
    {
        float temp = axisRanges[0];
        axisRanges[0] = axisRanges[1];
        axisRanges[1] = temp;
    }

    private void flipAxes(float[][] axisRanges)
    {
        float[] temp = axisRanges[0];
        axisRanges[0] = axisRanges[1];
        axisRanges[1] = temp;
    }

    private void flipLabels(String[] axisLabels)
    {
        String temp = axisLabels[0];
        axisLabels[0] = axisLabels[1];
        axisLabels[1] = temp;
    }

    public float getGridInterval()
    {
        if(glPanel3d.getCurrentDirNo()  != StsCursor3d.ZDIR)
	        return model.getProject().getZGridFrequency();
        else
            return model.getProject().getGridFrequency();
    }



    /*
        private float[] getVerticalDisplayRange(float[][]axisDisplayRanges, int currentDirNo)
        {
            if(currentDirNo != ZDIR) return axisDisplayRanges[2];
            else                     return axisDisplayRanges[0];
        }


        private float[] getHorizontalDisplayRange(float[][]axisDisplayRanges, int currentDirNo)
        {
            if(currentDirNo != YDIR) return axisDisplayRanges[1];
            else                     return axisDisplayRanges[0];
        }
    */
    /**
     * Key Release event handling
     *
     *
     *
     * @param mouse mouse object
     * @param e     key event
     * @return true if successful
     */
    public void keyReleased(StsMouse mouse, KeyEvent e)
    {
		int releasedKeyCode = e.getKeyCode();
        super.keyReleased(mouse, e);
        if(mouse.isButtonDown(StsMouse.VIEW))
        {

            glPanel3d.restoreCursor();
            return;
        }

        if(releasedKeyCode == KeyEvent.VK_LEFT)
            glPanel3d.window.displayPreviousObject(StsSeismicVolume.class);
        else if(releasedKeyCode == KeyEvent.VK_RIGHT)
            glPanel3d.window.displayNextObject(StsSeismicVolume.class);
        else if(releasedKeyCode == KeyEvent.VK_UP)
        {
            if(glPanel3d.checkAddView(StsViewXP.class) == null)
                glPanel3d.checkAddView(StsView3d.class);
        }
        else if(releasedKeyCode == KeyEvent.VK_DOWN)
        {
            glPanel3d.checkAddView(StsView3d.class);
        }
        else
            return;
        viewChangedRepaint();
    }

    /**
     * Set or Reset the view to the default view
     *
     * @see StsViewCursor#computePixelScaling()
     * @see StsCursor3d#resetAxisRanges()
     */
    //TODO fix this so if we haven't changed project bounds, we use the same position alorithm
    public void setDefaultView()
    {
        initializeCursor();
        setLabels();
        setViewOrientation();
        computePixelScaling();
    }

    /**
     * Clear the crossplot overlay from the current display and force a redraw.
     */
    /*
        public void clearCrossplotDisplay()
        {
            if(currentDirNo == -1) return;
            cursor3d.cursorSectionDisplayables[currentDirNo].clearTextureDisplay();
        }
    */

    /**
     * Clear the current display and force a redraw. Clear the textures as well.
     * @param andCursors specify whether to clear cursors as well
     */
    /*
        public void clearDisplay(boolean andCursors)
        {
           if(andCursors)
               cursor3d.cursorSectionDisplayables[currentDirNo].setDeleteTextures();
           clearDisplay();
        }
    */
    /**
     * Compute the mouse pick point which will be applied to the cursor and crossplot views
     *
     * @param mouse the mouse object
     * @return the picked point
     */
    public StsPoint computePickPoint3d(StsMouse mouse)
    {
        StsCursor3d cursor3d = getCursor3d();
        float planeCoordinate = cursor3d.getCurrentDirCoordinate(cursor3d.selectedDirection);

        StsPoint point2d = computePickPoint(mouse);
        switch(cursor3d.selectedDirection)
        {
            case XDIR:
                return new StsPoint(planeCoordinate, point2d.v[0], point2d.v[1]);
            case YDIR:
                return new StsPoint(point2d.v[0], planeCoordinate, point2d.v[1]);
            case ZDIR:
                if(!axesFlipped)
                    return new StsPoint(point2d.v[0], point2d.v[1], planeCoordinate);
                else
                    return new StsPoint(point2d.v[1], point2d.v[0], planeCoordinate);
            default:
                return null;
        }
    }

    /**
     * Output the mouse tracking readout to the information panel on the main screen
     *
     * @param mouse mouse object
     */
    public void cursor3dLogReadout2d(StsMouse mouse)
    {
        StsCursorPoint cursorPoint = getCursor3d().getCursorPoint2d(glPanel3d, mouse, axesFlipped);
        if(cursorPoint == null) return;
        setCursorXOR(glPanel3d, mouse, cursorPoint.point.v);
        if(mouse.getButtonStateCheckClear(StsMouse.LEFT) == StsMouse.RELEASED)
        {
        	// don't move with readout if user tracking is off.
     		StsWin3dBase window = glPanel3d.window;
    		StsWin3dBase parentWindow = window.parentWindow;
    		if(parentWindow == null)
    			parentWindow = window;
    		if(parentWindow != null && parentWindow.isLocked)
    			moveCursor3d(mouse, glPanel3d); // Still want to reposition cursor even if mouse tracking
            getCursor3d().logReadout2d(glPanel3d, cursorPoint, axesFlipped);
        }
    }

    public StsCursorPoint getCursorPoint(StsMouse mouse)
    {
        return getCursor3d().getCursorPoint2d(glPanel3d, mouse, axesFlipped);
    }

    /**
     * Move the other two cursor planes to the point selected on this cursor plane
     *
     * @param mouse mouse object
     */
    public boolean moveCursor3d(StsMouse mouse, StsGLPanel3d glPanel3d)
    {
        if(!mouse.isButtonStateReleased(StsMouse.LEFT)) return false;
        StsCursorPoint cursorPoint = getCursor3d().getNearestPointInCursorPlane(glPanel3d, mouse);
        StsPoint p = computePickPoint3d(mouse);
        if(p != null)
            cursorPoint = new StsCursorPoint(getCursor3d().selectedDirection, p);
        if(cursorPoint == null) return false;

        int dirNo = getCursor3d().getCurrentDirNo();
        getCursor3d().moveCursor3d(glPanel3d, cursorPoint);
        getCursor3d().setCurrentDirNo(dirNo);
        return true;
    }

    public StsCursorPoint logReadout(StsGLPanel3d glPanel3d, StsMouse mouse)
    {
        StsCursorPoint cursorPoint = getCursor3d().getNearestPointInCursorPlane(glPanel3d, mouse);
        StsPoint p = computePickPoint3d(mouse);
        if(p != null)
            cursorPoint = new StsCursorPoint(getCursor3d().selectedDirection, p);
        if(cursorPoint == null) return null;

        setCursorXOR(glPanel3d, mouse, p.v);
        glPanel3d.getCursor3d().logReadout(glPanel3d, cursorPoint);
		return cursorPoint;
        /*
        if (mouse.getButtonStateCheckClear(StsMouse.LEFT) == StsMouse.RELEASED)
        {
            moveCursor3d(mouse, glPanel3d); // Still want to reposition cursor even if mouse tracking
            glPanel3d.getCursor3d().logReadout(glPanel3d, cursorPoint);
        }
        */
    }
    /*
      * custom serialization requires versioning to prevent old persisted files from barfing.
      * if you add/change fields, you need to bump the serialVersionUID and fix the
      * reader to handle both old & new
      */
    public boolean viewObjectChanged(Object source, Object object)
    {
        if(model.isStsClassCursor3dTextureDisplayable(object))
        {
            // model.clearDisplayTextured3dCursors(object);
            if(glPanel3d.window.getCursor3d().viewObjectChanged(source, object))
                return true;
        }
       // if(getCursor3d().objectChanged(object)) return true;
       if(object instanceof StsMultiAttributeVector)
       {
           // ?? Call repaint using viewObjectRepaint
//           glPanel3d.repaint();
           return true;
       }
       return false;
    }

    public void setDefaultView(byte axis)
    {
//	   retainZoomLevel(currentDirNo);
        super.setDefaultView();
    }

    public void adjustCursor(int dirNo, float dirCoor)
    {
        if(glPanel3d == null) return;
        if(glPanel3d.window.isLocked) centerViewOnCursors();
    }

    protected void centerViewOnCursors()
    {
        boolean canRecenter = canRecenter(axisRanges[0], totalAxisRanges[0]) ||
                              canRecenter(axisRanges[1], totalAxisRanges[1]);
        if(!canRecenter) return;

        float xAxisPos = 0.0f;
        float yAxisPos = 0.0f;
        // Range is the same, but min and max change.
        switch(getCursor3d().selectedDirection)
        {
            case StsCursor3d.XDIR:
                xAxisPos = getCursor3d().getCurrentDirCoordinate(StsCursor3d.YDIR);
                yAxisPos = getCursor3d().getCurrentDirCoordinate(StsCursor3d.ZDIR);
                break;
            case StsCursor3d.YDIR:
                xAxisPos = getCursor3d().getCurrentDirCoordinate(StsCursor3d.XDIR);
                yAxisPos = getCursor3d().getCurrentDirCoordinate(StsCursor3d.ZDIR);
                break;
            case StsCursor3d.ZDIR:
                xAxisPos = getCursor3d().getCurrentDirCoordinate(StsCursor3d.XDIR);
                yAxisPos = getCursor3d().getCurrentDirCoordinate(StsCursor3d.YDIR);
                break;
        }

        if(axesFlipped)
        {
            float temp = xAxisPos;
            xAxisPos = yAxisPos;
            yAxisPos = temp;
        }

        // center on screen x and y position specified by two other cursor planes,
        // but don't allow the range outside of the total range limits

        centerViewOnCursor(xAxisPos, axisRanges[0], totalAxisRanges[0]);
        centerViewOnCursor(yAxisPos, axisRanges[1], totalAxisRanges[1]);
        computePixelScaling();
        viewChangedRepaint();
    }

    private boolean canRecenter(float[] axisRange, float[] totalAxisRange)
    {
        float range = Math.abs(axisRange[1] - axisRange[0]);
        float totalRange = Math.abs(totalAxisRange[1] - totalAxisRange[0]);
        return totalRange <= range;
    }

    private void centerViewOnCursor(float axisPos, float[] axisRanges, float[] totalAxisRanges)
    {
        float range = axisRanges[1] - axisRanges[0];
        float min = axisPos - range / 2;
        float max = axisPos + range / 2;
        float minLimit = totalAxisRanges[0];
        float maxLimit = totalAxisRanges[1];
        boolean minBetween = StsMath.betweenInclusive(min, minLimit, maxLimit);
        boolean maxBetween = StsMath.betweenInclusive(max, minLimit, maxLimit);
        if(minBetween && maxBetween)
        {
            axisRanges[0] = min;
            axisRanges[1] = max;
        }
        else if(!minBetween && !maxBetween)
        {
            axisRanges[0] = minLimit;
            axisRanges[1] = maxLimit;
        }
        else if(!minBetween)
        {
            axisRanges[0] = minLimit;
            axisRanges[1] = minLimit + range;
        }
        else if(!maxBetween)
        {
            axisRanges[1] = maxLimit;
            axisRanges[0] = maxLimit - range;
        }
    }

    /*
       public void recomputeZoomedScaling(float[][] oldRanges, byte axis)
        {
            if(prestackLineSet3d == null) return;

            float inches, axisRange;
            if((axis == BOTH) || (axis == VERTICAL))
            {
                axisRange = (axisRanges[1][0] - axisRanges[1][1])/1000.0f;  //seconds
                inches = (float)glPanel3d.glc.getHeight() / (float)getPixelsPerInch();  // inches
                prestackLineSet3d.getWiggleDisplayProperties().setInchesPerSecond(inches/axisRange);
            }

            if((axis == BOTH) || (axis == HORIZONTAL))
            {
                float tracesPerInch = prestackLineSet3d.getWiggleDisplayProperties().getTracesPerInch();
                inches = (float)glPanel3d.glc.getWidth() / (float)getPixelsPerInch();  // inches
                float oldRange = (oldRanges[0][1] - oldRanges[0][0]);
                float unitsPerTrace = (oldRange/inches)/tracesPerInch;
                axisRange = (axisRanges[0][1] - axisRanges[0][0]);  // distance
                float nTraces = axisRange/unitsPerTrace;
                prestackLineSet3d.getWiggleDisplayProperties().setTracesPerInch(nTraces/inches);
            }
    //        retainZoomLevel(currentDirNo);
            // Not optimal but without it a zoom operation will not be persisted......SAJ
            // scale parameters now moved to StsClass level and are persisted with propertiesPersistManager. TJL 3/5/07
    //        lineSet.getWiggleDisplayProperties().saveState();
    //        lineSet.getWiggleDisplayProperties().commitChanges();
        }
    */
    /* Calculate the horizontal axis limits based on user input
    */
    public void resetToOrigin()
    {
        float xRange = axisRanges[0][1] - axisRanges[0][0];
        float yRange = axisRanges[1][0] - axisRanges[1][1];
        axisRanges[0][0] = totalAxisRanges[0][0];
        axisRanges[0][1] = totalAxisRanges[0][0] + xRange;
        axisRanges[1][0] = totalAxisRanges[1][1] + yRange;
        axisRanges[1][1] = totalAxisRanges[1][1];
        computePixelScaling();
        viewChangedRepaint();
//        retainZoomLevel(currentDirNo);
    }

    /*
       public void rescaleInchesPerSecond()
        {
            if(prestackLineSet3d == null) return;
            // with our current size, how many samples can we fit ?
            float inches = (float)glPanel3d.glc.getHeight() / (float)getPixelsPerInch();
            float nSeconds = inches / prestackLineSet3d.getWiggleDisplayProperties().getInchesPerSecond();

    //        axisRanges[1][1] = axisRanges[1][1];   do nothing statement
            axisRanges[1][0] = axisRanges[1][1] + (nSeconds * 1000.0f);
            if(axisRanges[1][0] > totalAxisRanges[1][0])
                axisRanges[1][0] = totalAxisRanges[1][0];
    //        retainZoomLevel(currentDirNo);
        }

       public void rescaleTracesPerInch()
        {
            if(prestackLineSet3d == null) return;
            // with our current size, how many traces can we fit ?
            float inches = (float)glPanel3d.glc.getWidth() / (float)getPixelsPerInch();
            float nTraces = inches * prestackLineSet3d.getWiggleDisplayProperties().getTracesPerInch();

            axisRanges[0][0] = axisRanges[0][0];
            axisRanges[0][1] = axisRanges[0][0] + nTraces;
            if(axisRanges[0][1] > totalAxisRanges[0][1])
                axisRanges[0][1] = totalAxisRanges[0][1];
    //        retainZoomLevel(currentDirNo);
        }
    */
    public void retainZoomLevel(int currentDirNo)
    {
        System.out.println("Retain Zoom Levels");
//   		cursorAxisRanges[currentDirNo] = axisRanges;
    }

    /*
       public void restoreZoomLevel(int newDirNo)
       {
           System.out.println("Restore Zoom Levels");
              axisRanges = cursorAxisRanges[newDirNo];
           totalAxisRanges = cursorTotalAxisRanges[newDirNo];
           viewObjectRepaint(this);
       }
    */
    public boolean viewObjectRepaint(Object source, Object object)
    {
        if(object instanceof StsProject)
        {
            glPanel3d.repaint();
            return true;
        }
        else if(glPanel3d.window.getCursor3d().isDisplayableObject(object))
        {
            glPanel3d.repaint();
            return true;
        }
        else if(getCursor3d().isDisplayingObject(object))
        {
            glPanel3d.repaint();
            return true;
        }
        if(object == StsCursor3dPanel.isGridCoordinateChanged)
        {
            getCursor3d().setAxisLabels();
            glPanel3d.repaint();
            return true;
        }
        else if(object instanceof StsSuperGather)
        {
            centerViewOnCursors();
            glPanel3d.repaint();
            return true;
        }
        else if(object instanceof StsMultiAttributeVector)
        {
            glPanel3d.repaint();
            return true;
        }
        else if(object instanceof StsSensor)
        {
            glPanel3d.repaint();
            return true;
        }
        else if(object instanceof StsPanelProperties)  // Catch all property changes.
        {
            glPanel3d.repaint();
            return true;
        }
        else if(object instanceof StsWell)
        {
            glPanel3d.repaint();
            return true;
        }
        return false;
    }

    public void resetToTopLeft() { }

    public byte getHorizontalAxisType()
    {
        return AXIS_TYPE_DISTANCE;
    }

    public byte getVerticalAxisType()
    {
        int selectedDirection = getCursor3d().selectedDirection;
        if(selectedDirection == ZDIR)
            return AXIS_TYPE_TIME;
        else
            return AXIS_TYPE_DISTANCE;
    }

    public boolean initializeDefaultAction()
    {
        return setDefaultAction3d();
    }

    public boolean setDefaultAction3d()
    {
        StsActionManager actionManager = glPanel3d.window.actionManager;
        if (actionManager == null) return true;
        StsAction defaultAction = getDefaultAction();
        if (defaultAction instanceof StsDefaultAction3d) return false;
        setDefaultAction(new StsDefaultAction3d(glPanel3d, actionManager));
        return true;
    }
}
