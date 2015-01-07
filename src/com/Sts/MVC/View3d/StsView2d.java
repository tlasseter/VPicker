package com.Sts.MVC.View3d;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: Base class used to present two-dimensional data.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.1
 */

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;
import com.magician.fonts.*;

import javax.media.opengl.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.text.*;

abstract public class StsView2d extends StsView implements StsSerializable
{
    boolean displayAxes = true;
    /** horizontal and vertical axes are flipped for this display */
    public boolean axesFlipped = false;
    /** range on project x axis is reversed */
    public boolean xAxisReversed = false;
    /** range on project y axis is reversed */
    public boolean yAxisReversed = false;
    /** Axis ranges, including zoom */
    public float[][] axisRanges = null;
    /** Total axis ranges, excluding zoom */
    public float[][] totalAxisRanges;
    /** Decimal format for labels */
    public String labelFormatString = "###0.#";
    /** Fractional digits format */
    public String decimalFormatString = "###0.0###";
    /** Axis Labels */
    public String[] axisLabels = new String[2];
    /** classInitialize a weighpoint */
    StsWeighPoint wp = null;

    /** Insets are to accommodate labeling and decorations. noInsets used when axes are off */
    transient public Insets insetsOff = new Insets(0, 0, 0, 0);
    /** Insets are to accommodate labeling and decorations. noInsets used when axes are on */
    transient public Insets insetsOn = new Insets(0, 0, 0, 0);
    /** Current insets */
    transient public Insets insets;
    /** Font dimensions */
    transient public Dimension horizontalFontDimension, verticalFontDimension;
    /** pixels per unit in X horizontal direction */
    transient public float pixelsPerXunit;
    /** pixels per unit in Y vertical direction */
    transient public float pixelsPerYunit;
    /** Decimal format for labels */
    transient public DecimalFormat labelFormat = null;
    /** Fractional digits format */
    transient public DecimalFormat decimalFormat = null;
    /** point in range where mouse was pressed; used in zooming about this point */
    transient public float[] mousePressedLocation = null;
    /** Font used on horizontal axis */
    transient public GLBitmapFont horizontalFont = null;
    /** title font */
    transient public GLBitmapFont titleFont = null;
    /** Vertical axis font */
    transient public StsVerticalFont verticalFont;
    /** limits the ability to move window beyond totalAxisRanges */
    transient public boolean limitPan = true;

    /** width of axis lines */
    static public final int lineWidth = 2;
    /** half width of axis lines */
    static public final int halfWidth = lineWidth / 2;
    /** space between tickLabel line and axisLabel */
    static public final int fontLineSpace = 2;
    static public final int majorTickLength = 6;
    static public final int minorTickLength = 3;
    /** minimum allowed space between tick labels */
    static public final int labelGap = 5;
    /** used in computing "nice" major/minor ticks */
    static public final int minMinorTickSpacing = 25;
    /** factor used in vertical exaggeration */
    static final float stretchFactor = 4.0f / 3.0f;

    static public final byte VERTICAL = 0;
    static public final byte HORIZONTAL = 1;
    static public final byte BOTH = 2;

    private transient StsMenuItem aspectBtn = new StsMenuItem();

    static public boolean debug = false;

    abstract public void resetToOrigin();

    abstract public void display(GLAutoDrawable component);

    abstract public void setAxisRanges();

    /** override in subclasses if 2D view scale factors (units per length) should be applied when range is changed */
    public void resetRangeWithScale()
    { }

    public StsView2d()
    {
    }

    /**
     * StsView2D constructor
     *
     * @param glPanel3d - the graphics context
     */
    public StsView2d(StsGLPanel3d glPanel3d)
    {
        super(glPanel3d);
    }

    /**
     * Initialize the 2d view plot by setting the graphics materials and view matrices
     *
     * @param drawable
     * @params component the drawable component
     */
    public void init(GLAutoDrawable drawable)
    {
        if (isViewGLInitialized) return;
        super.init(drawable);
        //		initialize(gl); // jbw
        isViewGLInitialized = true;
    }

    protected void initializeView()
    {
        gl = glPanel.getGL();
        glu = glPanel.getGLU();
        gl.glShadeModel(GL.GL_FLAT);
        gl.glEnable(GL.GL_ALPHA_TEST);
        gl.glAlphaFunc(GL.GL_GREATER, 0.1f);

        initializeFont(gl);
        initializeInsets();
        computePixelScaling();
        computeProjectionMatrix();
        computeModelViewMatrix();
    }

    public void initializeTransients(StsGLPanel3d glPanel)
    {
        super.initializeTransients(glPanel);
        if(xBtn.getActionListeners().length == 0)
        {
            xBtn.setMenuActionListener("X / Crossline", this, "changeActiveSlice", xBtn);
            yBtn.setMenuActionListener("Y / Inline", this, "changeActiveSlice", yBtn);
            zBtn.setMenuActionListener("Time / Depth", this, "changeActiveSlice", zBtn);
            aspectBtn.setMenuActionListener("1:1 Aspect Ratio", this, "aspectRatio", null);
        }
        //        computePixelScaling();
    }

    /**
     * Initialize the 2D view
     * Cannot override (final).
     * Subclass must implement doInitialize().
     */
    public void initialize()
    {
        setAxisRanges();
    }

    public void initialize(GL gl)
    {
        if (isInitialized) return;
        initializeFont(gl);
        initializeInsets();
        computePixelScaling();
        isInitialized = true;
    }

    /** Initialize the graphics fonts */
    protected void initializeFont(GL gl)
    {
        horizontalFont = GLHelvetica12BitmapFont.getInstance(gl);
        // NOTE: 18 pixel font is too big. Pre-stack windows are usually narrow and title
        // disappears.
        titleFont = GLHelvetica12BitmapFont.getInstance(gl);
        verticalFont = new StsVerticalFont(horizontalFont, gl);

        GLBitmapChar fontChar;
        fontChar = horizontalFont.getCharData('O');
        horizontalFontDimension = new Dimension(fontChar.width, fontChar.height);
        fontChar = verticalFont.getCharData('O');
        verticalFontDimension = new Dimension(fontChar.width, fontChar.height);
    }

    public void setLabelFormatString(String string)
    {
        labelFormatString = new String(string);
        labelFormat = new DecimalFormat(labelFormatString);
    }

    public DecimalFormat getLabelFormat()
    {
        if (labelFormat == null)
            labelFormat = new DecimalFormat(labelFormatString);
        return labelFormat;
    }

    public void setDecimalFormatString(String string)
    {
        decimalFormatString = new String(string);
        decimalFormat = new DecimalFormat(decimalFormatString);
    }

    public DecimalFormat getDecimalFormat()
    {
        if (decimalFormat == null)
            decimalFormat = new DecimalFormat(decimalFormatString);
        return decimalFormat;
    }

    //	abstract public void doInitialize();

    /**
     * The total axis ranges regardless of zoom
     *
     * @return X and Y data ranges [2][2] (xMin,xMax,yMin,yMax)
     */
    public float[][] getTotalAxisRanges()
    {
        return totalAxisRanges;
    }

    /**
     * The axis range in view, including zoom
     *
     * @return X and Y data ranges [2][2] (xMin,xMax,yMin,yMax)
     */
    public float[][] getAxisRanges()
    {
        return axisRanges;
    }

    public float getVerticalScale()
    {
        return (axisRanges[1][1] - axisRanges[1][0]) / getInsetHeight();
    }

    public float getHorizontalScale()
    {
        return (axisRanges[0][1] - axisRanges[0][0]) / getInsetWidth();
    }

    /** Inset the viewport to acommodate the labeling */
    public void insetViewPort()
    {
        if (insets == null) return;
        //        int oldHeight = glPanel3d.getViewPortHeight();
        //        int oldWidth = glPanel3d.getViewPortWidth();
        int width = getInsetWidth();
        int height = getInsetHeight();
        int x = glPanel.viewPort[0];
        int y = glPanel.viewPort[1];
        gl.glViewport(x + insets.left, y + insets.bottom, width, height);
        //       glPanel3d.debugPrintViewport("StsView2d.insetViewPort()");
        //        glPanel3d.setViewPort(x + insets.left, y + insets.bottom, width, height);
        // Need to avoid resetting viewPort on every draw operation. TJL 2/21/07
        //		glPanel3d.viewPortChanged = true;
        //System.out.println("StsVewXP.insetViewPort() StsViewXP: " + toString() + " glPanel3d: " + glPanel3d.toString());
    }

    public boolean isInsideInsetViewPort(int x, int y)
    {
        if (x < insets.left || x > getWidth() - insets.right) return false;
        if (y < insets.bottom || y > getHeight() - insets.top) return false;
        return true;
    }

    final protected int getInsetWidth()
    {
        int width = getWidth();
        if (!displayAxes)
            return width;
        else
            return width - insets.left - insets.right;
    }


    final protected int getInsetHeight()
    {
        int height = getHeight();
        if (!displayAxes || insets == null)
            return height;
        else
            return height - insets.top - insets.bottom;
    }

    /** initialize insets after font has been set. */
    public void initializeInsets()
    {
        int leftInset = halfWidth + majorTickLength + 2 * verticalFontDimension.width + 2 * fontLineSpace;
        int bottomInset = halfWidth + majorTickLength + 2 * horizontalFontDimension.height + 2 * fontLineSpace;
        int topInset = halfWidth + majorTickLength + 2 * horizontalFontDimension.height + 2 * fontLineSpace;
        int rightInset = majorTickLength + lineWidth * 2;
        insetsOn = new Insets(topInset, leftInset, bottomInset, rightInset);

        setInsets(displayAxes);
    }

    /**
     * Set the insets based on whether the axis labeling in on or off.
     *
     * @param axisOn Is labeling on, true mean on
     */
    public void setInsets(boolean axisOn)
    {
        insets = axisOn ? insetsOn : insetsOff;
    }

    /**
     * Get the current data diaplayable area excluding any insets for labeling
     *
     * @return a Rectangle describing the data area only
     */
    public Rectangle getInsetRectangle()
    {
        return new Rectangle(insets.left, insets.bottom, getInsetWidth(), getInsetHeight());
    }

    /**
     * Get the inset values
     *
     * @return the left, bottom, right, and top inset values
     */
    public int[] getInsetViewPort()
    {
        return new int[]
            {insets.left, insets.bottom, getInsetWidth(), getInsetHeight()};
    }

    /**
     * Get the Rectangle of the drawable area excluding insets
     *
     * @return a Rectangle describing the entire graphics area
     */
    public Rectangle getDrawableRectangle()
    {
        return new Rectangle(0, 0, getWidth(), getHeight());
    }

    public float[] getInsetCoorPoint(StsMouse mouse)
    {
        // mousePoint is measured from upper left of view; remove left and top offsets from x & y
        // and convert to coordinates using axisRanges
        float[] coorXY = new float[2];
        StsMousePoint mousePoint = mouse.getMousePoint();
        int x = mousePoint.x - insets.left;
        float f = (float) x / getInsetWidth();
        coorXY[0] = axisRanges[0][0] + f * (axisRanges[0][1] - axisRanges[0][0]);
        int y = mousePoint.y - insets.top;
        f = (float) y / getInsetHeight();
        coorXY[1] = axisRanges[1][1] + f * (axisRanges[1][0] - axisRanges[1][1]);
        return coorXY;
    }

    public int getWidth(){ return glPanel.getWidth(); }

    public int getHeight(){ return glPanel.getHeight(); }

    /**
     * Draw the vertical like horizontal axis only rotated CCW to vertical
     *
     * @param axisLabel   the axis label
     * @param range       the data range[2]
     * @param timingLines timing lines on/off
     */
    public void drawVerticalAxis(String axisLabel, float[] range, boolean timingLines)
    {
        if (axisLabel == null) return;
        drawVerticalAxis(axisLabel, range, timingLines, gl);
    }

    public void drawVerticalAxis(String axisLabel, float[] range, boolean timingLines, GL gl)
    {
        int y;
        if (range[0] == range[1]) return;
        gl.glDisable(GL.GL_LIGHTING);
        getGridColor().setGLColor(gl);
        gl.glLineWidth(lineWidth);

        setFractionalFormat(range);

        int leftAxisX = insets.left - halfWidth;
        int rightAxisX = getWidth() - insets.right + halfWidth;
        int tickX = insets.left - lineWidth;
        int tickLabelX = tickX - majorTickLength;
        int axisLabelX = tickLabelX - verticalFontDimension.width - fontLineSpace;

        int axisMinX = insets.left;
        int axisMaxX = getWidth() - insets.right;
        axisMaxX = Math.max(axisMaxX, axisMinX + 1);
        int axisMinY = insets.bottom;
        int axisMaxY = getHeight() - insets.top;
        axisMaxY = Math.max(axisMaxY, axisMinY + 1);

        int maxNumMinorTicks = (axisMaxY - axisMinY) / minMinorTickSpacing;
        //        System.out.println("Vertical axis scale: " + range[0] + " " +  range[1] + " " + maxNumMinorTicks);
        if (range[0] == range[1]) return;
        double[] niceScale = StsMath.niceScaleBold(range[0], range[1], maxNumMinorTicks, true);
        //        System.out.println("Vertical axis niceScale: " + niceScale[0] + " " + niceScale[1] + " " + niceScale[2] + " " + niceScale[3]);
        float majorTickStartValue = (float) niceScale[0];
        float majorTickEndValue = (float) niceScale[1];
        float minorTickInterval = (float) niceScale[2];
        float majorTickInterval = (float) niceScale[3];

        float gridInterval = getGridInterval();

        if (gridInterval != StsProject.AUTO)
        {
            if ((majorTickInterval < 0) && (gridInterval > 0))
                gridInterval = -gridInterval;
            else if ((majorTickInterval > 0) && (gridInterval < 0))
                gridInterval = -gridInterval;
            minorTickInterval = gridInterval;
            majorTickInterval = minorTickInterval * 5;
        }

        /*
                if(model.getProject().getLabelFrequency() != StsProject.AUTO)
                {
                    if((majorTickInterval < 0) && (model.getProject().getLabelFrequency() > 0))
                        model.getProject().setLabelFrequency( -model.getProject().getLabelFrequency());
                    else if((majorTickInterval > 0) && (model.getProject().getLabelFrequency() < 0))
                        model.getProject().setLabelFrequency( -model.getProject().getLabelFrequency());
                    majorTickInterval = model.getProject().getLabelFrequency();
                    minorTickInterval = majorTickInterval / 5.0f;
                }
        */
        float pixelsPerUnit = (axisMaxY - axisMinY) / (range[1] - range[0]);
        float tickStartY = (majorTickStartValue - range[0]) * pixelsPerUnit + axisMinY - halfWidth;
        float majorTickIncY = majorTickInterval * pixelsPerUnit;
        float minorTickIncY = minorTickInterval * pixelsPerUnit;

        // draw vertical line along left-side
        gl.glBegin(GL.GL_LINE_STRIP);
        gl.glVertex2i(leftAxisX, axisMinY - lineWidth);
        gl.glVertex2i(leftAxisX, axisMaxY + lineWidth);
        gl.glEnd();
        gl.glBegin(GL.GL_LINE_STRIP);
        gl.glVertex2i(rightAxisX, axisMinY - lineWidth);
        gl.glVertex2i(rightAxisX, axisMaxY + lineWidth);
        gl.glEnd();

        // draw minRangeLabel from bottom of line up

        String minRangeLabel = getLabelFormat().format(range[0]);
        int minLabelLength = StsGLDraw.getFontStringLength(verticalFont, minRangeLabel);
        StsGLDraw.verticalFontOutput(gl, tickLabelX, axisMinY, minRangeLabel, verticalFont);

        // draw maxRangeLabel from top of line down
        String maxRangeLabel = getLabelFormat().format(range[1]);
        int maxLabelLength = StsGLDraw.getFontStringLength(verticalFont, maxRangeLabel);
        y = axisMaxY - maxLabelLength;
        StsGLDraw.verticalFontOutput(gl, tickLabelX, y, maxRangeLabel, verticalFont);

        float tickValue;
        float tickYF;
        int tickY;

        // draw major ticks
        tickValue = majorTickStartValue;
        tickYF = tickStartY;
        int openSpaceYMin = axisMinY + minLabelLength;
        int openSpaceYMax = axisMaxY - maxLabelLength;
        while (true)
        {
            if (beyondRange(tickValue, range)) break;
            if (StsMath.betweenInclusive(tickValue, range[0], range[1]))
            {
                getGridColor().setGLColor(gl);
                gl.glLineWidth(lineWidth);
                tickY = Math.round(tickYF);

                gl.glBegin(GL.GL_LINE_STRIP);
                gl.glVertex2i(tickX, tickY);
                gl.glVertex2i(tickX - majorTickLength, tickY);
                gl.glEnd();
                gl.glBegin(GL.GL_LINE_STRIP);
                gl.glVertex2i(axisMaxX + lineWidth, tickY);
                gl.glVertex2i(axisMaxX + lineWidth + majorTickLength, tickY);
                gl.glEnd();

                String tickLabel = getLabelFormat().format(tickValue);
                int tickLabelLength = StsGLDraw.getFontStringLength(verticalFont, tickLabel);
                y = tickY - tickLabelLength / 2;
                if (y - labelGap > openSpaceYMin && y + tickLabelLength + labelGap < openSpaceYMax)
                {
                    StsGLDraw.verticalFontOutput(gl, tickLabelX, y, tickLabel, verticalFont);
                    openSpaceYMin = y + tickLabelLength;
                }

                if (timingLines)
                {
                    model.getProject().getStsTimingColor().setGLColor(gl);
                    gl.glLineWidth(2);
                    gl.glBegin(GL.GL_LINE_STRIP);
                    gl.glVertex2i(tickX - majorTickLength, tickY);
                    gl.glVertex2i(axisMaxX, tickY);
                    gl.glEnd();
                }
            }
            tickValue += majorTickInterval;
            tickYF += majorTickIncY;
        }

        // draw minor ticks
        tickValue = majorTickStartValue;
        tickYF = tickStartY;
        getGridColor().setGLColor(gl);
        gl.glLineWidth(lineWidth);
        while (true)
        {
            if (beyondRange(tickValue, range)) break;
            if (StsMath.betweenInclusive(tickValue, range[0], range[1]))
            {
                tickY = Math.round(tickYF);
                gl.glBegin(GL.GL_LINE_STRIP);
                gl.glVertex2i(tickX, tickY);
                gl.glVertex2i(tickX - minorTickLength, tickY);
                gl.glEnd();
                gl.glBegin(GL.GL_LINE_STRIP);
                gl.glVertex2i(axisMaxX + lineWidth, tickY);
                gl.glVertex2i(axisMaxX + lineWidth + minorTickLength, tickY);
                gl.glEnd();
                if (timingLines)
                {
                    model.getProject().getStsTimingColor().setGLColor(gl);
                    gl.glLineWidth(1);
                    gl.glBegin(GL.GL_LINE_STRIP);
                    gl.glVertex2i(tickX - minorTickLength, tickY);
                    gl.glVertex2i(axisMaxX, tickY);
                    gl.glEnd();
                }
            }
            tickValue += minorTickInterval;
            tickYF += minorTickIncY;
        }

        int axisLabelLength = StsGLDraw.getFontStringLength(verticalFont, axisLabel);
        y = (axisMaxY - axisLabelLength) / 2;
        StsGLDraw.verticalFontOutput(gl, axisLabelX, y, axisLabel, verticalFont);

        /** Don't forget to finally turn lighting back on when axes are drawn
         *  as this is the default.
         */
    }

    public float getGridInterval()
    {
        return model.getProject().getGridFrequency();
    }

    protected StsColor getGridColor()
    {
        return model.getProject().getGridColor();
    }

    /**
     * Draw horizontal axis with this layout:
     * -------------------------------------------------------------
     * lineSpacing
     * AXIS LABEL                                 horizontalFontHeight
     * lineSpacing
     * 0    10    20     30    40    50    60    70   80    90   100  horizontalFontHeight
     * |.....|.....|.....|.....|.....|.....|.....|.....|.....|.....|  majorTickLength
     * <p/>
     * where "|" is a major tick and "." is a minor tick.
     * Only major ticks are labeled.  First and last tick label are
     * aligned with end, and interior tick labels are centered on tick.
     * Any tick label which crowds a neighbor is not drawn.
     * Parameters controlling layout are specified as static finals
     * (see class members for definitions).
     *
     * @param axisLabel   the axis label
     * @param range       the data range[2]
     * @param timingLines timing lines on/off
     */
    public void drawHorizontalAxis(String title, String axisLabel, float[] range, boolean timingLines)
    {
        if (axisLabel == null) return;
        drawHorizontalAxis(title, axisLabel, range, timingLines, false, gl);
    }

    public void drawHorizontalAxis(String title, String axisLabel, float[] range, boolean timingLines, boolean onTop, GL gl)
    {
        drawHorizontalAxis(title, axisLabel, range, timingLines, false, gl, false, null);
    }


    private float findAux(float ix, double[][] auxValues)
    {
        float ret = 0;
        for (int i = 0; i < auxValues[0].length; i++)
        {
            if (((int) ix) == ((int) auxValues[0][i]))
                return (float) auxValues[1][i];
        }
        return 0;
    }

    public void drawHorizontalAxis(String title, String axisLabel, float[] range, boolean timingLines, boolean onTop, GL gl, boolean round, double[][] auxValues)
    {
        int x;

        if (range[0] == range[1]) return;

        gl.glDisable(GL.GL_LIGHTING);

        getGridColor().setGLColor(gl);
        gl.glLineWidth(lineWidth);

        setFractionalFormat(range);

        int btmAxisY = insets.bottom - halfWidth;
        int topAxisY = getHeight() - insets.top + halfWidth;
        int tickY = insets.bottom - lineWidth;
        int tickLabelY, axisLabelY;

        if (!onTop)
        {

            tickLabelY = tickY - majorTickLength - horizontalFontDimension.height;
            axisLabelY = tickLabelY - fontLineSpace - horizontalFontDimension.height;

        }
        else
        {
            int t = getHeight() - insets.top + lineWidth;
            tickLabelY = t + majorTickLength;
            axisLabelY = tickLabelY + (fontLineSpace + horizontalFontDimension.height);
        }

        int axisMinX = insets.left;
        int axisMaxX = getWidth() - insets.right;
        int axisMinY = insets.bottom;
        int axisMaxY = getHeight() - insets.top;
        int axisMiddleX = (getWidth() - insets.right) / 2;

        int maxNumMinorTicks = (axisMaxX - axisMinX) / minMinorTickSpacing;
        double[] niceScale = StsMath.niceScaleBold(range[0], range[1], maxNumMinorTicks, true);
        //        System.out.println("Horizontal axis niceScale: " + niceScale[0] + " " +  + niceScale[1] + " " + niceScale[2] + " " + niceScale[3]);
        float majorTickStartValue = (float) niceScale[0];
        float majorTickEndValue = (float) niceScale[1];
        float minorTickInterval = (float) niceScale[2];
        float majorTickInterval = (float) niceScale[3];
        float gridInterval = model.getProject().getGridFrequency();
        if (model.getProject().getGridFrequency() != StsProject.AUTO)
        {
            if ((majorTickInterval < 0) && (gridInterval > 0))
                gridInterval = -gridInterval;
            else if ((majorTickInterval > 0) && (gridInterval < 0))
                gridInterval = -gridInterval;
            minorTickInterval = gridInterval;
            majorTickInterval = minorTickInterval * 5.0f;
        }
        if (majorTickInterval == 0)
        {
            System.out.println("axis tick = 0 !");
            return;
        }
        /*
                if(model.getProject().getLabelFrequency() != StsProject.AUTO)
                {
                    if((majorTickInterval < 0) && (model.getProject().getLabelFrequency() > 0))
                        model.getProject().setLabelFrequency( -model.getProject().getLabelFrequency());
                    else if((majorTickInterval > 0) && (model.getProject().getLabelFrequency() < 0))
                        model.getProject().setLabelFrequency( -model.getProject().getLabelFrequency());
                    majorTickInterval = model.getProject().getLabelFrequency();
                    minorTickInterval = majorTickInterval / 5.0f;
                }
        */
        float pixelsPerUnit = (axisMaxX - axisMinX) / (range[1] - range[0]);
        float tickStartX = (majorTickStartValue - range[0]) * pixelsPerUnit + axisMinX; // jbw - halfWidth;
        float majorTickIncX = majorTickInterval * pixelsPerUnit;
        float minorTickIncX = minorTickInterval * pixelsPerUnit;

        // draw horizontal line along bottom
        gl.glBegin(GL.GL_LINE_STRIP);
        gl.glVertex2i(axisMinX - lineWidth, btmAxisY);
        gl.glVertex2i(axisMaxX + lineWidth, btmAxisY);
        gl.glEnd();
        gl.glBegin(GL.GL_LINE_STRIP);
        gl.glVertex2i(axisMinX - lineWidth, topAxisY);
        gl.glVertex2i(axisMaxX + lineWidth, topAxisY);
        gl.glEnd();

        int titlePosition = axisMiddleX - StsGLDraw.getFontStringLength(titleFont, title) / 2;
        if (titlePosition < 1) titlePosition = 1;

        if (!onTop)
            StsGLDraw.fontOutput(gl, titlePosition, topAxisY + insets.top / 2, title, titleFont);
        else
            StsGLDraw.fontOutput(gl, titlePosition, axisLabelY + (fontLineSpace + horizontalFontDimension.height + fontLineSpace)
                , title, titleFont);

        // draw minRangeLabel from left of line to the right
        int minLabelLength, maxLabelLength;
        if (auxValues == null)
        {
            if (!round)
            {
                String minRangeLabel = getLabelFormat().format(range[0]);
                minLabelLength = StsGLDraw.getFontStringLength(horizontalFont, minRangeLabel);
                StsGLDraw.fontOutput(gl, axisMinX, tickLabelY, minRangeLabel, horizontalFont);
            }
            else
            {
                float val = (float) Math.ceil((double) range[0]);
                if (val == -0.0f) val = 0.0f;
                int locX = (int) ((val - range[0]) * pixelsPerUnit) + axisMinX;
                String minRangeLabel = getLabelFormat().format(val);
                minLabelLength = StsGLDraw.getFontStringLength(horizontalFont, minRangeLabel);
                locX -= minLabelLength / 2;
                StsGLDraw.fontOutput(gl, locX, tickLabelY, minRangeLabel, horizontalFont);
            }
            // draw maxRangeLabel from right of line to the left
            if (!round)
            {
                String maxRangeLabel = getLabelFormat().format(range[1]);
                maxLabelLength = StsGLDraw.getFontStringLength(horizontalFont, maxRangeLabel);
                x = axisMaxX - maxLabelLength;
                StsGLDraw.fontOutput(gl, x, tickLabelY, maxRangeLabel, horizontalFont);
            }
            else
            {
                float val = (float) Math.floor((double) range[1]);
                String maxRangeLabel = getLabelFormat().format(val);
                maxLabelLength = StsGLDraw.getFontStringLength(horizontalFont, maxRangeLabel);
                int locX = axisMaxX - maxLabelLength;
                StsGLDraw.fontOutput(gl, locX, tickLabelY, maxRangeLabel, horizontalFont);
            }
        }
        else
        {
            {
                String minRangeLabel = getLabelFormat().format(findAux(range[0], auxValues));
                minLabelLength = StsGLDraw.getFontStringLength(horizontalFont, minRangeLabel);
                StsGLDraw.fontOutput(gl, axisMinX - minLabelLength / 2, tickLabelY, minRangeLabel, horizontalFont);

                String maxRangeLabel = getLabelFormat().format(findAux(range[1], auxValues));
                maxLabelLength = StsGLDraw.getFontStringLength(horizontalFont, maxRangeLabel);
                x = axisMaxX - maxLabelLength;
                StsGLDraw.fontOutput(gl, x, tickLabelY, maxRangeLabel, horizontalFont);

            }

        }
        float tickValue;
        float tickXF;
        int tickX;

        // draw  major ticks
        tickValue = majorTickStartValue;
        tickXF = tickStartX;
        int openSpaceXMin = axisMinX + minLabelLength;
        int openSpaceXMax = axisMaxX - maxLabelLength;
        while (true)
        {
            if (beyondRange(tickValue, range)) break;
            if (StsMath.betweenInclusive(tickValue, range[0], range[1]))
            {
                getGridColor().setGLColor(gl);
                gl.glLineWidth(lineWidth);
                tickX = Math.round(tickXF);
                gl.glBegin(GL.GL_LINE_STRIP);
                gl.glVertex2i(tickX, tickY);
                gl.glVertex2i(tickX, tickY - majorTickLength);
                gl.glEnd();
                gl.glBegin(GL.GL_LINE_STRIP);
                gl.glVertex2i(tickX, axisMaxY + halfWidth);
                gl.glVertex2i(tickX, axisMaxY + halfWidth + majorTickLength);
                gl.glEnd();

                String tickLabel = "";
                if (auxValues == null)
                    tickLabel = getLabelFormat().format(tickValue);
                else
                    tickLabel = getLabelFormat().format(findAux(tickValue, auxValues));
                int tickLabelLength = StsGLDraw.getFontStringLength(horizontalFont, tickLabel);
                x = tickX - tickLabelLength / 2;
                if (x - labelGap > openSpaceXMin && x + tickLabelLength + labelGap < openSpaceXMax)
                {
                    StsGLDraw.fontOutput(gl, x, tickLabelY, tickLabel, horizontalFont);
                    openSpaceXMin = x + tickLabelLength;
                }

                if (timingLines)
                {
                    model.getProject().getStsTimingColor().setGLColor(gl);
                    gl.glLineWidth(2);
                    gl.glBegin(GL.GL_LINE_STRIP);
                    gl.glVertex2i(tickX, tickY - majorTickLength);
                    gl.glVertex2i(tickX, axisMaxY);
                    gl.glEnd();
                }

            }
            tickValue += majorTickInterval;
            tickXF += majorTickIncX;
        }

        // draw  minor ticks
        tickValue = majorTickStartValue;
        tickXF = tickStartX;
        getGridColor().setGLColor(gl);
        gl.glLineWidth(lineWidth);
        while (true)
        {
            if (beyondRange(tickValue, range)) break;
            if (StsMath.betweenInclusive(tickValue, range[0], range[1]))
            {
                tickX = Math.round(tickXF);
                gl.glBegin(GL.GL_LINE_STRIP);
                gl.glVertex2i(tickX, tickY);
                gl.glVertex2i(tickX, tickY - minorTickLength);
                gl.glEnd();
                gl.glBegin(GL.GL_LINE_STRIP);
                gl.glVertex2i(tickX, axisMaxY + halfWidth);
                gl.glVertex2i(tickX, axisMaxY + halfWidth + minorTickLength);
                gl.glEnd();

                if (timingLines)
                {
                    model.getProject().getStsTimingColor().setGLColor(gl);
                    gl.glLineWidth(1);
                    gl.glBegin(GL.GL_LINE_STRIP);
                    gl.glVertex2i(tickX, tickY - minorTickLength);
                    gl.glVertex2i(tickX, axisMaxY);
                    gl.glEnd();
                }
            }
            tickValue += minorTickInterval;
            tickXF += minorTickIncX;
        }

        int axisLabelLength = StsGLDraw.getFontStringLength(horizontalFont, axisLabel);
        x = (axisMaxX - axisLabelLength) / 2;

        StsGLDraw.fontOutput(gl, x, axisLabelY, axisLabel, horizontalFont);

        /*
          if(auxValues != null)
          {
              int locY;

              locY = tickLabelY;
              int lastX = -999;
              for(int i = 0; i < auxValues[0].length; i++)
              {
                  int locX = (int)((auxValues[0][i] - range[0]) * pixelsPerUnit) + axisMinX;
                  String label = labelFormat.format(auxValues[1][i]);
                  int labelLength = StsGLDraw.getFontStringLength(horizontalFont, label);
                  locX -= labelLength / 2;
                  if(locX > lastX + (labelLength/2))
                  {
                      StsGLDraw.fontOutput(gl, locX, locY, label, horizontalFont);
                      lastX = locX + labelLength;
                  }
              }
          }
        */
        /** Don't forget to finally turn lighting back on when axes are drawn
         *  as this is the default.
         */
    }

    private boolean beyondRange(float value, float[] range)
    {
        if (range[1] >= range[0]) return value > range[1];
        else return value < range[1];
    }

    private void setFractionalFormat(float[] range)
    {
        float maxAbsValue = Math.max(Math.abs(range[0]), Math.abs(range[1]));
        if (maxAbsValue > 1000.0f)
            getDecimalFormat().setMaximumFractionDigits(0);
        else
            getDecimalFormat().setMaximumFractionDigits(10);
    }

    /**
     * Compute the relative point picked by the mouse position
     *
     * @param mouse the StsMouse object
     * @return X,Y point represented by current mouse position
     * @see #computeMousePoint(StsPoint2D)
     */
    public StsPoint computePickPoint(StsMouse mouse)
    {
        StsMousePoint mousePoint = mouse.getMousePoint();
        Rectangle viewRectangle = getInsetRectangle();
        float mx = (float) (mousePoint.x - viewRectangle.x);
        float my = (float) (viewRectangle.height - mousePoint.y + insets.top);
        float fx = mx / viewRectangle.width;
        float fy = my / viewRectangle.height;
        float x = axisRanges[0][0] + fx * (axisRanges[0][1] - axisRanges[0][0]);
        float y = axisRanges[1][0] + fy * (axisRanges[1][1] - axisRanges[1][0]);
        return new StsPoint(x, y, 0.0f);
    }

    public void setDefaultView()
    {
        axisRanges = StsMath.copyFloatArray(totalAxisRanges);
        resetRangeWithScale();
        repaint();
    }

    private void moveUp()
    {
        float min = axisRanges[1][0];
        float max = axisRanges[1][1];
        float pageSize = (axisRanges[1][1] - axisRanges[1][0]) * 0.75f;
        axisRanges[1][0] = min + pageSize;
        axisRanges[1][1] = max + pageSize;
        if (limitPan) limitPan(true);
        repaint();
    }

    private void moveDown()
    {
        float min = axisRanges[1][0];
        float max = axisRanges[1][1];
        float pageSize = (axisRanges[1][1] - axisRanges[1][0]) * 0.75f;
        axisRanges[1][0] = min - pageSize;
        axisRanges[1][1] = max - pageSize;
        if (limitPan) limitPan(true);
        repaint();
    }

    private void moveLeft()
    {
        float min = axisRanges[0][0];
        float max = axisRanges[0][1];
        float pageSize = (axisRanges[0][1] - axisRanges[0][0]) * 0.75f;
        axisRanges[0][0] = min - pageSize;
        axisRanges[0][1] = max - pageSize;
        if (limitPan) limitPan(true);
        repaint();
    }

    private void moveRight()
    {
        float min = axisRanges[0][0];
        float max = axisRanges[0][1];
        float pageSize = (axisRanges[0][1] - axisRanges[0][0]) * 0.75f;
        axisRanges[0][0] = min + pageSize;
        axisRanges[0][1] = max + pageSize;
        if (limitPan) limitPan(true);
        repaint();
    }

    public StsPoint computePickPoint(StsMousePoint mousePoint)
    {
        Rectangle viewRectangle = getInsetRectangle();
        float mx = (float) (mousePoint.x - viewRectangle.x);
        float my = (float) (viewRectangle.height - mousePoint.y + insets.top);
        float fx = mx / viewRectangle.width;
        float fy = my / viewRectangle.height;
        float x = axisRanges[0][0] + fx * (axisRanges[0][1] - axisRanges[0][0]);
        float y = axisRanges[1][0] + fy * (axisRanges[1][1] - axisRanges[1][0]);
        //if (insets != null)
        //   return new StsPoint(x-insets.left, y+insets.top, 0.0f);
        //else
        return new StsPoint(x, y, 0.0f);
    }

    /**
     * Key Pressed event handling.
     * Up key or A key pressed, right button mouse released increases vertical scale by 2.0.
     * Down key or S key pressed, right button mouse released increases vertical scale by 0.5.
     * If Z key pressed, mouse pressed motion results in zooming.
     * If X key pressed, mouse pressed motion results in panning
     * If no key pressed, mouse pressed motion results in rotation.
     *
     *
     *
     * @param mouse mouse object
     * @param e     key event
     * @return true if successful
     */
    public boolean keyPressed(StsMouse mouse, KeyEvent e)
    {
        super.keyPressed(mouse, e);
        if (getKeyCode() != KeyEvent.VK_R && glPanel3d != null) // cancel any rectangle zoom
        {
            StsMouseActionToolbar toolbar = glPanel3d.window.getMouseActionToolbar();
            if (toolbar != null) toolbar.resetMouseToggle();
        }
        switch (getKeyCode())
        {
            case KeyEvent.VK_UP:
                moveUp();
                break;
            case KeyEvent.VK_DOWN:
                moveDown();
                break;
            case KeyEvent.VK_LEFT:
                moveLeft();
                break;
            case KeyEvent.VK_RIGHT:
                moveRight();
                break;
            case KeyEvent.VK_PAGE_UP:
                if (this instanceof StsViewPreStack)
                    ((StsViewPreStack) this).next();
                break;
            case KeyEvent.VK_PAGE_DOWN:
                if (this instanceof StsViewPreStack)
                    ((StsViewPreStack) this).next();
                break;
            case KeyEvent.VK_D:
                if (this instanceof StsViewPreStack)
                    ((StsViewPreStack) this).displayMode();
                break;
            case KeyEvent.VK_E:
                if (this instanceof StsViewPreStack)
                    ((StsViewPreStack) this).editMode();
                break;
            default:
                break;
        }
        super.keyPressed(mouse, e);

        if (mouse.isButtonDown(StsMouse.RIGHT))
        {
            if (e.isShiftDown())
            {
                if (getKeyCode() == KeyEvent.VK_A || getKeyCode() == KeyEvent.VK_UP) horizontalStretch(mouse);
                else if (getKeyCode() == KeyEvent.VK_S || getKeyCode() == KeyEvent.VK_DOWN) horizontalShrink(mouse);
            }
            else
            {
                if (getKeyCode() == KeyEvent.VK_A || getKeyCode() == KeyEvent.VK_UP) verticalStretch(mouse);
                else if (getKeyCode() == KeyEvent.VK_S || getKeyCode() == KeyEvent.VK_DOWN) verticalShrink(mouse);
                else if (getKeyCode() == KeyEvent.VK_Z) setMotionCursor(StsCursor.ZOOM);
                else if (getKeyCode() == KeyEvent.VK_X) setMotionCursor(StsCursor.PAN);
            }
            return true;
        }
        return false;
    }

    /**
     * Compute the mouse position from the relative point
     *
     * @return mouse position
     * @see #computePickPoint(StsMouse)
     */
    public StsMousePoint computeMousePoint(StsPoint2D point)
    {
        Rectangle viewRectangle = getInsetRectangle();
        float fx = (point.x - axisRanges[0][0]) / (axisRanges[0][1] - axisRanges[0][0]);
        float fy = (point.y - axisRanges[1][0]) / (axisRanges[1][1] - axisRanges[1][0]);
        //        int mx = (int)(viewRectangle.width*fx);
        int mx = (int) ((viewRectangle.width - insets.right) * fx);
        //        int my = (int)(viewRectangle.height*(1.0f - fy));
        int my = (int) ((viewRectangle.height - insets.top) * (1.0f - fy));
        return new StsMousePoint(mx, my);
    }

    /**
     * Are the two pick points near each other. Check if the distance between picks is less than pickSize
     *
     * @param pickA    first pick point
     * @param pickB    second pick point
     * @param pickSize distance check
     * @return true if the distance between points is less than pickSize
     */
    public boolean mousePicksNear(StsPoint pickA, StsPoint pickB, int pickSize)
    {
        if (pickA == null || pickB == null) return false;
        Rectangle viewRectangle = getInsetRectangle();
        int dx = (int) (viewRectangle.width * Math.abs((pickA.v[0] - pickB.v[0]) / (axisRanges[0][1] - axisRanges[0][0])));
        if (dx > pickSize) return false;
        int dy = (int) (viewRectangle.height * Math.abs((pickA.v[1] - pickB.v[1]) / (axisRanges[1][1] - axisRanges[1][0])));
        return dy <= pickSize;
    }

    /**
     * Draws in foreground XORed against current view.  Must be called
     * again and draw same objects to erase.
     */
    public void drawXORRectangle()
    {

        gl.glDisable(GL.GL_BLEND);
        gl.glDrawBuffer(GL.GL_FRONT);
        gl.glEnable(GL.GL_COLOR_LOGIC_OP);
        gl.glLogicOp(GL.GL_XOR);
        gl.glDepthMask(false);
        gl.glDepthFunc(GL.GL_ALWAYS);

        gl.glMatrixMode(GL.GL_PROJECTION);
        gl.glPushMatrix();
        gl.glLoadIdentity();
        Rectangle viewRectangle = getInsetRectangle();
        glu.gluOrtho2D(0, viewRectangle.width, 0, viewRectangle.height);
        gl.glMatrixMode(GL.GL_MODELVIEW);
        gl.glPushMatrix();
        gl.glLoadIdentity();

        gl.glDisable(GL.GL_LIGHTING);
        gl.glDisable(GL.GL_DITHER);

        float bx = mousePressedPoint.x;
        float by = getHeight() - mousePressedPoint.y;

        insetViewPort();
        if (insets != null)
        {
            bx -= insets.left;
            by -= insets.top;
        }
        if (paintOldRectangle)
        {

            gl.glColor3f(0.3f, 0.5f, 1.f);
            gl.glBegin(GL.GL_POLYGON);
            gl.glVertex2f(bx, by);
            gl.glVertex2f(bx + XORoldMouseDelta.x, by);
            gl.glVertex2f(bx + XORoldMouseDelta.x, by - XORoldMouseDelta.y);
            gl.glVertex2f(bx, by - XORoldMouseDelta.y);
            gl.glEnd();

        }
        if (paintNewRectangle)
        {

            gl.glColor3f(0.3f, 0.5f, 1.f);
            gl.glBegin(GL.GL_POLYGON);
            gl.glVertex2f(bx, by);
            gl.glVertex2f(bx + mouseDragDelta.x, by);
            gl.glVertex2f(bx + mouseDragDelta.x, by - mouseDragDelta.y);
            gl.glVertex2f(bx, by - mouseDragDelta.y);
            gl.glEnd();

        }
        XORoldMouseDelta = mouseDragDelta;

        gl.glMatrixMode(GL.GL_PROJECTION);
        gl.glPopMatrix();
        gl.glMatrixMode(GL.GL_MODELVIEW);
        gl.glPopMatrix();

        gl.glLogicOp(GL.GL_COPY);
        gl.glDisable(GL.GL_COLOR_LOGIC_OP);
        gl.glDepthMask(true);
        gl.glDepthFunc(GL.GL_LESS);
        gl.glDrawBuffer(GL.GL_BACK);

        gl.glFlush();
        resetViewPort();
    }

    public void computeRectangleZoom()
    {
        StsPoint corner1, corner2;
        int dx, dy;

        corner1 = computePickPoint(mousePressedPoint);

        dx = mousePressedPoint.x + mouseDragDelta.x;
        dy = mousePressedPoint.y + mouseDragDelta.y;
        if ((Math.abs(mouseDragDelta.x) < 10) || (Math.abs(mouseDragDelta.y) < 10))
        {
            StsToolkit.beep();
            return;
        }

        corner2 = computePickPoint(new StsMousePoint(dx, dy));

        //System.out.println("xMin: " + axisRanges[0][0] + " xMax: " + axisRanges[0][1] + " yMin: " + axisRanges[1][0] + " yMax: " + axisRanges[1][1]);
        if (axisRanges[0][0] < axisRanges[0][1])
        {
            axisRanges[0][0] = Math.min(corner1.getX(), corner2.getX());
            axisRanges[0][1] = Math.max(corner1.getX(), corner2.getX());
        }
        else
        {
            axisRanges[0][0] = Math.max(corner1.getX(), corner2.getX());
            axisRanges[0][1] = Math.min(corner1.getX(), corner2.getX());
        }
        if (axisRanges[1][0] < axisRanges[1][1])
        {
            axisRanges[1][0] = Math.min(corner1.getY(), corner2.getY());
            axisRanges[1][1] = Math.max(corner1.getY(), corner2.getY());
        }
        else
        {
            axisRanges[1][0] = Math.max(corner1.getY(), corner2.getY());
            axisRanges[1][1] = Math.min(corner1.getY(), corner2.getY());
        }
        //System.out.println("out "+axisRanges[0][0]+" "+axisRanges[1][0]+"      "+axisRanges[0][1]+" "+axisRanges[1][1]);
        //System.out.println("xMin: " + axisRanges[0][0] + " xMax: " + axisRanges[0][1] + " yMin: " + axisRanges[1][0] + " yMax: " + axisRanges[1][1]);
        computePixelScaling();
        repaint();
        moveLockedWindows();
        StsMouseActionToolbar toolbar = glPanel3d.window.getMouseActionToolbar();
        if (toolbar != null) toolbar.resetMouseToggle();
    }

    /**
     * Pans, Zooms and Rotates the model view.
     *
     * @param mouse the mouse object
     */
    public void moveWindow(StsMouse mouse)
    {
        StsMousePoint loc = mouse.getMousePoint();

        StsMousePoint mouseDelta = mouse.getMouseDelta();
        int rightButtonState = mouse.getButtonStateCheckClear(StsMouse.RIGHT);
        //System.out.println("keyCode = "+keyCode);
        // right mouse down; if no movement (yet) classInitialize some mouse actions
        if (rightButtonState == StsMouse.RELEASED)
        {
            if (getKeyCode() == KeyEvent.VK_R || (mouseMode == StsCursor.RECTZOOM && getKeyCode() == KeyEvent.VK_UNDEFINED))
            {
                paintNewRectangle = false;
                mouseDragDelta = new StsMousePoint(loc.x - mousePressedPoint.x, loc.y - mousePressedPoint.y);
                setXORRectangle(true);

                glPanel3d.gld.setAutoSwapBufferMode(false);
                glPanel3d.gld.display();
                glPanel3d.gld.setAutoSwapBufferMode(true);

                setXORRectangle(false);

                computeRectangleZoom();
            }
            mousePressedLocation = null;
            return;
        }

        if (mouseDelta.x == 0 && mouseDelta.y == 0)
        {
            // rectangle zoom: R-key & right-button or rectangle-zoom button selected
            if (getKeyCode() == KeyEvent.VK_R || (mouseMode == StsCursor.RECTZOOM && getKeyCode() == KeyEvent.VK_UNDEFINED))
            {

                if (rightButtonState == mouse.PRESSED)
                {
                    mousePressedPoint = new StsMousePoint(mouse.getMousePoint());
                    mouseDragDelta = new StsMousePoint(loc.x - mousePressedPoint.x, loc.y - mousePressedPoint.y);

                    paintNewRectangle = true;
                    paintOldRectangle = false;
                    setXORRectangle(true);

                    glPanel3d.gld.setAutoSwapBufferMode(false);
                    glPanel3d.gld.display();
                    glPanel3d.gld.setAutoSwapBufferMode(true);

                    paintOldRectangle = true;
                    setXORRectangle(false);
                }
            }
            // mouse zoom: Z key & rightbutton or Zoom button selected
            /*
            else
            if (keyCode == KeyEvent.VK_Z || (glPanel3d.getMouseMode() == StsCursor.ZOOM && keyCode == KeyEvent.VK_UNDEFINED))
            {
                mouseZoom(mouse);
            }
            */
            return;
        }

        if (getKeyCode() == KeyEvent.VK_Z || (mouseMode == StsCursor.ZOOM && getKeyCode() == KeyEvent.VK_UNDEFINED))
        {
            mouseZoom(mouse);
        }
        else if (getKeyCode() == KeyEvent.VK_X || (mouseMode == StsCursor.PAN && getKeyCode() == KeyEvent.VK_UNDEFINED))
        {
            mousePan(mouseDelta);
        }
        else if (getKeyCode() == KeyEvent.VK_R || (mouseMode == StsCursor.RECTZOOM && getKeyCode() == KeyEvent.VK_UNDEFINED))
        {

            setXORRectangle(true);
            mouseDragDelta = new StsMousePoint(loc.x - mousePressedPoint.x, loc.y - mousePressedPoint.y);
            glPanel3d.gld.setAutoSwapBufferMode(false);
            glPanel3d.gld.display();
            glPanel3d.gld.setAutoSwapBufferMode(true);
            setXORRectangle(false);
        }
        else
        {
            mouseZoom(mouse); // for wheel
        }
        viewChangedRepaint();
        moveLockedWindows();
    }
      /**
       * Display the popup menu)
       *
       * @param mouse
       */
      public void showPopupMenu(StsMouse mouse)
      {
          JPopupMenu popupMenu = new JPopupMenu("2D View Popup");
          glPanel3d.add(popupMenu);

          JMenu sliceMenu = new JMenu("Active Slice");
          sliceMenu.add(xBtn);
          sliceMenu.add(yBtn);
          sliceMenu.add(zBtn);

          popupMenu.add(aspectBtn);
          popupMenu.add(sliceMenu);

          popupMenu.show(glPanel3d, mouse.getX(), mouse.getY());
          clearCurrentKeyPressed();
      }

    public StsCursor3d getCursor3d()
    {
        return glPanel3d.window.getCursor3d();
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

    public void mouseZoom(StsMouse mouse)
    {
        StsMousePoint mouseDelta = mouse.getMouseDelta();
        int rightButtonState = mouse.getButtonStateCheckClear(StsMouse.RIGHT);

        setMotionCursor(StsCursor.ZOOM);

        if (mousePressedLocation == null)
        {
            mousePressedLocation = getInsetCoorPoint(mouse);
            return;
        }
        else if (rightButtonState == StsMouse.RELEASED)
        {
            mousePressedLocation = null;
            return;
        }
        float zoom;
        //        System.out.println("zoom: " + zoom + " mouseDelta.y: " + mouseDelta.y + "height: " + glPanel3d.winRectGL.height);
        for (int n = 0; n < 2; n++)
        {
            float min = axisRanges[n][0];
            float max = axisRanges[n][1];

            if (n == 0)
                zoom = 1.0f + (float) mouseDelta.x / getInsetHeight();
            else
                zoom = 1.0f + (float) mouseDelta.y / getInsetHeight();

            float center;
            if (mousePressedLocation != null)
                center = mousePressedLocation[n];
            else
                center = (max + min) / 2;

            //		System.out.println(" n, zoom, center, min, max " + n + " " + zoom + " " + center + " " + min + " " + max);
            axisRanges[n][0] = zoom * (min - center) + center;
            axisRanges[n][1] = zoom * (max - center) + center;
            min = axisRanges[n][0];
            max = axisRanges[n][1];
            //			System.out.println(" new min, max " + min + " " + max);
        }
        boolean zoomedOut = mouseDelta.y > 0;
        if (limitPan) limitPan(zoomedOut);
        computePixelScaling();
    }

    protected void limitPan(boolean zoomedOut)
    {
        for (int n = 0; n < 2; n++)
        {
            float min = axisRanges[n][0];
            float max = axisRanges[n][1];
            float dRange = max - min;
            float dTotalRange = totalAxisRanges[n][1] - totalAxisRanges[n][0];
            if (Math.abs(dRange) > Math.abs(dTotalRange))
            {
                axisRanges[n][0] = totalAxisRanges[n][0];
                axisRanges[n][1] = totalAxisRanges[n][1];
            }
            else if (zoomedOut) // we have zoomed out; limit range to totalRanges
            {
                if (max > min)
                {
                    if (totalAxisRanges[n][0] > axisRanges[n][0])
                    {
                        axisRanges[n][0] = totalAxisRanges[n][0];
                        axisRanges[n][1] = axisRanges[n][0] + dRange;
                    }
                    else if (totalAxisRanges[n][1] < axisRanges[n][1])
                    {
                        axisRanges[n][1] = totalAxisRanges[n][1];
                        axisRanges[n][0] = axisRanges[n][1] - dRange;
                    }
                }
                if (max < min)
                {
                    if (totalAxisRanges[n][0] < axisRanges[n][0])
                    {
                        axisRanges[n][0] = totalAxisRanges[n][0];
                        axisRanges[n][1] = axisRanges[n][0] + dRange;
                    }
                    else if (totalAxisRanges[n][1] > axisRanges[n][1])
                    {
                        axisRanges[n][1] = totalAxisRanges[n][1];
                        axisRanges[n][0] = axisRanges[n][1] - dRange;
                    }
                }
            }
        }
    }

    private void mousePan(StsMousePoint mouseDelta)
    {
        setMotionCursor(StsCursor.PAN);

        if (mouseDelta.x != 0)
        {
            float mouseDx = mouseDelta.x;
            float xUnitsPerPixel = (axisRanges[0][1] - axisRanges[0][0]) / getInsetWidth();
            if (xUnitsPerPixel != 0.0f)
            {
                if (limitPan)
                {
                    if (mouseDx < 0)
                    {
                        float limitDx = -(totalAxisRanges[0][1] - axisRanges[0][1]) / xUnitsPerPixel;
                        if (limitDx > mouseDx) mouseDx = limitDx;
                    }
                    else if (mouseDx > 0)
                    {
                        float limitDx = -(totalAxisRanges[0][0] - axisRanges[0][0]) / xUnitsPerPixel;
                        if (limitDx < mouseDx) mouseDx = limitDx;
                    }
                }
                float dx = -mouseDx * xUnitsPerPixel;
                axisRanges[0][0] += dx;
                axisRanges[0][1] += dx;
            }
        }
        if (mouseDelta.y != 0)
        {
            float mouseDy = mouseDelta.y;
            float yUnitsPerPixel = (axisRanges[1][1] - axisRanges[1][0]) / getInsetHeight();
            if (yUnitsPerPixel != 0.0f)
            {
                if (limitPan)
                {
                    if (mouseDy < 0)
                    {
                        float limitDy = -(totalAxisRanges[1][1] - axisRanges[1][1]) / yUnitsPerPixel;
                        if (limitDy > mouseDy) mouseDy = limitDy;
                    }
                    else if (mouseDy > 0)
                    {
                        float limitDy = -(totalAxisRanges[1][0] - axisRanges[1][0]) / yUnitsPerPixel;
                        if (limitDy < mouseDy) mouseDy = limitDy;
                    }
                }
                float dy = -mouseDy * yUnitsPerPixel;
                axisRanges[1][0] += dy;
                axisRanges[1][1] += dy;
            }
        }
        repaint();
    }

    public void reCenterOnPoint(float[] xy)
    {
        float halfX = Math.abs(axisRanges[0][1] - axisRanges[0][0]) / 2.0f;
        float halfY = Math.abs(axisRanges[1][1] - axisRanges[1][0]) / 2.0f;
        axisRanges[0][0] = xy[0] + halfX;
        axisRanges[0][1] = xy[0] - halfX;
        axisRanges[1][0] = xy[1] + halfY;
        axisRanges[1][1] = xy[1] - halfY;
        viewChangedRepaint();
        repaint();
    }

    /** Increase the vertical stretch by stretchFactor times */
    public void verticalStretch(StsMouse mouse)
    {
        stretch(mouse, 1, 1.0f / stretchFactor);
    }

    /**
     * stretches or shrinks horizontal or vertical scale by factor
     *
     * @param mouse  StsMouse
     * @param index  0 for horizontal stretch, 1 for vertical stretch
     * @param factor factor for stretch or shrink
     */
    public void stretch(StsMouse mouse, int index, float factor)
    {
        float min = axisRanges[index][0];
        float max = axisRanges[index][1];
        float center;
        if (mouse == null)
            center = (min + max) / 2;
        else
        {
            StsPoint mousePoint = computePickPoint(mouse); // returns xy in world coordinates
            center = mousePoint.v[index];
        }
        axisRanges[index][0] = center + (min - center) * factor;
        axisRanges[index][1] = center + (max - center) * factor;
        if (limitPan) limitPan(factor > 1.0f);
        computePixelScaling();
    }

    /** Decrease the vertical stretch by 1/stretchFactor times */
    public void verticalShrink(StsMouse mouse)
    {
        stretch(mouse, 1, stretchFactor);
    }

    /** Increase the vertical stretch by stretchFactor times */
    public void horizontalStretch(StsMouse mouse)
    {
        stretch(mouse, 0, 1.0f / stretchFactor);
    }

    /** Decrease the vertical stretch by 1/stretchFactor times */
    public void horizontalShrink(StsMouse mouse)
    {
        stretch(mouse, 0, stretchFactor);
    }

    public void computePixelScaling()
    {
        computeHorizontalPixelScaling();
        computeVerticalPixelScaling();
    }

    public void computeHorizontalPixelScaling()
    {
        if (insets == null) return;
        int axisMinX = insets.left;
        int axisMaxX = getWidth() - insets.right;
        pixelsPerXunit = (axisMaxX - axisMinX) / (axisRanges[0][1] - axisRanges[0][0]);
    }

    public void computeVerticalPixelScaling()
    {
        if (insets == null) return;
        int axisMinY = insets.bottom;
        int axisMaxY = getHeight() - insets.top;
        pixelsPerYunit = (axisMaxY - axisMinY) / (axisRanges[1][1] - axisRanges[1][0]);
    }

    /**
     * Key Released event handling.
     * Up key or A key pressed, right button mouse released increases vertical scale by 2.0.
     * Down key or S key pressed, right button mouse released increases vertical scale by 0.5.
     * Q key pressed, right mouse button released increase view shift factor by 2.0.
     * W key pressed, right mouse button released increase view shift factor by 2.0.
     * Left key changes to crossplot view
     * Right key changes to cursor view
     *
     *
     *
     * @param mouse mouse object
     * @param e     key event
     * @return true if successful
     */
    public void keyReleased(StsMouse mouse, KeyEvent e)
    {
        glPanel3d.restoreCursor();
		super.keyReleased(mouse, e);
    }

    protected void computeModelViewMatrix()
    {
        gl.glMatrixMode(GL.GL_MODELVIEW);
        gl.glLoadIdentity();
        gl.glGetDoublev(GL.GL_MODELVIEW_MATRIX, glPanel.modelViewMatrix, 0);
    }

    /**
     * Set the view parameters
     *
     * @param ranges 4 floats [2][2] (xMin, xMax, yMin, yMax)
     */
    public void setAxesRange(float[][] ranges)
    {
        setHorizontalAxisRange(ranges[0]);
        setVerticalAxisRange(ranges[1]);
    }

    public void setVerticalAxisRange(float[] range)
    {
        axisRanges[1][0] = range[0];
        axisRanges[1][1] = range[1];
    }

    public void setHorizontalAxisRange(float[] range)
    {
        axisRanges[0][0] = range[0];
        axisRanges[0][1] = range[1];
    }

    public float[] getVerticalAxisRange(){ return axisRanges[1]; }

    public float[] getHorizontalAxisRange(){ return axisRanges[0]; }

    public void viewChanged()
    {
        if (glPanel3d == null || glPanel3d.window == null) return;
        StsMouseActionToolbar toolbar = glPanel3d.window.getMouseActionToolbar();
        if (toolbar != null) toolbar.viewChanged(false);
        glPanel3d.viewChanged();
    }

    /*
    * custom serialization requires versioning to prevent old persisted files from barfing.
    * if you add/change fields, you need to bump the serialVersionUID and fix the
    * reader to handle both old & new
    */
    static final long serialVersionUID = 1l;

    public boolean checkSetHorizontalRangeWithScale(StsView2d movedView)
    {
        if (!sameHorizontalAxisType(movedView)) return false;
        pixelsPerXunit = movedView.pixelsPerXunit;
        if (StsMath.sameAs(getHorizontalAxisRange(), movedView.getHorizontalAxisRange())) return false;
        axisRanges[0][0] = movedView.axisRanges[0][0];
        axisRanges[0][1] = axisRanges[0][0] + getInsetWidth() / pixelsPerXunit;
        return true;
    }

    public boolean checkSetVerticalRangeWithScale(StsView2d movedView)
    {
        if (!sameVerticalAxisType(movedView)) return false;
        if (StsMath.sameAs(getVerticalAxisRange(), movedView.getVerticalAxisRange())) return false;
        axisRanges[1][0] = movedView.axisRanges[1][0];
        pixelsPerYunit = movedView.pixelsPerYunit;
        if (pixelsPerYunit == 0.0f) return false;
        axisRanges[1][1] = axisRanges[1][0] + getInsetHeight() / pixelsPerYunit;
        return true;
    }

    public void resetLimits()
    {
        //       computePixelScaling();
        axisRanges[0][1] = axisRanges[0][0] + getInsetWidth() / pixelsPerXunit;
        axisRanges[1][1] = axisRanges[1][0] + getInsetHeight() / pixelsPerYunit;
    }

    public boolean moveWithView(StsView movedView)
    {
        if (movedView instanceof StsView3d) return false;
        boolean moved;
        moved = checkSetVerticalRangeWithScale((StsView2d) movedView);
        moved = moved | checkSetHorizontalRangeWithScale((StsView2d) movedView);
        if (!moved) return false;
        viewChangedRepaint();
        return true;
    }

    public void setMouseMode(int mouseMode)
    {
        if (mouseMode == StsCursor.ROTATE) return;
        super.setMouseMode(mouseMode);
    }

    /**
     * Set the model view parameters
     *
     * @param axisRange 4 floats (min, max display values of both axes)
     * @return true if successful
     */
    public boolean changeModelView2d(float[][] axisRange)
    {
        setDefaultView();
        setAxesRange(axisRange);
        computePixelScaling();
        viewChangedRepaint();
        repaint();
        return true;
    }

    /**
     *
     */
    public void addWeighPoint()
    {
        wp = new StsWeighPoint(model, model.getWindowFamily(this.getWindow()));
        //		clearCurrentKeyPressed();
    }
}
