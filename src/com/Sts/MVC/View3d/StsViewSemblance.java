package com.Sts.MVC.View3d;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: Class to draw and manage the 2D view of semblanceBytes plotted data</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.1
 */

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Types.PreStack.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;
import java.awt.event.*;

abstract public class StsViewSemblance extends StsViewPreStack implements StsTextureSurfaceFace, ItemListener, StsSerializable
{
    public float[][] dataRanges = new float[2][2];
    /** Scaling factors when zooming */
    public float[][] rangeScaling;
    boolean isPixelMode = false;

    /** User horizontal scaling */
    float velocityMin = 5.f;
    float velocityMax = 15.f;

    /** Current prestack volume being viewed */
//	transient public StsPreStackLineSet lineSet = null;
    /** Model prestack volume class which maintains a list of all prestack volume objects */
    //	transient public StsPreStackLineSetClass lineSetClass = null;
    /** Active cursor plane XDIR, YDIR or ZDIR */
    transient public int currentPlaneIndex = -1;
    transient public int nTextureRows, nTextureCols;
    transient int nBackgroundRows, nBackgroundCols;
    //	transient public byte semblanceType = SEMBLANCE_STANDARD; // standard semblanceBytes type
    //    transient public int windowWidth = 0;
    transient public boolean flatten;
    transient StsSemblanceDisplayProperties currentProperties;
    transient int texture = 0;
    //    transient boolean deleteTexture = true;
    /** Tiles on which texture is generated */
    transient public StsTextureTiles textureTiles;
    transient boolean textureChanged = true;
    transient boolean geometryChanged = true;
    /** Display lists should be used (controlled by View:Display Options) */
    transient boolean useDisplayLists = false;
    /** Display lists currently being used for surface geometry */
    transient boolean usingDisplayLists = false;

    transient double[][] semblanceTimes;
    //    transient double[] semblanceVelocities;
    transient float velocityInc;
    transient float zDisplayInc;

    transient int colorDisplayListNum = 0;

    transient boolean runTimer = false;
    transient StsTimer timer;

    transient byte[][] axisData;
    transient int nCursorRows, nCursorCols;

    transient byte[] xData, yData;
    transient StsSpectrum spectrum;

    transient StsPoint pickedPoint;

    static final StsColor muteColor = StsColor.DARK_GREEN;

//	public final String entireVolumeLabel = "Entire PostStack3d";

    static public final int SEMBLANCE_STANDARD = StsSemblanceComputeProperties.SEMBLANCE_STANDARD;

    /** Default constructor */
    public StsViewSemblance()
    {
        limitPan = true;
    }

    public StsViewSemblance(StsGLPanel3d glPanel3d)
    {
        super(glPanel3d);
    }

    public void initialize()
    {
        limitPan = true;
        superGather = lineSet.getSuperGather(glPanel3d.window);
//        currentRow = superGather.superGatherRow;
//        currentCol = superGather.superGatherCol;
        axisRanges = null;
        flatten = this.lineSetClass.getFlatten();
        initializeAxisLabels();
        setAxisRanges();
//        StsMouseActionToolbar toolbar = glPanel3d.window.getMouseActionToolbar();
//       if(toolbar != null) toolbar.zoom();
//		setDefaultView();
    }
    
    public void init(GLAutoDrawable drawable)
    {
        if (isViewGLInitialized) return;
        super.init(drawable);
        viewChanged();
        //isViewGLInitialized = true;
        glPanel3d.glc.addMouseMotionListener(StsVelocityCursor.createVelocityCursorMouseMotionListener(this));
    }

    public void initializeTransients(StsGLPanel3d glPanel3d)
    {
        super.initializeTransients(glPanel3d);
//		  glPanel3d.mouseMode = StsCursor.ZOOM;
        superGather = lineSet.getSuperGather(glPanel3d.window);
//        superGather.initializeSuperGather(currentRow, currentCol);
        isPixelMode = lineSetClass.getIsPixelMode();
        initializeAxisIncs();
        initializeAxisLabels();
        initializeColorList();
//        StsMouseActionToolbar toolbar = glPanel3d.window.getMouseActionToolbar();
//        if(toolbar != null) toolbar.zoom();
        //	semblanceType = line2d.getSemblanceDisplayProperties().semblanceType;

        //resurrect(this);
    }

//    public int getDefaultShader() { return StsJOGLShader.NONE; }
    //public int getDefaultShader() { return StsJOGLShader.ARB_TLUT_NO_SPECULAR_LIGHTS; }
     public int getDefaultShader()
     {
         byte semblanceType = lineSet.semblanceComputeProperties.semblanceType;
         if (semblanceType == StsSemblanceComputeProperties.SEMBLANCE_ENHANCED && StsSuperGather.phaseDisplay)
             return StsJOGLShader.ARB_TLUT_NO_SPECULAR_LIGHTS_ALPHA;
         else
             return StsJOGLShader.ARB_TLUT_NO_SPECULAR_LIGHTS;
     }

    public boolean getUseShader() { return lineSetClass.getContourColors(); }

    private void initializeAxisIncs()
    {
        velocityInc = lineSet.semblanceRangeProperties.velocityStep;
        zDisplayInc = lineSet.semblanceRangeProperties.zInc;
    }

    protected void initializeAxisLabels()
    {
        axisLabels = new String[]{"Velocity", "Time"};
    }

    private void initializeColorList()
    {
//		StsColorList colorList = lineSetClass.semblanceColorList;
//		colorList.addItemListener(this);
    }
/*
	private void initializeGather()
	{
		currentRow = line2d.currentRow;
		currentCol = line2d.currentCol;
		gather.classInitialize(currentRow, currentCol);
		dataRanges[1][0] = line2d.zMax;
		dataRanges[1][1] = line2d.zMin;
		setRange();
		// colorscale & colorscaleDisplayList are static: there is only one for all semblanceBytes displays
		initializeColorList();
//		StsColorscale colorscale = lineSet.getSemblanceColorscale();
//		colorscale.setName("Semblance");
//		colorscale.addItemListener(this);
		//		line2d.addColorscale(colorscale);

		isPixelMode = lineSet.getIsPixelMode();
		isShader = lineSet.getContourColors();

		axisLabels = new String[]
			{"Velocity", "Time"};
	}
*/

    public void setAxisRanges()
    {
        //    StsSemblanceDisplayProperties properties = line2d.getSemblanceDisplayProperties();
        //    semblanceType = properties.semblanceType;

        StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
        totalAxisRanges = new float[][]{{rangeProperties.velocityMin, rangeProperties.velocityMax}, {rangeProperties.zMax, rangeProperties.zMin}};
        axisRanges = StsMath.copyFloatArray(totalAxisRanges);
        dataRanges = StsMath.copyFloatArray(totalAxisRanges);
        velocityInc = rangeProperties.velocityStep;
        zDisplayInc = rangeProperties.zInc;
//        rescaleInchesPerSecond();
//        matchLockedWindows();
//        glPanel3d.viewChanged = true;
    }

    // Tom's version
    public boolean adjustAxisRanges()
    {
        StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
        StsWiggleDisplayProperties wiggleProperties = lineSet.getWiggleDisplayProperties();
        boolean changed = false;
        float velMin = rangeProperties.velocityMin;
        float velMax = rangeProperties.velocityMax;
        boolean velocityRangeChanged = velMin != totalAxisRanges[0][0] || velMax != totalAxisRanges[0][1];
        changed = velocityRangeChanged;
        if(velocityRangeChanged)
        {
            totalAxisRanges[0][0] = velMin;
            totalAxisRanges[0][1] = velMax;
            dataRanges[0][0] = velMin;
            dataRanges[0][1] = velMax;

            float velWidth = axisRanges[0][1] - axisRanges[0][0];
            // amount current left-side must move right to match new display min
            float leftMove = Math.max(0, velMin - axisRanges[0][0]);
            // amount current right-side must move left to match new display max
            float rightMove = Math.max(0, axisRanges[0][1] - velMax);
            //  move the minimum amount and apply the same scaling (same width) unless it goes beyond limit values (velMin and velMax)
            if (leftMove <= rightMove)
            {
                axisRanges[0][0] = velMin;
                axisRanges[0][1] = Math.min(velMin + velWidth, velMax);
                changed = true;
            }
            else
            {
                axisRanges[0][1] = velMax;
                axisRanges[0][0] = Math.max(velMax - velWidth, velMin);
                changed= true;
            }
        }
        if (velocityInc != rangeProperties.velocityStep)
        {
            changed = true;
            velocityInc = rangeProperties.velocityStep;
        }

        float zMin = rangeProperties.zMin;
        float zMax = rangeProperties.zMax;
        boolean zRangeChanged = zMin != totalAxisRanges[1][1] || zMax != totalAxisRanges[1][0];
        changed = changed | zRangeChanged;
        if(zRangeChanged)
        {
            totalAxisRanges[1][1] = zMin;
            totalAxisRanges[1][0] = zMax;
            dataRanges[1][1] = zMin;
            dataRanges[1][0] = zMax;

            float zHeight = zMax - zMin;
            // amount current top must move down to match new display min
            float topMove = Math.min(0, zMin - axisRanges[1][1]);
            // amount current rbottom must move up to match new display max
            float botMove = Math.min(0, axisRanges[1][0] - zMax);
            //  move the minimum amount and apply the same scaling (same width) unless it goes beyond limit values (velMin and velMax)
            if(topMove != 0.0f || botMove != 0.0f)
            {
                if (topMove <= botMove)
                {
                    axisRanges[1][1] = zMin;
                    axisRanges[1][0] = Math.min(zMin + zHeight, zMax);
                    changed = true;
                }
                else
                {
                    axisRanges[1][0] = zMax;
                    axisRanges[1][1] = Math.max(zMax - zHeight, zMin);
                    changed = true;
                }
            }
        }
        if (zDisplayInc != rangeProperties.zInc)
        {
            zDisplayInc = rangeProperties.zInc;
            changed = true;
        }
        if (wiggleProperties.isRangeChanged())
        {
//            if(mainDebug) System.out.println("InchesPerSecond=" + inchesPerSecond + " New inchesPerSecond=" + rangeProperties.inchesPerSecond);
//            inchesPerSecond = wiggleProperties.inchesPerSecond;
            resetRangeWithScale();
            changed = true;
        }
        if (changed) rangeChanged();
        return changed;
    }

    public void resetRangeWithScale()
    {
        rescaleInchesPerSecond();
    }

    protected void rangeChanged()
    {
        textureChanged();
        viewDataChanged();
    }

    public boolean textureChanged()
    {
        textureChanged = true;
        if(glPanel3d != null)
            glPanel3d.viewChanged = true;
        return true;
    }

    public boolean dataChanged()
    {
        if(superGather != null) superGather.semblanceBytes = null;
        return textureChanged();
    }

    public void geometryChanged()
    {
        geometryChanged = true;
    }

    protected void viewDataChanged()
    {
        if (superGather != null) superGather.semblanceBytes = null;
    }
/*
    public boolean adjustAxisRangesStuart()
    {
        StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
        StsWiggleDisplayProperties wiggleProperties = lineSet.getWiggleDisplayProperties();
        boolean changed = false;

        float velMin = rangeProperties.velocityMin;
        if (velMin != totalAxisRanges[0][0])
            changed = true;
        totalAxisRanges[0][0] = velMin;
        axisRanges[0][0] = velMin;
        dataRanges[0][0] = velMin;

        float velMax = rangeProperties.velocityMax;
        if (velMax != totalAxisRanges[0][1])
            changed = true;
        totalAxisRanges[0][1] = velMax;
        axisRanges[0][1] = velMax;
        dataRanges[0][1] = velMax;
        if (velocityInc != rangeProperties.velocityStep)
        {
            changed = true;
            velocityInc = rangeProperties.velocityStep;
        }

        float zMax = rangeProperties.zMax;
        if (zMax != totalAxisRanges[1][0]) changed = true;
        axisRanges[1][0] = zMax;
        totalAxisRanges[1][0] = zMax;
        dataRanges[1][0] = zMax;
        float zMin = rangeProperties.zMin;
        if (zMin != totalAxisRanges[1][1])
            changed = true;
        axisRanges[1][1] = zMin;
        totalAxisRanges[1][1] = zMin;
        dataRanges[1][1] = zMin;
        if (zDisplayInc != rangeProperties.zInc)
        {
            zDisplayInc = rangeProperties.zInc;
            changed = true;
        }
        if (wiggleProperties.isRangeChanged())
        {
//            if(mainDebug) System.out.println("InchesPerSecond=" + inchesPerSecond + " New inchesPerSecond=" + rangeProperties.inchesPerSecond);
//            inchesPerSecond = wiggleProperties.inchesPerSecond;
            rescaleVertical();
            moveLockedWindows();
            glPanel3d.viewPortChanged = true;
            changed = true;
        }
        this.setAxisRanges();
        if (changed)
            glPanel3d.viewChanged = true;

        return changed;
    }
*/
    public void computeProjectionMatrix()
    {
        if (axisRanges == null) return;
        gl.glMatrixMode(GL.GL_PROJECTION);
        gl.glLoadIdentity();
        glu.gluOrtho2D(axisRanges[0][0], axisRanges[0][1], axisRanges[1][0], axisRanges[1][1]);
        gl.glGetDoublev(GL.GL_PROJECTION_MATRIX, glPanel3d.projectionMatrix, 0);
        gl.glMatrixMode(GL.GL_MODELVIEW);
    }

    private boolean isSubVolumeApplied()
    {
        StsSubVolumeClass sv = (StsSubVolumeClass) model.getStsClass(StsSubVolume.class);
        return sv.getIsApplied();
    }

    /**
     * Is the crossplot viewable
     *
     * @return true if viewable
     */
    public boolean isViewable()
    {
        return lineSet != null;
    }

    public boolean initializeDefaultAction()
    {
        StsActionManager actionManager = glPanel3d.actionManager;
        if (actionManager == null) return true;
        setDefaultAction(null);
        return true;
    }

    protected boolean isDisplayData()
    {
        return lineSet.getSemblanceDisplayProperties().displaySemblance;
    }

    protected synchronized boolean checkComputeData()
    {
        if (isCursor3dDragging()) return false;
        if (!isDisplayData()) return false;

        if (superGather.computingSemblance) return true;
        if (!isDataChanged()) return false;
        if (!isDataNull()) return false;
//        if(!textureDataChanged()) return false;
        if (debug) debugThread("computeSemblanceProcess started");
        if (!computeSemblanceData()) return false;
//        glPanel3d.viewPortChanged = true;
        return true;
    }

    protected void checkTextureSizeAndRange()
    {
        if (!textureChanged) return;
        if (debug) debugThread("setting texture range");
        StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
        int nVel = 1 + StsMath.ceiling((rangeProperties.velocityMax - rangeProperties.velocityMin) / rangeProperties.velocityStep);
        nTextureRows = nVel;
        byte semblanceType = lineSet.semblanceComputeProperties.semblanceType;
        if (semblanceType != StsSemblanceComputeProperties.SEMBLANCE_ENHANCED)
        {
            nTextureCols = superGather.nSemblanceSamples;
        }
        else
        {
            semblanceTimes = superGather.semblanceTimes;
            nTextureCols = superGather.nSemblanceSamples;
        }
    }

    protected boolean isDataNull()
    {
        return superGather.semblanceBytes == null;
    }

    /**
     * If semblanceBytes are null, they haven't been computed, so data potentially is available.
     * If we have computed it and it failed, semblanceBytes will be set to a single null byte as a flag;
     * if this is the case, we return true, indicating that compute has failed and data unavailable so
     * that we don't try the compute again.
     *
     * @return true if data is unavailable or the compute thereof has failed.
     */
    protected boolean isDataUnavailable()
    {
        if (super.isDataUnavailable()) return true;
        return superGather.semblanceBytes != null && superGather.semblanceBytes.length == 1;
    }

    protected boolean isDataChanged()
    {
        return textureChanged;
    }

    protected byte[] getData()
    {
        /*
        if (debug)
        {
            int n = 0;
            float scale = 255.0f/nTextureRows;
            for (int row = 0; row < nTextureRows; row++)
            {
                byte debugValue = (byte) (row*scale);
                for (int col = 0; col < nTextureCols; col++, n++)
                    superGather.semblanceBytes[n] = debugValue;
            }
        }
        */
        return superGather.semblanceBytes;
    }

    protected boolean computeSemblanceData()
    {
        return superGather.computeSemblanceProcess(glPanel3d.window);
    }

    public float[] getScaledHorizontalAxisRange()
    {
        // Change velocity to feet or meters per second
        return new float[]{axisRanges[0][0] * 1000, axisRanges[0][1] * 1000};
    }

    protected void exceptionCleanup(GL gl)
    {
        textureTiles.deleteTextures(gl);
    }

    protected void debugThread(String message)
    {
        System.out.println(Thread.currentThread().getName() + " window " + glPanel3d.window.getTitle() + " view class " + StsToolkit.getSimpleClassname(this) + " " + message +
            " textureDataOK " + !isDataNull() + " textureChanged " + textureChanged + " computingSemblance " + superGather.computingSemblance + " currentRow " + superGather.superGatherRow + " currentCol " + superGather.superGatherCol);
        if (superGather.semblanceBytes == null) return;
        int min = Math.min(100, superGather.semblanceBytes.length);
        byte value = 0;
        for (int n = 0; n < min; n++)
            if ((value = superGather.semblanceBytes[n]) != 0) break;
//        System.out.println("    first non-zero semblance byte is " + value);
    }

    /**
     * Output the mouse tracking readout to the information panel on the main screen
     *
     * @param glPanel3d
     * @param mouse mouse object
     */
    public StsCursorPoint logReadout(StsGLPanel3d glPanel3d, StsMouse mouse)
    {
        float[] xyz = new float[3];

        StsCursorPoint cursorPoint = glPanel3d.getCursor3d().getNearestPointInCursorPlane(glPanel3d, mouse);
        if (cursorPoint == null) return null;

        StsPoint point = glPanel3d.getPointInPlaneAtMouse(mouse);
        xyz[0] = glPanel3d.getCursor3d().getCurrentDirCoordinate(StsCursor3d.XDIR);
        xyz[1] = glPanel3d.getCursor3d().getCurrentDirCoordinate(StsCursor3d.YDIR);
        xyz[2] = point.v[1];

        setCursorXOR(glPanel3d, mouse, xyz);

        int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
        if (leftButtonState == StsMouse.PRESSED || leftButtonState == StsMouse.DRAGGED)
            StsMessageFiles.infoMessage(" Velocity= " + "= " + point.v[0] + " " + model.getProject().getVelocityUnits() +
                " " + model.getProject().getZDomainString() + "= " + point.v[1]);

        return cursorPoint;
    }

    private void displayLimitIntervalVelocities(StsGLPanel3d glPanel3d)
    {
        if (currentVelocityProfile == null) return;

        double[] limitVInt = superGather.centerGather.velocities.checkComputeLimitIntervalVelocities(currentVelocityProfile);
        if (limitVInt == null) return;
        try
        {
            GL gl = glPanel3d.getGL();
            gl.glDisable(GL.GL_LIGHTING);
            gl.glLineWidth(4.0f);
            StsColor.BLUE.setGLColor(gl);
            gl.glBegin(GL.GL_LINES);
            int nSlices = lineSet.nSlices;
            double tInc = lineSet.zInc;
            double t1 = lineSet.zMin + tInc;
            for (int n = 0; n < nSlices - 1; n++)
            {
                double t0 = t1;
                t1 += tInc;
                gl.glVertex2d(limitVInt[n], t0);
                gl.glVertex2d(limitVInt[n], t1);
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsViewSemblance.displayIntervalVelocities() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glLineWidth(1.0f);
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    private void displayRMSVelocities(StsGLPanel3d glPanel3d)
    {
        if (currentVelocityProfile == null) return;
        StsPoint[] points = currentVelocityProfile.getProfilePoints();
        if (superGather.centerGather.velocities.vInt == null) return;
        try
        {
            GL gl = glPanel3d.getGL();
            gl.glDisable(GL.GL_LIGHTING);
            gl.glLineWidth(8.0f);
            StsColor.BLUE.setGLColor(gl);
            gl.glBegin(GL.GL_LINE_STRIP);
            double tMin = superGather.centerGather.velocities.tMin;
            double tInc = superGather.centerGather.velocities.tInc;
            double t0 = tMin;
            double[] vi = superGather.centerGather.velocities.vInt;
            int nIntervals = vi.length;
            double vrms = vi[0];
            gl.glVertex2d(vrms, tMin);
            gl.glVertex2d(vrms, t0);
            double vsum = vrms * vrms * (t0 - tMin);
            double t = t0 - tMin;
            for (int n = 1; n < nIntervals; n++)
            {
                vsum += vi[n] * vi[n] * tInc;
                t += tInc;
                vrms = Math.sqrt(vsum / (t - tMin));
                gl.glVertex2d(vrms, t);
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsViewSemblance.displayIntervalVelocities() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glLineWidth(1.0f);
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    private void displayStraightIntervalVelocities(StsGLPanel3d glPanel3d)
    {
        if (currentVelocityProfile == null) return;
        StsPoint[] points = currentVelocityProfile.getProfilePoints();
        int nPoints = points.length;
        if (nPoints < 2) return;
        float[] intervalVelocities = new float[nPoints - 1];
        try
        {
            GL gl = glPanel3d.getGL();
            gl.glDisable(GL.GL_LIGHTING);
            gl.glLineWidth(2.0f);
            double v1 = points[0].v[0];
            double t1 = points[0].v[1];
            for (int n = 0; n < nPoints - 1; n++)
            {
                double v0 = v1;
                double t0 = t1;
                v1 = points[n + 1].v[0];
                t1 = points[n + 1].v[1];
                intervalVelocities[n] = (float) Math.sqrt((v1 * v1 * t1 - v0 * v0 * t0) / (t1 - t0));

            }
            StsColor.BLACK.setGLColor(gl);
            gl.glBegin(GL.GL_LINES);
            for (int n = 0; n < nPoints - 1; n++)
            {
                gl.glVertex2f(intervalVelocities[n], points[n].v[1]);
                gl.glVertex2f(intervalVelocities[n], points[n + 1].v[1]);
            }
            gl.glEnd();
            float vv1 = intervalVelocities[0];
            for (int n = 0; n < nPoints - 2; n++)
            {
                float vv0 = vv1;
                vv1 = intervalVelocities[n + 1];
                float t = points[n + 1].v[1];
                StsGLDraw.drawDottedLine2d(gl, StsColor.WHITE, StsColor.BLACK, 1.0f, vv0, t, vv1, t);
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawLine2d() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glLineWidth(1.0f);
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    private void drawMutes(StsVelocityProfile velocityProfile, GL gl)
    {
        StsPoint mute = velocityProfile.getTopMute();
        if (mute != null) StsGLDraw.drawPoint2d(mute, muteColor, gl, 6);
        mute = velocityProfile.getBottomMute();
        if (mute != null) StsGLDraw.drawPoint2d(mute, muteColor, gl, 6);

    }

    public void setLineSet(StsPreStackLineSet lineSet)
    {
        super.setLineSet(lineSet);
        if (superGather == null) return; // indicates view is not really setup completely
        // superGather.semblanceChanged();
        superGather = null;
        textureChanged();
    }

    protected void semblanceChanged()
    {
        if (superGather == null) superGather = lineSet.getSuperGather(glPanel3d.window);

        superGather.semblanceChanged();
        textureChanged();
    }

    protected boolean semblanceChanged(StsVelocityProfile velocityProfile)
    {
        if (superGather.semblanceChanged(velocityProfile))
        {
            textureChanged();
            return true;
        }
        return false;
    }

    /** filters on gather traces have changed */
    protected void gatherTracesChanged()
    {
        superGather.gatherTracesChanged();
        textureChanged();
    }

    public boolean viewObjectRepaint(Object source, Object object)
    {
        if (debug) System.out.println("Object repaint backbone or semblance");
        if (object instanceof StsPreStackLineSet)
        {
//            currentRow = superGather.superGatherRow;
//            currentCol = superGather.superGatherCol;
//			textureChanged = true;
//			setAxisRanges();  // Causes panned display to maintain zoom but reset to 0,0.
            glPanel3d.repaint();
            return true;
        }
        else if (object instanceof StsVelocityProfile)
        {
            glPanel3d.repaint();
            return true;

        }
        if (object instanceof StsSuperGather)
        {
            glPanel3d.repaint();
            return true;
        }
        else if (object instanceof StsSeismicVelocityModel)
        {
            glPanel3d.repaint();
            return true;

        }
        else if (object == lineSet.semblanceColorList)
        {
            glPanel3d.repaint();
            return true;
        }
        else if (object instanceof StsPanelProperties)  // Catch all property changes.
        {
            glPanel3d.repaint();
            return true;
        }
        return false;
    }

    public void viewChanged()
    {
        textureChanged();
        if(glPanel3d != null)
            glPanel3d.repaint();
    }

    public boolean viewObjectChanged(Object source, Object object)
    {
        if (debug) System.out.println("Object changed semblance");
        if (StsPreStackLineSet.class.isInstance(object))
        {
            if (object != lineSet)
            {
                setLineSet((StsPreStackLineSet) object);
                return true;
            }
            else
                return false;
        }
        if (object instanceof StsSuperGather)
        {
            semblanceChanged();
            return true;
        }
        else if (object instanceof StsVelocityProfile)
        {
            return semblanceChanged((StsVelocityProfile) object);
        }
        else if (object instanceof StsPreStackVelocityModel)
        {
            return true;
        }
        else if (object == lineSet.semblanceColorList)
        {
            dataChanged(); // because we have possibly rescaled, we need to recompute semblance
            return true;
        }
        else if (object instanceof StsWiggleDisplayProperties)
        {
            StsWiggleDisplayProperties wiggleProperties = (StsWiggleDisplayProperties) object;
            if (wiggleProperties.isRangeChanged())
            {
                resetRangeWithScale();
            }
            if (wiggleProperties.stretchMuteChanged)
            {
                // TODO should be able to only adjust the mutes
                semblanceChanged();
            }
            return true;
        }
        else if (object instanceof StsAGCProperties || object instanceof StsFilterProperties)
        {
            gatherTracesChanged();
            semblanceChanged();
            return true;
        }
        else if (object instanceof StsSuperGatherProperties)
        {
            superGather.clearGathers();
            semblanceChanged();

            /*
                if(superGather.superGatherChanged((StsSuperGatherProperties)object))
                {
                    semblanceChanged();
                    return false;
                }
                return true;
            */
            return true;
        }
        else if (object instanceof StsSemblanceRangeProperties)
        {
            superGather.setSemblanceRange();
            adjustAxisRanges();
            textureChanged();
//            deleteTexture = true;
            glPanel3d.viewChanged = true;
            return true;
        }
        else if (object instanceof StsSemblanceComputeProperties)
        {
            // StsViewResidualSemblance overrides this method
            if (checkSemblanceComputeProperties((StsSemblanceComputeProperties) object))
            {
                semblanceChanged();
            }
        }
        else if (object instanceof StsDatumProperties)
        {
            if (superGather.checkSetDatumShift())
            {
                gatherTracesChanged();
                semblanceChanged(currentVelocityProfile);
                return true;
            }
            return false;
        }
        return false;
    }

    //TODO check whether return should be conditional on computeProperties
    protected boolean checkSemblanceComputeProperties(StsSemblanceComputeProperties semblanceComputeProperties)
    {
        return semblanceComputeProperties.isChanged();
    }

    protected void displayData(GL gl, GLU glu)
    {
        StsSemblanceDisplayProperties semblanceDisplayProperties = lineSet.getSemblanceDisplayProperties();
        if (!semblanceDisplayProperties.displaySemblance) return;
//        System.out.println("StsViewSemblance.isDisplayData() called. semblanceBytes " + gather.semblanceBytes + " textureChanged " + textureChanged);
        superGather = lineSet.getSuperGather(glPanel3d.window);

        if (isPixelMode != lineSetClass.getIsPixelMode())
        {
            isPixelMode = !isPixelMode;
            textureChanged();
        }

        // size of textureTiles may have changed above, so check and rebuild if necessary
        if(!checkTextureTiles(gl)) return;

        if(textureTiles.shaderChanged()) textureChanged();
        if(textureTiles.checkBuildDisplayLists(gl, false)) geometryChanged();

        if(textureChanged) textureTiles.deleteTextures(gl);
        if(geometryChanged) textureTiles.deleteDisplayLists(gl);

        try
        {
            gl.glDisable(GL.GL_LIGHTING);
            gl.glEnable(GL.GL_TEXTURE_2D);
            gl.glEnable(GL.GL_BLEND);
//             gl.glTexEnvf( GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE );
            gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
            gl.glShadeModel(GL.GL_FLAT);
//			gl.glDisable(GL.GL_DEPTH_TEST);

            if (lineSet == null || lineSet.semblanceColorList == null) return;
//            if (!lineSet.semblanceColorList.setGLColorList(gl, true, StsJOGLShader.NONE))
            if (!lineSet.semblanceColorList.setGLColorList(gl, true, textureTiles.shader)) return;

            if (debug) System.out.println("draw semblance texture");
       //     if(textureTiles.shader != StsJOGLShader.NONE)
       //         StsJOGLShader.enableARBShader(gl, textureTiles.shader);
            if (textureChanged && !isDataNull())
            {
                if (debug) debugThread("rebuilding new texture");
                textureTiles.displayTiles2d(this, gl, false, isPixelMode, getData(), nullByte);
                textureChanged = false;
//				gather.semblanceBytes = null;
            }
            else
            {
                if (debug)
                    debugThread("displaying current texture");
                textureTiles.displayTiles2d(this, gl, false, isPixelMode, (byte[]) null, nullByte);
            }

//			if(shader != StsJOGLShader.NONE) StsJOGLShader.disableARBShader(gl);
            if (runTimer)
                timer.stopPrint("display semblanceBytes");
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
        finally
        {
            gl.glDisable(GL.GL_TEXTURE_2D);
            gl.glDisable(GL.GL_BLEND);
//			gl.glEnable(GL.GL_DEPTH_TEST);
            if (textureTiles.shader != StsJOGLShader.NONE)
            {
//                System.out.println("SHADER DISABLED.");
                StsJOGLShader.disableARBShader(gl);
            }
            gl.glDisable(GL.GL_TEXTURE_2D);
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    private float[][] getTextureRanges()
    {
        float[][] textureRanges = new float[2][2];
        textureRanges[0][0] = dataRanges[1][1];
        textureRanges[0][1] = dataRanges[1][0];
        textureRanges[1][0] = dataRanges[0][0];
        textureRanges[1][1] = dataRanges[0][1];
        return textureRanges;
    }

    private boolean checkTextureTiles(GL gl)
    {
        if(lineSet == null) return false;
        if(textureTiles != null && !tilesChanged()) return true;

        float[][] textureRanges = new float[2][2];
        textureRanges[0][0] = dataRanges[1][1];
        textureRanges[0][1] = dataRanges[1][0];
        textureRanges[1][0] = dataRanges[0][0];
        textureRanges[1][1] = dataRanges[0][1];
        if(textureTiles != null) deleteTexturesAndDisplayLists(gl);
        textureTiles = StsTextureTiles.constructor(model, this, nTextureRows, nTextureCols, isPixelMode, textureRanges);
        if(textureTiles == null) return false;
        if(debug) System.out.println("TextureTiles created with shader mode " + textureTiles.shader);
        textureChanged();
        geometryChanged();
        return true;
    }

    /**
     * Check if number of texture rows and columns or texture range have changed.
     * Tile rows start in upper left and go down.  Data rows start in lower left and go right.
     * If only the range has changed, we could reset the ranges on existing tiles rather than rebuild,
     * but this hasn't been implemented yet.  TJL 1/25/07.
     */
    private boolean tilesChanged()
    {
        if (textureTiles == null) return true;
        if (textureTiles.nTotalRows != nTextureRows) return true;
        if (textureTiles.nTotalCols != nTextureCols) return true;
        float[][] textureRanges = textureTiles.axisRanges;
        if (textureRanges[0][0] != dataRanges[1][1]) return true;
        if (textureRanges[0][1] != dataRanges[1][0]) return true;
        if (textureRanges[1][0] != dataRanges[0][0]) return true;
        if (textureRanges[1][1] != dataRanges[0][1]) return true;
        return false;
    }

    public StsTextureTiles getTextureTiles() { return textureTiles; }

    public void setTextureTiles(StsTextureTiles textureTiles) {this.textureTiles = textureTiles; }

   /**
     * This puts texture display on delete list.  Operation is performed
     * at beginning of next draw operation.
     */
/*
    public void addTextureToDeleteList()
    {
        if(textureTiles != null)
        {
            StsTextureList.addTextureToDeleteList(this);
        }
//		if(gather != null)gather.semblanceBytes = null;
        textureChanged();
    }
*/
    /** Called to actually delete the displayables on the delete list. */
    public void deleteTexturesAndDisplayLists(GL gl)
    {
        if(textureTiles == null) return;
        textureTiles.deleteTextures(gl);
        deleteDisplayLists(gl);
        textureChanged();
    }
    
    public void drawTextureTileSurface(StsTextureTile tile, GL gl, boolean is3d)
    {
        drawTextureTileSurface2d(tile, gl);
    }

    public void drawTextureTileSurface2d(StsTextureTile tile, GL gl)
    {
        byte semblanceType = lineSet.semblanceComputeProperties.semblanceType;
        tile.drawQuadSurface2d(gl);
    }

    public void deleteDisplayLists(GL gl)
    {
        if (!usingDisplayLists) return;
        usingDisplayLists = false;
        if (textureTiles != null)
            textureTiles.deleteDisplayLists(gl);
    }

    public void itemStateChanged(ItemEvent e)
    {
        if (e.getItem() instanceof StsColorscale) updateColors();
    }

    protected void updateColors()
    {
        textureChanged();
        repaint();
    }

    /**
     * Draws in foreground XORed against current view.  Must be called
     * again and draw same objects to erase.
     */
    public void drawForeground(GL gl)
    {

        gl.glDisable(GL.GL_BLEND);
        gl.glDrawBuffer(GL.GL_FRONT);
        gl.glEnable(GL.GL_COLOR_LOGIC_OP);
        gl.glLogicOp(GL.GL_XOR);
        gl.glDepthMask(false);
        gl.glDepthFunc(GL.GL_ALWAYS);

        drawForegroundCursor();
        gl.glFlush();
        if (cursorButtonState == StsMouse.RELEASED)
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
        if (previousXYZ != null)
            drawCursorPoint(previousXYZ, this.glPanel3d);
        if (currentXYZ != null)
            drawCursorPoint(currentXYZ, glPanel3d);
        previousXYZ = currentXYZ;
        currentXYZ = null;
    }

    /** Draw  horizontal lines thru the cursor point. */
    private void drawCursorPoint(float[] xyz, StsGLPanel3d glPanel3d)
    {
        StsColor.GRAY.setGLColor(gl);
        gl.glLineWidth((float) 3.f);
        gl.glBegin(GL.GL_LINES);
        gl.glVertex2f((float) axisRanges[0][0], xyz[2]);
        gl.glVertex2f((float) axisRanges[0][1], xyz[2]);
        gl.glEnd();

    }

    public void showPopupMenu(StsMouse mouse)
    {
        lineSet.semblancePropertiesDialog(this.glPanel3d.window, mouse);
    }

    /*
      * custom serialization requires versioning to prevent old persisted files from barfing.
      * if you add/change fields, you need to bump the serialVersionUID and fix the
      * reader to handle both old & new
      */
    static final long serialVersionUID = 1l;

    public void rescaleVertical()
    {
        rescaleInchesPerSecond();
        glPanel3d.viewChanged = true;
    }
    
    public void stretch(StsMouse mouse, int index, float factor)
    {
        super.stretch(mouse, index, factor, VERTICAL);
    }

    public void setDefaultView()
    {
        super.setDefaultView(VERTICAL);
        rescaleInchesPerSecond();
    }

    public void adjustCursor(int dir, float dirCoor)
    {
        if(lineSet == null)return;
        if(superGather != null)
        superGather.resetReprocessSemblanceBlockFlag();
        viewChanged();
    }

    public byte getHorizontalAxisType() { return AXIS_TYPE_VELOCITY; }
    public byte getVerticalAxisType() { return AXIS_TYPE_TIME; }
}