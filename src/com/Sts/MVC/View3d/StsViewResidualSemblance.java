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
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.event.*;
import java.text.*;

abstract public class StsViewResidualSemblance extends StsViewSemblance implements StsTextureSurfaceFace, ItemListener, StsSerializable
{
	transient float percentRangeInc;
    transient boolean computingSemblance = false;
    //transient Runnable computeSemblance = new Runnable() { public void run() { runComputeSemblance(); } };

    static final boolean debug = false;

    /** Default constructor */
	public StsViewResidualSemblance()
	{
		limitPan = true;
	}

	/**
	 * Semblance constructor
	 * @param glPanel3d the graphics context to use for the plot
	 */
	public StsViewResidualSemblance(StsGLPanel3d glPanel3d)
	{
		super(glPanel3d);
    }

    protected void initializeAxisLabels()
	{
		axisLabels = new String[] {"Delta Velocity Percent", "Time"};
	}

	public void initializeTransients(StsGLPanel3d glPanel3d)
	{
		super.initializeTransients(glPanel3d);
		initializeAxisIncs();
	}

	private void initializeAxisIncs()
	{
		percentRangeInc = lineSet.semblanceRangeProperties.percentRangeInc;
		zDisplayInc = lineSet.semblanceRangeProperties.zInc;
	}

	public void setAxisRanges()
	{
		StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
		float percentRange = rangeProperties.percentRange;
		percentRangeInc = rangeProperties.percentRangeInc;
		int nVel = 1 + 2*Math.round(percentRange/percentRangeInc);
		this.nTextureRows = nVel;
		float zMin = rangeProperties.zMin;
		float zMax = rangeProperties.zMax;
		totalAxisRanges = new float[][] { {-percentRange, percentRange}, {zMax, zMin} };
		axisRanges = StsMath.copyFloatArray(totalAxisRanges);
		dataRanges = StsMath.copyFloatArray(totalAxisRanges);
//        rescaleInchesPerSecond();
//		glPanel3d.viewChanged = true;
	}

    public boolean adjustAxisRanges()
	{
		StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
        StsWiggleDisplayProperties wiggleProperties = lineSet.getWiggleDisplayProperties();
        boolean changed = false;
        float percentRange = rangeProperties.percentRange;
        boolean percentRangeChanged = percentRange != -totalAxisRanges[0][0] || percentRange != totalAxisRanges[0][1];
        changed = percentRangeChanged;
        if(percentRangeChanged)
        {
            glPanel3d.viewChanged = true;
            totalAxisRanges[0][0] = -percentRange;
            totalAxisRanges[0][1] = percentRange;
            dataRanges[0][0] = -percentRange;
            dataRanges[0][1] = percentRange;

            float velWidth = axisRanges[0][1] - axisRanges[0][0];
            // amount current left-side must move right to match new display min
            float leftMove = Math.min(0, -percentRange - axisRanges[0][0]);
            // amount current right-side must move left to match new display min
            float rightMove = Math.min(0, axisRanges[0][1] - percentRange);
            //  move the minimum amount and apply the same scaling (same width) unless it goes beyond limit values (velMin and velMax)

            if(leftMove <= rightMove)
            {
                axisRanges[0][0] = -percentRange;
                axisRanges[0][1] = Math.min(-percentRange + velWidth, percentRange);
                changed = true;
            }
            else
            {
                axisRanges[0][1] = percentRange;
                axisRanges[0][0] = Math.max(percentRange - velWidth, -percentRange);
                changed = true;
            }
        }
        if(percentRangeInc != rangeProperties.percentRangeInc)
         {
             changed = true;
             percentRangeInc = rangeProperties.velocityStep;
         }

        float zMin = rangeProperties.zMin;
        float zMax = rangeProperties.zMax;
        boolean zRangeChanged = zMin != totalAxisRanges[1][1] || zMax != totalAxisRanges[1][0];
        changed = changed | zRangeChanged;
        if(zRangeChanged)
        {
            glPanel3d.viewChanged = true;

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
                if(topMove <= botMove)
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
        if(zDisplayInc != rangeProperties.zInc)
		{
			zDisplayInc = rangeProperties.zInc;
			changed = true;
		}
        if(wiggleProperties.isRangeChanged())
        {
//            if(mainDebug) System.out.println("InchesPerSecond=" + inchesPerSecond + " New inchesPerSecond=" + rangeProperties.inchesPerSecond);
//            inchesPerSecond = wiggleProperties.inchesPerSecond;
            resetRangeWithScale();
            changed = true;
		}
        return changed;
	}

    public float[] getScaledHorizontalAxisRange()
    {
        return new float[]{axisRanges[0][0], axisRanges[0][1]};
    }


    protected boolean isDataNull()
    {
        return superGather.residualSemblance == null;
    }

    public void resetRangeWithScale()
    {
        rescaleInchesPerSecond();
    }

   /** If semblanceBytes are null, they haven't been computed, so data potentially is available.
     *  If we have computed it and it failed, semblanceBytes will be set to a single null byte as a flag;
     *  if this is the case, we return true, indicating that compute has failed and data unavailable so
     *  that we don't try the compute again.
     * @return true if data is unavailable or the compute thereof has failed.
     */
    protected boolean isDataUnavailable()
    {
        if(super.isDataUnavailable()) return true;
        return superGather.residualSemblance != null && superGather.residualSemblance.length == 1;
    }

    protected byte[] getData()
    {
        return superGather.residualSemblance;
    }

    protected boolean computeSemblanceData()
    {
        return superGather.computeResidualSemblanceProcess(glPanel3d.window);
    }

    protected void displayProfiles()
    {
        StsSemblanceDisplayProperties semblanceDisplayProperties = lineSet.getSemblanceDisplayProperties();
        currentVelocityProfile = getVelocityProfile(superGather.superGatherRow, superGather.superGatherCol);

        // display picked velocity profile
        if (currentVelocityProfile != null)
        {
            StsPoint[] displayPoints = getDisplayedProfilePoints(currentVelocityProfile);
            StsVelocityProfile.displayOnSemblance(glPanel3d, StsColor.RED, true, false, semblanceDisplayProperties.showLabels, displayPoints, pickedIndex);
        }
        else
            currentVelocityProfile = null;
    }

    public StsCursorPoint logReadout(StsGLPanel3d glPanel3d, StsMouse mouse)
    {
        float[] xyz = new float[3];

        StsCursorPoint cursorPoint = glPanel3d.getCursor3d().getNearestPointInCursorPlane(glPanel3d, mouse);
        if(cursorPoint == null) return null;

        StsPoint point = glPanel3d.getPointInPlaneAtMouse(mouse);
        xyz[0] = glPanel3d.getCursor3d().getCurrentDirCoordinate(StsCursor3d.XDIR);
        xyz[1] = glPanel3d.getCursor3d().getCurrentDirCoordinate(StsCursor3d.YDIR);
        xyz[2] = point.v[1];

        setCursorXOR(glPanel3d, mouse, xyz);

        int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
        if(leftButtonState == StsMouse.PRESSED || leftButtonState == StsMouse.DRAGGED)
            StsMessageFiles.infoMessage(" Velocity= " + "= " + point.v[0] + " " + model.getProject().getVelocityUnits() +
                                        " " + model.getProject().getZDomainString() + "= " + point.v[1]);
		return cursorPoint;
    }
/*
	public void displayOnSemblance(StsVelocityProfile velocityProfile, StsGLPanel3d glPanel3d, StsColor stsColor, boolean drawVertices, boolean active)
	{
		StsPoint[] displayPoints = computeSemblancePoints(velocityProfile);
		displayOnSemblance(glPanel3d, stsColor, displayPoints, drawVertices, active);
    }
*/
 	public void drawTextureTileSurface2d(StsTextureTile tile, GL gl)
	{
        byte semblanceType = lineSet.semblanceComputeProperties.semblanceType;
        if(semblanceType == SEMBLANCE_STANDARD)
			tile.drawQuadSurface2d(gl);
		else
        {
            StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
            float percentMin = -rangeProperties.percentRange;
            float percentInc = rangeProperties.percentRangeInc;
            tile.drawQuadStripSurface2d(gl, percentMin, percentInc, semblanceTimes);
        }
    }

    void drawLabels(StsPoint[] displayPoints, GL gl)
	{
		if(currentVelocityProfile == null) return;
		StsPoint[] semblancePoints = currentVelocityProfile.getProfilePoints();
		DecimalFormat labelFormat = new DecimalFormat("#,###.#");
		int nPoints = semblancePoints.length;
		for (int n = 0; n < nPoints; n++)
		{
			String text = new String(labelFormat.format(semblancePoints[n].v[0] * 1000) + ", " + labelFormat.format(semblancePoints[n].v[1]));
			StsGLDraw.fontHelvetica12(gl, displayPoints[n].getXYZorT(), text);
		}
	}

	public StsPoint[] getDisplayedProfilePoints(StsVelocityProfile velocityProfile)
	{
		StsPoint[] profilePoints = velocityProfile.getProfilePoints();
		int nPoints = profilePoints.length;
        StsPoint[] displayPoints = new StsPoint[nPoints];
        StsPoint[] initialProfilePoints = velocityProfile.getInitialProfilePoints();
        if(initialProfilePoints == null)
        {
            for(int n = 0; n < nPoints; n++)
            {
                float time = profilePoints[n].v[1];
                displayPoints[n] = new StsPoint(0.0f, time);
            }
        }
        else
        {
            for(int n = 0; n < nPoints; n++)
            {
                float velocity = profilePoints[n].v[0];
                float time = profilePoints[n].v[1];
                float backboneVelocity = StsMath.interpolateValue(initialProfilePoints, time, 1, 0);
                float percent = 100.0f * (velocity / backboneVelocity - 1.0f);
                displayPoints[n] = new StsPoint(percent, time);
            }
        }
        return displayPoints;
	}
 
    protected void semblanceChanged()
    {
        superGather.residualSemblanceChanged();
        textureChanged();
    }

    protected boolean semblanceChanged(StsVelocityProfile velocityProfile)
    {
        if(superGather.residualSemblanceChanged(velocityProfile))
        {
            textureChanged();
            return true;
        }
        return false;
    }


    protected boolean checkSemblanceComputeProperties(StsSemblanceComputeProperties semblanceComputeProperties)
    {
        byte order = semblanceComputeProperties.order;
        textureChanged();
        return superGather.checkSetOrder(order);
	}
/*
	protected void isDisplayData(GL gl, GLU glu)
	{
		if (isPixelMode != lineSetClass.getIsPixelMode())
		{
			isPixelMode = !isPixelMode;
			deleteTexture = true;
			textureChanged = true;
		}
		// size of textureTiles may have changed above, so check and rebuild if necessary
		if (!checkTextureTiles(gl))
			return;

		try
		{
			gl.glDisable(GL.GL_LIGHTING);
			gl.glEnable(GL.GL_TEXTURE_2D);
			gl.glEnable(GL.GL_BLEND);
			//        gl.glTexEnvf( GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE );
			gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
			gl.glShadeModel(GL.GL_FLAT);
//			gl.glDisable(GL.GL_DEPTH_TEST);

			if (lineSet == null || lineSet.semblanceColorList == null)
				return;
			if (!lineSet.semblanceColorList.setGLColorList(gl, true, shader))
				return;

			useDisplayLists = StsObject.getCurrentModel().getBooleanProperty("Use Display Lists");
			if (!useDisplayLists && usingDisplayLists)
			{
				deleteDisplayLists(gl);
			}
			else if (useDisplayLists && !usingDisplayLists)
			{
				if (textureTiles == null)
					StsException.systemError("StsSurface.displaySurfaceFill() failed. textureTiles should not be null.");
				textureTiles.constructSurface(this, gl, useDisplayLists, false);
				usingDisplayLists = true;
			}

			if (debug) System.out.println("draw semblanceBytes texture");
			if (shader != StsJOGLShader.NONE)
				StsJOGLShader.enableARBShader(gl, shader);
			if (textureChanged && gather.semblanceBackbone != null)
			{
				textureTiles.displayTiles2d(this, gl, false, isPixelMode, gather.semblanceBackbone);
				textureChanged = false;
				//gather.semblanceBackbone = null;
				//debugThread("isDisplayData() set textureChange false");
			}
			else
				textureTiles.displayTiles2d(this, gl, false, isPixelMode, (byte[])null);

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
			if (shader != StsJOGLShader.NONE)
				StsJOGLShader.disableARBShader(gl);
			gl.glDisable(GL.GL_TEXTURE_2D);
			gl.glEnable(GL.GL_LIGHTING);
		}
	}
*/
    /*
    private void runComputeSemblance()
	{
		try
		{
            if(debug) debugThread("runComputeSemblance() called");
            if (gather.computeBackboneSemblance(currentVelocityProfile, lineSet, lineSetClass, glPanel3d.window))
			{
				initializeTextureSizeAndRange();
                gather.setComputingGather(false);
                model.viewObjectRepaintFamily(glPanel3d, lineSet);
                gather.setComputingGather(false);
            }
		}
		catch(Exception e)
		{
            StsException.outputWarningException(this, "runComputeSemblance", e);
            gather.setComputingGather(false);
		}
	}*/
/*
    private void runComputeSemblance()
	{
        debugThread("runComputeSemblance() called");
        if(debug) System.out.println("runComputeSemblance() called");
        computingGather = false;
		textureChanged = false;
		if(!gather.computeBackboneSemblance(currentVelocityProfile, lineSet, lineSetClass, glPanel3d.window)) return;
		textureChanged = true;
		initializeTextureSizeAndRange();
        computingGather = false;
    }
*/
	protected void checkTextureSizeAndRange()
	{
        if(!textureChanged) return;
        if(debug) debugThread("setting texture range");
        StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
		float percentRange = rangeProperties.percentRange;
		percentRangeInc = rangeProperties.percentRangeInc;
		float zMin = rangeProperties.zMin;
		float zInc = rangeProperties.zInc;
		int nVel = 1 + 2 * Math.round(percentRange / percentRangeInc);
		nTextureRows = nVel;
        byte semblanceType = lineSet.semblanceComputeProperties.semblanceType;
        if(semblanceType == SEMBLANCE_STANDARD)
		{
			nTextureCols = superGather.nSemblanceSamples;
			dataRanges[1][0] = (float)superGather.semblanceZMax;
			dataRanges[1][1] = (float)superGather.semblanceZMin;
		}
		else
		{
            semblanceTimes = superGather.semblanceTimes;
            nTextureCols = superGather.nEnhancedSemblancePoints;
		}
	}

	/** This puts texture display on delete list.  Operation is performed
	 *  at beginning of next draw operation.
	 */
/*
    public void clearTextureTileSurface(StsGLPanel glPanel)
	{
		if(textureTiles != null)
		{
			glPanel.deleteTextureTileSurface(this);
		}
		if(superGather != null) superGather.semblanceBackbone = null;
		textureChanged = true;
	}
 */

	public StsPoint computePickPoint(StsMouse mouse)
	{
		StsPoint pick = super.computePickPoint(mouse);
		if(pick == null) return null;
//		if(currentVelocityProfile == null) return null;
		float time = pick.v[1];
		float percent = pick.v[0];
		StsPoint[] initialProfilePoints = currentVelocityProfile.getInitialProfilePoints();
		float velocity = StsMath.interpolateValue(initialProfilePoints, time, 1, 0);
		pick.v[0] = velocity*(1.0f + percent/100);
		return pick;
	}
	/*
	 * custom serialization requires versioning to prevent old persisted files from barfing.
	 * if you add/change fields, you need to bump the serialVersionUID and fix the
	 * reader to handle both old & new
	 */
	static final long serialVersionUID = 1l;

    public void adjustCursor(int dir, float dirCoor)
    {
        if(lineSet == null)return;
        superGather.resetReprocessResidualBlockFlag();
        textureChanged();
    }

    public byte getHorizontalAxisType() { return AXIS_TYPE_VELOCITY_PERCENT; }
    public byte getVerticalAxisType() { return AXIS_TYPE_TIME; }
}
