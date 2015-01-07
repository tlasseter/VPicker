package com.Sts.MVC.View3d;

import com.Sts.DBTypes.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Apr 8, 2008
 * Time: 1:25:56 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsViewVVStacks extends StsViewStacks
{
    static public final String viewNameVVS = "VVS View";
    static public final String shortViewNameVVS = "VVS";

    public StsViewVVStacks()
	{
	}

    public StsViewVVStacks(StsGLPanel3d glPanel3d)
	{
		super(glPanel3d);
    }

    public void init(GLAutoDrawable drawable)
    {
        if (isViewGLInitialized) return;
        super.init(drawable);
        StsVvsUpdateToolbar.checkAddToolbar(glPanel3d.window);
        viewChanged();
        //isViewGLInitialized = true;
        
        glPanel3d.glc.addMouseMotionListener(StsVelocityCursor.createVelocityCursorMouseMotionListener(this));
    }

    protected boolean isDataUnavailable()
    {
        if(super.isDataUnavailable()) return true;
        if(lineSet.velocityModel == null) return true;
        if(lineSet.velocityModel.getCurrentVelocityProfile(getWindow()) == null) return true;
        return !lineSet.hasVelocityProfiles() && lineSet.velocityModel.getCurrentVelocityProfile(getWindow()).getProfilePoints().length > 1;  //Tom??? 1 point isn't really velocity profile - this causes NPE for 1st point in VVS on new project;
    }

    protected boolean computeStacksData()
    {
        return superGather.computeVVSPanelsProcess();
    }

    protected void stackPanelsChanged()
    {
        superGather.vvsTraces = null;
    }

    protected boolean rangeChanged()
    {
        StsCVSProperties cvsProperties = lineSet.cvsProperties;
        if(nPanels != cvsProperties.numberPanels) return true;
        if(nPanelTraces != cvsProperties.tracesPerPanel) return true;
        StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
        if(totalAxisRanges[1][1] != rangeProperties.zMin) return true;
        if(totalAxisRanges[1][0] != rangeProperties.zMax) return true;
        if(velocityStep != rangeProperties.vvsPercentRangeInc) return true;
        return false;
    }

    //TODO this is called multiple times on start up: figure out why and fix
    protected void recomputeAxisRanges()
    {
        StsCVSProperties cvsProperties = lineSet.cvsProperties;
        nPanels = cvsProperties.getNumberPanels();
        nPanelTraces = cvsProperties.getTracesPerPanel();
        StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
        float zMin = rangeProperties.zMin;
        float zMax = rangeProperties.zMax;

        velocityStep = rangeProperties.vvsPercentRangeInc;
        velocityMin = -(nPanels/2)*velocityStep;
        velocityMax = (nPanels/2)*velocityStep;

        horizontalTraceSpacing = velocityStep /(nPanelTraces+1);

        float halfPanelPlusMargin = horizontalTraceSpacing*(nPanelTraces+1)/2;
        float xMin =  velocityMin - halfPanelPlusMargin;
        float xMax =  velocityMax + halfPanelPlusMargin;
        totalAxisRanges = new float[][] { { xMin, xMax }, { zMax, zMin } };
        axisRanges = new float[][] { { xMin, xMax }, { zMax, zMin } };
        dataRanges = new float[][] { { xMin, xMax }, { zMax, zMin } };
        originalAxisRanges = new float[][] { { xMin, xMax }, { zMax, zMin } };
    }

	public StsPoint[] getDisplayedProfilePoints(StsVelocityProfile velocityProfile)
	{
		StsPoint[] profilePoints = velocityProfile.getProfilePoints();
		int nPoints = profilePoints.length;
        StsPoint[] displayPoints = new StsPoint[nPoints];
        StsPoint[] vvsInitialProfilePoints = velocityProfile.getSetVvsInitialProfilePoints();
        if(vvsInitialProfilePoints == null)
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
                float backboneVelocity = 0;
                if (time < vvsInitialProfilePoints[0].v[1]) 
                    backboneVelocity = vvsInitialProfilePoints[0].v[0];
                else
                    backboneVelocity = StsMath.interpolateValue(vvsInitialProfilePoints, (float)time, 1, 0);
                float percent = 100.0f * (velocity / backboneVelocity - 1.0f);
                displayPoints[n] = new StsPoint(percent, time);
            }
        }
        return displayPoints;
	}

	public StsPoint computePickPoint(StsMouse mouse)
	{
		StsPoint pick = super.computePickPoint(mouse);
		if(pick == null) return null;
		float time = pick.v[1];
		float percent = pick.v[0];
		if (currentVelocityProfile == null) currentVelocityProfile = superGather.velocityProfile;  //allows StsVelocityCursor to work before first pick is made (uses interpolated velocity)
		if (currentVelocityProfile == null) return null;
		StsPoint[] vvsInitialProfilePoints = currentVelocityProfile.getSetVvsInitialProfilePoints();
		if (vvsInitialProfilePoints == null || vvsInitialProfilePoints.length == 0) return null;
		float velocity = 0;
        if (time < vvsInitialProfilePoints[0].v[1]) 
            velocity = vvsInitialProfilePoints[0].v[0];
        else
            velocity = StsMath.interpolateValue(vvsInitialProfilePoints, (float)time, 1, 0);
		pick.v[0] = velocity*(1.0f + percent/100);
		return pick;
	}

    /** Superclass will make changes to appropriate objects, but will not delete the VVS view.
     *  If however we have edited a profile (object == StsVelocityProfile) on an StsViewVVStacks view or
     *  we have pushed the vvsUpdate button, then delete and recompute the view.
     *  we want don't want to update the object until the user pushes the VVS Update button.
     *  If any other object has been changed, e.g. AGCProperties, we do update the object.
     *  If the "VVS Update" button has been pushed (source == StsVvsUpdateToolbar, object == StsVvsUpdateToolbar),
     *  update the view as if the object was StsVelocityProfile.
     *
     * @param source
     * @param object
     * @return
     */

    public boolean viewObjectRepaint(Object source, Object object)
	{
        if(source instanceof StsVvsUpdateToolbar)
        {
			glPanel3d.repaint();
			return true;
        }
        return super.viewObjectRepaint(source, object);
    }

    public void viewChanged()
    {
        if(superGather == null) return;
        currentVelocityProfile = getVelocityProfile(superGather.superGatherRow,  superGather.superGatherCol);
        if(currentVelocityProfile == null) return;
        currentVelocityProfile.clearVvsInitialProfilePoints();
        if(superGather != null) superGather.vvsStackChanged();
        textureChanged();
    }

    static public String getStaticViewName()
    {
        return viewNameVVS;
    }

    static public String getStaticShortViewName()
    {
        return shortViewNameVVS;
    }

    public byte getHorizontalAxisType() { return AXIS_TYPE_VELOCITY_PERCENT_STACKS; }
    public byte getVerticalAxisType() { return AXIS_TYPE_TIME; }

    public boolean viewObjectChanged(Object source, Object object)
    {
        if(source instanceof StsVvsUpdateToolbar)
        // if((source instanceof StsViewVVStacks && object instanceof StsVelocityProfile) || source instanceof StsVvsUpdateToolbar)
        {
            viewChanged();
            return true;
        }
        else if(object instanceof StsVelocityProfile)
        {
            StsVelocityProfile velocityProfile = (StsVelocityProfile)object;
            superGather.velocityProfileChanged(velocityProfile);
            if(velocityProfile.changeType == StsVelocityProfile.CHANGE_MUTE)
                gatherTracesChanged();
            if(lineSet.cvsProperties.getAutoStackUpdate()) viewChanged();
            return true;
        }
        return super.viewObjectChanged(source, object);
    }

    protected boolean  isDataNull() { return superGather.vvsTraces == null; }
    protected float[][][] getData() { return superGather.vvsTraces; }
}
