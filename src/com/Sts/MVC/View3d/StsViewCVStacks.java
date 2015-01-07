package com.Sts.MVC.View3d;

import com.Sts.DBTypes.*;
import com.Sts.Types.*;
import com.Sts.UI.*;

import javax.media.opengl.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Apr 8, 2008
 * Time: 1:25:56 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsViewCVStacks extends StsViewStacks
{
    static public final String viewNameCVS = "CVS View";
    static public final String shortViewNameCVS = "CVS";

    public StsViewCVStacks()
	{
	}

    public StsViewCVStacks(StsGLPanel3d glPanel3d)
	{
		super(glPanel3d);
	}

    protected boolean computeStacksData()
    {
        return superGather.computeCVSPanelsProcess(glPanel3d.window);
    }

    protected void stackPanelsChanged()
    {
        superGather.cvsTraces = null;
    }

    protected boolean rangeChanged()
    {
        StsCVSProperties cvsProperties = lineSet.cvsProperties;
        if(nPanels != cvsProperties.numberPanels) return true;
        if(nPanelTraces != cvsProperties.tracesPerPanel) return true;
        StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
        if(totalAxisRanges[1][1] != rangeProperties.zMin) return true;
        if(totalAxisRanges[1][0] != rangeProperties.zMax) return true;

        if(velocityMin != rangeProperties.velocityMin) return true;
        if(velocityMax != rangeProperties.velocityMax) return true;
        
        
        //if(ignoreSuperGather != cvsProperties.ignoreSuperGather) return true;
        return false;
    }
    
    public void init(GLAutoDrawable drawable)
    {
        if (isViewGLInitialized) return;
        super.init(drawable);
        viewChanged();
        //isViewGLInitialized = true;
        glPanel3d.glc.addMouseMotionListener(StsVelocityCursor.createVelocityCursorMouseMotionListener(this));
    }

    protected void recomputeAxisRanges()
    {
        StsCVSProperties cvsProperties = lineSet.cvsProperties;
        nPanels = cvsProperties.numberPanels;
        nPanelTraces = cvsProperties.tracesPerPanel;
        StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
        float zMin = rangeProperties.zMin;
        float zMax = rangeProperties.zMax;

        velocityMin = rangeProperties.velocityMin;
        velocityMax = rangeProperties.velocityMax;
        velocityStep = (velocityMax - velocityMin)/(nPanels-1);

        horizontalTraceSpacing = velocityStep/(nPanelTraces+1);

        float halfPanelPlusMargin = horizontalTraceSpacing*(nPanelTraces+1)/2;
        float xMin =  velocityMin - halfPanelPlusMargin;
        float xMax =  velocityMax + halfPanelPlusMargin;
        totalAxisRanges = new float[][] { { xMin, xMax }, { zMax, zMin } };
        axisRanges = new float[][] { { xMin, xMax }, { zMax, zMin } };
        dataRanges = new float[][] { { xMin, xMax }, { zMax, zMin } };
        originalAxisRanges = new float[][] { { xMin, xMax }, { zMax, zMin } };
        //ignoreSuperGather = lineSet.cvsProperties.ignoreSuperGather;
    }
/*
    protected void displayProfiles()
     {
         StsSemblanceDisplayProperties semblanceDisplayProperties = lineSet.getSemblanceDisplayProperties();
         if(lineSet.velocityModel != null)
            currentVelocityProfile = lineSet.velocityModel.getExistingVelocityProfile(superGather.superGatherRow, superGather.superGatherCol);
         else
            currentVelocityProfile = null;

         // display picked velocity profile
         if (currentVelocityProfile != null)
         {
             StsPoint[] displayPoints = getComputeProfilePoints(currentVelocityProfile);
             StsVelocityProfile.displayOnSemblance(glPanel3d, StsColor.RED, true, false, semblanceDisplayProperties.showLabels, displayPoints, pickedIndex);
         }
     }
*/
	public StsPoint[] getDisplayedProfilePoints(StsVelocityProfile velocityProfile)
	{
		return velocityProfile.getProfilePoints();
	}

    static  public String getStaticViewName()
    {
        return viewNameCVS;
    }

    static public String getStaticShortViewName()
    {
        return shortViewNameCVS;
    }    

    public float[] getScaledHorizontalAxisRange()
    {
        // Change velocity to feet or meters per second
        return new float[]{axisRanges[0][0] * 1000, axisRanges[0][1] * 1000};
    }

    public byte getHorizontalAxisType() { return AXIS_TYPE_VELOCITY_STACKS; }
    public byte getVerticalAxisType() { return AXIS_TYPE_TIME; }

    public Class getDisplayableClass() { return StsPreStackLineSet.class; }

    public boolean viewObjectChanged(Object source, Object object)
    {
        if(object instanceof StsVelocityProfile)
        {
            StsVelocityProfile velocityProfile = (StsVelocityProfile)object;
            if(velocityProfile.changeType == StsVelocityProfile.CHANGE_MUTE)
                gatherTracesChanged();
            return true;
        }
        return super.viewObjectChanged(source, object);
    }
    
    protected boolean  isDataNull() { return superGather.cvsTraces == null; }
    protected float[][][] getData() { return superGather.cvsTraces; }

    public void viewChanged()
    {
        textureChanged();
//        superGather.cvsStackChanged();
    }
}
