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
import com.Sts.Types.PreStack.*;
import com.Sts.Types.PreStack.StsMicroseismicGatherFile;
import com.Sts.Types.*;

public class StsViewGatherMicroseismic extends StsViewGather implements StsSerializable
{

    static public final String viewClassnameMicroGather = "View Gather Microseismic";
    /**
	 * Default constructor
	 */
	public StsViewGatherMicroseismic()
	{
	}

	/**
	 * StsPreStackView2d constructor
	 * @param glPanel3d the graphics context
	 */
    public StsViewGatherMicroseismic(StsGLPanel3d glPanel3d)
    {
        super(glPanel3d);
        //TJL setLineSet(StsPreStackMicroseismicSetClass.currentProjectPreStackMicroseismicSet);
		initialize();
    }
    
    public void setLineSet(StsPreStackMicroseismicSet lineSet)
    {
        //TJL this.lineSet = lineSet;
        //TJL lineSetClass = (StsPreStackMicroseismicSetClass)model.getCreateStsClass(lineSet);
		//TJL superGather = lineSet.getSuperGather(glPanel3d);
    }
    
    public StsViewGatherMicroseismic(StsGLPanel3d glPanel3d, StsPreStackMicroseismicSet seismicVolume)
	{
		super(glPanel3d);
        setLineSet(seismicVolume);
        initialize();
    }

    public String getViewClassname()
    {
        return viewClassnameMicroGather;
    }
    
    public void initializeTransients(StsGLPanel3d glPanel3d)
	{
		super.initializeTransients(glPanel3d);
		//TJL resurrect(this);
        //TJL setLineSet(StsPreStackMicroseismicSetClass.currentProjectPreStackMicroseismicSet);
	}
    /**
	 * Make a copy of the current cursor view
	 * @param glPanel3d the graphics context
	 * @return in new instance of StsPreStackView2d
	 */
	public StsView copy(StsGLPanel3d glPanel3d)
	{
		StsPreStackMicroseismicSet seismicVolume = (StsPreStackMicroseismicSet)model.getCurrentObject(StsPreStackMicroseismicSet.class);
        StsViewGather newView = new StsViewGatherMicroseismic(glPanel3d, seismicVolume);
		newView.copy(this);
		return newView;
	}
	
    public void setLineSet(StsPreStackMicroseismicSet lineSet, StsMicroseismicGatherFile microGatherFile)
    {
        //TJL this.lineSet = lineSet;
        //TJL lineSetClass = (StsPreStackMicroseismicSetClass)model.getCreateStsClass(lineSet);
		//TJL superGather = lineSet.getSuperGather(glPanel3d, microGatherFile);
    }	
}