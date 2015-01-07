package com.Sts.MVC.View3d;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: Class defining two-dimensional cursor view. The cursor is a reference to the
 * three planes that are isVisible in the 3D view. An object based on this class wuold result in a 2D
 * view of any one of the three cursor planes.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DB.*;
import com.Sts.DBTypes.*;

public class StsViewGather2d extends StsViewGather implements StsSerializable
{
    static public final String viewNameGather2d = "View Gather 2d";
    static public final String shortViewNameGather2d = "G2d";

    /**
	 * Default constructor
	 */
	public StsViewGather2d()
	{
    }

	/**
	 * StsPreStackView2d constructor
	 * @param glPanel3d the graphics context
	 */
	public StsViewGather2d(StsGLPanel3d glPanel3d)
	{
		super(glPanel3d, StsPreStackLineSetClass.currentProjectPreStackLineSet);
	}

	public StsViewGather2d(StsGLPanel3d glPanel3d, StsPreStackLineSet2d seismicVolume)
	{
		super(glPanel3d, seismicVolume);
	}

    static public String getStaticViewName()
    {
        return viewNameGather2d;
    }

    static public String getStaticShortViewName()
    {
        return shortViewNameGather2d;
    }

//    public void adjustCursor(int dir, float dirCoor)
//    {
//        //((StsPreStackLineSet2d)lineSet).setDirCoordinate(dir, dirCoor, glPanel3d.window);
//    	System.out.println("StsViewGather2d: adjustCursor doesn't do anything!");
//    }
}
