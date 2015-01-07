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

public class StsViewGather3d extends StsViewGather implements StsSerializable
{
    static public final String viewNameGather3d = "View Gather 3d";
    static public final String shortViewNameGather3d = "G3d";

    /**
	 * Default constructor
	 */
	public StsViewGather3d()
	{
    }

	/**
	 * StsPreStackView2d constructor
	 * @param glPanel3d the graphics context
	 */
	public StsViewGather3d(StsGLPanel3d glPanel3d)
	{
		super(glPanel3d, StsPreStackLineSetClass.currentProjectPreStackLineSet);
	}

	public StsViewGather3d(StsGLPanel3d glPanel3d, StsPreStackLineSet3d seismicVolume)
	{
		super(glPanel3d, seismicVolume);
	}

    static public String getStaticViewName()
    {
        return viewNameGather3d;
    }

    static public String getStaticShortViewName()
    {
        return shortViewNameGather3d;
    }

    public void adjustCursor(int dir, float dirCoor)
    {
        if (glPanel3d == null || lineSet == null) return;
        ((StsPreStackLineSet3d)lineSet).setDirCoordinate(dir, dirCoor, glPanel3d.window);
    }
}
