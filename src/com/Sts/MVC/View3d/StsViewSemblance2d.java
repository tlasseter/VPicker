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

import java.awt.event.*;

public class StsViewSemblance2d extends StsViewSemblance implements StsTextureSurfaceFace, ItemListener, StsSerializable
{
    static public final String viewNameSemblance2d = "Semblance View 2D";
    static public final String shortViewNameSemblance2d = "S2d";

    /** Default constructor */
	public StsViewSemblance2d()
	{
        super();
	}

    public StsViewSemblance2d(StsGLPanel3d glPanel3d)
    {
        super(glPanel3d);
        setLineSet((StsPreStackLineSet)model.getCurrentObject(StsPreStackLineSet2d.class));
		initialize();
    }

    public StsViewSemblance2d(StsGLPanel3d glPanel3d, StsPreStackLineSet2d seismicVolume)
	{
		super(glPanel3d);
        setLineSet(seismicVolume);
        initialize();
    }

	public void initializeTransients(StsGLPanel3d glPanel3d_)
	{
		this.model = glPanel3d_.model;
		this.glPanel3d = glPanel3d_;
        setLineSet((StsPreStackLineSet)model.getCurrentObject(StsPreStackLineSet2d.class));
		super.initializeTransients(glPanel3d_);
	}

    static public String getStaticViewName()
    {
        return viewNameSemblance2d;
    }

    static public String getStaticShortViewName()
    {
        return shortViewNameSemblance2d;
    }

}
