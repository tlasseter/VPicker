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

public class StsViewSemblance3d extends StsViewSemblance implements StsTextureSurfaceFace, ItemListener, StsSerializable
{
    static public final String viewNameSemblance3d = "Semblance View 3d";
    static public final String shortViewNameSemblance3d = "S3d";

    /** Default constructor */
	public StsViewSemblance3d()
	{
        super();
	}

    public StsViewSemblance3d(StsGLPanel3d glPanel3d)
    {
        super(glPanel3d);
        setLineSet((StsPreStackLineSet3d)model.getCurrentObject(StsPreStackLineSet3d.class));
		initialize();
    }

	public void initializeTransients(StsGLPanel3d glPanel3d)
	{
//        super.initializeTransients(glPanel3d);
        this.model = glPanel3d.model;
        this.glPanel3d = glPanel3d;
        setLineSet((StsPreStackLineSet3d)model.getCurrentObject(StsPreStackLineSet3d.class));
		super.initializeTransients(glPanel3d);
	}

    static public String getStaticViewName()
    {
        return viewNameSemblance3d;
    }

    static public String getStaticShortViewName()
    {
        return shortViewNameSemblance3d;
    }
}
