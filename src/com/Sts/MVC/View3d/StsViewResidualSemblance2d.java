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

public class StsViewResidualSemblance2d extends StsViewResidualSemblance implements StsTextureSurfaceFace, ItemListener, StsSerializable
{
    static public final String viewNameResidual2d = "Residual Semblance View 2d";
    static public final String shortViewNameResidual2d = "SR2d";

    public StsViewResidualSemblance2d()
	{
    }

    public StsViewResidualSemblance2d(StsGLPanel3d glPanel3d, StsPreStackLineSet2d seismicVolume)
	{
		super(glPanel3d);
        setLineSet(seismicVolume);
		initialize();
    }

	public StsViewResidualSemblance2d(StsGLPanel3d glPanel3d)
	{
		super(glPanel3d);
        setLineSet(StsPreStackLineSetClass.currentProjectPreStackLineSet);
		initialize();
	}


	public void initializeTransients(StsGLPanel3d glPanel3d)
	{
        this.glPanel3d = glPanel3d;
        this.model = glPanel3d.model;
        setLineSet((StsPreStackLineSet)model.getCurrentObject(StsPreStackLineSet2d.class));
		super.initializeTransients(glPanel3d);
	}

    static public String getStaticViewName()
    {
        return viewNameResidual2d;
    }

    static public String getStaticShortViewName()
    {
        return shortViewNameResidual2d;
    }
    /*
	 * custom serialization requires versioning to prevent old persisted files from barfing.
	 * if you add/change fields, you need to bump the serialVersionUID and fix the
	 * reader to handle both old & new
	 */
	static final long serialVersionUID = 1l;
}
