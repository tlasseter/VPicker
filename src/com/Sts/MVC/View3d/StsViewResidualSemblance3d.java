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

public class StsViewResidualSemblance3d extends StsViewResidualSemblance implements StsTextureSurfaceFace, ItemListener, StsSerializable
{
    static public final String viewNameResidual3d = "Residual Semblance View 3d";
    static public final String shortViewNameResidual3d = "SR3d";

    public StsViewResidualSemblance3d()
	{
    }

    public StsViewResidualSemblance3d(StsGLPanel3d glPanel3d, StsPreStackLineSet3d seismicVolume)
	{
		super(glPanel3d);
        setLineSet(seismicVolume);
		initialize();
    }

	public StsViewResidualSemblance3d(StsGLPanel3d glPanel3d)
	{
		super(glPanel3d);
        setLineSet(StsPreStackLineSetClass.currentProjectPreStackLineSet);
		initialize();
	}

    static public String getStaticViewName()
    {
        return viewNameResidual3d;
    }

    static public String getStaticShortViewName()
    {
        return shortViewNameResidual3d;
    }

    public void initializeTransients(StsGLPanel3d glPanel3d)
	{
        this.glPanel3d = glPanel3d;
        this.model = glPanel3d.model;
        setLineSet((StsPreStackLineSet3d)model.getCurrentObject(StsPreStackLineSet3d.class));
		super.initializeTransients(glPanel3d);
	}

    /*
      * custom serialization requires versioning to prevent old persisted files from barfing.
      * if you add/change fields, you need to bump the serialVersionUID and fix the
      * reader to handle both old & new
      */
	static final long serialVersionUID = 1l;
}
