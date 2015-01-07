package com.Sts.Interfaces;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

/** subclasses of StsClass which can be displayed implement this interface */

import com.Sts.MVC.View3d.*;

public interface StsClassDisplayable
{
	public void displayClass(StsGLPanel3d glPanel3d);
}