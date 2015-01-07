package com.Sts.Interfaces;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.MVC.View3d.*;

/** subClasses of StsClass which are timeDisplayable (i.e., they may or may not be displayed depending on project time),
 *  implement this interface.
 */
public interface StsClassTimeDisplayable extends StsClassDisplayable
{
    public void displayTimeClass(StsGLPanel3d glPanel3d, long time);
}
