package com.Sts.Interfaces;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

/** subclasses of StsClass which can be isVisible implement this interface */

import com.Sts.MVC.StsGLPanel;
import com.Sts.MVC.View3d.StsGLPanel3d;
import com.Sts.Types.StsMouse;

public interface StsViewSelectable
{
    public void mouseSelectedEdit(StsMouse mouse);
    public void showPopupMenu(StsGLPanel glPanel, StsMouse mouse);

}