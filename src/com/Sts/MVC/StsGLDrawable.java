package com.Sts.MVC;

/**
 * Title:        S2S Well Viewer
 * Description:  Well Model-Viewer
 * Copyright:    Copyright (c) 2001
 * Company:      4D Systems LLC
 * @author T.J.Lasseter
 * @version 1.0
 */

//import com.Sts.MVC.*;

import com.Sts.Actions.*;
import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;
import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.nio.*;
import java.util.*;

//import javax.media.opengl.*; import javax.media.opengl.glu.*;

/**
 * An abstract class subclassed from JPanel which contains the GL drawing canvas GLComponent.

 */


public interface StsGLDrawable
{
	abstract public GL getGL();
	abstract public GLU getGLU();
	abstract public float getScreenZ(float [] xyz);
	//abstract public StsActionManager getActionManager();
	abstract public void setViewShift(GL gl, double a);
	abstract public void resetViewShift(GL gl);
	public boolean panelViewChanged = true;
	//abstract public StsView getStsView();
	abstract public int getWidth();
	abstract public int getHeight();
	abstract public void initAndStartGL(GraphicsDevice g);
}
