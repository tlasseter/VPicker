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
 * <p/>
 * After subclassing from StsGLPanel3d, implement the abstract methods:
 * <ul>
 * <li>{@link #computeProjectionMatrix()}</li>
 * <li>{@link #keyPressed(KeyEvent,StsMouse)}</li>
 * <li>{@link #keyReleased(KeyEvent,StsMouse)}</li>
 * <li>{@link #moveWindow(StsMouse)}</li>
 * <li>{@link #showPopupMenu(com.Sts.Types.StsMouse)}</li>
 * </ul>
 * <p/>
 * GLPanel construction and display <em> MUST </em> be done on the AWT Event thread.
 * This is handled for you, but you should be aware of it.
 * <p/>
 * Mouse and key events will come to the GLComponent.  The actionManager is the listener for these
 * events and handles them accordingly.
 * <p/>
 * Starting and running GL must follow this procedure:
 * <ul>
 * <li>construct glPanel</li>
 * <li>complete all window construction</li>
 * <li>call glPanel.startGL() </li>
 * </ul>
 * <p/>
 * For details on GL commands, see {@linkplain "OpenGL Programming Guide, Addison-Wesley"}
 */


public class StsGLJPanel extends StsGLPanel implements StsSerializable, GLEventListener, KeyListener, MouseListener, MouseMotionListener, FocusListener, MouseWheelListener
// abstract public class StsGLPanel extends GLJPanel implements GLEventListener, KeyListener, StsSerializable
{
	//transient public GLJPanel glc = null;

   /** default constructor: not normally used */
    public StsGLJPanel()
    {
		super();

		initAndStartGL();
		setLayout(null);

    }

    public StsGLJPanel(StsModel model, StsActionManager actionManager, int width, int height, StsView view)
    {
        this(model, actionManager);
        this.view = view;
		setPanelSize(width, height);
    }

    public StsGLJPanel(StsModel model, StsActionManager actionManager)
    {
		super();
        initializeTransients(model, actionManager);
		initAndStartGL();
		setLayout(null);

    }



    public void initAndStartGL()
    {
         if(model == null)
        {
            StsException.systemError(this, "initAndStartGL", "model is null.");
            return;
        }
        if(state == GL_STATE_STARTED)
        {
            //paint(getGraphics());
            return;
        }
        else if(state != GL_STATE_NONE)
            return;

        try
        {
			GLContext shareableContext = Main.sharedContext;
			glc = new GLJPanel(glCapabilities, null, shareableContext); //g);
			gld = (GLAutoDrawable) glc;
			//shareableContext.makeCurrent();
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "initGL", e);
        }
        // sets GL which might be debug or trace GL
        //setGL();
		// start GL by adding GL listener
		//System.out.println("gljpanel add "+this);
        glu = new GLU();
		addGLlisteners();
        // add panel and initialize panel size
		add(glc);
		initializeSize();



        state = GL_STATE_INITIALIZED;
       // if(glStateDebug) System.out.println("GL initialized and started.");
       // if(Main.isGLDebug || Main.isGLTrace) System.out.println("GL STATE: STARTED");
		//paint(getGraphics()); // jbw jdk 1.6.11 and above -- force initial paint
    }

    public void addNotify()
	{
		//System.out.println("add notify glj \n"+this);
		super.addNotify();
		validate();
		if (glc != null)
		{
			Dimension D = getSize();
			glc.setSize(D);
			//System.out.println("set size glj addnot"+D);
			glc.paint(getGraphics());
		}
    }


    /** seems to be required for proper splitPane behavior when a glPanel is in a splitPane */
    protected void initializeSize()
    {
        setMinimumSize(size00);
        Dimension size = getSize();
        //setPreferredSize(size);
		//System.out.println("initializesize "+size);
		//if (glc != null)
        //glc.setPreferredSize(size);
    }


    public void init(GLAutoDrawable drawable)
    {
		setGL();
		String renderer = gl.glGetString(GL.GL_RENDERER);
		if(StsStringUtils.stringContainsString(renderer, "Mesa"))
			throw new RuntimeException("Missing Linux Accelerated 3D OpenGL drivers");
		if(StsStringUtils.stringContainsString(renderer, "GDI Gen"))
			throw new RuntimeException("Missing Windows Accelerated 3D OpenGL drivers");
        //System.out.println("Debug: glJPanel "+renderer);
        if(view != null) view.init(drawable);
		state = GL_STATE_STARTED;
		initAndStartGL();
		initializeSize();
		repaint();
    }


    public void setPanelSize(int width, int height)
    {
        setPanelSize(new Dimension(width, height));
    }

    public void setPanelSize(Dimension size)
    {
        this.size = size;
        setMinimumSize(size00);
        setSize(size);
        setPreferredSize(size);
		//System.out.println("setpanelsize "+size);
    }

    public void setBounds(int x, int y, int width, int height)
    {
       size.width = width;
       size.height = height;
       super.setBounds(x, y, width, height);
	   //glc.setBounds(x,y,width,height);
	   //System.out.println("setbounds glj glc="+glc+" "+width+" "+height);

	   //repaint();
    }


    // Component initialization from DB
    public void initializeTransients(StsModel model, StsActionManager actionManager)
    {
        try
        {
            this.model = model;
            this.actionManager = actionManager;
            cursor = new StsCursor(this);
            setLayout(new BorderLayout());
            setBorder(null);
            if(model != null) setClearColor(model.getProject().getBackgroundStsColor());
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "initializeTransients", e);
        }
    }


    public void paint(Graphics g)
	{
	   // System.out.println(" GLPanel paint "+glc);
	   super.paint(g);
	   //super.paintChildren(g);
	   if (glc!=null)
       {

		   glc.setVisible(true); // jbw hack for maximize bug

		   glc.paint(g);

		   int count=0;
		   while (glc.getWidth() <=0 && count < 10) // peer not realized
		   {
			   //System.out.println("glj paintwidth 0"+this.getSize());
			   glc.setSize(this.getSize());
			   glc.paint(g);
			   count++;
		   }
	   }
       else
           StsException.systemDebug(this, "paint", "glc is null");
    }

      /** called when a window setBounds event has occurred */
    public void reshape(GLAutoDrawable drawable, int x, int y, int width, int height)
    {
        int oldWidth = viewPort[2];
        int oldHeight = viewPort[3];
        viewPortChanged = oldWidth != width || oldHeight != height;
        if(Main.isDrawDebug) System.out.println("viewPort changed: " + viewPortChanged);
        viewPort[0] = x;
        viewPort[1] = y;
        viewPort[2] = width;
        viewPort[3] = height;
        if(!viewPortChanged)
        {
            if(!model.displayOK) model.enableDisplay();
            return;
        }
        if(model.displayOK) model.disableDisplay();

        setViewPort(x, y, width, height);
        //glc.reshape(x, y, width, height);

        if(view != null)
        {
            view.reshape(x, y, width, height);
            view.computeProjectionMatrix();
        }
        viewPortChanged = false;
        model.enableDisplay();
		glc.setVisible(false); // jbw hack! for maximize bug
		glc.reshape(x, y, width, height);
        //glc.setVisible(true); // jbw no go
		repaint();
    }



    private void setGL()
    {
        gl = gld.getGL();
        if(Main.isGLDebug)
        {
            gl = new DebugGL(gl);
            gld.setGL(gl);
        }

        else if(Main.isGLTrace)
        {
            gl = new TraceGL(gl, System.err);
            gld.setGL(gl);
        }
    }



/* using these messes geometry up

	public Dimension getPreferredSize()
	{
		Dimension d = super.getPreferredSize();
		System.out.println(" pref sizes "+d+" "+glc.getPreferredSize());
		return new Dimension(0,0);
	}

    public void setPreferredSize(Dimension size)
    {
        System.out.println( "glj setPreferredSize " + size.width + ", " + size.height);
        super.setPreferredSize(size);
		glc.setPreferredSize(size);
    }

	public void reshape(int x, int y, int width, int height)
	{
		System.out.println("glj reshape glc"+glc+" "+x+" "+y+" "+width+" "+height);
		super.reshape(x,y,width,height);
		if (glc == null) return;
		glc.setSize(width,height);
		glc.setPreferredSize(new Dimension(width,height));
		repaint();
	}
*/
}
