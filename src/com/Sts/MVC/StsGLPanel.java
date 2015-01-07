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


public class StsGLPanel extends JPanel implements StsGLDrawable, StsSerializable, GLEventListener, KeyListener, MouseListener, MouseMotionListener, FocusListener, MouseWheelListener
// abstract public class StsGLPanel extends GLJPanel implements GLEventListener, KeyListener, StsSerializable
{
    /** window coordinates with origin in lower left */
    //    public Rectangle winRectGL = null;
   /** GL window pixel coordinates: same as JPanel.size, but repeated here so it can be easily persisted */
    public Dimension size = new Dimension();
    /** viewport parameters */
    public int[] viewPort = new int[4];
    /** current view being displayed */
    public StsView view = null;

    /** current background color */
    transient public StsColor clearColor = StsColor.BLACK;
    /** current foreground color for text and ticks */
    transient public StsColor foregroundColor = StsColor.WHITE;
    /** Convenience copy of the model. */
    transient public StsModel model = null;
    /** Convenience copy of the controller. */
    transient public StsActionManager actionManager = null;
    /** The OpenGL component to render onto */
    transient public Component glc = null;
	transient public GLAutoDrawable gld = null;
    /** The OpenGL pipeline interface */
    transient protected GL gl;
    /** The GL utility library interface */
    transient protected GLU glu;
    /** GL window pixel coordinates used for picking when winRectGL is not set correctly */
    transient public Rectangle insetWinRectGL = null;
    /** 4x4 modelView matrix (local copy) */
    transient public double[] modelViewMatrix = new double[16];
    /** 4x4 matrix which shifts view around in 3D screen coors */
    transient public double[] shiftProjMatrix =
            {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
    /** 4x4 perspective matrix (local copy) */
    transient public double[] projectionMatrix = new double[16];
    /** view shift adjust factor ("O" increase, "o" decrease) */
    transient protected double viewShiftFactor = 1.0;

    /** use view shifts if true */
    transient protected boolean viewShiftOn = true;
    /** keyboard char currently pressed */
    transient protected char currentKeyPressed = ' ';

    /** listeners to change in view; used by sliders controlling view for example */
    transient protected ChangeListener viewChangeListener = null;
    /**
     * Indicates view has changed; generally calls for redoing the projection matrix.
     * Because this needs to be done in the appropriate GLContext, the easiest thing
     * to do is set this flag and then when a redraw occurs, redo the matrix before
     * starting the draw routines.
     */
    transient public boolean viewChanged = true;
    /** indicates the view should be reinitialized as the view on this glPanel has changed */
    transient public boolean panelViewChanged = true;
    /**
     * Indicates the viewPort size has changed.  Before starting the redraw,
     * respond to this flag and make any necessary changes.
     */
    transient public boolean viewPortChanged = true;
    /** indicates projection changed and projectionMatrix should be recomputed */
    transient public boolean projectionChanged = true;
    /** Intended to handle printing of window.  Printing is not currently implemented. */
    transient protected boolean print = false;

    /**
     * View mode for the mouse; the 3 current views are ROTATE, ZOOM, and PAN.
     * ROTATE only applies in 3D windows but ZOOM and PAN can be used in 2d or 3d.
     * RECTZOOM also applies to 2d and 3d as a fourth mode.
     * mouseMode: zoom {@see com.Sts.MVC.StsCursor}
     */
//    transient private int mouseMode = StsCursor.ZOOM;

    /** Indicates view information for mouse should be displayed in info panel. */
    transient public boolean mouseInfo = false;
    /** The displayed mouse cursor. */
    transient public StsCursor cursor = null;
    /** indicates GL is initialized and ready to go. */
    // transient public boolean initialized = false;
    /** state of the GLPanel construction processed: NONE, STARTED, INITIALIZED or SUSPENDED */
    transient protected byte state = GL_STATE_NONE;

    transient protected Vector changeListeners = new Vector();

    transient public boolean hasFocus = false;

    /** max texture tile dimension to be used */
    static public int maxTextureSize = 0;

    /** has all mouse location and button information; singleton */
    static public StsMouse mouse = StsMouse.getInstance();

    static public final boolean textureDebug = false;

    /** state of GLPanel construction process; {@link #state} */
    static public final byte GL_STATE_NONE = StsGraphicParameters.GL_STATE_NONE; // glEventListener setup
    static public final byte GL_STATE_INITIALIZED = StsGraphicParameters.GL_STATE_INITIALIZED; // GLCanvas or GLJPanel has been initialized
    static public final byte GL_STATE_STARTED = StsGraphicParameters.GL_STATE_STARTED; // glEventListener has been started (GLEVentListener attached)
    static public final byte GL_STATE_SUSPENDED = StsGraphicParameters.GL_STATE_SUSPENDED;

    static public boolean debugPicker = false;
    static public boolean debugPicking = false;
    static public boolean debugGLmessages = false;

    static public final boolean debugProjectionMatrix = false;
    static public final boolean debugModelViewMatrix = false;
    static public final boolean glStateDebug = false;
    static public final double[] debugMatrix = new double[16];

    transient ByteBuffer backingBuffer = null;

    protected static GLCapabilities glCapabilities;

	// transient Color backgroundColor = Color.BLACK;

    static
    {
        glCapabilities = new GLCapabilities();
        glCapabilities.setDepthBits(24);
        glCapabilities.setBlueBits(8);
        glCapabilities.setRedBits(8);
        glCapabilities.setGreenBits(8);
        glCapabilities.setAlphaBits(8);
        glCapabilities.setDoubleBuffered(true);
    }

    static public final Dimension size00 = new Dimension(0, 0);

    /** default constructor: not normally used */
    public StsGLPanel()
    {
	       super();

    }

    public StsGLPanel(StsModel model, StsActionManager actionManager, int width, int height, StsView view)
    {
        this(model, actionManager);
        this.view = view;
        setPanelSize(width, height);

    }

    public StsGLPanel(StsModel model, StsActionManager actionManager)
    {
        initializeTransients(model, actionManager);

    }

    /** All operations are handled by the GLCanvas.  Adding the GLEventListener triggers the start of the GLEvent operations */
    protected void addGLlisteners()
    {
        glc.setFocusable(true);
        glc.addKeyListener(this);
        glc.addMouseListener(this);
        glc.addMouseMotionListener(this);
        glc.addFocusListener(this);
        glc.addMouseWheelListener(this);
        ((GLAutoDrawable)glc).addGLEventListener(this);
    }

    /** version for GLJPanel */
    /*
    public StsGLPanel(StsModel model, StsWindowActionManager actionManager, Dimension size, GraphicsDevice g)
    {
        super();
//        super(glCapabilities, null, null);
        printGLCapabilities();
        startGL();
        this.model = model;
        this.actionManager = actionManager;
        if (model.project == null) return;
        try
        {
            setSize(size);
            cursor = new StsCursor(this);
            setBorder(null);
            setLayout(new BorderLayout());
            //     initGL(g);
            setFocusable(true);
            addKeyListener(this);
            if (actionManager != null) actionManager.addListeners(this);
            if (model != null) setClearColor(model.project.getBackgroundStsColor());
//            startGL();
        }
        catch (Exception e)
        {
            e.printStackTrace();
            StsException.outputFatalException(this, "constructor", e);

        }
    }
    */
    /** Initialize GL.  Adding the GLEventListener actually starts it.
     *  GL then calls init(GLAutoDrawable) where additional
     */

    public void initAndStartGL(GraphicsDevice g)
    {
         if(model == null)
        {
            StsException.systemError(this, "initAndStartGL", "model is null.");
            return;
        }
        if(state == GL_STATE_STARTED)
        {
            paint(getGraphics());
            return;
        }
        else if(state != GL_STATE_NONE)
            return;

        try
        {
			GLContext shareableContext = Main.sharedContext;
			glc = new GLCanvas(glCapabilities, null, shareableContext, g);
			gld = (GLAutoDrawable) glc;
			//shareableContext.makeCurrent();
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "initGL", e);
        }
        // sets GL which might be debug or trace GL

        setGL();
        glu = new GLU();

        // add panel and initialize panel size
        add(glc);
        initializeSize();

        String renderer = gl.glGetString(GL.GL_RENDERER);
        if(StsStringUtils.stringContainsString(renderer, "Mesa"))
            throw new RuntimeException("Missing Linux Accelerated 3D OpenGL drivers");
        if(StsStringUtils.stringContainsString(renderer, "GDI Gen"))
            throw new RuntimeException("Missing Windows Accelerated 3D OpenGL drivers");

        // start GL by adding GL listener
        addGLlisteners();
        state = GL_STATE_STARTED;
        if(glStateDebug) System.out.println("GL initialized and started.");
        if(Main.isGLDebug || Main.isGLTrace) System.out.println("GL STATE: STARTED");
		if (getGraphics() != null)
		   paint(getGraphics()); // jbw jdk 1.6.11 and above -- force initial paint
    }



    /** seems to be required for proper splitPane behavior when a glPanel is in a splitPane */
    protected void initializeSize()
    {
        setMinimumSize(size00);
        Dimension size = getSize();
        setPreferredSize(size);
        glc.setPreferredSize(size);
    }
/*
    public void startGL()
    {
        if(state != GL_STATE_INITIALIZED) return;
        glc.addGLEventListener(this);
        state = GL_STATE_STARTED;
        if(glStateDebug) System.out.println("GL started.");
        if(Main.isGLDebug || Main.isGLTrace) System.out.println("GL STATE: STARTED");
    }
*/

    public void init(GLAutoDrawable drawable)
    {
        if(view != null) view.init(drawable);
    }

    public void display(GLAutoDrawable drawable)
    {
        if(view != null) view.display(drawable);
    }
    /*public void setSize(int width, int height)
      {
          setSize(new Dimension(size.width,size.height));
      }
    */

    // Component initialization from DB
    /*
        public void initializeTransients(StsModel model, StsWindowActionManager actionManager)
        {
            try
            {
                this.model = model;
                this.actionManager = actionManager;
                //setSize(size);
                cursor = new StsCursor(this);
                initGL();
                glc.setFocusable(true);
                glc.addKeyListener(this);
                if(actionManager != null)actionManager.addListeners(this);
                if(model != null)setClearColor(model.getProject().getBackgroundStsColor());
            }
            catch(Exception e)
            {
                StsException.outputException("StsGLPanel.constructor() failed.", e, StsException.WARNING);
            }
        }
    */

    protected void printGLCapabilities()
    {
        System.out.println(((GLAutoDrawable)glc).getChosenGLCapabilities());
        // System.out.println(getChosenGLCapabilities());
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
		//System.out.println("setpanelsize "+size);
        setPreferredSize(size);
    }

    public void setBounds(int x, int y, int width, int height)
    {
       size.width = width;
       size.height = height;
       super.setBounds(x, y, width, height);
	   //System.out.println("setbounds glc="+glc);
	   repaint();
    }

    /**
     * GL window has origin in lower-left while a Java window has the origin in the upper left.
     * So GLXOrigin = WindowXOrigin && GLYOrigin = WindowYOrigin - WindowHeight
     */
    public int geGLX()
    { return getX(); }

    /**
     * GL window has origin in lower-left while a Java window has the origin in the upper left.
     * So GLXOrigin = WindowXOrigin && GLYOrigin = WindowYOrigin - WindowHeight
     */
    public int getGLY()
    { return getHeight() - getY() - 1; }

    /**
     * GL window has origin in lower-left while a Java window has the origin in the upper left.
     * So GLX = WindowX && GLY = WindowYOrigin - WindowY
     */
    public int getGLX(int windowX)
    { return windowX; }

    public int getGLY(int windowY) { return getHeight() - windowY - 1; }

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
//            if(model != null) setClearColor(model.getProject().getBackgroundStsColor());
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "initializeTransients", e);
        }
    }

    /** standard constructor: use this for preferred size */
    /*
        public StsGLPanel(StsModel model, StsWindowActionManager actionManager, Dimension size)
        {
       this(model, actionManager, size, false);
        }
    */
    /** standard constructor: use this for fixed size */
    /*
        public StsGLPanel(StsModel model, StsWindowActionManager actionManager, Dimension size, boolean isFixed)
        {
       try
       {
           this.model = model;
           this.actionManager = actionManager;

           if (isFixed)
           {
               setMinimumSize(size);
               setPreferredSize(size);
           }
           else
               setPreferredSize(size);

     winRectGL = new Rectangle(size);
     cursor = new StsCursor(this);

//            initializeSize(size);
//            initializeGL();
     initGL();
           glc.setFocusable(true);
           glc.addKeyListener(this);
           if(actionManager != null) actionManager.addListeners(this);
//            createGLRunnable();
           if(model != null) setClearColor(model.getProject().getBackgroundStsColor());
       }
       catch (Exception e)
       {
           StsException.outputException("StsGLPanel.constructor() failed.", e, StsException.WARNING);
       }
        }
    */
 /*
    public void setBackgroundColor(Color c)
	{
		super.setBackground(c);
		backgroundColor = c;
		setClearColor(new StsColor(c));
	}
*/
    public void repaint()
    {
        if(glc != null) glc.repaint();
    }

    public void paint(Graphics g)
	{
	   // System.out.println(" GLPanel paint "+glc);
//	   setBackground(Color.BLACK);
	   super.paint(g);
	   //super.paintChildren(g);
	   if (glc!=null)
       {
	       //System.out.println("glc paint parent"+this);

		   glc.setVisible(true); // jbw hack for maximize bug

		   glc.paint(g);

		   int count=0;
		   while (glc.getWidth() <=0 && count < 1000) // peer not realized
		   {
			//   System.out.println("width 0");
			   glc.setSize(this.getSize());
			   glc.paint(g);
			   count++;
		   }
	   }
       else
           StsException.systemDebug(this, "paint", "glc is null");
    }

    private boolean isZeroSize(Dimension size)
    {
        return size.width == 0 || size.height == 0;
    }


    public double[] getProjectionMatrix()
    {
        return projectionMatrix;
    }

    public double[] getModelViewMatrix()
    {
        return modelViewMatrix;
    }

    /**
     * The action manager sends key events to both this glPanel and to the currentAction.
     * So put any key event handling here in coordination with current actions.
     */
    public void keyPressed(KeyEvent e, StsMouse mouse)
    {
       	if(view != null) view.keyPressed(mouse, e);
    }
    /**
     * The action manager sends key events to both this glPanel and to the currentAction.
     * So put any key event handling here in coordination with current actions.
     */
    public void keyReleased(KeyEvent e, StsMouse mouse)
    {
        if(view != null) view.keyReleased(mouse, e);
    }

    /** Implement this method to handle window motion in response to mouse events. */
    public void moveWindow(StsMouse mouse)
    {
        if(view != null) view.moveWindow(mouse);
    }

    /**
     * A middle mouse or shift-right mouse button can trigger a popup menu.
     * Implement this event handling here.
     * @param mouse
     */
    public void showPopupMenu(StsMouse mouse)
     {
         if(view == null) return;
         view.showPopupMenu(mouse);
     }

    /** each view is responsible for the projectionMatrix being used */
    public void computeProjectionMatrix()
    {
        if(view == null) return;
        view.computeProjectionMatrix();
    }

    /**
     * A left mouse pressed asks for a default edit of an object under the mouse.
     * Implement this event handling here.
     */
    public StsObject defaultMousePressed(StsMouse mouse)
    {
        return null;
    }

    /**
     * A left mouse motion asks for a default edit of an object under the mouse.
     * Implement this event handling here.
     */
    public void defaultMouseMotion(StsMouse mouse)
    {
        return;
    }

	public void cancelPopupMenu(StsMouse mouse)
	{
       if(view == null) return;
	      view.cancelPopupMenu(mouse);
     }

    /*
      public void initializeSize(Dimension size)
      {
          setPreferredSize(size);
          winRectGL = new Rectangle(size);
          cursor = new StsCursor(this);
      }

     /** insures that we are on the eventDispatchThread and initializes GL */
    /*
       private void initializeGL()
       {
           if(state >= INITIALIZED)
           {
               System.err.println("System error.  StsGLPanel already initialized.");
               return;
           }

           try
           {
               if(SwingUtilities.isEventDispatchThread())
               {
                   if(debugGLmessages) printDebugGL("initializeGL: already on eventDispatchThread");
                   initGL();
               }
               else
               {
                   if(debugGLmessages) printDebugGL("initializeGL: switching to eventDispatchThread");
                   Runnable initGLRunnable = new Runnable()
                   {
                       public void run()
                       {
                           initGL();
                       }
                   };
                   SwingUtilities.invokeAndWait(initGLRunnable);
               }
           }
           catch (Exception e)
           {
               StsException.outputException("StsGLPanel.repaint() failed.", e, StsException.WARNING);
           }
       }
    */

    /** classInitialize the GL capabilities */
    /*
        private void initGL()
        {
          initGL(null);
        }

        private void initGL(GraphicsDevice g)
        {
            setLayout(new BorderLayout());
            setSize(winRectGL.width, winRectGL.height);
        }
    */
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

        // jbw 1/23/2010 dump the reshapes to hopefully let swing do it better.

        setViewPort(x, y, width, height);
        glc.reshape(x, y, width, height);

        if(view != null)
        {
            view.reshape(x, y, width, height);
            view.computeProjectionMatrix();
        }
        viewPortChanged = false;
        model.enableDisplay();
		glc.setVisible(false); // jbw hack! for maximize bug
		//glc.reshape(x, y, width, height);

		repaint();
    }

    public void setViewPort(int x, int y, int width, int height)
    {
//         StsException.systemDebug(this, "setViewPort", " x: " + x + " y: " + y + " width: " + width + " height: " + height);
        gl.glViewport(x, y, width, height);
        gl.glGetIntegerv(GL.GL_VIEWPORT, viewPort, 0);
        // debugPrintViewport("setViewPort(x,y,width.height)");
    }

 //   public int getWidth() { return viewPort[2]; }
 //   public int getHeight() { return viewPort[3]; }
    public double getAspectRatio() { return viewPort[2]/viewPort[3]; }

    public void displayChanged(GLAutoDrawable drawable, boolean modeChanged, boolean deviceChanged)
    {
    }

    /*
       private void setGL()
       {
           coregl = new CoreGL();
           gl = (GL)coregl;
           coreglu = new CoreGLU();
           glu = (GLU)coreglu;
       }
    */
    /*
       protected void setGLdebugLockAndTrace(boolean mainDebug)
       {
           if (mainDebug) setGLdebug(GL_LOCK_TRACE);
           else
               setGLdebug(GL_NO_DEBUG);
       }
    */
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

    public GL getGL()
    {
        return gl;
    }

    public GLU getGLU()
    {
        return glu;
    }

    /**
     * Needed only when GLcontext must be manually set (normally handled by Magician).
     * Currently only called in picking operations.
     */
    /*
       public void setGLContext()
       {
           if (debugGLmessages) System.out.println("glPanel3d.setGLContext() called.");
           GLContext context = glc.getContext();
           try
           {
               if (context == null)
                   System.out.println("Null GLContext in setGLContext()!");
               else
                   context.makeCurrent();
           }
           catch (Exception e)
           {
               System.out.println("Exception in setGLContext()! \n" + e);
               if (context != null)
               {
                   System.out.print("Retrying...");
                   context.unlock();
                   try { context.makeCurrent(); System.out.println("successful."); }
                   catch (Exception ex) { System.out.println("failed."); }
               }
           }
       }

       public void freeGLContext()
       {
           GLContext context = glc.getContext();
           try
           {
               if (context == null)
                   System.out.println("Null GLContext in setGLContext()!");
               else
                   context.unlock();
           }
           catch (Exception e)
           {
               System.out.println("Exception in freeGLContext()! \n" + e);
           }
       }
    */
    public void destroy()
    {
        // StsTextureList.deleteAllTextures(gl);
        // deleteAllDisplayLists(gl);
        if(glc != null) gld.removeGLEventListener(this);
    }

    public void resetViewPort()
    {
        gl.glViewport(viewPort[0], viewPort[1], viewPort[2], viewPort[3]);
//        viewPortChanged = true;
//        gl.glGetIntegerv(GL.GL_VIEWPORT, viewPort, 0);
//        debugPrintViewport("resetViewPort()");
    }

    public void debugPrintViewport(String string)
    {
        int[] viewPort = new int[4];
        gl.glGetIntegerv(GL.GL_VIEWPORT, viewPort, 0);
        System.out.println(string + " GL     viewport. x: " + viewPort[0] + " y: " + viewPort[1] + " width: " + viewPort[2] + " height: " + viewPort[3]);
        System.out.println(string + " window viewport. x: " + geGLX() + " y: " + getGLY() + " width: " + getWidth() + " height: " + getHeight());
    }

    public StsMousePoint getMousePoint() { return mouse.getMousePoint(); }

    static public StsMouse getMouse() { return mouse; }

    static public boolean isLeftMouseDown()
    {
        return mouse.isButtonStatePressed(StsMouse.LEFT);
    }

    static public boolean isLeftMouseDragging()
    {
        return mouse.isButtonStateDragging(StsMouse.LEFT);
    }
    /**
     * mouse coordinates in Java are from the upper left corner.  GL mouse
     * coordinates are from the lower left corner.  This routine returns
     * a new mousePoint with the origin shifted to the lower left.
     */
    public StsMousePoint getMousePointGL()
    {
        if(actionManager != null)
        {
            StsMousePoint mousePoint = getMousePoint();
            return mousePoint.lowerLeftCoor(viewPort);
        }
        return null;
    }

    public int getGLMouseY(int mouseY)
    {
        return viewPort[3] - mouseY - 1;
    }

    public void changeModelView()
    {
        viewChanged = true;
        repaint();
    }

    public void projectionChanged()
    {
        projectionChanged = true;
        repaint();
    }

    public void addViewChangeListener(ChangeListener listener)
    {
        viewChangeListener = listener;
    }

    public void fireViewChangeEvent()
    {
        if(this.viewChangeListener != null)
            viewChangeListener.stateChanged(new ChangeEvent(this));
    }

    public void printDebugGL(String message)
    {
        System.out.println("GL mainDebug: " + message);
    }

    public void statusBarMessage(String msg)
    {
        StsMessageFiles.logMessage(msg);
    }

    public void printWindow()
    {
        print = true;
        repaint();
        print = false;
    }

    /** utility routine used to return screen coordinates given world coordinates. */
    public double[] getScreenCoordinates(StsPoint worldPoint)
    {
        double[] screenXYZ = new double[3];
        float[] xyz = worldPoint.getPointXYZ();
        glu.gluProject(xyz[0], xyz[1], xyz[2], modelViewMatrix, 0, projectionMatrix, 0, viewPort, 0, screenXYZ, 0);
        return screenXYZ;
    }

    public double[] getScreenCoordinates(double x, double y, double z)
    {
        double[] screenXYZ = new double[3];
        glu.gluProject(x, y, z, modelViewMatrix, 0, projectionMatrix, 0, viewPort, 0, screenXYZ, 0);
        return screenXYZ;
    }

    /** utility routine used to return screen coordinates given world coordinates. */
    public double[] getScreenCoordinates(float[] v)
    {
        double[] screenXYZ = new double[3];
        glu.gluProject((double)v[0], (double)v[1], (double)v[2], modelViewMatrix, 0, projectionMatrix, 0, viewPort, 0, screenXYZ, 0);
        return screenXYZ;
    }

    public double[] getScreenCoordinates(double[] v)
    {
        double[] screenXYZ = new double[3];
        glu.gluProject(v[0], v[1], v[2], modelViewMatrix, 0, projectionMatrix, 0, viewPort, 0, screenXYZ, 0);
        return screenXYZ;
    }

    /** utility routine used to return screen coordinates given 2d world coordinates. */
    public double[] getScreenCoordinates(float x, float y)
    {
        double[] screenXYZ = new double[3];
        glu.gluProject((double)x, (double)y, 0.0, modelViewMatrix, 0, projectionMatrix, 0, viewPort, 0, screenXYZ, 0);
        return screenXYZ;
    }

    public double getScreenCoordinateDistance(float[] p1, float[] p2)
    {
        double[] xyz1 = getScreenCoordinates(p1);
        double[] xyz2 = getScreenCoordinates(p2);
        return StsMath.distance(xyz1, xyz2, 3);
    }

    public void loadProjectionMatrix(GL gl)
    {
        gl.glMatrixMode(GL.GL_PROJECTION);
        gl.glLoadMatrixd(projectionMatrix, 0);
        gl.glMatrixMode(GL.GL_MODELVIEW);
    }

    /** utility routine used to return screen z coordinate given world coordinates. */
    public float getScreenZ(float[] v)
    {
        double[] screenXYZ = new double[3];
        glu.gluProject(v[0], v[1], v[2], modelViewMatrix, 0, projectionMatrix, 0, viewPort, 0, screenXYZ, 0);
        return (float)screenXYZ[0];
    }

    /** utility routine used to return world coordinates given screen coordinates. */
    public StsPoint getWorldCoordinatesPoint(double[] screenPoint)
    {
        double[] worldXYZ = new double[3];
        glu.gluUnProject(screenPoint[0], screenPoint[1], screenPoint[2],
                modelViewMatrix, 0, projectionMatrix, 0, viewPort, 0,
                worldXYZ, 0);
        return new StsPoint(worldXYZ);
    }

    /** utility routine used to return world coordinates given screen coordinates. */
    public double[] getWorldCoordinates(double[] screenPoint)
    {
        double[] worldXYZ = new double[3];
        glu.gluUnProject(screenPoint[0], screenPoint[1], screenPoint[2],
                modelViewMatrix, 0, projectionMatrix, 0, viewPort, 0,
                worldXYZ, 0);
        return worldXYZ;
    }

    public float[] getWorldCoordinateFloats(double[] screenPoint)
    {
        double[] xyz = getWorldCoordinates(screenPoint);
        return StsMath.convertDoubleToFloatArray(xyz);
    }

    /** Shift the view towards (+) or away (-) from the viewer */
    private double currentShift = 0.0;
    private final double viewShift = 0.00025; /** This is a hard-wired offset value */

    /**
     * Shift the view towards (+) or away (-) from the viewer
     * Note that the absolute amount of shift can be scaled with viewShiftFactor
     * which is adjusted with the "O" key (increase by 2.0) and the "o" key
     * (decrease by 2.0).
     * Since several viewShifts might be applied in sequence, we rely on the existence of the "original" projectionMatrix
     * which is used in the multiply here and in the restore when resetViewShift is called.
     */
    public void setViewShift(GL gl, double shift)
    {
        if(shift == 0.0) return;
        if(!viewShiftOn || shift == currentShift) return;
        currentShift = shift;
    /*
        if(StsGLPanel.debugProjectionMatrix)
        {
            double totalShift = shift * viewShift * viewShiftFactor;
            debugPrintProjectionMatrix("StsGLPanel.setViewShift(). shift: " + shift + " totalShift: " + totalShift + " proj matrix before shift: ");
        }
   */
        gl.glMatrixMode(GL.GL_PROJECTION);
        shiftProjMatrix[14] = -shift * viewShift * viewShiftFactor;
        gl.glLoadMatrixd(shiftProjMatrix, 0);
        gl.glMultMatrixd(projectionMatrix, 0);
        gl.glMatrixMode(GL.GL_MODELVIEW);
    }

    public double getViewShift()
    {
        return viewShift;
    }

    public void resetViewShift(GL gl)
    {
        if(!viewShiftOn) return;
        if(currentShift == 0.0) return;
        currentShift = 0.0;
        gl.glMatrixMode(GL.GL_PROJECTION);
        gl.glLoadMatrixd(projectionMatrix, 0);
        gl.glMatrixMode(GL.GL_MODELVIEW);
//        if(StsGLPanel.debugProjectionMatrix)
//            debugPrintProjectionMatrix("StsGLPanel.resetViewShift(). proj matrix after reset: ");
    }

    public void toggleViewShift(ItemEvent e)
    {
        boolean viewShiftOn = (e.getStateChange() == ItemEvent.SELECTED);
        toggleViewShift(viewShiftOn);
    }

    public void changeViewShiftFactor(double factor)
    {
        viewShiftFactor *= factor;
        //        System.out.println("viewShiftFactor: " + viewShiftFactor);
        changeModelView();
    }

    public void toggleViewShift(boolean viewShiftOn)
    {
        this.viewShiftOn = viewShiftOn;
        changeModelView();
    }

/*
    public int getMouseMode()
    {
        return mouseMode;
    }

    public void setMouseMode(int mouseMode)
    {
        this.mouseMode = mouseMode;
    }
*/
    public void restoreCursor()
    {
        cursor.restoreCursor();
    }
/*
    public void setMotionCursor()
    {
        cursor.setStsCursor(getMouseMode());
    }
*/
    public void setMotionCursor(int mode)
    {
        cursor.setStsCursor(mode);
    }

    public void setClearColor(StsColor color)
    {
        clearColor = color;
        setClearColor();
    }

    public void setClearColor()
    {
        Color javaColor = clearColor.getColor();
        super.setBackground(javaColor);
		if (glc != null)
        glc.setBackground(javaColor);
    }

    public void applyClearColor(int clearBits)
    {
		glc.setBackground(clearColor.getColor());
        gl.glClearColor(clearColor.red, clearColor.green, clearColor.blue, clearColor.alpha);
        gl.glClear(clearBits);
    }
/*
    public void setClearColor(Color color)
    {
        clearColor = new StsColor(color);
    }
*/
/*
     static public void toggleGLDrawDebug(ItemEvent e)
     {
         debugGLmessages = (e.getStateChange() == ItemEvent.SELECTED);
     }
*/
/*

               if(SwingUtilities.isEventDispatchThread())
               {
                   if(debugGLmessages) printDebugGL("repaint: already on eventDispatchThread");
                   glc.repaint();
               }
               else
               {
                   if(debugGLmessages)
                   {
                       System.err.println("repaint: switching to eventDispatchThread");
                       debugPrintGLEventListener();
                   }

                   Runnable glRunnable = new Runnable()
                   {
                       public void run()
                       {
                           if(debugGLmessages) printDebugGL("glRunnable called.");
                           // setGLContext();
                           glc.repaint();
                           // freeGLContext();
                       }
                   };

                   SwingUtilities.invokeAndWait(glRunnable);
                   if(debugGLmessages) printDebugGL("StsGLPanel.repaint().invokeAndWait() completed.");
               }
           }
           catch (Exception e)
           {
               StsException.outputException("StsGLPanel.repaint() failed.", e, StsException.WARNING);
           }
       }
    */

    /*
        public void addTextureToDeleteList(StsTextureTiles deleteTextureTiles)
        {
            deleteTextureTiles.textureChanged();
        }

        public void addTextureToList(StsTextureTiles textureTiles)
        {
            if (textureList.contains(textureTiles)) return;
            if (textureDebug) System.out.println("Added " + textureTiles.toString() + " to list.");
            textureList.add(textureTiles);
        }
    */


    /*
        public void checkDeleteDisplayLists(GL gl)
        {
            StsTextureTiles textureTiles;
            Iterator iterator = textureList.iterator();
            while (iterator.hasNext())
            {
                textureTiles = (StsTextureTiles) iterator.next();
                textureTiles.checkDeleteDisplayLists(gl);
            }
        }

        public void deleteAllDisplayLists(GL gl)
        {
            StsTextureTiles textureTiles;
            Iterator iterator = textureList.iterator();
            while (iterator.hasNext())
            {
                textureTiles = (StsTextureTiles) iterator.next();
                textureTiles.deleteDisplayLists(gl);
            }
        }
    */

/*
    public void clearAllDisplayLists()
    {
        StsTextureTiles textureTiles;
        Iterator iterator = textureList.iterator();
        while (iterator.hasNext())
        {
            textureTiles = (StsTextureTiles) iterator.next();
            textureTiles.geometryChanged();
        }
    }
*/

    /*
        static public void toggleGLDrawDebug(ItemEvent e)
        {
            debugGLmessages = (e.getStateChange() == ItemEvent.SELECTED);
        }
    */
    static public void toggleDebugPicker(ItemEvent e)
    {
        debugPicker = (e.getStateChange() == ItemEvent.SELECTED);
    }

    public double[] gluUnProject(double screenX, double screenY, double screenZ,
                                 double[] modelViewMatrix, double[] projectionMatrix, int[] viewPort)
    {
        double[] worldXYZ = new double[3];
        glu.gluUnProject(screenX, screenY, screenZ, modelViewMatrix, 0, projectionMatrix, 0, viewPort, 0,
                worldXYZ, 0);
        return worldXYZ;
    }

    public double[] gluProject(double worldX, double worldY, double worldZ,
                               double[] modelViewMatrix, double[] projectionMatrix, int[] viewPort)
    {
        double[] screenXYZ = new double[3];

        glu.gluProject(worldX, worldY, worldZ, modelViewMatrix, 0, projectionMatrix, 0, viewPort, 0,
                screenXYZ, 0);
        return screenXYZ;
    }

    public int getGLWidth() { return getWidth(); }

    public int getGLHeight() { return getHeight();}

    transient boolean pixelsSaved = false;
    transient boolean doPixelsSaved = false;

    public void savePixels(boolean saveem)
    {
        if(saveem)
        {
            doPixelsSaved = true;
            pixelsSaved = false;
        }
        else
        {
            doPixelsSaved = pixelsSaved = false;
            deleteBacking();
        }
    }

    public void doSavePixels()
    {
        //System.out.println("save");
        pixelsSaved = saveBacking();
        doPixelsSaved = false;
    }

    public void doRestorePixels()
    {
        //System.out.println("restore");
        restoreBacking();

    }

    public boolean saveBacking()
    {
        int numBytes = getGLWidth() * getGLHeight() * 4;
        if(backingBuffer == null || backingBuffer.capacity() < numBytes)
        {
            try
            {
                backingBuffer = ByteBuffer.allocateDirect(numBytes);
            }
            catch(Exception e)
            {
                return false;
            }
            catch(java.lang.OutOfMemoryError e2)
            {
                return false;
            }
        }
        backingBuffer.rewind();
        gl.glViewport(0, 0, getGLWidth(), getGLHeight());
        gl.glPixelStorei(GL.GL_PACK_ALIGNMENT, 1);
        gl.glReadPixels(0, 0, getGLWidth(), getGLHeight(), GL.GL_BGRA, GL.GL_UNSIGNED_BYTE, backingBuffer);
        return true;
    }

    public boolean restoreBacking()
    {
        int numBytes = getGLWidth() * getGLHeight() * 4;
        if(backingBuffer == null || backingBuffer.capacity() < numBytes)
            return false;

        backingBuffer.rewind();
        gl.glViewport(0, 0, getGLWidth(), getGLHeight());
        gl.glRasterPos2i(0, 0);
        gl.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1);
        gl.glDrawPixels(getGLWidth(), getGLHeight(), GL.GL_BGRA, GL.GL_UNSIGNED_BYTE, backingBuffer);
        return true;
    }

    public void deleteBacking()
    {
        if(backingBuffer != null)
            backingBuffer = null;
    }

    public void debugPrintMatrixMode(String message)
    {
        int[] mm = new int[1];
        gl.glGetIntegerv(GL.GL_MATRIX_MODE, mm, 0);
        int matrixMode = mm[0];
        String modeString = "";
        if(matrixMode == GL.GL_PROJECTION)
            modeString = "Projection";
        else if(matrixMode == GL.GL_MODELVIEW)
            modeString = "ModelView";
        System.out.println("matrix mode: " + modeString + " " + message);
    }

    public void debugPrintProjectionMatrix(String message)
    {
        double[] matrix = new double[16];
        gl.glGetDoublev(GL.GL_PROJECTION_MATRIX, matrix, 0);
        System.out.println("projection matrix " + message + " [0][0]: " + matrix[0] + " [1][1]: " + matrix[5] + " [2][2]: " + matrix[10] + " [3][3]: " + matrix[15]);
    }

    public void debugPrintModelViewMatrix(String message)
    {
        double[] matrix = new double[16];
        gl.glGetDoublev(GL.GL_MODELVIEW_MATRIX, matrix, 0);
        System.out.println("modelView matrix " + message + " [0][0]: " + matrix[0] + " [1][1]: " + matrix[5] + " [2][2]: " + matrix[10] + " [3][3]: " + matrix[15]);
    }

    public void debugPrintModelViewMatrixChanged(String message)
    {
        debugPrintCompareMatrices("MODELVIEW MATRIX " + message, debugMatrix, modelViewMatrix);
    }

    public void debugPrintProjectionMatrixChanged(String message)
    {
        debugPrintCompareMatrices("PROJECTION MATRIX " + message, debugMatrix, projectionMatrix);
    }

    public void debugPrintCompareMatrices(String message, double[] oldMatrix, double[] newMatrix)
    {
        int n = 0;
        boolean changed = false;
        for(int i = 0; i < 4; i++)
        {
            for(int j = 0; j < 4; j++)
            {
                if(oldMatrix[n] != newMatrix[n])
                changed = true;
                n++;
            }
        }
        if(!changed)
        {
            System.out.println(message + " NOT CHANGED.");
            return;
        }
        System.out.println(message + " CHANGED.");
        System.out.print("    ");
        n = 0;
        for(int i = 0; i < 4; i++)
        {
            for(int j = 0; j < 4; j++)
            {
                if(oldMatrix[n] != newMatrix[n])
                System.out.print( " [" + i + "][" + j + "]: " + oldMatrix[n] + " to " + newMatrix[n]);
                n++;
            }
        }
        System.out.println();
     }
    /** ********** custom (shallow) serialization ********* */
    /*
    * custom serialization requires versioning to prevent old persisted files from barfing.
    * if you add/change fields, you need to bump the serialVersionUID and fix the
    * reader to handle both old & new
    */
    static final long serialVersionUID = 1l;

    public void resurrect(StsGLPanel o)
    {

    }

    public boolean isInside(StsMousePoint mousePoint)
    {
        if(mousePoint.x < 0 || mousePoint.x > getWidth()) return false;
        if(mousePoint.y < 0 || mousePoint.y > getHeight()) return false;
        return true;
    }

    public boolean hasFocus()
    {
        if(glc == null) return false;
        return glc.hasFocus();
    }

    public void performMouseAction(MouseEvent e)
    {
        performMouseAction();
    }

    public void performMouseAction()
    {
        if(view != null) view.performMouseAction(actionManager, mouse);
    }

    /** test routine */
    protected void showLine(String s)
    {
//		System.out.println(s);
//		System.out.flush();
    }

    /** test routine */
    protected String mousemods(MouseEvent e)
    {
        int mods = e.getModifiers();
        String s = "";
        if(e.isShiftDown()) s += "Shift ";
        if(e.isControlDown()) s += "Ctrl ";
        if((mods & InputEvent.BUTTON1_MASK) != 0) s += "Button 1 ";
        if((mods & InputEvent.BUTTON2_MASK) != 0) s += "Button 2 ";
        if((mods & InputEvent.BUTTON3_MASK) != 0) s += "Button 3 ";
        return s;
    }

    /**
     * gets key pressed event from StsGLPanel:
     * passes it back with the mouse object for processing
     * by concrete subclass of StsGLPanel.
     */
    public void keyPressed(KeyEvent e)
    {
//        currentKeyPressed = e.getKeyChar();
        keyPressed(e, mouse);
        if(actionManager == null) return;
        StsAction currentAction = actionManager.getCurrentAction();
        if(currentAction != null) currentAction.keyPressed(e, mouse);
        // following code is for debugging key event
        /*
		String eventtype, modifiers, code, character;
		eventtype = "KEY_PRESSED";

		// Convert the list of modifier keys to a string
		modifiers = KeyEvent.getKeyModifiersText(e.getModifiers());

		// Get string and numeric versions of the key code, if any.
		if (e.getID() == KeyEvent.KEY_TYPED) code = "";
		else code = "Code=" + KeyEvent.getKeyText(e.getKeyCode()) + " (" + e.getKeyCode() + ")";

		// Get string and numeric versions of the Unicode character, if any.
		if (e.isActionKey()) character = "";
		else character = "Character=" + e.getKeyChar() + " (Unicode=" + ((int)e.getKeyChar()) + ")";

		// Display it all.
		showLine(eventtype + ": " + modifiers + " " + code + " " + character);
        */
    }

// Method: keyTyped

    public void keyTyped(KeyEvent e)
    {
        // following code is for debugging key event
        /*
		int keyCode;
		String eventtype, modifiers, code, character;

		keyCode = e.getKeyCode();
	//	if(keyCode == KeyEvent.VK_ALT) middleMouseDown = true;

		eventtype = "KEY_TYPED";

		// Convert the list of modifier keys to a string
		modifiers = KeyEvent.getKeyModifiersText(e.getModifiers());

		// Get string and numeric versions of the key code, if any.
		if (e.getID() == KeyEvent.KEY_TYPED) code = "";
		else code = "Code=" + KeyEvent.getKeyText(e.getKeyCode()) + " (" + e.getKeyCode() + ")";

		// Get string and numeric versions of the Unicode character, if any.
		if (e.isActionKey()) character = "";
		else character = "Character=" + e.getKeyChar() + " (Unicode=" + ((int)e.getKeyChar()) + ")";

		// Display it all.
		showLine(eventtype + ": " + modifiers + " " + code + " " + character);
        */
    }

    public void keyReleased(KeyEvent e)
    {
        keyReleased(e, mouse);
        StsAction currentAction = actionManager.getCurrentAction();
        if(currentAction != null) currentAction.keyReleased(e, mouse, this);

//        currentKeyPressed = ' ';
        // following code is for debugging key event
        /*
		int keyCode;
		String modifiers, code, character;

		// Convert the list of modifier keys to a string
		modifiers = KeyEvent.getKeyModifiersText(e.getModifiers());

		// Get string and numeric versions of the key code, if any.
		if (e.getID() == KeyEvent.KEY_TYPED) code = "";
		else code = "Code=" + KeyEvent.getKeyText(e.getKeyCode()) + " (" + e.getKeyCode() + ")";

		// Get string and numeric versions of the Unicode character, if any.
		if (e.isActionKey()) character = "";
		else character = "Character=" + e.getKeyChar() + " (Unicode=" + ((int)e.getKeyChar()) + ")";

		// Display it all.
		showLine("Key released: " + modifiers + " " + code + " " + character);
        */
    }

//    public char getCurrentKeyPressed() { return currentKeyPressed; }

    public void addChangeListener(ChangeListener listener)
    {
        changeListeners.addElement(listener);
//        fireChangeEvent();
    }

    public void removeChangeListener(ChangeListener listener)
    {
        changeListeners.removeElement(listener);
    }

    /*
        public void fireChangeEvent()
        {
            ChangeEvent e = new ChangeEvent(this);
            Enumeration enum = changeListeners.elements();
            while( enum.hasMoreElements() )
            {
                ChangeListener listener = (ChangeListener) enum.nextElement();
                listener.stateChanged(e);
            }
        }
    */
    public boolean isButtonDown(int button)
    {
        return mouse.isButtonDown(button);
    }

    public void focusGained(FocusEvent e)
    {
//        System.out.println("actionManager focus gained.");
//        Component component = e.getComponent();
//        System.out.println(component.toString());
//        glPanel.glc.requestFocus();
        hasFocus = true;
        // jbw repaint();
    }

    public void focusLost(FocusEvent e)
    {
//        System.out.println("actionManager focus lost.");
        hasFocus = false;
        // jbw repaint();
    }

    public void mousePressed(MouseEvent e)
    {
        glc.requestFocus();
        mouse.setState(e, StsMouse.PRESSED);
        performMouseAction(e);
    }

    public void mouseReleased(MouseEvent e)
    {
        mouse.setState(e, StsMouse.RELEASED);
        performMouseAction(e);
    }

    public void mouseDragged(MouseEvent e)
    {
        mouse.setState(e, StsMouse.DRAGGED);
        performMouseAction(e);
    }

    public void mouseClicked(MouseEvent e)
	{
		// glc.requestFocus();
        // mouse.setState(e, StsMouse.CLICKED);
        // performMouseAction(e);
	}

    public void mouseMoved(MouseEvent e)
    {
    }

    public void mouseWheelMoved(MouseWheelEvent e)
    {
        mouse.setState(e, StsMouse.WHEEL);
        performMouseAction(e);
    }

    public void mouseEntered(MouseEvent e)
    {
//        glPanel.glc.requestFocus();
//		System.out.println("Entered " + glPanel.getName());
    }

    public void mouseExited(MouseEvent e)
    {
//		System.out.println("Exited " + glPanel.getName());
    }

    public StsView getView() { return view; }

    public String getCurrentViewName()
    {
        if(view == null)
            return "null";
        else
            return view.getViewName();
    }

    public void setView(StsView view)
    {
        view.glPanel = this;
        this.view = view;
        if(this.view != null) this.view.viewChanged();
        // currentView.isViewGLInitialized = false;
//        viewPortChanged = true;
    }

    public void setPreferredSize(Dimension size)
    {
 //       if(size.width != 0 || size.height != 0)
 //           StsException.systemDebug(this, "setPreferredSize", "size not 0, 0: " + size.width + ", " + size.height);
        super.setPreferredSize(size);

    }

	public GraphicsDevice getGraphicsDevice()
	{
	   GraphicsConfiguration gConfig = getGraphicsConfiguration();
	   if(gConfig == null) return null;
	   return gConfig.getDevice();
    }

	public void reshape(int x, int y, int width, int height)
	{
		//System.out.println("reshape"+x+" "+y+" "+width+" "+height);
		super.reshape(x,y,width,height);
		if (glc == null) return;
		glc.setSize(width,height);
		repaint();
	}
}
