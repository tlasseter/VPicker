package com.Sts.MVC.View3d;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: An abstract view subclassed by StsView3d and StsView2d defining views in StsGLPanel3d</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.Actions.*;
import com.Sts.Actions.Wizards.SensorXplot.*;
import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.Interfaces.StsViewSelectable;
import com.Sts.MVC.*;
import com.Sts.Reflect.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;
import java.awt.*;
import java.awt.event.*;
import java.lang.reflect.*;
import java.util.*;

abstract public class StsView implements /* GLEventListener, */ /* Externalizable, */ StsSerializable /*, Serializable */
{
    /** able to view: set in subclasses as needed */
    transient public boolean isViewable = true;
    /** 3d point where XOR cursor is being drawn */
    transient public float[] currentXYZ = null;
    /** previous 3d point which must be XOR redrawn */
    transient public float[] previousXYZ = null;
    /** Indicates mouse is being dragged in this window. Used to handle focus. */
    transient public boolean isCursorWindow = false;
    /** intialize(GLDrawable) has been called and executed. */
    transient public boolean isViewGLInitialized = false;
    /** intialize() has been called and executed. */
    transient public boolean isInitialized = false;
    /** button state of leftButton, used for cursor dsplay. */
    transient public int cursorButtonState = StsMouse.CLEARED;
    transient protected int keyCode = KeyEvent.VK_UNDEFINED;
	transient public int releasedKeyCode;
    /** current mouse mode for this view */
    transient protected int mouseMode = StsCursor.ZOOM;
    /** action performed if no current active action */
    transient public StsAction defaultAction = null;
    transient public StsGLPanel3d glPanel3d;
    transient public StsGLPanel glPanel; // for 2d
	transient public StsGLJPanel glJPanel; //
    transient protected StsModel model;
    transient protected GL gl;
    transient protected GLU glu;
	transient double fovy = 70;
    /** XOR rectangle mode */
    transient public boolean XORrectangle = false;
    transient public boolean paintOldRectangle = false;
    transient public boolean paintNewRectangle = false;
    transient public StsMousePoint XORoldMouseDelta = null;
    transient public StsMousePoint mouseDragDelta = null;
    transient public StsMousePoint mousePressedPoint = null;

    static public final int maxAllowableTextureSize = 1024;

    static public final byte AXIS_TYPE_NONE = 0;
    static public final byte AXIS_TYPE_DEPTH = 1;
    static public final byte AXIS_TYPE_TIME = 2;
    static public final byte AXIS_TYPE_VELOCITY = 3;
    static public final byte AXIS_TYPE_VELOCITY_PERCENT = 4;
    static public final byte AXIS_TYPE_MEASURED_DEPTH = 5;
    static public final byte AXIS_TYPE_VALUE = 6;
    static public final byte AXIS_TYPE_DISTANCE = 7;
    static public final byte AXIS_TYPE_TRACES = 8;
    static public final byte AXIS_TYPE_VELOCITY_STACKS = 9;
    static public final byte AXIS_TYPE_VELOCITY_PERCENT_STACKS = 10;

    static public int pixelsPerInch;

    public transient StsMenuItem yBtn = new StsMenuItem();
    public transient StsMenuItem xBtn = new StsMenuItem();
    public transient StsMenuItem zBtn = new StsMenuItem();
    static
    {
        pixelsPerInch = Toolkit.getDefaultToolkit().getScreenResolution();
    }

	/** provides debug for mouse and key interactions controlling the cursor motion (ROTATE, PAN, etc) */
	static public final boolean mouseKeyDebug = false;

    /**
     * object being viewed is changed. Returns true if view is interested.
     * subclass should prepare for a redraw by clearing texture, etc.
     * Implement as needed in concrete subclasses.
     */
    abstract public boolean viewObjectChanged(Object source, Object object); // { return false; }

    /**
     * each subclass must define the axis type so windows can be moved in lock-step; i.e., if two views in a family have
     * the same (or a convertible) axis type, moving one moves the other to the same new range.
     *
     * @return axisType
     */
    abstract public byte getHorizontalAxisType();

    /**
     * each subclass must define the axis type so windows can be moved in lock-step; i.e., if two views in a family have
     * the same (or a convertible) axis type, moving one moves the other to the same new range.
     *
     * @return axisType
     */
    abstract public byte getVerticalAxisType();

    /** compute pixels/unit in X and Y directions. Override in subclasses as needed.  Not used in vie3d */
    public void computePixelScaling()
    {

    }

    /**
     * Every view must be able to classInitialize a defaultAction (which could be null) on the requesting actionManager.
     * Override in subclasses as needed.
     */
    public boolean initializeDefaultAction()
    {
        return false;
    }

//    abstract public boolean initializeDefaultAction(StsActionManager actionManager);

    /** Default constructor */
    public StsView()
    {
    }

    /**
     * StsView constructor
     *
     * @param glPanel3d graphics context
     */
    public StsView(StsGLPanel3d glPanel3d)
    {
        if(glPanel3d == null) return;
        this.glPanel3d = glPanel3d;
        this.glPanel = glPanel3d;
        this.model = glPanel3d.model;
//        initGL();
//		this.name = getViewName();
//        gl = glPanel3d.getGL();
        //       glu = glPanel3d.getGLU();
    }

    /** Initialize the view */
    //public abstract void classInitialize();
    public void displayChanged(GLAutoDrawable drawable, boolean modeChanged, boolean displayChanged)
    { }

    /** Compute the projection matrix */
    public void computeProjectionMatrix()
    { }

    public void reshape(int x, int y, int width, int height) { }

    public StsActionManager getActionManager()
    {
        if(glPanel == null) return null;
        return glPanel.actionManager;
    }
    /** Compute the model view matrix */
    protected void computeModelViewMatrix() { }

    /** Clear all graphics */
    public void clearToBackground(int clearBits)
    {
        gl.glDrawBuffer(GL.GL_BACK);
        StsColor backgroundColor = getBackgroundColor();
        glPanel.setClearColor(backgroundColor);
        glPanel.applyClearColor(clearBits);
    }

    public StsColor getBackgroundColor()
    {
        return model.getProject().getBackgroundColor();
    }

    public StsColor getForegroundColor()
    {
        return model.getProject().getForegroundColor();
    }

    /** Set view port changed */
    protected void viewPortChanged() { }

    /** Display graphics */
    public void display(GLAutoDrawable component) { }

    /** XOR mode */
    public void drawXORRectangle() {}

    /** Set default view */
    public void setDefaultView() { }

    /** Move graphics */
    public void moveWindow(StsMouse mouse) { }

    /** Stretch vertically */
    public void verticalStretch(StsMouse mouse) { }

    /** Shrink vertically */
    public void verticalShrink(StsMouse mouse) { }

    /** Stretch horizontally */
    public void horizontalStretch(StsMouse mouse) { }

    /** Shrink horizontally */
    public void horizontalShrink(StsMouse mouse) { }

    /** Redraw whole window */
    protected void drawCompleteWindow() throws StsException
    {
        throw new StsException();
    }
    /** Only draw foreground objects */
    //    protected void drawForeground() { }
    /** Change view parameters */
    public boolean changeModelView3d(float[] parameters) { return true; }
    public boolean changeModelView2d(float[][] axisRange) { return true; }

    /** Key pressed handler */
    public boolean keyPressed(StsMouse mouse, KeyEvent e)
    {
		int keyCode = e.getKeyCode();
		if(this.keyCode == keyCode) return false;
		if(mouseKeyDebug) StsException.systemDebug(this, "keyPressed", "keyCode changed from " + getKeyText(this.keyCode) + " to " + getKeyText(keyCode));
		this.keyCode = keyCode;
        return false;
    }

    /** Key released handler */
    public void keyReleased(StsMouse mouse, KeyEvent e)
    {
        int keyCode = KeyEvent.VK_UNDEFINED;
   		if(this.keyCode == keyCode) return;
        if(mouseKeyDebug) StsException.systemDebug(this, "keyReleased", "keyCode changed from " + getKeyText(this.keyCode) + " to " + getKeyText(keyCode));
		this.keyCode = keyCode;
    }

    /**
     * Show the popup menu
     *
     * @param mouse
     */
    public void showPopupMenu(StsMouse mouse)
    {}

    /**
     * Show popup specific to a selected object menu
     *
     * @param mouse
     * @param object
     */
    public void showObjectPopupMenu(StsMouse mouse, StsObject object)
    {}
	public void cancelPopupMenu(StsMouse mouse) {}
    /** Edit object under mouse */
    public StsObject defaultMousePressed(StsMouse mouse) { return null; }

    /** Edit object under mouse */
    public void defaultMouseMotion(StsMouse mouse) { }

    /** Is the view viewable */
    public boolean isViewable() { return true; }

    /** Clear data display */
    public void clearDataDisplay() { }

    /** Clear texture display on next redraw cycle for this view */
    public void clearTextureDisplay() { }

    /** remove this class from display in this view */
    public void removeDisplayableClass(StsObject object) { }

    /** View Changed? */
    public void viewChanged() { }

    /**
     * object being viewed is changed. Repaint this view if affected.
     * Implement as needed in concrete subclasses.
     */
    public boolean viewObjectRepaint(Object source, Object object) { return false; }

    /** Cursor3d changed; notify view */
    public void adjustCursor(int dir, float dirCoor) { }

    /** Set view to object */
    public boolean setViewObject(Object object) { return false; }

    /** Log the mouse position to information panel */
    public StsCursorPoint logReadout(StsGLPanel3d glPanel3d, StsMouse mouse) { return null; }

    public int getMouseMode() { return mouseMode; }

    public void setMouseModeZoom() { setMouseMode(StsCursor.ZOOM); }

    public void setMouseMode(int mouseMode)
	{
		if(mouseKeyDebug)
		{
			if(this.mouseMode != mouseMode)
				StsException.systemDebug(this, "setMouseMode", " mouseMode changed from " + StsCursor.mouseModeStrings[this.mouseMode] + " to " + StsCursor.mouseModeStrings[mouseMode]);
			else
				StsException.systemDebug(this, "setMouseMode", " mouseMode reset to same " + StsCursor.mouseModeStrings[this.mouseMode]);
		}
		this.mouseMode = mouseMode;
	}

    public void setMouseModeFromToolbar(int mouseMode)
	{
		if(mouseKeyDebug)
		{
			if(this.mouseMode != mouseMode)
				StsException.systemDebug("MouseActionToolbar.setMouseMode() mouseMode changed from " + StsCursor.mouseModeStrings[this.mouseMode] + " to " + StsCursor.mouseModeStrings[mouseMode]);
			else
				StsException.systemDebug("MouseActionToolbar.setMouseMode() mouseMode reset to same " + StsCursor.mouseModeStrings[this.mouseMode]);
		}
		this.mouseMode = mouseMode;
	}
    /** classInitialize from DB restored state * */
    public void initializeTransients(StsGLPanel3d glPanel3d)
    {
        this.glPanel3d = glPanel3d;
        this.glPanel = glPanel3d;
        this.model = glPanel3d.model;
//        gl = glPanel3d.getGL();
//        glu = glPanel3d.getGLU();
        // reinitializeGL();
        initializeDefaultAction();
    }

    public void init(GLAutoDrawable drawable)
    {
        if(isViewGLInitialized) return;
        if(glPanel == null) return;
        gl = glPanel.getGL();
        glu = glPanel.getGLU();
        if(Main.isGLDebug)
        {
            GLContext glContext = drawable.getContext();
            System.out.println("GLContext: " + glContext.toString());
            System.out.println("GLAutoDrawable: " + drawable.toString());
        }
        if(StsGLPanel.glStateDebug) System.out.println("init(GLAutoDrawable) called.");
//        glu = new GLU();
        StsJOGLShader.initialize(gl);
        // glPanel.setBackgroundColor(glPanel.clearColor.getColor()); // jbw
        glPanel.maxTextureSize = getMaxTextureSize() / 2;
        isViewGLInitialized = true;
//        StsTextureList.deleteAllTextures(gl);
    }

    public int getMaxTextureSize()
    {
        int[] sizeArray = new int[1];
        gl.glGetIntegerv(GL.GL_MAX_TEXTURE_SIZE, sizeArray, 0);
        return Math.min(sizeArray[0], maxAllowableTextureSize);
    }

    public void reinitializeGL()
    {
        isViewGLInitialized = false;
    }

    /** Reset the insets to the entire graphics screen. Eliminating any labeling and filling with graphics. */
    public void resetViewPort()
    {
        glPanel.resetViewPort();
    }

    /**
     * Get the key currently being pressed
     *
     * @return key
     */
    public int getCurrentKeyCodePressed()
    { return getKeyCode(); }

    /** Set the current key to null */
    public void clearCurrentKeyPressed()
    {
        setKeyCode(KeyEvent.VK_UNDEFINED);
    }

    static public StsView constructSubclass(Class viewClass, StsGLPanel3d glPanel3d)
    {
        try
        {
            Class[] argTypes = new Class[]{StsGLPanel3d.class};
            Constructor constructor = viewClass.getDeclaredConstructor(argTypes);
            Object[] args = new Object[]{glPanel3d};
            return (StsView)constructor.newInstance(args);
        }
        catch(Exception e)
        {
            StsException.systemError("StsView.constructSubclass() " + viewClass.getName() + " failed.  Error: " + e.getMessage());
            return null;
        }
    }

    /**
     * Move the other two cursor planes to the point selected on this cursor plane
     *
     * @param mouse mouse object
     */
    public boolean moveCursor3d(StsMouse mouse, StsGLPanel3d glPanel3d)
    {
        return false;
    }

    /**
     * Make the current view look exactly like the other view
     *
     * @param otherView the other view
     */
    public StsView copy(StsView otherView)
    {

//        name = otherView.name;
//        glPanel3d.actionManager.checkStartAuxiliaryWindowAction(otherView);
        return this;
    }

    //   public StsGLPanel3d getGlPanel3d() { return glPanel3d; }

    public StsWin3dBase getWindow() { return glPanel3d.window; }

    public StsWin3dFull getParentWindow() { return glPanel3d.window.getParentWindow(); }

    /*
        public void compressColorscale()
        {
            if((this instanceof StsView3d) || (this instanceof StsViewCursor))
                glPanel3d.cursor3d.getCurrentSeismicCursorSection().
                    line2d.getColorscale().compressColorscale(5);
            else
                ((StsViewXP)this).getCrossplot().getColorscale().compressColorscale(5);
        }

        public void uncompressColorscale()
        {
            if((this instanceof StsView3d) || (this instanceof StsViewCursor))
                glPanel3d.cursor3d.getCurrentSeismicCursorSection().
                    line2d.getColorscale().uncompressColorscale(5);
            else
                ((StsViewXP)this).getCrossplot().getColorscale().compressColorscale(5);
        }
    */
    public boolean moveLockedWindows()
    {
        StsWin3dBase window = glPanel3d.window;
        StsWin3dBase parentWindow = window.parentWindow;
        // If it is the main window parent = window
//        if (parentWindow == null) parentWindow = window;
        boolean familyLocked = parentWindow.isLocked;
        if(!familyLocked) // if family not locked, move the views only in this window
            return moveLockedViews(window);
        else
        {
            boolean moved = false;
            StsWin3dBase[] windows = parentWindow.getFamilyWindows();
            for(StsWin3dBase w : windows)
                moved = moved | moveLockedViews(w);
            return moved;
        }
    }

    private boolean moveLockedViews(StsWin3dBase window)
    {
        boolean moved = false;
        StsView[] displayedViews = window.getDisplayedViews();
        for(StsView displayedView : displayedViews)
        {
            if(displayedView == this) continue;
            if(displayedView == null || displayedView.glPanel3d == null) continue;
            if(displayedView.moveWithView(this)) moved = true;
        }
        return moved;
    }

    /** we have created a new view; match it's range to existing views */
    protected boolean matchLockedWindows()
    {
        boolean moved = false;
        StsWin3dBase window = glPanel3d.window;
        StsWin3dBase parentWindow = window.parentWindow;
        // If it is the main window parent = window
        if(parentWindow == null)
            parentWindow = window;
        if(parentWindow == null || !parentWindow.isLocked) return moved;
        StsWin3dBase[] familyWindows = parentWindow.getFamilyWindows();
        int nWindows = familyWindows.length;
        for(int n = 0; n < nWindows; n++)
        {
            StsWin3dBase otherWindow = familyWindows[n];
            if(otherWindow == window) continue;
            StsView[] displayedViews = familyWindows[n].getDisplayedViews();
            for(StsView displayedView : displayedViews)
            {
                if(displayedView == this) continue;
                if(displayedView == null || displayedView.glPanel3d == null) continue;
                if(moveWithView(displayedView)) moved = true;
            }
        }
        return moved;
    }

    /** override in subclasses as needed */
    public boolean moveWithView(StsView movedView)
    {
        return false;
    }

    public void setDefaultMouseMode()
    {
        setMouseModeZoom();
        //   StsMouseActionToolbar toolbar = glPanel3d.window.getMouseActionToolbar();
        //    if (toolbar != null) toolbar.zoom();
    }

    public void setCursorXOR(StsGLPanel3d glPanel3d, StsMouse mouse, float[] xyz)
    {
        int currentCursorButtonState = mouse.getButtonState(StsMouse.LEFT);
        setCursorXOR(currentCursorButtonState, xyz, true);

        if(model.win3d.isMouseTracking)
        {
            StsWindowFamily windowFamily = model.getWindowFamily(glPanel3d.window);
            Iterator<StsView> viewIterator = windowFamily.getWindowViewIterator();
            while(viewIterator.hasNext())
            {
                StsView view = viewIterator.next();
                view.setCursorXOR(currentCursorButtonState, xyz, false);
            }
        }
    }

    public void setCursorXOR(int cursorButtonState, float[] xyz, boolean isCursorWindow)
    {
        this.cursorButtonState = cursorButtonState;
        this.isCursorWindow = isCursorWindow;
        // System.out.println("XOR Window: " + glPanel3d.window.getName() + " Mouse state: " + StsMouse.stateLabels[cursorButtonState]);
        if(cursorButtonState == StsMouse.PRESSED)
        {
            previousXYZ = null;
            currentXYZ = xyz;
            setXOR(true);
        }
        else if(cursorButtonState == StsMouse.DRAGGED)
            currentXYZ = xyz;
        else if(cursorButtonState == StsMouse.RELEASED)
            setXOR(false);
        glPanel3d.repaint();
    }

    /** If we want XOR (on = true), turn off clear and swap flags */
    private void setXOR(boolean on)
    {
//        glPanel3d.glc.setNeedClear(!on);
//        glPanel3d.glc.setNeedSwap(!on);
        glPanel3d.gld.setAutoSwapBufferMode(!on);
    }

    public void setXORRectangle(boolean on)
    {
        this.XORrectangle = on;
    }

    public boolean getXORRectangle()
    {
        return this.XORrectangle;
    }

    public void setMotionCursor()
    {
        glPanel.cursor.setStsCursor(mouseMode);
    }

    public boolean setMotionCursor(int mouseMode)
    {
		if(this.mouseMode == mouseMode) return false;
        setMouseMode(mouseMode);
        glPanel.cursor.setStsCursor(mouseMode);
		return true;
    }

    public boolean isLeftMouseDownInMainWindow()
    {
        StsWin3dBase mainWindow = model.getWindowFamilyParent(getWindow());
        if(mainWindow == null) return false;
        boolean downInWindow = mainWindow.hasFocus() && StsGLPanel.isLeftMouseDown();
        boolean downOnCursor3d = mainWindow.getCursor3d().getIsDragging();
        return downInWindow || downOnCursor3d;
    }

    public boolean isCursor3dDragging()
    {
        StsWin3dFull parentWindow = glPanel3d.window.getParentWindow();
        if(parentWindow == null) return false;
        if(parentWindow.getCursor3d().getIsDragging()) return true;
        if(!StsGLPanel.isLeftMouseDragging()) return false;
        // if left mouse is down in a 3d window, we are probably selecting a new view location which is equivalent to dragging the cursor3d selector
        StsWindowFamily windowFamily = getWindow().getWindowFamily();
        Iterator viewIterator = windowFamily.getWindowViewIteratorOfType(StsView3d.class);
        while(viewIterator.hasNext())
        {
            StsView view = (StsView)viewIterator.next();
            if(view.hasFocus()) return true;
        }
        return false;
    }

    public boolean hasFocus()
    {
        return glPanel.hasFocus();
    }

    /*
    * custom serialization requires versioning to prevent old persisted files from barfing.
    * if you add/change fields, you need to bump the serialVersionUID and fix the
    * reader to handle both old & new
    */
    static final long serialVersionUID = 1l;

    public String getViewName()
    {
        return (String)StsMethod.invokeStaticMethod(getClass(), "getStaticViewName", new Class[0]);

    }

    public String getName() { return getShortViewName(); }

    public String getShortViewName()
    {
        return (String)StsMethod.invokeStaticMethod(getClass(), "getStaticShortViewName", new Class[0]);
    }

    /** override in concrete subclassses */
    static public String getStaticViewName()
    { return ""; }

    /** override in concrete subclassses */
    static public String getStaticShortViewName()
    { return ""; }

    public void setGlPanel3d(StsGLPanel3d glPanel3d)
    {
        this.glPanel3d = glPanel3d;
    }

    static public StsViewItem[] constructViewList(StsModel model)
    {
        StsArrayList viewClasses = model.getDisplayableViewClasses();
        int nViewClasses = viewClasses.size();
        StsViewItem[] viewItems = new StsViewItem[nViewClasses + 2];
        viewItems[0] = new StsViewItem(StsView3d.class);
        viewItems[1] = new StsViewItem(StsViewCursor.class);
        for(int i = 0, n = 2; i < nViewClasses; i++, n++)
        {
            Class viewClass = (Class)viewClasses.get(i);
            viewItems[n] = new StsViewItem(viewClass);
        }
        return viewItems;
    }

    /** Repaint the graphics screen */
    public void viewChangedRepaint()
    {
        if(glPanel == null) return;
        glPanel.viewChanged = true;
        glPanel.repaint();
    }

    public void repaint()
    {
        if(glPanel == null) return;
        glPanel.repaint();
    }

    public void setDefaultAction(StsAction action)
    {
        if(StsActionManager.actionDebug) StsException.systemDebug(this, "setDefaultAction", action.getName());
        defaultAction = action;
    }

    public StsAction getDefaultAction() { return defaultAction; }

    public boolean sameVerticalAxisType(StsView view)
    {
        return getVerticalAxisType() == view.getVerticalAxisType();
    }

    public boolean sameHorizontalAxisType(StsView view)
    {
        return getHorizontalAxisType() == view.getHorizontalAxisType();
    }

    static public int getPixelsPerInch() { return pixelsPerInch; }


    public void initAndStartGL()
    {
        if(glPanel == null) return;
        glPanel.initAndStartGL(glPanel.getGraphicsDevice());
    }

    public double getFovy()
	{
		return fovy;
	}

    public void performMouseAction(StsActionManager actionManager, StsMouse mouse)
    {
        try
        {
            //	if (e.getID() == MouseEvent.MOUSE_WHEEL) {
            //	   glPanel.moveWindow(mouse);
//      //          System.out.println("Right Mouse - Module: " + Main.usageModule + " Message: " + Main.usageMessage);
            //       Main.logUsageTimer("com.Sts.MVC.StsGLPanel","Graphic Actions");
            //	}

            StsAction currentAction = actionManager.getCurrentAction();

            int but = mouse.getCurrentButton();
            {
                /** If there is an active listener for leftButton, notify it if any left mouse action	*/
                if (but == StsMouse.LEFT)
                {
                    if(currentAction == null || (currentAction instanceof StsSensorXplotWizard))
                    {
                        int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
                        if(leftButtonState == StsMouse.PRESSED)
                        {
                            defaultMousePressed(mouse);
                        }
                        else if(leftButtonState == StsMouse.RELEASED) // jbw
                        {
                            //StsViewSelectable stsObject = getSelectableObject(mouse);
                            //if(stsObject != null)
                            //    stsObject.showPopupMenu(glPanel, mouse);
                            //else
                              if(defaultAction != null)
                                defaultAction.performMouseAction(mouse, glPanel);
                        }
                        // jbw need to process a down-drag for well panel cursor;
                        else if(defaultAction != null && leftButtonState != StsMouse.DRAGGED) // jbw
                        {
                            defaultMouseMotion(mouse);
                            defaultAction.performMouseAction(mouse, glPanel);
                        }
                    }
                    else if(!currentAction.isStarted())
                        Toolkit.getDefaultToolkit().beep();
                    else if(!currentAction.performMouseAction(mouse, glPanel))
                    {
                        actionManager.abortCurrentAction();
                        actionManager.removeCurrentAction();
                    }
                    else
                        actionManager.logUsageTimer(currentAction);
                }


                /** If any right mouse action, move view */
                if (but == StsMouse.VIEW || (mouse.getButtonStateCheckClear(StsMouse.NONE) == StsMouse.WHEEL))
                {
                    int buttonState = mouse.getButtonStateCheckClear(StsMouse.VIEW);
                    if(buttonState == StsMouse.DRAGGED)
                        setMotionCursor();
                    else if(buttonState == StsMouse.RELEASED)
                    {
                        glPanel.restoreCursor();
                        repaint();  // needed if mouse change changes how objects are drawn
                    }
                    moveWindow(mouse);

                    if((buttonState == StsMouse.WHEEL || mouse.getButtonStateCheckClear(StsMouse.NONE) == StsMouse.WHEEL))
                    {
						glPanel.restoreCursor();
						mouse.clearButtonState(StsMouse.NONE);
						mouse.clearButtonState(StsMouse.VIEW);
					}
//                    System.out.println("Right Mouse - Module: " + Main.usageModule + " Message: " + Main.usageMessage);
                    Main.logUsageTimer("com.Sts.MVC.StsGLPanel", "Graphic Actions");
//                    Main.logUsageTimer();
                }


                /** If middle mouse button clicked, terminate any active function. 	*/
                /** If none active, trigger pop-up menu								*/
                if (but == StsMouse.POPUP)
                {
                    /*
                         if( mouse.getButtonState(StsMouse.MIDDLE) == StsMouse.RELEASED )
                             endCurrentAction();
                         else
                 */
                    int buttonState = mouse.getButtonStateCheckClear(StsMouse.POPUP);

                    if(buttonState == StsMouse.PRESSED)  // jbw to pressed
                    {
//                        newPopupObject = null;
//                        if(model != null)
//                            newPopupObject = model.objectPopup(mouse);
//                        if(newPopupObject == null)
                        // StsObject stsObject = mouseSelectedEdit(mouse);
                        showPopupMenu(mouse);
                        mouse.clearButtonState(StsMouse.POPUP);
                    }
					else if(buttonState == StsMouse.RELEASED)  // jbw to pressed
					{
						// StsObject stsObject = mouseSelectedEdit(mouse);
						// if(stsObject == null) cancelPopupMenu(mouse);
						mouse.clearButtonState(StsMouse.POPUP);
					}

                }

            }
        }
        catch(Exception ex)
        {
            StsException.outputWarningException(this, "performMouseAction", ex);
        }
        catch(OutOfMemoryError e2)
        {
            StsException.outputWarningException(this, "performMouseAction", "(out of memory)", e2);
        }
    }

     public StsViewSelectable getSelectableObject(StsMouse mouse)
     {
         ArrayList selectableClasses = model.selectable3dClasses;
          for(Object selectableClass : selectableClasses)
          {
            StsViewSelectable selectedObject = (StsViewSelectable)StsJOGLPick.pickVisibleClass3d(glPanel3d, (StsClass)selectableClass, StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_CLOSEST);
            if(selectedObject != null) return selectedObject;
          }
          return null;
     }

      public void changeActiveSlice(StsMenuItem btn)
      {
          StsWindowFamily family = model.getWindowFamily(this.getWindow());
          if (btn == xBtn)
              family.setSelectedDirection(StsCursor3d.XDIR, true);
          else if (btn == yBtn)
              family.setSelectedDirection(StsCursor3d.YDIR, true);
          else if (btn == zBtn)
              family.setSelectedDirection(StsCursor3d.ZDIR, true);
      }

      public int getGLMouseY(StsMouse mouse)
      {
        return glPanel.getGLMouseY(mouse.getY());
      }
/*
    public JMenuItem getViewEditMenuItem(StsWin3dBase win3d)
    {
        StsMenuItem viewEditBtn = new StsMenuItem();
        viewEditBtn.setMenuActionListener("Edit multiView panels...", this, "editMultiViewPanels", win3d);
        return viewEditBtn;
    }

    public void editMultiViewPanels(StsWin3d win3d)
    {
        editMultiViewPanels((StsWin3dBase)win3d);
    }

    public void editMultiViewPanels(StsWin3dBase win3d)
    {
        StsViewItem[] viewItems = StsView.constructViewList(model, win3d);
        Object[] selectedViewItems = StsObjectTransferDialog.getDialogSelectedObjects(win3d, "Select views", viewItems);
        win3d.constructViewPanel(selectedViewItems);
    }

    class EditMultiViewPanelsDialog extends JDialog
    {
        EditMultiViewPanelsDialog(Frame frame, String title, boolean modal)
        {
            super(frame, title, modal);
        }
    }
*/

	/** keyboard code currently pressed */
	public int getKeyCode()
	{
		return keyCode;
	}

	public void setKeyCode(int keyCode)
	{
		if(this.keyCode == keyCode) return;
		if(mouseKeyDebug)
				StsException.systemDebug(this, "setKeyCode", "keyCode changed from " + getKeyText(this.keyCode) + " to " + getKeyText(keyCode));
		this.keyCode = keyCode;
	}

	protected String getKeyText(int keyCode)
	{
		if(keyCode == 0) return "NONE";
		else return KeyEvent.getKeyText(keyCode);
	}
}
