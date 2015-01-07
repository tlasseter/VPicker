package com.Sts.MVC.View3d;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: Class defining three-dimensional view. This view includes all visible data objects and
 * the 3D cursor planes. </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.Actions.*;
import com.Sts.Collaboration.*;
import com.Sts.DB.DBCommand.*;
import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.Interfaces.StsView3dDisplayable;
import com.Sts.Interfaces.StsViewSelectable;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Types.PreStack.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;
import com.magician.fonts.GLBitmapFont;
import com.magician.fonts.GLHelvetica18BitmapFont;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;


public class StsView3d extends StsView implements StsSerializable
{
	/** lightedSurfaces: surfaces are to be lighted */
	public boolean lightedSurfaces;
	/** centerPoint location (for gluPerspective) */
	public StsPoint centerViewPoint = new StsPoint(3);
    /** current screen location of compass */
	public double[] compassScreenPoint = new double[] { 0.1, 0.9, 0.9 };
	/** Distance from eye view point to centerViewPoint */
	public float distance;
	/** view azimuth angle */
	public float azimuth;
	/** view inclination angle */
	public float elevation;
	/** Min distance from eye view point to centerViewPoint */
	public float minViewDistance;
	/** distances to near clipping plane */
	public double clipNear;
	/** distances to far clipping plane */
	public double clipFar;
	/** clipFar/clipNear ratio; */
	public double clipFarOverNear;
	/** perspective y-direction field of view */
	public double perspectiveFovy;
	/** perspective y-direction/x-direction aspect ratio */
	//    public double perspectiveAspect;
	/** horizontal deg/pixel used for window motion */
	public float degPerXPixel;
	/** vertical deg/pixel used for window motion */
	public float degPerYPixel;
	private float zScale = 1.0f;
	/** multiplier for MouseDelta.y to distance changes */
	public float viewDistanceMult;
	/** xyz point where eye is located */
	public StsPoint eyePoint;
	/** classInitialize a weighpoint */
	StsWeighPoint wp = null;

	/** Material properties */
	transient public int solidMaterial, transparentMaterial;
	/** Current material type is solid (not transparent) */
	transient public boolean isSolid = false;
	/** isRGB: current mode is RGB (not colorMap) */
	transient public boolean isRGB;
	/** is3dOverlay: drawing in 3d window is in overlay mode */
	transient public boolean is3dOverlay;
	/** drawing3dOverlay: true if 1st pass (background saved) */
	transient public boolean drawing3dOverlay = false;
	/** Light at view position (moves with viewer). */
	transient public float viewLightPosition[] =
			{0.0f, 0.15f, 0.3f, 0.0f};
	transient public float viewLightPosition2[] =
			{0.0f, 0.15f, -0.3f, 0.0f};
	/** Fixed light above origin. */
	transient public float fixedLightTopPosition[] =
			{0.0f, 0.15f, 1.0f, 0.0f};
	transient public float fixedLightTopPosition2[] =
			{0.0f, 0.15f, -1.0f, 0.0f};
	/** 4x4 modelView matrix (local copy) */
	transient public double[] unrotatedModelViewMatrix;
    	/** compass data */
	transient Compass compass;

	/** angle from unrotated to rotated view: used in panning window */
	//	transient public float rotationAngle = 0.0f;

	{
		/** default window mouse mode is ZOOM; for 3d windows make it ROTATE
		 *  mode for this view is maintained throughout the session so as we flip between windows,
		 *  the mode for that view will be restored.
		 */
		setMouseMode(StsCursor.ROTATE);
	}

	/** Timer for draw routines. */
	private transient StsTimer timer = new StsTimer();
	private transient boolean runTimer = false;
	private transient StsCursor cursor = null;

    private transient StsMenuItem propertyWindow = new StsMenuItem();
    private transient StsMenuItem classPropertyWindow = new StsMenuItem();
	private transient StsMenuItem addWeighPointBtn = new StsMenuItem();
	private transient StsMenuItem topDownViewBtn = new StsMenuItem();
	private transient StsMenuItem normalViewBtn = new StsMenuItem();
	private transient StsMenuItem eastViewBtn = new StsMenuItem();
    private transient StsMenuItem northViewBtn = new StsMenuItem();    
	private transient StsMenuItem deleteBtn = new StsMenuItem();
	private transient JPopupMenu popupMenu = null;
	/*
		 private transient Menu captureBtn = new Menu("Capture Image");
			private transient StsMenuItem captureDesktopBtn = new StsMenuItem();
			private transient StsMenuItem captureWindowBtn = new StsMenuItem();
			private transient StsMenuItem captureGraphicBtn = new StsMenuItem();
		  */
	private float[] viewParams; // for persistence

	static final boolean debug = false;

	static public final String viewName3d = "3D View";
	static public final String shortViewName3d = "3D";

	/** Default constructor */
	public StsView3d()
	{
	}

	/**
	 * StsView3d constructor
	 *
	 * @param glPanel3d the graphics context
	 */
	public StsView3d(StsGLPanel3d glPanel3d)
	{
		super(glPanel3d);
		initializeNew();
	}

	static public String getStaticViewName()
	{
		return viewName3d;
	}

	static public String getStaticShortViewName()
	{
		return shortViewName3d;
	}

	/**
	 * Initialize the 3d view plot by setting the viewport and initializing the graphics materials.
	 *
	 * @param drawable
	 * @params component the drawable component
	 */
	public void init(GLAutoDrawable drawable)
	{
		if(isViewGLInitialized) return;
		super.init(drawable);
		if(Main.isDrawDebug)
		{
			StsException.systemDebug(this, "init", "initialize view. GL: " + gl.toString());
		}
		//         initializeView();
	}

	protected void initializeView()
	{
		gl.glDrawBuffer(GL.GL_BACK);
		gl.glShadeModel(GL.GL_SMOOTH);
		gl.glEnable(GL.GL_NORMALIZE);
		gl.glEnable(GL.GL_DEPTH_TEST);
		gl.glDepthMask(true);
		gl.glEnable(GL.GL_ALPHA_TEST);
		gl.glAlphaFunc(GL.GL_GREATER, 0.1f);

		//         glPanel3d.setViewPort();

		is3dOverlay = false;

		initializeMaterials();
		initializeLighting();

		/** lightedSurfaces will be selected by user.	*/
		/** isRGB indicates light is enabled.			*/

		enableLighting();
		isRGB = true;
		lightedSurfaces = true;

		computeProjectionMatrix();
		setViewLights();
		computeModelViewMatrix();
		setFixedLights();

		viewParams = getCenterAndViewParameters();

		if(addWeighPointBtn.getActionListeners().length == 0)
		{
            if(Main.viewerOnly)
            {
                classPropertyWindow.setMenuActionListener("Group Properties...", model.getProject(), "popupClassPropertyPanel", null);
                propertyWindow.setMenuActionListener("Display Properties...", model.getProject(), "popupPropertyPanel", null);
            }
			addWeighPointBtn.setMenuActionListener("Add WayPoint...", this, "addWeighPoint", null);
			topDownViewBtn.setMenuActionListener("Top-Down", this, "prepackagedViews", topDownViewBtn);
			normalViewBtn.setMenuActionListener("Normal", this, "prepackagedViews", normalViewBtn);
			northViewBtn.setMenuActionListener("North", this, "prepackagedViews", northViewBtn);
            eastViewBtn.setMenuActionListener("East", this, "prepackagedViews", eastViewBtn);
			xBtn.setMenuActionListener("X / Crossline", this, "changeActiveSlice", xBtn);
			yBtn.setMenuActionListener("Y / Inline", this, "changeActiveSlice", yBtn);
			zBtn.setMenuActionListener("Time / Depth", this, "changeActiveSlice", zBtn);
		}
		/*
						 captureDesktopBtn.setMenuActionListener("Desktop...", glPanel3d.window, "outputImage", new Integer(StsWin3dBase.DESKTOP));
						 captureWindowBtn.setMenuActionListener("Window...", glPanel3d.window, "outputImage", new Integer(StsWin3dBase.WINDOW) );
						 captureGraphicBtn.setMenuActionListener("Graphic...", glPanel3d.window, "outputImage", new Integer(StsWin3dBase.GRAPHIC));

						 captureBtn.add(captureDesktopBtn);
						 captureBtn.add(captureWindowBtn);
						 captureBtn.add(captureGraphicBtn);
				  */
        compass = new Compass(compassScreenPoint);

		isViewGLInitialized = true;
	}

	/** Initialize the 3D cursor and set to a default view */
	public void initializeNew()
	{
		if(isInitialized) return;
		setDefaultView(false);
		cursor = new StsCursor(glPanel3d);
		//        initializeDefaultAction();
		isInitialized = true;
	}

	// (re) classInitialize from database restore
	public void initializeTransients(StsGLPanel3d glPanel3d)
	{
		super.initializeTransients(glPanel3d);
		//System.out.println("re-classInitialize "+this);
		cursor = new StsCursor(glPanel3d);

		float[] vp = new float[7];
		if(viewParams != null)
			System.arraycopy(viewParams, 0, vp, 0, 7);
		setDefaultView(false);
		//        initializeDefaultAction();
		isInitialized = true;

		if(viewParams != null)
		{
			setViewCenter(new StsPoint(vp[0], vp[1], vp[2]));

			setViewDistance(vp[3]);
			setViewAzimuth(vp[4]);
			setViewElevation(vp[5]);
			setZScale(vp[6]);
		}
		//setCenterAndViewParameters(viewParams);
	}

	/** Compute the graphics projection matrix */
	public void computeProjectionMatrix()
	{
		gl.glMatrixMode(GL.GL_PROJECTION);
		if(StsGLPanel3d.debugProjectionMatrix)
			gl.glGetDoublev(GL.GL_PROJECTION_MATRIX, glPanel3d.debugMatrix, 0);
		gl.glLoadIdentity();

		float width = glPanel3d.getWidth();
		float height = glPanel3d.getHeight();
		double aspectRatio = (double) width / height;
		if(model.getProject().getIsPerspective())
		{
			glu.gluPerspective(perspectiveFovy, aspectRatio, clipNear, clipFar);
		}
		else
		{
			double orthoWidth = aspectRatio*distance;
			double orthoHeight = distance;
			gl.glOrtho(-orthoWidth, orthoWidth, -orthoHeight, orthoHeight, clipNear, clipFar);
		}

		if(glPanel3d.debugPicker && glPanel3d.debugPicking)
		{
			System.out.println("StsView3d: debugPicking: true");
			StsJOGLPick.initializeDebugPickMatrix(glPanel3d, StsJOGLPick.PICKSIZE_MEDIUM);
		}

		degPerXPixel = (float) (perspectiveFovy / width);
		degPerYPixel = degPerXPixel;
		/** Save the projectionMatrix locally */
		gl.glGetDoublev(GL.GL_PROJECTION_MATRIX, glPanel3d.projectionMatrix, 0);
		if(StsGLPanel3d.debugProjectionMatrix)
			glPanel3d.debugPrintProjectionMatrixChanged("StsView3d.computeProjectionMatrix()");
		gl.glMatrixMode(GL.GL_MODELVIEW);
	}

	protected void setViewLights()
	{
		/** replace this with a list of lights to be processed */
		gl.glLightfv(GL.GL_LIGHT0, GL.GL_POSITION, viewLightPosition, 0);
		// JBW SIMULATE 2 SIDED LIGHTING
		gl.glLightfv(GL.GL_LIGHT2, GL.GL_POSITION, viewLightPosition2, 0);
	}

	protected void setFixedLights()
	{
		/* replace this with a list of lights to be processed */
		gl.glLightfv(GL.GL_LIGHT1, GL.GL_POSITION, fixedLightTopPosition, 0);
		// JBW SIMULATE 2 SIDED LIGHTING
		gl.glLightfv(GL.GL_LIGHT3, GL.GL_POSITION, fixedLightTopPosition2, 0);
	}

	/** Initialize the graphics material */
	protected void initializeMaterials()
	{
		float[] solidAmbient = {0.2f, 0.2f, 0.2f, 1.0f};
		float[] solidDiffuse = {0.8f, 0.8f, 0.8f, 1.0f};
		//float[] solidSpecular = {0.8f, 0.8f, 0.8f, 1.0f};
		float[] solidSpecular = {0.5f, 0.5f, 0.5f, 1.0f};
		float solidShininess = 50.0f;

		float alphaValue = 0.3f;
		float[] transparentAmbient = {0.2f, 0.2f, 0.2f, 1.0f};
		float[] transparentDiffuse = {0.8f, 0.8f, 0.8f, alphaValue};
		//float[] transparentSpecular = {0.8f, 0.8f, 0.8f,  1.0f};
		float[] transparentSpecular = {1.0f, 1.0f, 1.0f, 1.0f};
		float transparentShininess = 100.0f;

		/** Make a solid material display list */
		solidMaterial = gl.glGenLists(1);
		gl.glNewList(solidMaterial, GL.GL_COMPILE);
		gl.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_AMBIENT, solidAmbient, 0);
		gl.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, solidDiffuse, 0);
		gl.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_SPECULAR, solidSpecular, 0);
		gl.glMaterialf(GL.GL_FRONT_AND_BACK, GL.GL_SHININESS, solidShininess);
		gl.glEndList();

		/** Make a transparent material display list */
		transparentMaterial = gl.glGenLists(1);
		gl.glNewList(transparentMaterial, GL.GL_COMPILE);
		gl.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_AMBIENT, transparentAmbient, 0);
		gl.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, transparentDiffuse, 0);
		gl.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_SPECULAR, transparentSpecular, 0);
		gl.glMaterialf(GL.GL_FRONT_AND_BACK, GL.GL_SHININESS, transparentShininess);
		gl.glEndList();

		/** enable material properties */
		gl.glEnable(GL.GL_COLOR_MATERIAL);
		gl.glColorMaterial(GL.GL_FRONT_AND_BACK, GL.GL_AMBIENT_AND_DIFFUSE);

		/** set current material to solid	*/
		gl.glCallList(solidMaterial);
		isSolid = true;
	}

	/** Set the solid materials for the graphics */
	public void setSolidMaterial()
	{
		gl.glCallList(solidMaterial);
		isSolid = true;
	}

	/** Set the transparent materials for the graphics */
	public void setTransparentMaterial()
	{
		gl.glCallList(transparentMaterial);
		isSolid = true;
	}

	/** Initialize the graphics lighting */
	protected void initializeLighting()
	{
		float[] lightModelAmbient = {0.0f, 0.0f, 0.0f, 1.0f};
		//float[] ambient = {0.5f, 0.5f, 0.5f, 1.0f};
		// JBW SIMULATE 2 SIDED -- half as much ambient from 2 active lights
		float[] ambient = {0.15f, 0.15f, 0.15f, 1.0f};
		float[] diffuse = {0.7f, 0.7f, 0.7f, 1.0f};
		float[] specular = {0.8f, 0.8f, 0.8f, 1.0f};

		/* Specify global ambient lighting model */
		gl.glLightModelfv(GL.GL_LIGHT_MODEL_AMBIENT, lightModelAmbient, 0);
		gl.glLightModeli(GL.GL_LIGHT_MODEL_TWO_SIDE, GL.GL_FALSE);
		gl.glLightModeli(GL.GL_LIGHT_MODEL_LOCAL_VIEWER, GL.GL_TRUE);
		gl.glLightModeli(GL.GL_LIGHT_MODEL_COLOR_CONTROL,
				GL.GL_SEPARATE_SPECULAR_COLOR);
		gl.glEnable(GL.GL_NORMALIZE);

		/* Specify light properties for each light	*/
		/* First light is infinite distance behind viewer */
		//gl.glEnable(GL.GL_LIGHT0);
		gl.glLightfv(GL.GL_LIGHT0, GL.GL_AMBIENT, ambient, 0);
		gl.glLightfv(GL.GL_LIGHT0, GL.GL_DIFFUSE, diffuse, 0);
		gl.glLightfv(GL.GL_LIGHT0, GL.GL_SPECULAR, specular, 0);
		// JBW SIMULATE 2 SIDED LIGHTING
		gl.glLightfv(GL.GL_LIGHT2, GL.GL_AMBIENT, ambient, 0);
		gl.glLightfv(GL.GL_LIGHT2, GL.GL_DIFFUSE, diffuse, 0);
		gl.glLightfv(GL.GL_LIGHT2, GL.GL_SPECULAR, specular, 0);
		/* Second light is infinite distance above origin */
		//gl.glEnable(GL.GL_LIGHT1);
		gl.glLightfv(GL.GL_LIGHT1, GL.GL_AMBIENT, ambient, 0);
		gl.glLightfv(GL.GL_LIGHT1, GL.GL_DIFFUSE, diffuse, 0);
		gl.glLightfv(GL.GL_LIGHT1, GL.GL_SPECULAR, specular, 0);
		// JBW SIMULATE 2 SIDED LIGHTING
		gl.glLightfv(GL.GL_LIGHT3, GL.GL_AMBIENT, ambient, 0);
		gl.glLightfv(GL.GL_LIGHT3, GL.GL_DIFFUSE, diffuse, 0);
		gl.glLightfv(GL.GL_LIGHT3, GL.GL_SPECULAR, specular, 0);

		enableLighting();
	}

	/** Turn the lights on */
	public void enableLighting()
	{
		gl.glEnable(GL.GL_LIGHTING);
		gl.glEnable(GL.GL_LIGHT0);
		gl.glEnable(GL.GL_LIGHT2);
		//		gl.glEnable(GL.GL_LIGHT1);
	}

	/** Turn the lights off */
	public void disableLighting()
	{
		gl.glDisable(GL.GL_LIGHTING);
		gl.glDisable(GL.GL_LIGHT0);
		gl.glDisable(GL.GL_LIGHT1);
		gl.glDisable(GL.GL_LIGHT2);
		gl.glDisable(GL.GL_LIGHT3);
	}

	/** Set material transparency @param makeSolid true solid */
	public void setMaterialTransparency(boolean makeSolid)
	{
		if(makeSolid == isSolid) return;

		if(makeSolid)
			gl.glCallList(solidMaterial);
		else
			gl.glCallList(transparentMaterial);

		isSolid = makeSolid;
	}

	protected void computeModelViewMatrix()
	{
		gl.glMatrixMode(GL.GL_MODELVIEW);

		gl.glGetDoublev(GL.GL_MODELVIEW_MATRIX, glPanel3d.debugMatrix, 0);
		gl.glLoadIdentity();

		gl.glTranslatef(0.0f, 0.0f, -distance);
		gl.glRotatef(-elevation, 1.0f, 0.0f, 0.0f);
		gl.glRotatef(-azimuth, 0.0f, 0.0f, 1.0f);
		gl.glScalef(1.0f, 1.0f, -getZScale());
		gl.glTranslatef(-centerViewPoint.v[0], -centerViewPoint.v[1], -centerViewPoint.v[2]);
		// Save the modelViewMatrix locally
		gl.glGetDoublev(GL.GL_MODELVIEW_MATRIX, glPanel3d.modelViewMatrix, 0);

		if(StsGLPanel.debugModelViewMatrix)
		{
			glPanel3d.debugPrintModelViewMatrixChanged("StsView3d.display() after modelViewChanged");
		}
		/*
				  if (debug)
				  {
					  System.out.println("computeModelViewMatrix distance " + this.distance);
					  System.out.println("computeModelViewMatrix center " + centerViewPoint.toString());
					  computeEyePoint(glu);
				  }
				  */
		// save an unrotated modelView matrix
		float angle = model.getProject().getAngle();
		if(angle == 0.0f) return;
		if(unrotatedModelViewMatrix == null) unrotatedModelViewMatrix = new double[16];
		gl.glPushMatrix();
		gl.glRotatef(angle, 0.0f, 0.0f, -1.0f);
		gl.glGetDoublev(GL.GL_MODELVIEW_MATRIX, unrotatedModelViewMatrix, 0);
		gl.glPopMatrix();
		viewParams = getCenterAndViewParameters();

	}

	private void computeEyePoint(GLU glu)
	{
		double[] eyeDouble = glPanel3d.gluUnProject(0.0, 0.0, 0.0, getUnrotatedModelViewMatrix(), glPanel3d.projectionMatrix, glPanel3d.viewPort);
		StsProject project = model.getProject();
		float[] rotatedXY = project.getRotatedRelativeXYFromUnrotatedRelativeXY((float) eyeDouble[0], (float) eyeDouble[1]);
		eyeDouble[0] = rotatedXY[0];
		eyeDouble[1] = rotatedXY[1];
		eyePoint = new StsPoint(eyeDouble);
		System.out.println("computeEyePoint " + eyePoint.toString());
		viewParams = getCenterAndViewParameters();
	}

	public void setTopViewModelViewMatrix(GL gl, float elevation, float azimuth, float[] xyz)
	{
		gl.glMatrixMode(GL.GL_MODELVIEW);
		gl.glPushMatrix();
		gl.glLoadIdentity();

		gl.glTranslatef(0.0f, 0.0f, -distance);
		gl.glRotatef(-elevation, 1.0f, 0.0f, 0.0f);
		gl.glRotatef(-azimuth, 0.0f, 0.0f, 1.0f);
		gl.glScalef(1.0f, 1.0f, -getZScale());
		if(xyz != null) gl.glTranslatef(-xyz[0], -xyz[1], -xyz[2]);
		viewParams = getCenterAndViewParameters();
	}

	public void popModelViewMatrix(GL gl)
	{
		gl.glPopMatrix();
	}

	protected double[] getUnrotatedModelViewMatrix()
	{
		//        return glPanel3d.modelViewMatrix;

		float angle = model.getProject().getAngle();
		if(angle == 0.0f) return glPanel3d.modelViewMatrix;
		else return unrotatedModelViewMatrix;
	}

	/*
		  protected void clearToBackground()
		  {
			  glPanel3d.setClearColor(model.getProject().getBackgroundStsColor());
			  gl.glDrawBuffer(GL.GL_BACK);
			  glPanel3d.applyClearColor();
			  gl.glClear(GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);
		  }
	   */
	protected void viewPortChanged()
	{
		int width = glPanel3d.getWidth();
	}

	public void reshape(int x, int y, int width, int height)
	{
		degPerXPixel = (float) (perspectiveFovy / width);
		degPerYPixel = degPerXPixel;
	}

	/*
			   public void reshape(GLAutoDrawable drawable, int x, int y, int w, int h)
			   {
				   System.out.println("reshape called");
			   }
		  */

	/** Display all 3D object and cursors */
	/*
		  *  Procedure for 3d overlay drawing
		  *
		  *      for overlay on 1st pass:
		  * 		    draw background (no wireframe)
		  * 		    lock zbuffer
		  * 		    draw foreground
		  *
		  *  	for overlay on subsequent passes
		  *  		draw foreground
		  *
		  *  	when overlay terminated:
		  *  		clear overlay
		  *  		unlock zbuffer
		  */
	public void display(GLAutoDrawable component)
	{
		if(glPanel.panelViewChanged)
		{
			initializeView();
			glPanel.panelViewChanged = false;
		}
		// debugGLState(component);
		if(glPanel3d.viewPortChanged)
		{
			glPanel3d.resetViewPort();
			computeProjectionMatrix();
			glPanel3d.viewPortChanged = false;
			glPanel3d.viewDistanceChanged = false;
			glPanel3d.projectionChanged = false;

		}
		else if(glPanel3d.viewDistanceChanged)
		{
			computeProjectionMatrix();
			glPanel3d.viewDistanceChanged = false;
			glPanel3d.projectionChanged = false;
		}
		else if(glPanel3d.projectionChanged)
		{
			computeProjectionMatrix();
			glPanel3d.projectionChanged = false;
		}
		if(glPanel3d.viewChanged)
		{
			computeModelViewMatrix();
			setFixedLights();
			glPanel3d.viewChanged = false;
		}
		// if (StsGLPanel.debugProjectionMatrix)
		//     glPanel3d.debugPrintProjectionMatrix("StsView3d.display() after projectionChanged");
		if(cursorButtonState != StsMouse.CLEARED)
		{
			// if this is window where cursor is being dragged and we have focus, draw foreground cursor.
			// If not the window where cursor is being dragged, but we are displaying cursor here,
			// draw the windows;
			if(isCursorWindow && glPanel3d.hasFocus() || !isCursorWindow)
			{
				drawForeground();
				if(cursorButtonState != StsMouse.CLEARED) return;
			}
		}
		viewParams = getCenterAndViewParameters();
		if(!isRGB && lightedSurfaces)
		{
			setWindowToRGB();
			enableLighting();
			clearToBackground(GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);
			//            glPanel3d.swapBuffers();
		}
		else if(isRGB && !lightedSurfaces)
		{
			setWindowToColorMap();
			disableLighting();
			clearToBackground(GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);
			//            glPanel3d.swapBuffers();
		}
		// if (StsGLPanel.debugProjectionMatrix)
		//     glPanel3d.debugPrintProjectionMatrix("StsView3d.display() before drawing overlay");
		if(runTimer) timer.start();
		if(is3dOverlay)
		{
			if(!drawing3dOverlay)
			{
				drawing3dOverlay = true;
				copyColorBuffer(GL.GL_FRONT);
				gl.glDrawBuffer(GL.GL_FRONT);
			}
			else
			{
				copyColorBuffer(GL.GL_BACK);
			}
			gl.glDepthMask(false);
			drawForegroundOverlay();
			gl.glFlush();
			if(runTimer) timer.stopPrint(" total overlay draw time");
		}
		else
		{
			if(drawing3dOverlay)
			{
				drawing3dOverlay = false;

				gl.glDrawBuffer(GL.GL_BACK);
				gl.glEnable(GL.GL_DEPTH_TEST);
				gl.glDepthMask(true);
			}
			gl.glDrawBuffer(GL.GL_BACK);

			try
			{
				// if (StsGLPanel.debugProjectionMatrix)
				//     glPanel3d.debugPrintProjectionMatrix("StsView3d.display() before drawCompleteWindow");
				//                System.out.println("Drawing complete window.");
				drawCompleteWindow();
				//                glPanel3d.swapBuffers();
				// if (StsGLPanel.debugProjectionMatrix)
				//     glPanel3d.debugPrintProjectionMatrix("StsView3d.display() after drawCompleteWindow");
			}
			catch(StsException ex)
			{
				ex.fatalError();
			}
			if(runTimer) timer.stopPrint(" total non-overlay draw time");
		}
	}

	private void debugGLState(GLAutoDrawable component)
	{
		System.out.println("GL " + debugToString(gl));
		System.out.println("glPanel " + debugToString(glPanel));
		System.out.println("glCanvas " + debugToString(glPanel.glc));
		int[] depthTest = new int[1];
		gl.glGetIntegerv(GL.GL_DEPTH_TEST, depthTest, 0);
		boolean isDepthTest = (depthTest[0] == GL.GL_TRUE);
		System.out.println("gl.depthTest " + isDepthTest);
	}

	private String debugToString(Object object)
	{
		return object.getClass().getName() + "@" + Integer.toHexString(object.hashCode());
	}

	protected void setWindowToRGB() { isRGB = true; }

	protected void setWindowToColorMap() { isRGB = false; }

	/**
	 * Set flag indicating we want to draw in 3d overlay mode.
	 *
	 * @return value indicates status is changed
	 */
	public boolean set3dOverlay(boolean set3dOverlay)
	{
		// don't allow overlay mode for Windows (avoid bug in Compaq laptop)
		if(set3dOverlay && Main.OS.startsWith("Windows")) return false;

		if(set3dOverlay) return false; // jbw doesn't work for dragging seismic on Linux either


		if(set3dOverlay != is3dOverlay)
		{
			this.is3dOverlay = set3dOverlay;
			return true;
		}
		else
			return false;
	}

	/** Get the center of the view point @return XYZ point */
	public StsPoint getCenterViewPoint()
	{ return centerViewPoint; }

	protected void copyColorBuffer(int colorBuffer)
	{
		/* Set rectcopy source and destination	*/
		/* turn zbuffer off during rectcopy		*/

		if(colorBuffer == GL.GL_FRONT)
		{
			//			System.out.println("Copying front to back.");

			gl.glReadBuffer(GL.GL_FRONT);
			gl.glDrawBuffer(GL.GL_BACK);
		}
		else
		{
			//			System.out.println("Copying back to front.");

			gl.glReadBuffer(GL.GL_BACK);
			gl.glDrawBuffer(GL.GL_FRONT);
		}

		gl.glDisable(GL.GL_DEPTH_TEST);

		/* GL rectangle is lower-left origin, Java window is upper-left origin	*/

		int mm[] = new int[1];

		gl.glGetIntegerv(GL.GL_MATRIX_MODE, mm, 0);
		gl.glMatrixMode(GL.GL_PROJECTION);
		gl.glPushMatrix();
		gl.glLoadIdentity();
		gl.glOrtho(0.0, (double) glPanel3d.getWidth(), 0.0, (double) glPanel3d.getHeight(), -1.0, 1.0);
		gl.glMatrixMode(GL.GL_MODELVIEW);
		gl.glPushMatrix();
		gl.glLoadIdentity();

		gl.glRasterPos2i(0, 0);
		gl.glCopyPixels(0, 0, glPanel3d.getWidth(), glPanel3d.getHeight(), GL.GL_COLOR);

		gl.glPopMatrix();
		gl.glMatrixMode(GL.GL_PROJECTION);
		gl.glPopMatrix();
		gl.glMatrixMode(mm[0]);
		/*
				  int gly = winRectGL.y - winRectGL.height;
				  gl.glRasterPos2i(winRectGL.x, gly);
				  gl.glCopyPixels(winRectGL.x, gly, winRectGL.width, winRectGL.height, GL.GL_COLOR);
		  */
		gl.glEnable(GL.GL_DEPTH_TEST);
	}

	/**
	 * Draws in foreground XORed against current view.  Must be called
	 * again and draw same objects to erase.
	 */
	public void drawForeground()
	{
		gl.glDisable(GL.GL_BLEND);
		gl.glDrawBuffer(GL.GL_FRONT);
		gl.glEnable(GL.GL_COLOR_LOGIC_OP);
		gl.glLogicOp(GL.GL_XOR);
		gl.glDepthMask(false);
		gl.glDepthFunc(GL.GL_ALWAYS);
		drawForegroundCursor();
		gl.glFlush();
		if(cursorButtonState == StsMouse.RELEASED)
		{
			gl.glLogicOp(GL.GL_COPY);
			gl.glDisable(GL.GL_COLOR_LOGIC_OP);
			gl.glDepthMask(true);
			gl.glDepthFunc(GL.GL_LESS);
			gl.glDrawBuffer(GL.GL_BACK);
			cursorButtonState = StsMouse.CLEARED;
		}
	}

	/**
	 * Draw the cursor in the front buffer.  If previously drawn,
	 * draw previous one again to erase it; then drawn new one.
	 */
	private void drawForegroundCursor()
	{
		if(previousXYZ != null)
		{
			//           System.out.println("drawing previous cursor at x: " + previousXYZ[0] + " y: " + previousXYZ[1]);
			drawCursorPoint(previousXYZ, this.glPanel3d);
		}
		if(currentXYZ != null)
		{
			//           System.out.println("drawing current cursor at x: " + currentXYZ[0] + " y: " + currentXYZ[1]);
			drawCursorPoint(currentXYZ, glPanel3d);
		}
		previousXYZ = currentXYZ;
		currentXYZ = null;
	}

	/**
	 * Draw vertical and horizontal lines thru the cursor point.
	 * because the crossed area is drawn twice it will be clear.
	 */
	private void drawCursorPoint(float[] xyz, StsGLPanel3d glPanel3d)
	{
		StsGLDraw.drawLine(xyz, StsColor.GRAY, glPanel3d, 3, 31, StsParameters.NORTH, 0.0, 0.0, -15.0);
		StsGLDraw.drawLine(xyz, StsColor.GRAY, glPanel3d, 3, 31, StsParameters.EAST, 0.0, -15.0, 0.0);
	}

	private void drawViewCenter()
	{
		float vec[] = new float[3];
		StsProject project = model.getProject();
		// draw line from grid to centerViewPoint with cyan
		project.getStsCogColor().setGLColor(gl);
		//StsColor.CYAN.setGLColor(gl);

		short DOTTED_LINE = (short) 0xCCCC;
		gl.glLineStipple(1, DOTTED_LINE);
		gl.glEnable(GL.GL_LINE_STIPPLE);

		//		StsPoint centerViewPoint = glPanel3d.getCenterViewPoint();
		vec[0] = centerViewPoint.v[0];
		vec[1] = centerViewPoint.v[1];
		/*
					 float[] unrotatedXY = project.getUnrotatedRelativeXYFromRotatedXY(centerViewPoint.v[0], centerViewPoint.v[1]);
					 vec[0] = unrotatedXY[0];
					 vec[1] = unrotatedXY[1];
				  */
		float gridZ = project.getGridZ();
		gl.glLineWidth(4.0f);
		gl.glDisable(GL.GL_LIGHTING);
		gl.glBegin(GL.GL_LINE_STRIP);
		{
			vec[2] = gridZ;
			gl.glVertex3fv(vec, 0);
			vec[2] = centerViewPoint.v[2];
			gl.glVertex3fv(vec, 0);
		}
		gl.glEnd();

		gl.glDisable(GL.GL_LINE_STIPPLE);

		//		StsGLDraw.drawPoint(vec, StsColor.BLACK, glPanel3d, 6, -1.0);
		StsGLDraw.drawPoint(vec, project.getStsCogColor(), glPanel3d, 4, 0.0);
		gl.glEnable(GL.GL_LIGHTING);
	}

	/** Set the view as changed. This will force a redraw. */
	public void viewChanged()
	{
		if(glPanel3d == null || glPanel3d.window == null) return;
		StsMouseActionToolbar toolbar = glPanel3d.window.getMouseActionToolbar();
		if(toolbar != null) toolbar.viewChanged(true);
        if(compass != null) compass.computeAxisPoints();
		glPanel3d.viewChanged();
	}

	public void adjustCursor(int dir, float dirCoor)
	{
		if(getCursor3d() == null || getCursor3d().cursorSections == null) return;
		getCursor3d().cursorSections[dir].dirCoordinate = dirCoor; //this seems to be necessary for StsPreStackLineSet2D.moveToProfile() to work (ensures current x/y stored correctly) SWC 8/27/09
	}

	/**
	 * Clear and redraw the entire 3D view contents
	 *
	 * @throws StsException
	 */
	protected void drawCompleteWindow() throws StsException
	{
		clearToBackground(GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);
		drawViewCenter();
		if(model.getProject().getDisplayCompass())
            compass.draw();
		if(model.displayOK)
		{
			model.display(glPanel3d);
			glPanel3d.getCursor3d().display3d(glPanel3d, this, gl, glu);
			checkCreateCollaborationViewTransaction();
		}
	}

	private void checkCreateCollaborationViewTransaction()
	{

		StsCollaboration collaboration = StsCollaboration.getCollaboration();
		if(collaboration != null && collaboration.hasPeers())
		{
			float[] viewParameters = getCenterAndViewParameters();
			StsDBMethodCmd cmd = new StsDBMethodCmd(this, "changeModelView3d", new Object[]{viewParameters});
			model.addTransientTransactionCmd("changeView", cmd);
		}
	}

	/**
	 * Move the other two cursor planes to the point selected on this cursor plane
	 *
	 * @param mouse mouse object
	 */
	public boolean moveCursor3d(StsMouse mouse, StsGLPanel3d glPanel3d)
	{
		// if (!mouse.isButtonStateReleased(StsMouse.LEFT)) return false;
		StsCursor3d cursor3d = getCursor3d();
		StsCursorPoint cursorPoint = cursor3d.getNearestPointInCursorPlane(glPanel3d, mouse);
		if(cursorPoint == null) return false;

		int dirNo = cursor3d.getCurrentDirNo();
		cursor3d.moveCursor3d(glPanel3d, cursorPoint);
		cursor3d.setCurrentDirNo(dirNo);
		return true;
	}

	private StsCursor3d getCursor3d()
	{
		if(glPanel3d == null) return null;
		return glPanel3d.getCursor3d();
	}

	/** Draw the foreground display objects */
	protected void drawForegroundOverlay()
	{
		glPanel3d.getCursor3d().display3d(glPanel3d, this, gl, glu);
	}

	/** Set or reset to the default view. This is the persepctive that is seen when a project is first started. */
	public void setDefaultView(boolean doLock)
	{
		int i;

		StsProject project = model.getProject();
		centerViewPoint.v[0] = project.getXCenter();
		centerViewPoint.v[1] = project.getYCenter();
		centerViewPoint.v[2] = (project.getZorTMin() + project.getZorTMax()) / 2.0f;

		float projectSize = project.getMaxProjectDimension();
		minViewDistance = 0.0001f * projectSize;
		clipFarOverNear = 200.0; // jbw from 50

		setViewElevation(45.0f);
		float angle = project.getAngle();
		setViewAzimuth(-angle);
		setViewDistance(projectSize);

		perspectiveFovy = Main.fovy; //jbw
		//        float width = glPanel3d.getWidth();
		//        float height = glPanel3d.getHeight();
		//        perspectiveAspect = width / height;

		//        degPerXPixel = (float) (perspectiveFovy / width);
		//       degPerYPixel = degPerXPixel;

		setZScale(1.0f);
		viewDistanceMult = 0.002f * distance;

		glPanel3d.viewChanged = true;
		glPanel3d.viewDistanceChanged = true;
		set3dOverlay(false);
		repaint();
		if(doLock)
			moveLockedWindows();
	}

	public void adjustView()
	{
		//        float projectSize = model.project.getMaxProjectDimension();
		//		minViewDistance = 0.0001f * projectSize;
		//        setViewDistance(projectSize);
		setZScale(1.0f);
		//        viewDistanceMult = 0.002f * distance;
		glPanel3d.viewChanged = true;
		glPanel3d.viewDistanceChanged = true;
		glPanel3d.repaint();
	}

	public void setDefaultView()
	{
		setDefaultView(true);
	}

	/** Reset the view to be looking straight down on the data from above. */
	public void setTopView()
	{
		setViewElevation(0.0f);
		setViewAzimuth(0.0f);
		setViewDistance(1.0f * model.getProject().getMaxProjectDimension());
		moveLockedWindows();
	}

	/**
	 * Get the current Z scale
	 *
	 * @return scale factor in Z
	 */
	public float getViewZScale()
	{
		return this.getZScale();
	}

	/**
	 * Get the current view elevation
	 *
	 * @return elevation above Center
	 */
	public float getViewElevation()
	{
		return this.elevation;
	}

	/**
	 * Get the view azimuth relative to North
	 *
	 * @return azimuth
	 */
	public float getViewAzimuth()
	{
		return azimuth;
	}

	/**
	 * Get the view distance. The view distance is the distance in model units from the Center of the Model
	 *
	 * @return distance from center
	 */
	public float getViewDistance()
	{
		return distance;
	}

	/**
	 * Get the current model Center
	 *
	 * @return center as window coordinates xyz triplet
	 */
	private double[] getViewCenter()
	{
		StsProject project = model.getProject();
		float[] unrotatedXY = project.getUnrotatedRelativeXYFromRotatedXY(centerViewPoint.v[0], centerViewPoint.v[1]);
		return glPanel3d.gluProject((double) unrotatedXY[0], (double) unrotatedXY[1], (double) centerViewPoint.v[2],
				getUnrotatedModelViewMatrix(), glPanel3d.projectionMatrix, glPanel3d.viewPort);
	}

	/**
	 * Pans, Zooms and Rotates the model view.
	 *
	 * @param mouse the mouse object
	 */
	public void moveWindow(StsMouse mouse)
	{
		StsMousePoint mouseDelta = mouse.getMouseDelta();
		if(mouseDelta.x == 0 && mouseDelta.y == 0) return;

		if((mouseMode == StsCursor.ZOOM && getKeyCode() == KeyEvent.VK_UNDEFINED) || getKeyCode() == KeyEvent.VK_Z)
		{
			setMotionCursor(StsCursor.ZOOM);
			setViewAzimuth(-degPerXPixel * mouseDelta.x + azimuth);
			// jbw This formulation limits zoom to near the centerpoint,
			// instead of being able to zoom all the way to the back of the model.
			//setViewDistance(((float)mouseDelta.y/glPanel3d.winRectGL.height + 1.0f)*(distance));
			setViewDistance((distance + ((float) mouseDelta.y / glPanel3d.getHeight()) * (model.getProject().getMaxProjectDimension())));
			//			System.out.println("moveWindow mouseDelta.y " + mouseDelta.y);
			//			zoomViewCenter((float)mouseDelta.y/glPanel3d.winRectGL.height + 1.0f);
		}
		else if((mouseMode == StsCursor.PAN && getKeyCode() == KeyEvent.VK_UNDEFINED) || getKeyCode() == KeyEvent.VK_X)
		{
			setMotionCursor(StsCursor.PAN);
			double[] windowCenter = getViewCenter();
			windowCenter[0] -= (double) mouseDelta.x;
			windowCenter[1] -= (double) mouseDelta.y;
			setViewCenter(windowCenter);
		}
		else if((mouseMode == StsCursor.ROTATE && getKeyCode() == KeyEvent.VK_UNDEFINED) || getKeyCode() == KeyEvent.VK_C) // rotation
		{
			setMotionCursor(StsCursor.ROTATE);
			setViewAzimuth(-degPerXPixel * mouseDelta.x + azimuth);
			setViewElevation(degPerYPixel * mouseDelta.y + elevation);
		}
		else // wheel
		{
			setMotionCursor(StsCursor.ZOOM);
			setViewAzimuth(-degPerXPixel * mouseDelta.x + azimuth);
			// jbw This formulation limits zoom to near the centerpoint,
			// instead of being able to zoom all the way to the back of the model.
			//setViewDistance(((float)mouseDelta.y/glPanel3d.winRectGL.height + 1.0f)*(distance));
			setViewDistance((distance + ((float) mouseDelta.y / glPanel3d.getHeight()) * (model.getProject().getMaxProjectDimension())));
			//			System.out.println("moveWindow mouseDelta.y " + mouseDelta.y);
			//			zoomViewCenter((float)mouseDelta.y/glPanel3d.winRectGL.height + 1.0f);
		}
		glPanel3d.fireViewChangeEvent(); // notifies 3d sliders of changes
		if(glPanel3d.mouseInfo)
			StsMessageFiles.infoMessage("Azimuth= " + azimuth + " Elevation= " + elevation + " Distance= " + distance + " Z Scale= " + getZScale());
		viewChangedRepaint();
		moveLockedWindows();
	}

	/**
	 * Shift the view towards (+) or away (-) from the viewer
	 * Note that the absolute amount of shift can be scaled with viewShiftFactor
	 * which is adjusted with the "O" key (increase by 2.0) and the "o" key
	 * (decrease by 2.0).
	 */
	public void zoomViewCenter(float factor)
	{
		StsPoint eye = new StsPoint(eyePoint);
		float eyeCenterDistance = eye.distance(centerViewPoint);
		StsPoint eyeCenterVector = new StsPoint(3);
		eyeCenterVector.subPoints(centerViewPoint, eye);
		StsPoint newCenterViewPoint = StsPoint.multByConstantAddPointStatic(eyeCenterVector, factor, eye);
		float distance = centerViewPoint.distance(eye);
		//System.out.println("old distance " + this.distance + " new distance " + distance);
		//System.out.println("old center " + centerViewPoint.toString() + " new center " + newCenterViewPoint.toString());
		centerViewPoint = newCenterViewPoint;

		setViewDistance(distance);
	}

	public boolean moveWithView(StsView movedView)
	{
		if(!(movedView instanceof StsView3d)) return false;
		StsView3d moved3dView = (StsView3d) movedView;
		setViewElevation(moved3dView.elevation);
		setViewAzimuth(moved3dView.azimuth);
		setViewDistance(moved3dView.distance);
		setViewCenter(moved3dView.centerViewPoint);
		setVerticalStretch(moved3dView.getZScale());
		repaint();
		return true;
	}

	private void setViewCenter(float x, float y, float z)
	{
		setViewCenter(new double[]{x, y, z});
	}

	private void setViewCenter(double[] center)
	{
		//        viewChanged = true;

		double[] viewCenter = glPanel3d.gluUnProject((double) center[0], (double) center[1], (double) center[2],
				getUnrotatedModelViewMatrix(), glPanel3d.projectionMatrix, glPanel3d.viewPort);
		StsProject project = model.getProject();
		float[] rotatedXY = project.getRotatedRelativeXYFromUnrotatedRelativeXY((float) viewCenter[0], (float) viewCenter[1]);
		viewCenter[0] = rotatedXY[0];
		viewCenter[1] = rotatedXY[1];
		if(debug)
			System.out.println("centerViewPoint changed from " + this.centerViewPoint.toString() + " to " + StsMath.toString(viewCenter));
		centerViewPoint.setValues(viewCenter);
		glPanel3d.viewChanged = true;
	}

	protected double[] getRotatedPointFromUnrotatedScreenPoint(double[] screenPoint)
	{
		return getRotatedPointFromUnrotatedScreenPoint(screenPoint[0], screenPoint[1], screenPoint[2]);
	}

	public double[] getRotatedPointFromUnrotatedScreenPoint(double screenPoint0, double screenPoint1, double screenPoint2)
	{
		double[] point = glPanel3d.gluUnProject(screenPoint0, screenPoint1, screenPoint2,
				getUnrotatedModelViewMatrix(), glPanel3d.projectionMatrix, glPanel3d.viewPort);
		StsProject project = model.getProject();
		float[] rotatedXY = project.getRotatedRelativeXYFromUnrotatedRelativeXY((float) point[0], (float) point[1]);
		point[0] = rotatedXY[0];
		point[1] = rotatedXY[1];
		return point;
	}

	public double[] getUnrotatedPointFromUnrotatedScreenPoint(double screenPoint0, double screenPoint1, double screenPoint2)
	{
		return glPanel3d.gluUnProject(screenPoint0, screenPoint1, screenPoint2,
				getUnrotatedModelViewMatrix(), glPanel3d.projectionMatrix, glPanel3d.viewPort);
	}

	public void setViewCenter(StsPoint centerViewPoint)
	{
		if(centerViewPoint.equals(this.centerViewPoint)) return;
		if(debug)
			System.out.println("centerViewPoint changed from " + this.centerViewPoint.toString() + " to " + centerViewPoint.toString());
		this.centerViewPoint.copyFrom(centerViewPoint);
		glPanel3d.viewChanged = true;
	}

	/*
			 public void setViewCenter(double[] viewParams)
			 {
				 if(centerViewPoint.equals(viewParams)) return;
				 if(centerViewPoint.equals(this.centerViewPoint)) return;
				 this.centerViewPoint.copyFrom(centerViewPoint);
				 glPanel3d.viewChanged = true;
			 }
		  */

	/**
	 * Get the center ofthe model and current view parameters
	 *
	 * @return 7 floats (x,y,z,distance,azimuth,elevation,zScale)
	 */
	public float[] getCenterAndViewParameters()
	{
		return new float[]{centerViewPoint.v[0], centerViewPoint.v[1], centerViewPoint.v[2],
				distance, azimuth, elevation, getZScale()};
	}

	/**
	 * Get the center ofthe model and current view parameters
	 *
	 * @return 7 floats (x,y,z,distance,azimuth,elevation,zScale)
	 */
	public void setCenterAndViewParameters(float[] parameters)
	{
		setViewCenter(new StsPoint(parameters[0], parameters[1], parameters[2]));
		setViewDistance(parameters[3]);
		setViewAzimuth(parameters[4]);
		setViewElevation(parameters[5]);
		setZScale(parameters[6]);
	}

	/**
	 * Set the model view parameters
	 *
	 * @param parameters 7 floats (x,y,z,distance,azimuth,elevation,zScale)
	 * @return true if successful
	 */
	public boolean changeModelView3d(float[] parameters)
	{
		centerViewPoint.setValues(parameters);
		setViewDistance(parameters[3]);
		setViewAzimuth(parameters[4]);
		setViewElevation(parameters[5]);
		setZScale(parameters[6]);
		glPanel3d.viewChanged = true;
		return true;
	}

	/**
	 * Display the popup menu)
	 *
	 * @param mouse
	 */
	public void showPopupMenu(StsMouse mouse)
	{
		StsViewSelectable object = getSelectableObject(mouse);
		if(object != null)
		{
			object.showPopupMenu(glPanel3d, mouse);
			return;
		}
		// check if well was picked
		popupMenu = new JPopupMenu("3D View Popup");
        JLabel title = new JLabel("Project Popup Menu", JLabel.CENTER);
        title.setFont(new java.awt.Font("Dialog", 1, 14));
        popupMenu.add(title);
        popupMenu.addSeparator();
        
        JMenu groupListMenu = new JMenu("Group Lists");

        StsArrayList viewClasses = model.getDisplayableStsClassesWithObjects();
        for(int i=0; i<viewClasses.size(); i++)
        {
            StsClass viewClass = (StsClass)viewClasses.get(i);
            StsMenuItem button = new StsMenuItem();
            //button.setMenuActionListener(viewClass.getName().substring(viewClass.getName().indexOf("DBTypes.Sts")+11), this, "showList", viewClass.getName());
            button.setMenuActionListener(viewClass.getUserName(), this, "showList", viewClass.getName());
            groupListMenu.add(button);
        }

        if(Main.viewerOnly)
            popupMenu.add(propertyWindow);

		JMenu viewMenu = new JMenu("Views");
		viewMenu.add(topDownViewBtn);
		viewMenu.add(normalViewBtn);
		viewMenu.add(northViewBtn);
        viewMenu.add(eastViewBtn);

		JMenu weighpointMenu = new JMenu("WayPoints");
		weighpointMenu.add(addWeighPointBtn);
		// Add all weighpoints to selection menu
		StsObject[] wps = model.getObjectList(StsWeighPoint.class);
		StsMenuItem[] weighPointBtns = new StsMenuItem[wps.length];
		for(int i = 0; i < wps.length; i++)
		{
			StsWeighPoint wp = (StsWeighPoint) wps[i];
			weighPointBtns[i] = new StsMenuItem();
			weighPointBtns[i].setMenuActionListener(wp.getName(), wp, "treeObjectSelected", null);
			weighpointMenu.add(weighPointBtns[i]);
		}
        popupMenu.add(groupListMenu);        
		popupMenu.add(weighpointMenu);
		popupMenu.add(viewMenu);

		glPanel3d.add(popupMenu);
		popupMenu.show(glPanel3d, mouse.getX(), mouse.getY());
		clearCurrentKeyPressed();
	}

    public void showList(String className)
    {
        StsClass c = model.getStsClass(className);
        if(!c.hasObjects())
            new StsMessage(model.win3d, StsMessage.WARNING, "No Objects of this type in project.", true);
        else        
            c.popupClassListPanel();
    }
 
	/* Cancel the popup menu
		   *
		   * @param mouse
		  */
	public void cancelPopupMenu(StsMouse mouse)
	{
		if(popupMenu != null) popupMenu.setVisible(false);
	}

	/*
			  public void editMultiViewPanels(StsWin3dBase win3d)
			  {
				  super.editMultiViewPanels(win3d);
			  }
		  */
	public void prepackagedViews(StsMenuItem btn)
	{
		if(btn == topDownViewBtn)
		{
			setViewElevation(0.0f);
		}
		else if(btn == normalViewBtn)
		{
			setViewElevation(90.0f);
		}
		else if(btn == northViewBtn)
		{
			setViewAzimuth(360.0f - model.getProject().getAngle());
		}
        else if(btn == eastViewBtn)
		{
			setViewAzimuth(270.0f - model.getProject().getAngle());
		}
	}

	/* Not used for now since measuring is done by an action
			   public StsObject defaultMousePressed(StsMouse mouse)
			   {
				   int sensorIdx = 0;
				   StsMouseActionToolbar mouseToolbar = getWindow().getMouseActionToolbar();
				   if(mouseToolbar != null) mouseToolbar.clearHighlights();

				   StsSensor sensor = (StsSensor) StsJOGLPick.pickVisibleClass3d(glPanel3d, StsSensor.class, StsJOGLPick.PICKSIZE_LARGE,
						   StsJOGLPick.PICK_FIRST);
				   if (sensor != null)
				   {
					   StsPickItem items = StsJOGLPick.pickItems[0];
					   if(items.names.length < 2) return sensor;
					   sensorIdx = items.names[1];
					   if(mouseToolbar == null) return sensor;
					   if(mouseToolbar.isMeasuring())
					   {
						  float[] pt = sensor.getXYZ(sensorIdx);
						  if(pt != null)
						  {
							  measureFrom = new StsPoint(pt);
							  //System.out.println("Measure from:" + sensor.getName() + " index=" + sensorIdx);
							  sensor.highlight(sensorIdx);
							  model.viewObjectRepaint(sensor);
						  }
						  return sensor;
					   }
				  }
				   return null;
			   }
		  */
	public void defaultMouseMotion(StsMouse mouse)
	{
		//System.out.println("Motion");
	}

	public StsViewSelectable mouseSelectedEdit(StsMouse mouse)
	{
		StsViewSelectable selectedObject = getSelectableObject(mouse);
		if(selectedObject != null)
			selectedObject.mouseSelectedEdit(mouse);
		return selectedObject;
	}

	/**
	 * Key Pressed event handling.
	 * If Z key pressed, mouse pressed motion results in zooming.
	 * If X key pressed, mouse pressed motion results in panning
	 * If no key pressed, mouse pressed motion results in rotation.
	 *
	 * @param mouse mouse object
	 * @param e	 key event
	 * @return true if successful
	 */
	 public boolean keyPressed(StsMouse mouse, KeyEvent e)
      {
          super.keyPressed(mouse, e);
		  // this code is commented out for now because there are continuous keyPresssed events coming if the key is held down
		  // which is causing an endless stream of resetMouseToggle calls to reset the mouse state; if the right mouse is down
		  // it then immediately resets the state to the pressed key state
		  // Have not tested issues with the RECTZOOM whic currently seems to be dysfunctional for 2D views.     TJL  4/4/2011
	  /*
          if (getKeyCode() != KeyEvent.VK_R) // cancel any rectangle zoom
          {
              StsMouseActionToolbar toolbar = glPanel3d.window.getMouseActionToolbar();
              if (toolbar != null) toolbar.resetMouseToggle();
          }
      */
          if (mouse.isButtonDown(StsMouse.RIGHT))
          {
              if (getKeyCode() == KeyEvent.VK_Z) setMotionCursor(StsCursor.ZOOM);
              else if (getKeyCode() == KeyEvent.VK_X) setMotionCursor(StsCursor.PAN);
              else if (getKeyCode() == KeyEvent.VK_R) setMotionCursor(StsCursor.RECTZOOM);
              else setMotionCursor(StsCursor.ROTATE);
              return true;
          }
          super.keyPressed(mouse, e);
          return true;
      }

	/**
	 * Key Released event handling.
	 * Up key or A key pressed, right button mouse released increases vertical scale by 2.0.
	 * Down key or S key pressed, right button mouse released increases vertical scale by 0.5.
	 * Q key pressed, right mouse button released increase view shift factor by 2.0.
	 * W key pressed, right mouse button released increase view shift factor by 2.0.
	 * Left key changes to crossplot view
	 * Right key changes to cursor view
	 *
	 * @param mouse mouse object
	 * @param e	 key event
	 * @return true if successful
	 */
	public void keyReleased(StsMouse mouse, KeyEvent e)
	{
		releasedKeyCode = e.getKeyCode();
		super.keyReleased(mouse, e);

		if(mouse.isButtonDown(StsMouse.RIGHT))
		{
			if(releasedKeyCode == KeyEvent.VK_A || releasedKeyCode == KeyEvent.VK_UP) verticalStretch(mouse);
			else if(releasedKeyCode == KeyEvent.VK_S || releasedKeyCode == KeyEvent.VK_DOWN) verticalShrink(mouse);
			else if(releasedKeyCode == KeyEvent.VK_Q) glPanel3d.changeViewShiftFactor(2.0);
			else if(releasedKeyCode == KeyEvent.VK_W) glPanel3d.changeViewShiftFactor(0.5);
			//            else if(releasedKeyCode==KeyEvent.VK_C) compressColorscale();
			//            else if(releasedKeyCode==KeyEvent.VK_U) uncompressColorscale();
			if(setMotionCursor(StsCursor.ROTATE))
			{
				repaint();
				glPanel3d.restoreCursor();
			}
			return;
		}

		// if rightmouse is not down, and we are clearing either X or Z when these keys are released: return
		if(releasedKeyCode == KeyEvent.VK_Z || releasedKeyCode == KeyEvent.VK_X)
		{
			if(setMotionCursor(StsCursor.ROTATE))
			{
				repaint();
				glPanel3d.restoreCursor();
			}
		}

		else if(releasedKeyCode == KeyEvent.VK_W) addWeighPoint();
		else if(releasedKeyCode == KeyEvent.VK_L) glPanel3d.window.displayPreviousObject(StsSeismicVolume.class);
		else if(releasedKeyCode == KeyEvent.VK_R) glPanel3d.window.displayNextObject(StsSeismicVolume.class);
		else
		{
			StsWindowFamily family = null;
			//keyCode = e.getKeyCode();
			family = model.getWindowFamily(this.getWindow());
			switch(releasedKeyCode)
			{
				case KeyEvent.VK_1:
					family.setSelectedDirection(StsCursor3d.XDIR, true);
					break;
				case KeyEvent.VK_2:
					family.setSelectedDirection(StsCursor3d.YDIR, true);
					break;
				case KeyEvent.VK_3:
					family.setSelectedDirection(StsCursor3d.ZDIR, true);
					break;
				case KeyEvent.VK_LEFT:
				case KeyEvent.VK_UP:
					family.decrementCursor();
					break;
				case KeyEvent.VK_RIGHT:
				case KeyEvent.VK_DOWN:
					family.incrementCursor();
					break;
				case KeyEvent.VK_HOME:
					setDefaultView();
					break;
				default:
					break;
			}
		}
	}

	/** Run the microTimer for this view. This will output timing information to the information panel. */
	public void runTimer()
	{
		runTimer = true;
	}

	/**
	 *
	 */
	public void addWeighPoint()
	{
		wp = new StsWeighPoint(model, model.getWindowFamily(this.getWindow()));
		//		clearCurrentKeyPressed();
	}

	/**
	 *
	 */
	public void printScreen()
	{
		Toolkit t = Toolkit.getDefaultToolkit();
		PrintJob pj = t.getPrintJob(model.win3d, "3D Screen", null);
		Graphics pg = pj.getGraphics();
		this.glPanel3d.printAll(pg);
		pg.dispose();
		pj.end();
	}

	/**
	 * Get the view parameters
	 *
	 * @return 4 floats (azimuth, elevation, distance, zScale)
	 */
	public float[] getViewParameters()
	{
		return new float[]{azimuth, elevation, distance, getZScale()};
	}

	/**
	 * Set the view parameters
	 *
	 * @param parameters 4 floats (azimuth, elevation, distance, zScale)
	 */
	public void setViewParameters(float[] parameters)
	{
		setViewAzimuth(parameters[0]);
		setViewElevation(parameters[1]);
		setViewDistance(parameters[2]);
		setZScale(parameters[3]);
	}

	/**
	 * Set the new view elevation as angle between 0 and 360 degrees
	 *
	 * @param elevation new elevation
	 */
	public void setViewElevation(float elevation)
	{
		if(elevation == this.elevation) return;
		if(elevation < -180.f)
			elevation += 360.f;
		else if(elevation > 180.f)
			elevation = (elevation + 180.f) % 360.f - 180f;
		if(debug) System.out.println("Elevation changed from " + this.elevation + " to " + elevation);
		this.elevation = elevation;
		if(glPanel3d != null)
			glPanel3d.viewChanged = true;
	}

	/**
	 * Set the new view azimuth
	 *
	 * @param azimuth new azimuth
	 */
	public void setViewAzimuth(float azimuth)
	{
		if(azimuth == this.azimuth) return;
		if(azimuth < 0)
			azimuth += 360.f;
		else if(azimuth > 360.f)
			azimuth = azimuth % 360.f;
		if(debug) System.out.println("Azimuth changed from " + this.azimuth + " to " + azimuth);
		this.azimuth = azimuth;
		glPanel3d.viewChanged = true;
	}

	/**
	 * Set view distance. The distance is the distance in model units from the center of the model.
	 *
	 * @param distance new distance
	 */
	public void setViewDistance(float distance)
	{
		distance = Math.max(distance, -1.0f * model.getProject().getMaxProjectDimension());
		if(this.distance == distance) return;
		if(debug) System.out.println("distance changed from " + this.distance + " to " + distance);
		this.distance = distance;
		if((model.getProject().getZoomLimit() != model.getProject().NOLIMIT) && (distance < model.getProject().getZoomLimit()))
		{
			StsMessageFiles.infoMessage("Exceeding zoom limit...reseting to limit.");
			distance = model.getProject().getZoomLimit();
		}
		// this clipping isn't really right. This formulation limits zoom to near the centerpoint,
		// instead of being able to zoom all the way to the back of the model.
		clipNear = 0.05f * distance;
		//clipFar = clipNear*clipFarOverNear;
		// jbw don't make clipFar a funky ratio, rather, set to back of model (worst-case)
		// the far distance doesn't matter nearly as much, and it's weird
		// that as one zooms in, the back of the model gets clipped.
		//
		// worst case is center of rotation is in foremost part of model, so all of the model
		// (2 * max dimension) is behind the center.
		//
		clipFar = 2.0f * model.getProject().getMaxProjectDimension() + this.distance;
		//System.out.println("first clipnear is "+clipNear);
		//
		// jbw maintain a maximum far/near limit to prevent Z errors
		// do so by clipping off the near, not the far
		//
		clipNear = Math.max(clipFar / clipFarOverNear, clipNear);
		//System.out.println("second clipnear is "+clipNear);
		glPanel3d.viewChanged = true;
		glPanel3d.viewDistanceChanged = true;
	}

	/** Increase the vertical stretch by 1.5 times */
	public void verticalStretch(StsMouse mouse)
	{
		setZScale(1.5f * getZScale()); // jbw from 2.0
		glPanel3d.changeModelView();
		moveLockedWindows();
	}

	/** Decrease the vertical stretch by 0.75 times */
	public void verticalShrink(StsMouse mouse)
	{
		setZScale(0.75f * getZScale()); // jbw from 0.5
		glPanel3d.changeModelView();
		moveLockedWindows();
	}

	public void setVerticalStretch(float zScale)
	{
		if(zScale == this.getZScale()) return;
		this.setZScale(zScale);
		glPanel3d.viewChanged = true;
	}

	/**
	 * Make the current view look exactly like the other 3d view
	 *
	 * @param otherView the other 3d view
	 */
	public void copy(StsView3d otherView)
	{
		// copy((StsView)otherView);
		solidMaterial = otherView.solidMaterial;
		transparentMaterial = otherView.transparentMaterial;
		isSolid = otherView.isSolid;
		isRGB = otherView.isRGB;
		is3dOverlay = otherView.is3dOverlay;
		drawing3dOverlay = otherView.drawing3dOverlay;
		lightedSurfaces = otherView.lightedSurfaces;
		centerViewPoint = new StsPoint(otherView.centerViewPoint);
		distance = otherView.distance;
		azimuth = otherView.azimuth;
		elevation = otherView.elevation;
		minViewDistance = otherView.minViewDistance;
		clipNear = otherView.clipNear;
		clipFar = otherView.clipFar;
		clipFarOverNear = otherView.clipFarOverNear;
		perspectiveFovy = otherView.perspectiveFovy;
		//        perspectiveAspect = otherView.perspectiveAspect;
		degPerXPixel = otherView.degPerXPixel;
		degPerYPixel = otherView.degPerYPixel;
		setZScale(otherView.getZScale());
		viewDistanceMult = otherView.viewDistanceMult;
		glPanel3d.viewChanged = true;
		glPanel3d.viewDistanceChanged = true;
		for(int i = 0; i < 3; i++)
			glPanel3d.getCursor3d().setSelectedDirection(i, otherView.glPanel3d.getCursor3d().isDisplayEnable(i));
	}

	/**
	 * Output the mouse tracking on cursor readout to the information panel on the main screen
	 * Display cursor on window and other family windows if mouseTracking is on.
	 *
	 * @param glPanel3d
	 * @param mouse	 mouse object
	 */
	public StsCursorPoint logReadout(StsGLPanel3d glPanel3d, StsMouse mouse)
	{
		StsCursorPoint cursorPoint = this.glPanel3d.getCursor3d().getNearestPointInCursorPlane(glPanel3d, mouse);
		if(cursorPoint == null)
			setCursorXOR(glPanel3d, mouse, null);
		else
			setCursorXOR(glPanel3d, mouse, cursorPoint.point.v);
//          if(!moveCursor3d(mouse, glPanel3d)) return; // Still want to reposition cursor even if mouse tracking
		glPanel3d.getCursor3d().logReadout(glPanel3d, cursorPoint);
		return cursorPoint;
	}

	/**
	 * Output the mouse tracking on cursor readout to the information panel on the main screen
	 * Display cursor on window and other family windows if mouseTracking is on.
	 *
	 * @param mouse mouse object
	 */
	public void surfaceReadout(StsSurface surface, StsMouse mouse)
	{
		boolean display = !mouse.isButtonDown(StsMouse.LEFT);
		StsGridPoint gridPoint = surface.getSurfacePosition(mouse, display, glPanel3d);
		if(gridPoint == null) return;
		setCursorXOR(glPanel3d, mouse, gridPoint.getXYZorT());
	}


    public void sensorReadout(StsMouse mouse)
	{
		  int sensorIdx = 0;
		  StsDynamicSensor sensor = (StsDynamicSensor) StsJOGLPick.pickVisibleClass3d(glPanel3d, StsDynamicSensor.class, StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_FIRST);
		  if (sensor != null)
		  {
			  StsPickItem items = StsJOGLPick.pickItems[0];
              if (items.names.length < 2) return;
              sensorIdx = items.names[1];
              sensor.setPicked(sensorIdx);
			  if(mouse.getCurrentButton() == StsMouse.MIDDLE)
				  sensor.showPopupMenu(glPanel3d, mouse);
			  sensor.logMessage(sensor.toString(sensorIdx));
			  sensor.getDynamicSensorClass().displayGather(sensor, sensorIdx, model.win3d);
		  }
	}

	public void removeDisplayableClass(StsObject object)
	{
		if(glPanel3d == null || glPanel3d.getCursor3d() == null) return;
		glPanel3d.removeDisplayableClass(object);
	}

	public void setDefaultMouseMode()
	{
		setMouseMode(StsCursor.ROTATE);
		//        StsMouseActionToolbar toolbar = glPanel3d.window.getMouseActionToolbar();
		//        if (toolbar != null) toolbar.rotate();
	}

	/*
		  * custom serialization requires versioning to prevent old persisted files from barfing.
		  * if you add/change fields, you need to bump the serialVersionUID and fix the
		  * reader to handle both old & new
		  */
	static final long serialVersionUID = 1l;

	public boolean initializeDefaultAction()
	{
		StsObject object = model.getCurrentObject(StsPreStackLineSet.class);
		if(object instanceof StsPreStackLineSet2d)
			setDefaultActionPreStack2d(object);
		else
			setDefaultAction3d();
		return true;
		/*
					  else if(object instanceof StsPreStackLineSet3d)
					  {
						  setDefaultAction3d();
						  return true;
					  }
					  else
						  return false;
				  */
	}

	public boolean setViewObject(Object object)
	{
		return setDefaultAction(object);
	}

	public boolean setDefaultAction(Object object)
	{
		if(object == null) return false;
		if(object instanceof StsPreStackLineSet2d)
			return setDefaultActionPreStack2d(object);
		else
			return setDefaultAction3d();
	}

	private boolean setDefaultActionPreStack2d(Object object)
	{
		if(object == null) return false;

		if(defaultAction instanceof StsDefaultActionPreStack2d) return false;
		setDefaultAction(new StsDefaultActionPreStack2d(glPanel3d.actionManager, (StsPreStackLineSet2d) object));
		return true;
	}

	public boolean setDefaultAction3d()
	{
		StsActionManager actionManager = glPanel3d.window.actionManager;
		if(actionManager == null) return true;
		StsAction defaultAction = getDefaultAction();
		if(defaultAction instanceof StsDefaultAction3d) return false;
		setDefaultAction(new StsDefaultAction3d(glPanel3d, actionManager));
		return true;
	}

	//TODO user existing and new interfaces to make this more generic
	public boolean viewObjectChanged(Object source, Object object)
	{
		if(model.isStsClassCursor3dTextureDisplayable(object))
		{
			if(glPanel3d.window.getCursor3d().viewObjectChanged(source, object))
				return true;
		}
		else if(model.isStsClassTextureDisplayable(object))
		{
			return model.getStsClassTextureDisplayable(object).textureChanged((StsObject) object);
		}
		if(StsView3dDisplayable.class.isAssignableFrom(object.getClass()))
			return true;
		else if(object instanceof StsSubVolume)
		{
			glPanel3d.getCursor3d().subVolumeChanged();
			return true;
		}
		/*
				  else if (object instanceof StsSeismicLineSet)
				  {
					  ((StsSeismicLineSet) object).textureChanged();
					  return true;
				  }
				  else if (object instanceof StsPreStackLineSet2d)
				  {
					  ((StsPreStackLineSet2d) object).textureChanged();
					  setDefaultAction(object);
					  return true;
				  }
				  else if (object instanceof StsPreStackLine2d)
				  {
					  ((StsPreStackLine2d) object).textureChanged();
					  return true;
				  }
			  */
		else if(object instanceof StsPreStackVelocityModel)
		{
			glPanel3d.getCursor3d().clearTextureDisplays();
			return true;
		}
		else if(object instanceof StsSuperGather)
		{
			return true;
		}
		return false;
	}

	public boolean viewObjectRepaint(Object source, Object object)
	{
		glPanel3d.repaint();
		return true;
	}

	public byte getHorizontalAxisType() { return AXIS_TYPE_NONE; }

	public byte getVerticalAxisType() { return AXIS_TYPE_NONE; }

	/** vertical scale factor: negative to handle Z+ down */
	public float getZScale()
	{
		return zScale;
	}

	public void setZScale(float zScale)
	{
		if(this.zScale == zScale) return;
		System.out.println("zScale change from " + this.zScale + " to " + zScale);
		this.zScale = zScale;

	}


	class Compass
	{
		double size = 0.05;
		int direction;
		double[] screenOrigin;
		double[][] axisPoints = new double[4][];
		StsColor[] axisColors = new StsColor[] { StsColor.WHITE, StsColor.RED, StsColor.GREEN, StsColor.CYAN };
		String[] axisLabels = new String[] { "", "E", "N", "Z" };

		Compass(double[] screenOrigin)
		{
			this.screenOrigin = screenOrigin;
			computeAxisPoints();
		}

		public void computeAxisPoints()
		{
			int screenWidth = glPanel3d.getWidth();
			int screenHeight = glPanel3d.getHeight();
			double screenXOriginPixels = screenOrigin[0]*screenWidth;
			double screenYOriginPixels = screenOrigin[1]*screenHeight;
			axisPoints[0] = getUnrotatedPointFromUnrotatedScreenPoint(screenXOriginPixels, screenYOriginPixels, screenOrigin[2]);
			for(int n = 1; n <= 3; n++)
				axisPoints[n] = StsMath.copy(axisPoints[0]);

			double[] lengthPoint = getUnrotatedPointFromUnrotatedScreenPoint(screenXOriginPixels + size*screenWidth, screenYOriginPixels, screenOrigin[2]);
			double distance = StsMath.distance(axisPoints[0], lengthPoint, 3);
			axisPoints[1][0] += distance;
			axisPoints[2][1] += distance;
			axisPoints[3][2] += distance/zScale;
			StsProject project = model.getProject();
			float angle = project.getAngle();
			if(angle != 0.0f)
			{
				for(int n = 0; n < 4; n++)
					project.rotatePoint(axisPoints[n]);
			}
			// System.out.println("compass origin: " + axisPoints[0][0] + " " + axisPoints[0][1] + " " + axisPoints[0][2]);
		}

		public void draw()
		{
			computeAxisPoints();
			gl.glDisable(GL.GL_LIGHTING);
			StsGLDraw.enableLineAntiAliasing(gl);

			for(int n = 1; n <= 3; n++)
			{
				axisColors[n].setGLColor(gl);
				gl.glBegin(GL.GL_LINE_STRIP);
				{
					gl.glVertex3dv(axisPoints[0], 0);
					gl.glVertex3dv(axisPoints[n], 0);
				}
				gl.glEnd();
			}
			StsGLDraw.disableLineAntiAliasing(gl);
			GLBitmapFont horizontalFont = GLHelvetica18BitmapFont.getInstance(gl);
			for(int n = 1; n <= 3; n++)
			{
				axisColors[n].setGLColor(gl);
				StsGLDraw.fontOutput(gl, axisPoints[n], axisLabels[n], horizontalFont);
			}
			gl.glEnable(GL.GL_LIGHTING);
		}
	}    
}
