//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.MVC.View3d;

import com.Sts.DB.DBCommand.*;
import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.StsFile;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Sounds.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;

import javax.imageio.*;
import javax.imageio.stream.*;
import javax.media.opengl.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.io.*;
import java.util.*;

/**
 * StsGLPanel3d holds the GLComponent for 3d drawing and is inside
   * a panel of StsWin3d.
   * <p/>
   * Display is actually done by the currentView which is one of the views contained herein.
   * <p/>
   * glPanel3d has a cursor3d which consists of 3 orthogonal planes cutting thru the project.
   * {@link com.Sts.Interfaces.StsClassCursor3dTextureDisplayable} classes are displayed on these cursor planes.
   * <p/>
   * Notes:
   * <p/>
   * to redraw, ALWAYS call the method in the GLPanel super class: repaint().
   * This calls GLComponent.repaint(this) which handles GL context
   * locking/unlocking: VERY IMPORTANT.
   * <p/>
   * Lines are drawn with LIGHTING disabled, surfaces with LIGHTING enabled.
   * By default, lighting is on. When you draw lines, turn LIGHTING off
   * and then back on when finished. If you don't, lines are dull typically.
   */

  public class StsGLPanel3d extends StsGLPanel implements StsSerializable
  {
      //private StsView[] views = new StsView[0];
      /** indicates the center-position to eye-distance has changed requiring recomputing the projection matrix for a true 3d view */
      transient public boolean viewDistanceChanged = true;
      /** convenience copy of the window containing this panel */
      transient public StsWin3dBase window;
      /** limit the number of images captured to 2 frames per second
       */
      transient public long lastCapture = 0l;
      /**
       * List of picker operators.  When a pick operation is generated, window display is initiated on an invokeAndWait by the method
       * executing the pick; it can then act on the resulting picks. So when display is called in this class, instead of displaying,
       * pick operations will be executed in FIFO fashion.
       */
      transient public ArrayList pickQueue = new ArrayList(10);
      /** flag for pick queue */
      transient static boolean isPicking = false;
      /** image size for screen capture operations */
      // transient int mheight, mwidth;

      static final long serialVersionUID = 1l;

      //    transient boolean reshaped = false;
      public StsGLPanel3d()
      {
          // null constructor for serialization
      }

      /**
       * @param model         the model
       * @param actionManager the controller
       * @param window        window holding this glPanel
       */
  /*
      public StsGLPanel3d(StsModel model, StsWindowActionManager actionManager, Dimension size, StsWin3dBase window)
      {
          super(model, actionManager, size, null);
          window.setGlPanel3d(this);                                  s
          this.window = window;
      }
  */
      public StsGLPanel3d(StsModel model, StsActionManager actionManager, StsWin3dBase window)
      {
          super(model, actionManager);
          this.window = window;
      }

      public void initializeTransients(StsModel model, StsActionManager actionManager, StsWin3dBase window)
      {
          super.initializeTransients(model, actionManager);
          this.window = window;
          view.glPanel3d = this;
          view.glPanel = this;
      }
   //   public StsCursor3d getCursor3d() { return cursor3d; }

      /*
      private void addResizeListener()
      {
  //           Toolkit.getDefaultToolkit().setDynamicLayout(true);
         ComponentListener componentResizeListener = new ComponentAdapter()
         {
            public void componentResized(ComponentEvent e)
            {
                System.out.println("Window resized.");
  //                  boolean down = StsMouse.isLeftButtonDown();
  //                  System.out.println("Left button down: " + down);
            }
            public void componentMoved(ComponentEvent e)
            {
                System.out.println("Window moved.");
  //                  boolean down = StsMouse.isLeftButtonDown();
  //                  System.out.println("Left button down: " + down);
            }
         };
         addComponentListener(componentResizeListener);
      }
      */
      // Component initialization
      public void initializeTransients(StsModel model, StsWin3dBase window)
      {
          this.model = model;
          this.window = window;
          if(window == null) return;
          super.initializeTransients(model, window.getActionManager());
          if(size != null) setPanelSize(size);
          view.initializeTransients(this);
      }

      protected void initializeSize()
      {
          setMinimumSize(size00);
          Dimension size = getSize();
//          setPreferredSize(size);
          glc.setSize(size);
      }
  /*
      public boolean setView(String viewName, Class viewClass)
      {
          for(int n = 0; n < views.length; n++)
          {
              if(views[n].getViewName().equals(viewName))
              {
                  // required to reset from 2d with axis to 3d without.
                  currentView.getWindow().getMouseActionToolbar().reconfigure(this);
                  if(currentView == views[n])
                      return true;
                  else
                  {
                      viewChanged(views[n]);
                      repaint();
                      return true;
                  }
              }
          }
          return addView(viewClass);
      }
  */
    public StsView checkAddView(Class viewClass)
    {
        if(view != null && view.getClass().equals(viewClass)) return null;
        StsView view = StsView.constructSubclass(viewClass, this);
        if(view == null)
        {
            StsException.systemError(this, "checkAddView", "failed to construct view for " + StsToolkit.getSimpleClassname(view));
            return null;
        }
        setView(view);
        return view;
    }

    public void setView(StsView view)
    {
        super.setView(view);
        view.glPanel3d = this;
        view.initializeTransients(this);
    }

      /**
       * Add this view to the view list and make it the current view.
       * If a view exists which is the same class, though, ignore this
       * new view and set the existing view as current view.
       */
  /*
      public StsView setView(StsView view)
      {
          Class viewClass = view.getClass();
          for(int n = 0; n < views.length; n++)
          {
              if(views[n].getClass().equals(viewClass))
              {
                  StsException.systemError("StsGLPanel3d.setView(StsView) failed.  Already have a view instance: " + viewClass.getName());
                  if(currentView != views[n])
                  {
                      viewChanged(views[n]);
                      repaint();
                  }
                  return currentView;
              }
          }
          views = (StsView[])StsMath.arrayAddElement(views, view);
          viewChanged(view);
          repaint();
          return currentView;
      }
  */
      private void viewChanged(StsView view)
      {
          if(view == null) return;

          if(this.view == view) return;
          view.glPanel3d = this;

          view.initializeTransients(this);
          window.selectToolbarItem(StsViewSelectToolbar.NAME, view.getViewName(), true);

          // currentView.initGL();
          // notify the actionManager that the view has changed who will pass it to the active action(s)
          if(actionManager != null)
          {
              actionManager.viewChanged(view);
              view.initializeDefaultAction();
          }
          window.getMouseActionToolbar().reconfigure();
          window.selectToolbarItem(StsViewSelectToolbar.NAME, view.getViewName(), true);
      }

      public void viewChanged()
      {
          viewChanged = true;
          viewDistanceChanged = true;
      }

    /**
     * for an StsGLPanel3d, each view has its own mouseMode which it manages
     */
    /*
        public int getMouseMode()
        {
            if(currentView == null) return super.getMouseMode();
            else return currentView.getMouseMode();
        }
    */
    public void setMouseMode(int mouseMode)
    {
        view.setMouseMode(mouseMode);
    }

    /**

      /** A left mouse pressed asks for a default edit of an object under the mouse. */
      public StsObject defaultMousePressed(StsMouse mouse)
      {
          if(view == null) return null;
          StsObject stsObject = view.defaultMousePressed(mouse);
          return stsObject;
      }

      /** A left mouse motion asks for a default edit of an object under the mouse. */
      public void defaultMouseMotion(StsMouse mouse)
      {
          if(view == null) return;
          view.defaultMouseMotion(mouse);
          return;
      }

      public void moveWindow(StsMouse mouse) { view.moveWindow(mouse); }

      public void keyPressed(KeyEvent e, StsMouse mouse) { view.keyPressed(mouse, e); }

      public void keyReleased(KeyEvent e, StsMouse mouse) { view.keyReleased(mouse, e); }

      public boolean changeModelView3d(float[] parameters)
      {
          if(view == null) return false;
          return view.changeModelView3d(parameters);
      }

      public void computeProjectionMatrix()
      {
          if(view == null) return;
          //        viewPortChanged = true;
          view.computeProjectionMatrix();
      }

  //    public void setDisplay3dCursor(boolean display3dCursor) { cursor3d.setDisplay3dCursor(display3dCursor); }

  //    public boolean checkCursor3dIsOn() { return cursor3d.isDisplay3dCursor(); }
/*
      public void setBounds(int x, int y, int width, int height)
      {
          super.setBounds(x, y, width, height);
          if(actionManager == null) return;
          StsTransaction transaction = model.getCurrentTransaction();
          if(transaction != null)
          {
              transaction.add(new StsDBMethodCmd(model, "reshape3d", new Object[]{new int[]{x, y, width, height}}));
          }
      }
*/
      public void reshape3d(int[] parameters)
      {
          setBounds(parameters[0], parameters[1], parameters[2], parameters[3]);
      }

      /* Override GLComponent.setBounds( GL, GLU, ... ) */
      /*
          public void reshape( GLAutoDrawable drawable, int x, int y, int width, int height)
          {
              super.reshape(drawable,x,y,width,height);
              reshape(x, y, width, height );
              if(currentView == null) return;
      //		currentView.viewPortChanged();
          }
      */
      public void doPick(StsGLPicker picker)
      {
          pickQueue.add(picker);
          gld.display();
          //display();
      }
  /*
      public void reshape(GLAutoDrawable drawable, int x, int y, int width, int height)
      {
          setViewPort(x, y, width, height);

  //        StsException.systemDebug(this, "reshape", "width: " + width + " height: " + height);
          //if(Main.isDrawDebug) System.out.println("StsGLPanel.reshape() width: " + width + " height: " + height);
          if(getWidth() == width && getHeight() == height)
          {
              //           StsException.systemDebug(this, "reshape", "reshaped: false");
              //          if (model.displayOK) return;
              //           StsException.systemDebug(this, "reshape", "display enabled.");
              //           model.enableDisplay();
              //           model.win3dDisplayAll(window);
              return;
          }
  //        reshaped = true;
  //        StsException.systemDebug(this, "reshape", "reshaped: true");
          //       StsException.systemDebug(this, "reshape", "display disabled.");
          //       if (model.displayOK) model.disableDisplay();
          computeProjectionMatrix();
          StsTextureList.deleteAllTextures(gl);
          // deleteAllDisplayLists(gl);
      }
  */
      private String getDebugMessage()
      {
          String windowString = "null";
          if(window != null)
              windowString = window.getName();
          String viewString = " view: null";
          if(view != null)
              viewString = view.getViewName();
          return "window: " + windowString + " view: " + viewString;
      }

      /** display called by GLEvent loop. The current view display method is called. */
      public void display(GLAutoDrawable drawable)
      {
          /*
          if(viewPortChanged)
          {
              System.out.println("Resizing window: no display.");
              viewPortChanged = false;
              return;
          }
          */
          StsTextureList.checkDeleteTextures(gl);
          if(Main.isDrawDebug) StsException.systemDebug(this, "display", getDebugMessage() + "size: " + getWidth() + "," + getHeight());

          if(pickQueue.size() > 0 && !isPicking)
          {
              isPicking = true;
              while(pickQueue.size() > 0)
              {
                  try
                  {
                      ((StsGLPicker)pickQueue.get(0)).execute(this, drawable);
                  }
                  catch(Exception e)
                  {
                      StsException.outputWarningException(this, "display", "pickQueue.execute failed", e);
                  }
                  pickQueue.remove(0);
              }
              isPicking = false;
          }
          else
          {
              if(isPicking)
              {
                  ((StsGLPicker)pickQueue.get(0)).doPick(drawable);
              }
              else
              {
                  //				 System.out.println("StsGLPanel3d currentView.display().");
                  if(!view.getXORRectangle())
                      view.display(drawable);
                  else
                      view.drawXORRectangle();
              }
          }

          if(window.captureNextDraw() || window.captureMovie())
          {
              try
              {
                  int xorigin, yorigin;
                  int mwidth, mheight;
                  String typeString = null;
                  Toolkit toolkit = Toolkit.getDefaultToolkit();
                  Dimension screenSize = toolkit.getScreenSize();

                  if(window.getCaptureType() == window.DESKTOP)
                  {
                      // For full screen
                      xorigin = 0;
                      yorigin = 0;
                      mheight = screenSize.height;
                      mwidth = screenSize.width;
                      typeString = "Desktop";
                  }
                  else if(window.getCaptureType() == window.WINDOW)
                  {
                      // For window
                      xorigin = window.getX();
                      yorigin = window.getY();
                      mheight = window.getHeight();
                      mwidth = window.getWidth();
                      typeString = "Window";
                  }
                  else if(window.getCaptureType() == window.GRAPHIC)
                  {
                      // For glPanel
                      xorigin = window.getX() + window.currentViewPanel.getX() + 4;
                      yorigin = window.getY() + window.currentViewPanel.getY() + window.toolbarPanel.getHeight() - 10;
                      if(window instanceof StsWin3d)
                          yorigin = yorigin + 60;
                      mheight = window.currentViewPanel.getHeight();
                      mwidth = window.currentViewPanel.getWidth();
                      typeString = "Graphic";
                  }
                  else
                      return;

                  long time = System.currentTimeMillis();
                  if((time -lastCapture) > 250)
                  {
                    lastCapture = time;
                    String timeS = String.valueOf(lastCapture);
                    timeS = timeS.substring(0, timeS.length() - 2);
                    String filename = null;

                    Rectangle screenRect = new Rectangle(xorigin, yorigin, mwidth, mheight);
                    Robot robot = new Robot();
                    BufferedImage image = robot.createScreenCapture(screenRect);
                    boolean highResolution = model.getProject().getHighResolution();
                    Iterator iter = ImageIO.getImageWritersByFormatName("jpeg");
                    ImageWriter writer = (ImageWriter)iter.next();
                    ImageWriteParam iwp = writer.getDefaultWriteParam();
                    iwp.setCompressionMode(ImageWriteParam.MODE_EXPLICIT);
                    if(highResolution)
                        iwp.setCompressionQuality(1);

                    StsFile.checkCreateDirectory(model.getProject().getRootDirString() + model.getProject().getMediaDirString());
                    filename = model.getProject().getRootDirString() + model.getProject().getMediaDirString() + "S2S" + typeString + timeS + ".jpg";
                    File file = new File(filename);
                    FileImageOutputStream output = new FileImageOutputStream(file);
                    writer.setOutput(output);
                    IIOImage ioimage = new IIOImage(image, null, null);
                    writer.write(null, ioimage, iwp);

                    StsSound.play(StsSound.CAMERA_CLICK);
                    StsMessageFiles.logMessage("Captured image of " + typeString + " named:" + filename);
                    output.close();
                  }
              }
              catch(Exception e)
              {
                  StsException.outputException("Error capturing image, check disk space.", e, StsException.WARNING);
              }
              window.setCaptureNextDraw(false);
          }
      }

      // public int getFrameWidth() { return mwidth; }

      // public int getFrameHeight() { return mheight; }

      /**
       * A view shift is applied to an object when we wish it to be slghtly in front of
       * or slightly behind some other object, e.g., gridlines on a surface.  Without a
       * view shift, objects will occupy identical cells in the zbuffer and will be drawn
       * intermingled.
       * <p/>
       * This method toggles view shift on and off so we can see its effect.
       */
      public void toggleViewShift(ItemEvent e)
      {
          super.toggleViewShift(e);
      }

      public StsPoint getPointInPlaneAtMouse(StsMouse mouse)
      {
          double f, x, y, z;
          double[] point = new double[3];

          //        StsMousePoint mousePoint = mouse.getMousePointGL(winRectGL);
          //        double mx = mousePoint.x;
          //        double my = winRectGL.height - mousePoint.y - 1;

          int[] mouseViewPort = ((StsView2d)view).getInsetViewPort();


          double[] worldCoors = gluUnProject((double)getX(), (double)getY(), 0.0,
                  modelViewMatrix, projectionMatrix, mouseViewPort);

          return new StsPoint(worldCoors[0], worldCoors[1], worldCoors[2]);
      }

      public double[][] getViewLineAtMouse(StsMouse mouse)
      {
          double[] pointNear = new double[3];
          double[] pointFar = new double[3];

          StsMousePoint mousePoint = mouse.getMousePointGL(viewPort);
          //        double mx = mousePoint.x;
          //        double my = winRectGL.height - mousePoint.y - 1;

          double x = mousePoint.x;
          double y = mousePoint.y;
          pointNear = gluUnProject(x, y, 0.0, modelViewMatrix, projectionMatrix, viewPort);

          pointFar = gluUnProject(x, y, 1.0, modelViewMatrix, projectionMatrix, viewPort);

          return new double[][]{pointNear, pointFar};
      }

      public double[] getViewVectorAtMouse(StsMouse mouse)
      {
          double[][] line = getViewLineAtMouse(mouse);
          return StsMath.subtract(line[1], line[0]);
      }

      /** get near and far projection points */
      public double[] getNearProjectedPoint(StsMousePoint mousePoint)
      {
          return getProjectedPoint(mousePoint, 0.0);
      }

      public double[] getFarProjectedPoint(StsMousePoint mousePoint)
      {
          return getProjectedPoint(mousePoint, 1.0);
      }

      private double[] getProjectedPoint(StsMousePoint mousePoint, double mz)
      {
          //        StsMousePoint mousePointGL = mousePoint.lowerLeftCoor(winRectGL);

          return gluUnProject((double)getGLX(mousePoint.x), (double)getGLY(mousePoint.y), mz,
                  modelViewMatrix, projectionMatrix, viewPort);

          //        StsPoint pointF = new StsPoint((float)point[0], (float)point[1], (float)point[2]);
          //        return pointF;
      }

      public double[] getMouseTriangleIntersect(StsMousePoint mousePoint, StsPoint[] triPoints)
      {
          int n;
          StsPoint origin;
          StsPoint[] vectors = new StsPoint[2];
          StsPoint vectorM;

          //        StsMousePoint mousePointGL = mousePoint.lowerLeftCoor(winRectGL);

          double[][] screenCoor = new double[3][3];

          for(n = 0; n < 3; n++)
              screenCoor[n] = gluProject(triPoints[n].v[0], triPoints[n].v[1], triPoints[n].v[2],
                      modelViewMatrix, projectionMatrix, viewPort);

          origin = new StsPoint(screenCoor[2]);

          /** convert points to relative coordinates from origin */

          /** Vector from origin (point2) to points 0 and 1 */
          for(n = 0; n < 2; n++)
          {
              vectors[n] = new StsPoint(screenCoor[n]);
              vectors[n].v[0] -= origin.v[0];
              vectors[n].v[1] -= origin.v[1];
          }

          /** Vector from origin (point3) to mouse */
          float x = getX() - origin.v[0];
          float y = getY() - origin.v[1];

          /** area of triangle */
          float area = vectors[0].v[0] * vectors[1].v[1] - vectors[0].v[1] * vectors[1].v[0];

          if(Math.abs(area) < StsParameters.roundOff) return null;

          float[] w = new float[3];

          w[0] = -(vectors[1].v[0] * y - vectors[1].v[1] * x) / area;
          w[1] = (vectors[0].v[0] * y - vectors[0].v[1] * x) / area;
          w[2] = 1.0f - w[0] - w[1];

          float mouseZ = 0.0f;

          for(n = 0; n < 3; n++)
              mouseZ += w[n] * screenCoor[n][2];

          double[] point = getProjectedPoint(mousePoint, (double)mouseZ);
          //        System.out.println("triangle weights: " + w[0] + " " + w[1] + " " + w[2]);

          return getProjectedPoint(mousePoint, (double)mouseZ);
      }

      public StsPoint getPointOnLineNearestMouse(StsMouse mouse, StsPoint point0, StsPoint point1)
      {
          StsMousePoint mousePoint = mouse.getMousePoint();
          double[] screen0 = getScreenCoordinates(point0);
          double[] screen1 = getScreenCoordinates(point1);
          StsMousePoint mousePointGL = mousePoint.lowerLeftCoor(viewPort);
          double[] mouseXY = new double[]{(double)mousePointGL.x, (double)mousePointGL.y};
          double f = StsMath.nearestPointOnLineInterpolant(screen0, screen1, mouseXY, 2);
          //		System.out.println("StsGLPanel3d.getPointOnLineNearestMouse() f: " + f);
          return StsPoint.staticInterpolatePoints(point0, point1, (float)f);
      }

      public StsPoint getPointOnLineNearestViewLine(StsPoint point0, StsPoint point1)
      {
          double[] screen0 = getScreenCoordinates(point0);
          double[] screen1 = getScreenCoordinates(point1);
          double f = StsMath.nearestPointOnLineInterpolant(screen0, screen1, new double[]{0, 0}, 2);
          //		System.out.println("StsGLPanel3d.getPointOnLineNearestViewLine() f: " + f);
          return StsPoint.staticInterpolatePoints(point0, point1, (float)f);
      }

      public double[] getPickedViewPoint(StsMouse mouse)
      {

          //        StsMousePoint mousePoint = mouse.getMousePointGL(winRectGL);
          //        double mx = mousePoint.x;
          //        double my = winRectGL.height - mousePoint.y - 1;

          double[] point = gluUnProject((double)getX(), (double)getY(), 0.0,
                  modelViewMatrix, projectionMatrix, viewPort);
          return point;
      }

      public void verticalStretch(StsMouse mouse) { view.verticalStretch(mouse); }

      public void verticalShrink(StsMouse mouse) { view.verticalShrink(mouse); }

      public void horizontalStretch(StsMouse mouse) { view.horizontalStretch(mouse); }

      public void horizontalShrink(StsMouse mouse) { view.horizontalShrink(mouse); }

      public void setDefaultView()
      {
          view.setDefaultView();
      }

      public boolean isMainWindow()
      {
          return window.isMainWindow();
      }

      public void removeDisplayableClass(StsObject object)
      {
          view.removeDisplayableClass(object);
      }

      /*
          public void addView(StsView view)
          {
              if(hasView(view.getClass())) return;
              addView(view);
          }
      */
      public void initializeCurrentView()
      {
          if(view == null) return;
          view.initializeTransients(this);
          window.selectToolbarItem(StsViewSelectToolbar.NAME, view.getViewName(), true);
      }
  /*
      public StsView getNextView()
      {
          int nViews = views.length;
          if(nViews <= 1) return null;
          for(int n = 0; n < nViews; n++)
              if(views[n] == currentView) return views[(n + 1) % nViews];
          return null;
      }

      public StsView getPreviousView()
      {
          int nViews = views.length;
          if(nViews <= 1) return null;
          for(int n = 0; n < nViews; n++)
              if(views[n] == currentView) return views[(n + nViews - 1) % nViews];
          return null;
      }
  */
  /*
      public void setObject(Object object)
      {
          cursor3d.setObject(object);
          // because this sets the defaultAction for the actionManager,
          // we can't set the viewObject for all views as each might have a
          // different defaultAction causing conflicts for the actionManager.
          // TJL 1/28/07
          if(currentView != null) currentView.setViewObject(object);
          //        for(int n = 0; n < views.length; n++)
          //            views[n].setViewObject(object);
      }
  */
  /*
      public boolean setDefaultAction()
      {
          if(currentView == null) return false;
          return currentView.initializeDefaultAction();
      }
  */
      public void copy(StsGLPanel3d oldGLPanel3d)
      {
          /*
          boolean isGridCoordinates = !oldGLPanel3d.window.getCursorDisplayXY();
          isPerspective = oldGLPanel3d.isPerspective;
          cursor3d.setIsGridCoordinates(isGridCoordinates);
          cursor3d.copy(oldGLPanel3d.cursor3d);
          */
      }

      public float getZScale()
      {
          if(!(view instanceof StsView3d)) return 1.0f;
          return ((StsView3d)view).getZScale();
      }

      public float getAzimuth()
      {
          if(!(view instanceof StsView3d)) return 1.0f;
          return ((StsView3d)view).azimuth;
      }

      public float getElevation()
      {
          if(!(view instanceof StsView3d)) return 1.0f;
          return ((StsView3d)view).elevation;
      }

      public StsCursor3d getCursor3d() { return window.getCursor3d(); }

      public int getCurrentDirNo()
      {
          return window.getCursor3d().getCurrentDirNo();
      }
  /*
      public StsCursor3dTexture getDisplayableSection(int dirNo, Class displayableSectionClass)
      {
          if(displayableSections == null || displayableSections[dirNo] == null) return null;

          for(StsCursor3dTexture displayableSection : displayableSections[dirNo] )
              if(displayableSection.getDisplayableClass() == displayableSectionClass)
                  return displayableSection;
          return null;
      }

      public StsCursor3dTexture[] getDisplayableSections(int dirNo)
      {
          if(displayableSections == null || displayableSections[dirNo] == null)
              return new StsCursor3dTexture[0];
          else
              return displayableSections[dirNo];

      }

      public StsCursor3dTexture[] getVisibleDisplayableSections(int dirNo)
      {
          if(displayableSections == null)
              return new StsCursor3dTexture[0];
          StsCursor3dTexture[] cursorDisplayableSections = displayableSections[dirNo];
          if(cursorDisplayableSections == null)
              return new StsCursor3dTexture[0];
          StsCursor3dTexture[] visibleSections = new StsCursor3dTexture[cursorDisplayableSections.length];
          int nVisibleSections = 0;
          for(StsCursor3dTexture displayableSection : cursorDisplayableSections)
              if(displayableSection.isVisible) visibleSections[nVisibleSections++] = displayableSection;
          return (StsCursor3dTexture[])StsMath.trimArray(visibleSections, nVisibleSections);
      }

      StsCursor3dTexture checkAddDisplayableSection(StsClassCursor3dTextureDisplayable displayableClass, int dir, StsCursor3d cursor3d)
      {
          if(displayableClass.getCurrentObject() == null)return null;
          StsCursor3dTexture[] displayableSections = this.displayableSections[dir];
          Class objectClass = displayableClass.getObjectClass();
          for(StsCursor3dTexture displayableSection : displayableSections)
          {
              if(displayableSection.canDisplayClass(objectClass)) return displayableSection;
          }
          StsCursorSection cursorSection = cursor3d.cursorSections[dir];
          StsCursor3dTexture displayableSection = cursorSection.getDisplayableSection(objectClass);
          if(displayableSection != null)
          {
              addDisplayableSection(dir, displayableClass, displayableSection);
              return displayableSection;
          }

          displayableSection = displayableClass.constructDisplayableSection(model, cursor3d, dir);
          if(displayableSection == null) return null;
          cursorSection.addDisplayableSection(displayableSection);
          addDisplayableSection(dir, displayableClass, displayableSection);
          return displayableSection;
      }

      void addDisplayableSection(int dir, StsClassCursor3dTextureDisplayable displayableClass, StsCursor3dTexture displayableSection)
      {
          if(textureDebug) System.out.println("Adding displayableSection "+ displayableSection.toString() + " in direction " + dir);
          if(displayableClass.drawLast())
              displayableSections[dir] = (StsCursor3dTexture[])StsMath.arrayAddElement(displayableSections[dir], displayableSection);
          else
              displayableSections[dir] = (StsCursor3dTexture[])StsMath.arrayInsertElement(displayableSections[dir], displayableSection, 0);
      }
  */
      public void saveGeometry()
      {
          size = getSize();
      }
  /*
      public void adjustCursor(int dir, float dirCoor)
      {
         if(displayableSections == null || displayableSections[dir] == null) return;
          StsCursor3dTexture[] dirDisplayableSections = displayableSections[dir];
          for(StsCursor3dTexture displayableSection : dirDisplayableSections)
              displayableSection.setDirCoordinate(dir, dirCoor);
      }
  */
     public void focusGained(FocusEvent e)
      {
  //        window.getMouseActionToolbar().reconfigure(this);
          super.focusGained(e);
      }
  }
