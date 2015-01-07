//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Utilities;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Reflect.*;
import com.Sts.Types.*;
import com.sun.opengl.util.*;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;
import java.nio.*;

/**
 * Static methods for picking
  *
  * @author TJLasseter
  */

 public class StsJOGLPick implements StsGLPicker
 {
     private StsGLPanel glPanel;
     private StsObject[] objects;
     private int pickSize;
     private int pickType;
     private StsMethod pickMethod;
     public StsObject pickedObject = null;

     private GL gl;
     private GLU glu;

     static public int hits;
     /** number of pick hits */
     static public StsPickItem[] pickItems = null;
     /** processed pick items */
     static private final int pickBufferSize = 100000;
     /** size of raw pick buffer */
     static private IntBuffer pickBuffer = BufferUtil.newIntBuffer(pickBufferSize); /** pick buffer				*/

     /** pickSize in pixels */
     static public final int PICKSIZE_SMALL = 2;
     static public final int PICKSIZE_MEDIUM = 5;
     static public final int PICKSIZE_LARGE = 10;
     static public final int PICKSIZE_EXTRA_LARGE = 25;

     /** pickType */
     static public final int PICK_ALL = 1;
     static public final int PICK_CLOSEST = 2;
     static public final int PICK_FIRST = 3;
     static public final int PICK_LAST = 4;

     private static boolean debug = false;

     private StsJOGLPick(StsGLPanel glPanel, StsObject[] objects, int pickSize, int pickType)
     {
         this.glPanel = glPanel;
         this.objects = objects;
         this.pickSize = pickSize;
         this.pickType = pickType;
         gl = glPanel.getGL();
         glu = glPanel.getGLU();
     }

     public StsJOGLPick(StsGLPanel glPanel, StsMethod pickMethod, int pickSize, int pickType)
     {
         this.glPanel = glPanel;
         this.pickMethod = pickMethod;
         this.pickSize = pickSize;
         this.pickType = pickType;
         gl = glPanel.getGL();
         glu = glPanel.getGLU();
     }

     private StsJOGLPick(StsGLPanel glPanel, int pickSize)
     {
         this.glPanel = glPanel;
         this.pickSize = pickSize;
         gl = glPanel.getGL();
         glu = glPanel.getGLU();
     }

     static public StsJOGLPick initializeDebugPickMatrix(StsGLPanel glPanel, int pickSize)
     {
         StsJOGLPick picker = new StsJOGLPick(glPanel, pickSize);
         picker.applyPickMatrix(glPanel, pickSize);
         return picker;
     }

     static public StsObject pickClass3d(StsGLPanel glPanel, StsClass stsClass, int pickSize, int pickType)
     {
         if (stsClass == null || stsClass.getSize() == 0) return null;
         return pickClass3d(glPanel, stsClass.getElements(), pickSize, pickType);
     }

     static public StsObject pickVisibleClass3d(StsGLPanel glPanel, Class c, int pickSize, int pickType)
     {
         StsModel model = StsObject.getCurrentModel();
         StsObject[] objects = model.getVisibleObjectList(c);
         if (objects == null || objects.length == 0) return null;
         return pickClass3d(glPanel, objects, pickSize, pickType);
     }

     static public StsObject pickVisibleClass3d(StsGLPanel glPanel, StsClass stsClass, int pickSize, int pickType)
     {
         StsModel model = StsObject.getCurrentModel();
         StsObject[] objects = model.getVisibleObjectList(stsClass);
         if (objects == null || objects.length == 0) return null;
         return pickClass3d(glPanel, objects, pickSize, pickType);
     }

     static public StsObject pickClass3d(StsGLPanel glPanel, StsObjectRefList refList, int pickSize, int pickType)
     {
         if (refList == null || refList.getSize() == 0) return null;
         return pickClass3d(glPanel, refList.getElements(), pickSize, pickType);
     }

     static public StsObject pickClass3d(StsGLPanel glPanel, StsObject[] objects, int pickSize, int pickType)
     {
         if (objects == null || objects.length == 0) return null;
         StsJOGLPick picker = new StsJOGLPick(glPanel, objects, pickSize, pickType);
         StsGLPanel3d panel3d = (StsGLPanel3d) glPanel;
         //System.out.println("StsJOGLPick.  calling doPick() ");
         panel3d.doPick(picker);

         if (debug)
         {
             String classname = StsToolkit.getSimpleClassname(objects[0]);
             if (picker.pickedObject != null)
                 System.out.println("StsJOGLPick.  return from doPick(). Class: " + classname + " pickedObject: " + picker.pickedObject.getIndex());
             else
                 System.out.println("StsJOGLPick.  return from doPick(). Class: " + classname + " pickedObject: null");
         }
         return picker.pickedObject;
     }

     public boolean execute(GLEventListener caller, GLAutoDrawable glDrawable) throws Exception
     {
         //System.out.println("StsJOGLPick.execute() called.");

         if (glPanel == null && objects == null) return false;
         if (!parmsAreValid(pickSize, pickType)) return false;
         boolean pickHit;

         //  		glPanel.getGL().setGLContext();
         setupPickMatrix(glPanel, pickSize);

         try
         {
             // execute display on GLEventQueue to actually do picking
             //			System.out.println("JOGLPick Executing display");
             GLContext ctx = GLContext.getCurrent();
             //			System.out.println("context is "+ctx);
             boolean synced = false;
             if (ctx != null)
                 synced = ctx.isSynchronized();
             //			if (ctx != null)
             //				System.out.println("is sync "+synced);
             if (ctx != null)
                 ctx.setSynchronized(true);

             pickBuffer.rewind(); // jbw
             caller.display(glDrawable);
             if (Main.isGLDebug) System.out.println("StsJOGLPick.execute() called for glPanel " + glPanel);

             //glDrawable.display();
             //System.out.println("JOGLPick Executed display");
             //            if (!doPickClass3d(glPanel, objects)) return false;
             resetMatrix(gl);
             gl.glFlush();
             pickBuffer.rewind(); // tjl
             pickHit = processHits(pickBuffer, pickType);
             //System.out.println("JOGLPick Executed process");
             if (!pickHit) return false;

             /* convert the pick info to an StsObject */
             if (objects != null)
             {
                 StsPickItem pickItem = pickItems[0];
                 if (pickItem.names == null) return false;
                 pickedObject = objects[pickItem.names[0]];
             }
             return true;
         }
         catch (Exception e)
         {
             System.out.println("JOGLPick Execute failed");
             StsException.outputException("StsJOGLPick execute failed", e, StsException.WARNING);
             e.printStackTrace();
             // jbw make sure to get out of opengl pick mode
             gl.glRenderMode(GL.GL_RENDER);
             resetMatrix(gl);
             return false;
         }
         /*
         finally
         {
          // restore the perspective matrix
          gl.glMatrixMode(GL.GL_PROJECTION);
          gl.glPopMatrix();
          gl.glFlush();
         }
         */
     }

     public boolean doPick(GLAutoDrawable glDrawable)
     {
         try
         {
             if (pickMethod != null)
             {
                 pickMethod.getMethod().invoke(pickMethod.getInstance(), pickMethod.getMethodArgs());
                 return true;
             }
             else
                 return doPickClass3d();
         }
         catch (Exception e)
         {
             StsException.systemError("StsJOGLPick.doPick() failed.");
             return false;
         }
     }

     static private boolean parmsAreValid(int pickSize, int pickType)
     {
         switch (pickSize)
         {
             case PICKSIZE_SMALL:
             case PICKSIZE_MEDIUM:
             case PICKSIZE_LARGE:
             case PICKSIZE_EXTRA_LARGE:
                 break;
             default:
                 return false;
         }

         switch (pickType)
         {
             case PICK_ALL:
             case PICK_CLOSEST:
             case PICK_FIRST:
             case PICK_LAST:
                 break;
             default:
                 return false;
         }

         return true;
     }

     static public boolean pick3d(StsGLPanel glPanel, StsMethod pickMethod, int pickSize, int pickType)
     {
         StsJOGLPick picker = new StsJOGLPick(glPanel, pickMethod, pickSize, pickType);
         return picker.methodPick3d();
     }

     public boolean methodPick3d()
     {
         StsGLPanel3d panel3d = (StsGLPanel3d) glPanel;
         //System.out.println("StsJOGLPick.  calling doPick() ");
         panel3d.doPick(this);
         return hits > 0;
     }

     /*
     IntBuffer selectBuffer = BufferUtils.newIntBuffer(buffsize);
      int hits = 0;

      gl.glGetIntegerv(GL.GL_VIEWPORT, viewPort);
      gl.glSelectBuffer(buffsize, selectBuffer);
      gl.glRenderMode(GL.GL_SELECT);
      gl.glInitNames();
      gl.glMatrixMode(GL.GL_PROJECTION);
      gl.glPushMatrix();
      gl.glLoadIdentity();
      glu.gluPickMatrix(x, (double) viewPort[3] - y, 5.0d, 5.0d, viewPort);
      glu.gluOrtho2D(0.0d, 1.0d, 0.0d, 1.0d);
      drawable.display();
      gl.glMatrixMode(GL.GL_PROJECTION);
      gl.glPopMatrix();
      gl.glFlush();
      hits = gl.glRenderMode(GL.GL_RENDER);
     */

     public void setupPickMatrix(StsGLPanel glPanel, int pickSize)
     {
         int[] viewPort = new int[4];
         //double[] pickMatrix = new double[16];

         gl = glPanel.getGL();
         if (gl == null) return;

         glu = glPanel.getGLU();
         if (glu == null) return;

         gl.glGetIntegerv(GL.GL_VIEWPORT, viewPort, 0);
         pickBuffer.rewind(); // jbw
         gl.glSelectBuffer(pickBufferSize, pickBuffer);
         gl.glRenderMode(GL.GL_SELECT);

         gl.glInitNames();
         gl.glPushName(-1);

         gl.glMatrixMode(GL.GL_PROJECTION);

         /** Save the current projection matrix */
         gl.glPushMatrix();

         /** Compute the pick matrix */
         gl.glLoadIdentity();
         StsMousePoint mousePoint = glPanel.getMousePointGL();

         glu.gluPickMatrix(mousePoint.x, mousePoint.y, pickSize, pickSize, viewPort, 0);

         /** get pickMatrix for debugging */
         // gl.glGetDoublev(GL.GL_PROJECTION_MATRIX, pickMatrix);

         /** multiply the pick matrix (on stack) by the projection matrix */
         gl.glMultMatrixd(glPanel.getProjectionMatrix(), 0);
         gl.glMatrixMode(GL.GL_MODELVIEW);
     }
     public void applyPickMatrix(StsGLPanel glPanel, int pickSize)
      {
          double[] pickMatrix = new double[16];

          gl = glPanel.getGL();
          if (gl == null) return;

          glu = glPanel.getGLU();
          if (glu == null) return;

          int[] viewPort = new int[4];
          gl.glGetIntegerv(GL.GL_VIEWPORT, viewPort, 0);
          /** Save the current projection matrix */
          gl.glPushMatrix();

          /** Compute the pick matrix */
          gl.glLoadIdentity();
          StsMousePoint mousePoint = glPanel.getMousePointGL();

          glu.gluPickMatrix(mousePoint.x, mousePoint.y, pickSize, pickSize, viewPort, 0);
                  /** get pickMatrix for debugging */
          gl.glGetDoublev(GL.GL_PROJECTION_MATRIX, pickMatrix, 0);
          /** multiply the pick matrix (on stack) by the projection matrix */
          gl.glMultMatrixd(glPanel.getProjectionMatrix(), 0);
      }

     static public void resetMatrix(GL gl)
     {
         gl.glMatrixMode(GL.GL_PROJECTION);
         gl.glPopMatrix();
         gl.glMatrixMode(GL.GL_MODELVIEW);
         gl.glFlush();
     }

     public boolean doPickClass3d()
     {
         if (glPanel == null || objects == null) return false;

         for (int n = 0; n < objects.length; n++)
         {
             //System.out.println("Object index " + objects[n].index() + " drawn for pick.");
             if (objects[n] == null) break;
             gl.glInitNames();
             gl.glPushName(n);
             try
             {
                 objects[n].pick(gl, glPanel);
             }
             catch (Exception e)
             {
                 return false;
             }
         }
         return true;
     }

     /** Processes the hits from the selection process */
     private boolean processHits(IntBuffer intBuffer, int pickType)
     {
         int index = 0;
         int i, n;
         StsPickItem pickItem;
         int closestHit = 0;
         intBuffer.rewind();
         try
         {
             hits = gl.glRenderMode(GL.GL_RENDER);
             if (debug) StsException.systemDebug(this, "processHits", " hits = " + hits);
             if (hits < 0) StsMessageFiles.errorMessage("WARNING: pick buffer overflow.");
             if (hits <= 0) return false;
             if (pickType == PICK_ALL)
             {
                 pickItems = new StsPickItem[hits];
                 for (i = 0; i < hits; i++)
                 {
                     pickItem = pickItems[i] = new StsPickItem();

                     pickItem.noNames = intBuffer.get();
                     pickItem.zMin = computePickZ(intBuffer.get());
                     pickItem.zMax = computePickZ(intBuffer.get());

                     if (debug) System.out.print("    Pick: " + i +
                             " zMin: " + pickItem.zMin +
                             " zMax: " + pickItem.zMax +
                             " noNames: " + pickItem.noNames);

                     pickItem.names = new int[pickItem.noNames];
                     intBuffer.get(pickItem.names);
                     if (debug)
                     {
                         System.out.print(" names: ");
                         for (n = 0; n < pickItem.noNames; n++)
                             System.out.print(n + ": " + pickItem.names[n] + " ");
                         System.out.println(" ");
                     }
                 }
             }
             else if (pickType == PICK_CLOSEST)
             {
                 pickItems = new StsPickItem[1];
                 pickItem = pickItems[0] = new StsPickItem();

                 pickItem.zMin = StsParameters.largeFloat;

                 for (i = 0; i < hits; i++)
                 {
                     int noNames = intBuffer.get();
                     if (noNames == 0) continue;
                     float zMin = computePickZ(intBuffer.get());
                     float zMax = computePickZ(intBuffer.get());

                     if (debug) System.out.println("Pick: " + i +
                             " zMin: " + zMin +
                             " zMax: " + zMax +
                             " noNames: " + noNames);

                     int[] names = new int[noNames];
                     intBuffer.get(names);
                     if (debug)
                     {
                         System.out.print("   names: ");
                         for (n = 0; n < noNames; n++)
                             System.out.print(n + ": " + names[n] + " ");
                         System.out.println(" ");
                     }
                     if (zMin < pickItem.zMin)
                     {
                         closestHit = i;
                         pickItem.noNames = noNames;
                         pickItem.zMin = zMin;
                         pickItem.zMax = zMax;
                         pickItem.names = names;
                     }
                 }
                 if (debug) System.out.println("Closest Pick: " + closestHit +
                         " zMin: " + pickItem.zMin +
                         " zMax: " + pickItem.zMax +
                         " noNames: " + pickItem.noNames);
             }
             else if (pickType == PICK_FIRST)
             {
                 pickItems = new StsPickItem[1];
                 pickItem = pickItems[0] = new StsPickItem();

                 int noNames = intBuffer.get();
                 if (noNames == 0) return false;
                 float zMin = computePickZ(intBuffer.get());
                 float zMax = computePickZ(intBuffer.get());
                 int[] names = new int[noNames];
                 intBuffer.get(names);
                 if (debug)
                 {
                     System.out.print("   names: ");
                     for (n = 0; n < noNames; n++)
                         System.out.print(n + ": " + names[n] + " ");
                     System.out.println(" ");
                 }
                 pickItem.noNames = noNames;
                 pickItem.zMin = zMin;
                 pickItem.zMax = zMax;
                 pickItem.names = names;
                 if (debug) System.out.println("First Pick: " + " noNames: " + pickItem.noNames);
             }
             else if (pickType == PICK_LAST)
             {
                 pickItems = new StsPickItem[1];
                 pickItem = pickItems[0] = new StsPickItem();
                 int noNames = 0;
                 float zMin = 0.0f, zMax = 0.0f;
                 int[] names = null;
                 for (i = 0; i < hits; i++)
                 {
                     noNames = intBuffer.get();
                     if (noNames == 0) continue;
                     zMin = computePickZ(intBuffer.get());
                     zMax = computePickZ(intBuffer.get());
                     names = new int[noNames];
                     intBuffer.get(names);
                 }
                 pickItem.noNames = noNames;
                 pickItem.zMin = zMin;
                 pickItem.zMax = zMax;
                 pickItem.names = names;
                 if (debug) System.out.println("Last Pick: " + " noNames: " + pickItem.noNames);
             }
             return true;
         }
         catch (Exception e)
         {
             StsException.outputException("StsJOGLPick.processHits() failed.",
                     e, StsException.WARNING);
             return false;
         }
     }

     static public float computePickZ(int intZ)
     {
         // ~ is the bitwise inversion operator
         float z1 = (float) ((float) ~intZ / (float) 0x7fffffffl);
         float z2 = (float) ((float) intZ / (float) 0x7fffffffl);
         float z = (1.0f + z1) / 2;
         return z;
         //    	return 1.0f - (float)( (float)~intZ / (float)0x7fffffffl );
     }

     static public float getScreenZ()
     {
         if (pickItems == null || pickItems.length <= 0) return StsParameters.largeFloat;
         return pickItems[0].zMin;
     }

     static public boolean mousePicksPoint(StsMousePoint mousePoint, StsPoint point, StsGLPanel3d glPanel3d, int pickSize)
     {
         int dif;
         double[] screenCoor = glPanel3d.getScreenCoordinates(point);

         dif = Math.abs(mousePoint.x - (int) screenCoor[0]);
         if (dif > pickSize) return false;

         dif = Math.abs((int) (glPanel3d.viewPort[3] - mousePoint.y - screenCoor[1]));
         if (dif > pickSize) return false;

         return true;
     }
 }
