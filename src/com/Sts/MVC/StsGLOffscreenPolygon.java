package com.Sts.MVC;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 * @version 1.2 jbw no color-index mode so use scaled RGBA for 32 distinct levels
 */

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.Utilities.*;
import com.sun.opengl.util.*;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;
import javax.swing.*;
import java.awt.*;
import java.nio.*;

public class StsGLOffscreenPolygon implements GLEventListener
{

    /** The OpenGL state machine */
    protected GL gl = null;
    protected GLU glu = null;

    /** Offscreen rendering buffer */
    protected GLPbuffer offScreenBuffer = null;

    /** Optional display panel */
    public GLCanvas glc = null;

    protected int width, height;
    protected double left, right, bottom, top;
    protected boolean display = false;
    protected IPolygon[] polygons = null;
    protected int[] coorIndexes = null;
//    protected int[] buffer = null;
    protected ByteBuffer displayBuffer = null;
    protected int currentPolygon = -1;
    protected boolean drawingConcave;
    protected int nPoints = 0;
    protected byte[] displayData;
	protected boolean renderToTexture = false;


    static StsColor[] colors = StsColor.colors8;

    static float[] rgba = new float[4];

    static double near = 0.0;
    static double far = 1000.0;



/*
    // four square convex polygons
    static double[][][] polygonPoints = new double[][][]
    {
        { { 0.00, 0.00,  0.0 }, { 0.25, 0.00,  0.0 }, { 0.25, 0.25,  0.0 }, { 0.00, 0.25,  0.0 } },
        { { 0.25, 0.00, -1.0 }, { 0.50, 0.00, -1.0 }, { 0.50, 0.25, -1.0 }, { 0.25, 0.25, -1.0 } },
        { { 0.00, 0.25, -2.0 }, { 0.25, 0.25, -2.0 }, { 0.25, 0.50, -2.0 }, { 0.00, 0.50, -2.0 } },
        { { 0.25, 0.25, -3.0 }, { 0.50, 0.25, -3.0 }, { 0.50, 0.50, -3.0 }, { 0.25, 0.50, -3.0 } }
    };
*/
    static StsConcavePolygon concavePolygon = null;
    static final float[] normal = new float[] { 0, 0, 1 };

    public StsGLOffscreenPolygon(int width, int height, float left, float right, float bottom, float top, boolean display)
    {
        this.width = width;
        this.height = height;
        this.left = left;
        this.right = right;
        this.bottom = bottom;
        this.top = top;
        this.display = display;

        float border = 0.5f*(right - left)/(width-1);
        left -= border;
        right += border;

        border = 0.5f*(top - bottom)/(height-1);
        bottom -= border;
        top += border;

        nPoints = width*height;
        displayData = new byte[nPoints];

		GLCapabilities cap = new GLCapabilities();
//        offscreenCap.setDepthBits( 24 );
	   cap.setRedBits(8);
	   cap.setGreenBits(8);
	   cap.setBlueBits(8);
	   cap.setDoubleBuffered( false);
	   // no longer exists in jogl b3;
	   //cap.setOffscreenRenderToTextureRectangle(false);

	   try
       {
		   offScreenBuffer = GLDrawableFactory.getFactory().createGLPbuffer(cap, null, width, height, null);
	   }
       catch (Exception E)
       {
		   System.err.println("Failed to create offscreen, retry...");
		   E.printStackTrace();
		   renderToTexture=false;
		   //cap.setOffscreenRenderToTextureRectangle(false);
		   try
           {
		    offScreenBuffer = GLDrawableFactory.getFactory().createGLPbuffer(cap, null, width, height, null);
		   }
           catch (Exception E2)
           {
		      E2.printStackTrace();
			  System.err.println("Failed utterly to create offscreen");
		   }
	   }
	   if (offScreenBuffer == null)
       {
		   System.err.println("Failed utterly to create offscreen");
		   return;
	   }

	   gl = offScreenBuffer.getGL();
	   glu = new GLU();
       if(!display) return;

       JDialog dialog = new JDialog((Frame)null, "Debug display", false);

       /** Create the OpenGL Component */
	   //cap.setOffscreenRenderToTextureRectangle(false);
        glc = new GLCanvas(cap);

        dialog.getContentPane().add(glc);
		dialog.setPreferredSize(new Dimension(width,height));
        dialog.pack();
        dialog.setVisible(true);

   }

    /** Returns the width of the demo */
    public int getWidth() { return width; }

    /** Returns the height of the demo */
    public int getHeight() { return height; }

    /** Gets called exactly once when GLComponent.classInitialize() is invoked */
    public void init(GLAutoDrawable component)
    {
    }

    public void startGL()
    {
        /** Add the listener and classInitialize the offscreen buffer */
        offScreenBuffer.addGLEventListener( this );
//        offScreenBuffer.classInitialize();

        if(glc == null) return;

        /** Initialize the onscreen window */
        glc.addGLEventListener( this );
//        glc.classInitialize();
    }

    /** Handles viewport resizing */
    public void reshape( GLAutoDrawable component, int x, int y, int width, int height )
    {
        gl.glViewport(0, 0, width, height );
        gl.glMatrixMode( GL.GL_PROJECTION );
        gl.glLoadIdentity();
        glu.gluOrtho2D(left, right, bottom, top);
        gl.glMatrixMode( GL.GL_MODELVIEW );
        gl.glLoadIdentity();
    }

	public void displayChanged(GLAutoDrawable gLDrawable, boolean modeChanged, boolean deviceChanged) {}

    public GL getGL() { return gl; }

    public void setPolygons(IPolygon[] polygons)
    {
        this.polygons = polygons;
        coorIndexes = null;
    }

    public void setPolygons(IPolygon[] polygons, int[] coorIndexes)
    {
        this.polygons = polygons;
        this.coorIndexes = coorIndexes;
    }

    public void setPolygons(StsObjectRefList polygonList)
    {
        int nPolygons = polygonList.getSize();
        if(nPolygons == 0) return;

        polygons = new IPolygon[nPolygons];
        for(int n = 0; n < nPolygons; n++)
            polygons[n] = (IPolygon)polygonList.getElement(n);

        coorIndexes = null;
    }

    public void display()
    {
        display(offScreenBuffer);
    }
    /** Renders the scene */
    public void display(GLAutoDrawable component)
    {
        if(component == offScreenBuffer) displayOffScreen(component);
        else if(component == glc)  displayOnScreen(component);
    }

    public void displayOffScreen(GLAutoDrawable component)
    {
        System.out.println("displayOffScreen called.");
        gl.glMatrixMode( GL.GL_PROJECTION );
		gl.glLoadIdentity();
		glu.gluOrtho2D(left, right, bottom, top);
		gl.glMatrixMode( GL.GL_MODELVIEW );
		gl.glLoadIdentity();

        //gl.glClearIndex(0.0f);
		gl.glClearColor(0.f,0.f,0.f,0.f);
        gl.glClear(GL.GL_COLOR_BUFFER_BIT|GL.GL_DEPTH_BUFFER_BIT);
		gl.glDisable(GL.GL_DEPTH_TEST);
//        gl.glEnable(GL.GL_DEPTH_TEST);
//        gl.glDepthFunc(GL.GL_GREATER);
//        gl.glClearDepth(near);
//        gl.glClear(GL.GL_DEPTH_BUFFER_BIT);

        drawConcavePolygons();
		gl.glFinish();
		gl.glReadBuffer(GL.GL_FRONT);
        // we should be able to get bytes directly, but GL fails for some unknown reason,
        // so we get a float array and convert to bytes.
        ByteBuffer buffer = BufferUtil.newByteBuffer(nPoints*4);
        gl.glReadPixels(0, 0, width, height, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE, buffer);
		displayData = new byte[nPoints];
		buffer.rewind();
		for(int n = 0; n < nPoints; n++) {
			byte r = buffer.get();
			byte g = buffer.get();
			byte b = buffer.get();
			byte a = buffer.get();
			displayData[n]=(byte) (((r) / 8.f) + 0.5); // jbw 32 levels
			//if (isDisplayData[n] != 0) System.out.println("!= 0"+isDisplayData[n]+" "+n+" "+r);
		}

        if(display)
        {
			displayBuffer = BufferUtil.newByteBuffer(4*nPoints);
            gl.glReadPixels(0, 0, width, height, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE, displayBuffer);
			displayBuffer.rewind();

        }
    }

    private void displayOnScreen(GLAutoDrawable component)
    {
        System.out.println("displayOnScreen called.");
        /** Clear the framebuffer */
		gl.glMatrixMode( GL.GL_PROJECTION );
		gl.glLoadIdentity();
		glu.gluOrtho2D(left, right, bottom, top);
		gl.glMatrixMode( GL.GL_MODELVIEW );
		gl.glLoadIdentity();

        gl.glEnable(GL.GL_DEPTH_TEST);
        gl.glDepthFunc(GL.GL_GREATER);
        gl.glClearDepth(near);
        gl.glClear(GL.GL_DEPTH_BUFFER_BIT);

        /** Copy the offscreen buffer contents to this window */
        if(displayBuffer == null) return;
        gl.glRasterPos2i( 0, 0 );
		displayBuffer.rewind();
        gl.glDrawPixels( width, height, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE, displayBuffer );
    }

    private void drawPolygons()
    {
        int nPolygons = polygons.length;
        for(int n = 0; n < nPolygons; n++)
        {
            StsColor color = colors[n%colors.length];
            color.setGLColor(gl);

            gl.glBegin( GL.GL_POLYGON );

            double[][] points = polygons[n].getPoints();
            for(int i = 0; i < points.length; i++) gl.glVertex3dv( points[i], 0 );

            gl.glEnd();
        }
    }

    public void close()
    {
        if(offScreenBuffer != null)
        {
            offScreenBuffer.removeGLEventListener(this);

            offScreenBuffer.destroy();
        }
        if(glc != null)
        {
            glc.removeGLEventListener(this);
            //glc.destroy();
        }
    }

    private void drawConcavePolygons()
    {
        int nPolygons = polygons.length;
        for(int n = 0; n < nPolygons; n++)
        {
//            Color color = colors[n%colors.length];
//            gl.glColor3fv(color.getComponents(rgba));
            boolean ok = drawConcavePolygon(polygons[n], n);
            if(!ok) System.out.println(polygons[n].getLabel() + "failed.");
        }
    }

    private boolean drawConcavePolygon(IPolygon polygon, int ID)
    {
        double[][] points = polygon.getPoints();
        if(points == null || points.length < 3) return false;
        int xCoor = 0, yCoor = 1;
//        double z = getDepthFromID(polygon.getID());
//        int ID = polygon.getID();
//        double z = far - 1 - (double)ID;
        double[] point;

        try
        {
            //gl.glIndexi(ID + 1);

			gl.glColor3f((ID+1)/32.f,(ID+1)/32.f,(ID+1)/32.f); // jbw
			//System.out.println("Index = "+(ID+1));
			if(concavePolygon == null) concavePolygon = new StsConcavePolygon(gl, glu);
			GLUtessellator tesselator = glu.gluNewTess();
            String name = polygon.getLabel();
            concavePolygon.initialize(tesselator, polygon, false, true, normal, name);
            currentPolygon = ID;
            drawingConcave = true;

            glu.gluTessBeginPolygon(tesselator,null);
			glu.gluTessBeginContour(tesselator);
            glu.gluTessNormal(tesselator, 0, 0, 1);

            if(coorIndexes != null)
            {
                xCoor = coorIndexes[0];
                yCoor = coorIndexes[1];
            }

            for(int n = 0; n < points.length; n++)
            {
                point = new double[] { points[n][xCoor], points[n][yCoor], 0 };
                glu.gluTessVertex(tesselator, point, 0, point);
            }
			glu.gluTessEndContour(tesselator);
            glu.gluTessEndPolygon(tesselator);
            return drawingConcave;
        }
        catch(Exception e)
        {
            StsException.outputException("StsGLOffscreenPolygon.drawConcavePolygon() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    public static void main( String argv[] )
    {
        StsGLOffscreenPolygon offScreenBuffer = new StsGLOffscreenPolygon(255, 255, 0.0f, 1.0f, 0.0f, 1.0f, false);
        OffscreenPolygon polygon = new OffscreenPolygon();
        IPolygon[] polygons = new IPolygon[] { polygon };
        offScreenBuffer.setPolygons(polygons);
        offScreenBuffer.startGL();
    }

    public String getLabel()
    {
        return new String("polygon: " + currentPolygon);
    }

    public void drawConcaveFailed(int error)
    {
        drawingConcave = false;
    }

	public void drawConcaveFailed(String errorString)
	{
		drawingConcave = false;
	}

    public void repaint()
    {
       offScreenBuffer.display();
	   if (glc!= null)
		   glc.display();
    }

    public byte[] getData() { return displayData; }
}

class OffscreenPolygon implements IPolygon
{
    double[][] points = new double[][]
        { { 0.00, 0.00, 0.0 }, { 0.50, 0.00, 0.0 }, { 0.50, 0.25f, 0.0 },
          { 0.25f, 0.25f, 0.0 }, { 0.25f, 0.50, 0.0 }, { 0.00, 0.50, 0.0 } };

    OffscreenPolygon()
    {
    }

    public double[][] getPoints() { return points; }
    public int getNPoints() { return points.length; }
    public String getLabel() { return "polygon"; }
    public int getID() { return 0; }
    public void drawConcaveFailed(int error) { }
	public void drawConcaveFailed(String error) { }
}
