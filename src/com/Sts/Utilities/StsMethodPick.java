
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Utilities;

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
  * @author TJLasseter
  */

public class StsMethodPick implements StsGLPicker
{
  	public int hits;               						/** number of pick hits 	*/
    public StsPickItem[] pickItems = null;				/** processed pick items 	*/
  	private int pickBufferSize = 1000;                   /** size of raw pick buffer */
	private IntBuffer pickBuffer = BufferUtil.newIntBuffer(pickBufferSize);	/** pick buffer				*/

	private GL gl;
    private GLU glu;

    /** pickSize in pixels */
    static public final int PICKSIZE_SMALL  = StsJOGLPick.PICKSIZE_SMALL;
   	static public final int PICKSIZE_MEDIUM = StsJOGLPick.PICKSIZE_MEDIUM;
   	static public final int PICKSIZE_LARGE  = StsJOGLPick.PICKSIZE_LARGE;
   	static public final int PICKSIZE_EXTRA_LARGE  = StsJOGLPick.PICKSIZE_EXTRA_LARGE;

    /** pickType */
    static public final int PICK_ALL     = StsJOGLPick.PICK_ALL;
    static public final int PICK_CLOSEST = StsJOGLPick.PICK_CLOSEST;
	static public final int PICK_FIRST = StsJOGLPick.PICK_FIRST;
	static public final int PICK_LAST = StsJOGLPick.PICK_LAST;

	private StsGLPanel glPanel;
	private StsMethod pickMethod;
	private int pickSize = PICKSIZE_EXTRA_LARGE;
	private int pickType = PICK_CLOSEST;

	static final boolean debug = false;

    public StsMethodPick( StsGLPanel glPanel, StsMethod pickMethod, int pickSize, int pickType)
    {
    	this.glPanel = glPanel;
    	this.pickMethod = pickMethod;
    	this.pickSize = pickSize;
    	this.pickType = pickType;
		gl = glPanel.getGL();
        glu = glPanel.getGLU();
    }

	static public boolean pick(StsGLPanel glPanel, StsMethod pickMethod, int pickSize, int pickType)
	{
		StsMethodPick picker = new StsMethodPick(glPanel, pickMethod, pickSize, pickType);
		return picker.methodPick();
	}

	public boolean methodPick()
	{
		StsGLPanel3d panel3d = (StsGLPanel3d)glPanel;
		//System.out.println("StsJOGLPick.  calling doPick() ");
		panel3d.doPick(this);
		return hits > 0;
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

	public boolean execute(GLEventListener listener, GLAutoDrawable glDrawable) throws Exception
	{
        setupPickMatrix(glPanel, pickSize, pickType);

        try
        {
            listener.display(glDrawable);
            return processHits(pickBuffer, pickType);
        }
        catch(Exception e)
        {
            StsException.outputException("StsMethodPick.pick() failed. Pick method: " + pickMethod.getMethod() +
                                         " instance: " + pickMethod.getInstance() + " number of args: " +
                                         pickMethod.getMethodArgs().length, e, StsException.WARNING);
            return false;
        }
        finally
        {
            /** restore the perspective matrix */
            gl.glMatrixMode(GL.GL_PROJECTION);
            gl.glPopMatrix();
            gl.glMatrixMode(GL.GL_MODELVIEW);
        }
	}


	public boolean doPick(GLAutoDrawable glDrawable)
	{
		try
		{
//        System.out.println("StsJOGLPick.doPick() called.");
			pickMethod.getMethod().invoke(pickMethod.getInstance(), pickMethod.getMethodArgs());
			return true;
		}
		catch(Exception e)
		{
			StsException.systemError("StsMethodPick.doPick() failed.");
			return false;
        }
	}

    public void setupPickMatrix(StsGLPanel glPanel, int pickSize, int pickType)
    {
  		int[] viewPort = new int[4];
        //double[] pickMatrix = new double[16];

        GL gl = glPanel.getGL();
        if(gl == null) return;

        GLU glu = glPanel.getGLU();
        if(glu == null) return;

        gl.glGetIntegerv(GL.GL_VIEWPORT, viewPort, 0);

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

   /** Processes the hits from the selection process */
    private boolean processHits(IntBuffer intBuffer, int pickType)
    {
        int index = 0;
        int i, n;
        StsPickItem pickItem;
        int closestHit = 0;

        try
        {
        	GL gl = glPanel.getGL();
            hits = gl.glRenderMode(GL.GL_RENDER);

         //   System.out.println( "Hits = " + hits );

            if(hits < 0) StsMessageFiles.errorMessage("WARNING: pick buffer overflow.");

            if(hits <= 0) return false;

            if(pickType == PICK_ALL)
            {
                pickItems = new StsPickItem[hits];

                for(i = 0; i < hits; i++)
                {
                    pickItem = pickItems[i] = new StsPickItem();

                    pickItem.noNames = intBuffer.get();
                    pickItem.zMin = computePickZ(intBuffer.get());
                    pickItem.zMax = computePickZ(intBuffer.get());
                /*
                    System.out.println("Pick: " + i +
                                       " zMin: " + pickItem.zMin +
                                       " zMax: " + pickItem.zMax +
                                       " noNames: " + pickItem.noNames);
                */
                    pickItem.names = new int[pickItem.noNames];

                    // System.out.print("   names: ");

                    for(n = 0; n < pickItem.noNames; n++)
                    {
                        pickItem.names[n] = intBuffer.get();
                    //    System.out.print(n + ": " + pickItem.names[n] + " ");
                    }

                    // System.out.println(" ");
                }
            }
            else if(pickType == PICK_CLOSEST)
            {
                pickItems = new StsPickItem[1];
                pickItem = pickItems[0] = new StsPickItem();

                pickItem.zMin = StsParameters.largeFloat;

                for(i = 0 ; i < hits ; i++ )
                {
                    int noNames = intBuffer.get();
                    float zMin = computePickZ(intBuffer.get());
                    float zMax = computePickZ(intBuffer.get());
                /*
                    System.out.println("Pick: " + i +
                                       " zMin: " + zMin +
                                       " zMax: " + zMax +
                                       " noNames: " + noNames);
                */
                    if(zMin < pickItem.zMin)
                    {
                        closestHit = i;
                        pickItem.noNames = noNames;
                        pickItem.zMin = zMin;
                        pickItem.zMax = zMax;

                        pickItem.names = new int[pickItem.noNames];

                  // 		System.out.print("   names: ");

                        for(n = 0 ; n < pickItem.noNames; n++)
                        {
                            pickItem.names[n] = intBuffer.get();
                  //      	System.out.print(n + ": " + pickItem.names[n] + " ");
                        }

                  //      System.out.println(" ");
                    }
                    else
                        index += noNames;
                }
             /*
                System.out.println("Closest Pick: " + closestHit +
                                   " zMin: " + pickItem.zMin +
                                   " zMax: " + pickItem.zMax +
                                   " noNames: " + pickItem.noNames);

                for(n = 0 ; n < pickItem.noNames; n++)
                    System.out.print(n + ": " + pickItem.names[n] + " ");
             */
            }
			else if(pickType == PICK_FIRST)
			{
				pickItems = new StsPickItem[1];
				pickItem = pickItems[0] = new StsPickItem();

				int noNames = intBuffer.get();
				if(noNames == 0) return false;
				float zMin = computePickZ(intBuffer.get());
				float zMax = computePickZ(intBuffer.get());
				int[] names = new int[noNames];
				intBuffer.get(names);
				if(debug)
				{
					System.out.print("   names: ");
					for(n = 0; n < noNames; n++)
						System.out.print(n + ": " + names[n] + " ");
					System.out.println(" ");
				}
				pickItem.noNames = noNames;
				pickItem.zMin = zMin;
				pickItem.zMax = zMax;
				pickItem.names = names;
				if(debug) System.out.println("First Pick: " +  " noNames: " + pickItem.noNames);
			}
			else if(pickType == PICK_LAST)
			{
				pickItems = new StsPickItem[1];
				pickItem = pickItems[0] = new StsPickItem();
				int noNames = 0;
				float zMin = 0.0f, zMax = 0.0f;
				int[] names = null;
				for(i = 0; i < hits; i++)
				{
					noNames = intBuffer.get();
					if(noNames == 0) continue;
					zMin = computePickZ(intBuffer.get());
					zMax = computePickZ(intBuffer.get());
					names = new int[noNames];
					intBuffer.get(names);
				}
				pickItem.noNames = noNames;
				pickItem.zMin = zMin;
				pickItem.zMax = zMax;
				pickItem.names = names;
				if(debug) System.out.println("Last Pick: " + " noNames: " + pickItem.noNames);
			}
           return true;
         }
         catch(Exception e)
         {
            StsException.outputException("StsJOGLPick.processHits() failed.",
                e, StsException.WARNING);
            return false;
         }
    }

    static public float computePickZ(int intZ)
    {
		// ~ is the bitwise inversion operator
		float z1 = (float)( (float)~intZ / (float)0x7fffffffl );
		float z2 = (float)( (float)intZ / (float)0x7fffffffl );
		float z = (1.0f + z1)/2;
		return z;
//    	return 1.0f - (float)( (float)~intZ / (float)0x7fffffffl );
    }

	public float getScreenZ()
	{
		if(pickItems == null || pickItems.length <= 0) return StsParameters.largeFloat;
		return pickItems[0].zMin;
	}

    static public boolean mousePicksPoint(StsMousePoint mousePoint, StsPoint point, StsGLPanel3d glPanel3d, int pickSize)
    {
        int dif;
        double[] screenCoor = glPanel3d.getScreenCoordinates(point);

        dif = Math.abs(mousePoint.x - (int)screenCoor[0]);
        if(dif > pickSize) return false;

        dif = Math.abs( (int)(glPanel3d.viewPort[3] - mousePoint.y - screenCoor[1]) );
        if(dif > pickSize) return false;

        return true;
    }
}
