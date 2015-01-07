package com.Sts.DBTypes;

import com.Sts.MVC.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.*;
import java.awt.event.*;
import java.nio.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */
public class StsColorList implements ItemListener
{
	/** color list number */
	public int colorListNum = 0;
	/** if colorscale has been edited, rebuild colorList displayList */
	private boolean colorListChanged = true;
	/** colorscale for this colorList */
	public StsColorscale colorscale;
	StsModel model = null;
	StsObject displayObject = null;
	StsItemListeners itemListeners = null;

	static final boolean debug = false;


	public StsColorList(StsColorscale colorscale)
	{
		this.colorscale = colorscale;
	}

	public StsColorList(StsColorscale colorscale, ItemListener listener)
	{
		this.colorscale = colorscale;
		addItemListener(listener);
		colorscale.addItemListener(this);
	}
/*
	public StsColorList(StsColorscale colorscale)
	{
		this.colorscale = colorscale;
		colorscale.addItemListener(this);
	}
*/
	public synchronized void addItemListener(ItemListener listener)
	{
		if(itemListeners == null) itemListeners = new StsItemListeners();
		itemListeners.add(listener);
	}

	public synchronized void removeItemListener(ItemListener listener)
	{
		if(itemListeners == null) return;
		itemListeners.remove(listener);
	}

	protected void fireItemStateChanged(ItemEvent e)
	{
		if (itemListeners == null) return;
		itemListeners.fireItemStateChanged(e);
	}

	public void setColorscale(StsColorscale colorscale)
	{
		this.colorscale = colorscale;
		setColorListChanged(true);
	}

	public void setColorListChanged(boolean changed)
	{
		colorListChanged = changed;
        if(changed) colorscale.colormapBuffer = null;
    }
/*
    public boolean setGLColorList(GL gl, boolean nullsFilled, int shader)
	{
		if(shader != StsJOGLShader.NONE)
		{
			if(debug) System.out.println("Shader used and available.  velocity display");
			// TODO: we need to restore colorListChanged as we are constantly rebuilding the arrayRGBA and JOGLShader
//			if(colorListChanged)
			{
				float[][] arrayRGBA = colorscale.computeRGBAArray(nullsFilled);
				if(arrayRGBA == null)return false;
				StsJOGLShader.setARBColormap(gl, arrayRGBA);
//				colorListChanged = false;
			}
			StsJOGLShader.enableARBShader(gl, shader);
			return true;
		}
		else if(colorListChanged || colorListNum == 0)
		{
		    if(colorListNum > 0)
				gl.glDeleteLists(colorListNum, 1);
			colorListNum = gl.glGenLists(1);
			if(colorListNum == 0)
			{
				StsMessageFiles.logMessage("System error in StsColorList for colorscale " + colorscale.getName() + ": Failed to allocate a display list");
				return false;
			}

			gl.glNewList(colorListNum, GL.GL_COMPILE);
			createColorList(gl, nullsFilled);
			gl.glEndList();
			colorListChanged = false;
			if(debug)System.out.println("Shader not used or not available.  velocity display");
		}
		gl.glCallList(colorListNum);
		return true;
	}

	private void createColorList(GL gl, boolean nullsFilled)
	{
		if(colorscale == null)
		{
			StsException.systemError("Can't create colorList for line2d " + colorscale.getName() + ". Colorscale is null.");
			return;
		}

		float[][] arrayRGBA = colorscale.computeRGBAArray(nullsFilled);
		if(arrayRGBA == null)
		{
			return;
		}
		int nColors = arrayRGBA[0].length;
		gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_R, nColors, arrayRGBA[0], 0);
		gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_G, nColors, arrayRGBA[1], 0);
		gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_B, nColors, arrayRGBA[2], 0);
		gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_A, nColors, arrayRGBA[3], 0);
		gl.glPixelTransferf(GL.GL_MAP_COLOR, 1);
//        if(mainDebug) System.out.println("Color 0: " + arrayRGBA[0][0] + " "  + arrayRGBA[1][0] + " "+ arrayRGBA[2][0] + " "+ arrayRGBA[3][0]);
		arrayRGBA = null;
	}
*/
	public boolean setGLColorList(GL gl, boolean nullsFilled, int shader)
	{
        return setGLColorList(gl, nullsFilled, shader, null);
    }

    public boolean setGLColorList(GL gl, boolean nullsFilled, int shader, Color nullColor)
	{
        if(colorscale.nullColorChanged(nullsFilled, nullColor)) colorListChanged = true;
        // colorListChanged = colorListChanged | (colorscale.nullColorChanged(nullsFilled, nullColor));
        // if colorListChanged, remove colormap and displaylist
        if(colorListChanged) clearColorList(gl);

        // if we are using shader, checkCreate colormap
        if(shader != StsJOGLShader.NONE)
		{
            FloatBuffer colormapBuffer = colorscale.getComputeColormapBuffer(nullsFilled, nullColor);
            if(colormapBuffer == null)return false;
		    StsJOGLShader.loadEnableARBColormap(gl, colormapBuffer, shader);
			return true;
		}
        // otherwise we are using displayList; rebuild if needed
        else
        {
            if(colorListNum == 0)
            {
                colorListNum = gl.glGenLists(1);
                if(colorListNum == 0)
                {
                    StsMessageFiles.logMessage("System error in StsColorList for colorscale " + colorscale.getName() + ": Failed to allocate a display list");
                    return false;
                }

			    gl.glNewList(colorListNum, GL.GL_COMPILE);
			    createColorList(gl, nullsFilled, nullColor);
			    gl.glEndList();
		    }
            gl.glCallList(colorListNum);
		    return true;
        }
	}

    private void clearColorList(GL gl)
    {
        if(colorListNum > 0)
        {
            gl.glDeleteLists(colorListNum, 1);
            colorListNum = 0;
        }
        colorListChanged = false;
    }


    private void createColorList(GL gl, boolean nullsFilled, Color nullColor)
	{
		if(colorscale == null)
		{
			StsException.systemError("Can't create colorList for line2d " + colorscale.getName() + ". Colorscale is null.");
			return;
		}

		float[][] arrayRGBA = colorscale.getComputeRGBAArray(nullsFilled, nullColor);
		if(arrayRGBA == null)
		{
			return;
		}
		int nColors = arrayRGBA[0].length;
		gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_R, nColors, arrayRGBA[0], 0);
		gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_G, nColors, arrayRGBA[1], 0);
		gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_B, nColors, arrayRGBA[2], 0);
		gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_A, nColors, arrayRGBA[3], 0);
		gl.glPixelTransferf(GL.GL_MAP_COLOR, 1);
//        if(mainDebug) System.out.println("Color 0: " + arrayRGBA[0][0] + " "  + arrayRGBA[1][0] + " "+ arrayRGBA[2][0] + " "+ arrayRGBA[3][0]);
		arrayRGBA = null;
	}

	public void itemStateChanged(ItemEvent e)
	{
		setColorListChanged(true);
		fireItemStateChanged(e);
	}

    public boolean getColorListChanged()
    {
        return colorListChanged;
    }
}
