package com.Sts.MVC.ViewWell;

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.nio.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Feb 26, 2009
 * Time: 3:48:56 PM
 * To change this template use File | Settings | File Templates.
 */
abstract public class StsWellTextureView extends StsWellView implements ActionListener, StsTextureSurfaceFace
{
    /** number of traces in Seismic */
    public int nTraces = -1;
    /** number of samples in each trace */
    public int nSamples = -1;
    /** x coordinates (left to right) for texture */
    public double[] xCoordinates;
    /** y coordinates (top down) for texture */
    public double[] yCoordinates;

    /** Texture tiles for variable density display */
    transient StsTextureTiles textureTiles = null;
    /** indicates texture is new or changed */
    transient boolean textureChanged = true;
    transient boolean useDisplayLists = true;
    transient boolean usingDisplayLists = false;
    transient boolean displayWiggles = false;
    transient boolean isPixelMode = false;
    transient boolean displayValues = false;
    transient boolean pixelsSaved = false;
    transient boolean doPixelsSaved = false;
    transient ByteBuffer backingBuffer = null;
	private transient JPopupMenu tp = null;

    static final byte nullByte = StsParameters.nullByte;

    abstract protected void initializeRange();
    abstract public void adjustColorscale();
    abstract protected boolean getIsPixelMode();
    abstract protected void setGLColorList(GL gl);
    abstract protected void clearShader();
    abstract protected byte[] getData();
    abstract protected String getValueLabel(double xCoordinate);
    abstract protected float compute2dValue(GL gl, double x, double y);

    public StsWellTextureView()
	{
	}

    public StsWellTextureView(StsWellViewModel wellViewModel)
    {
        this.wellViewModel = wellViewModel;
        this.well = wellViewModel.well;
    }

    public void savePixels(boolean b)
	{
		glPanel.savePixels(b);
    }

    public void rebuild(int nSubWindow)
    {
        GridBagLayout g = new GridBagLayout();
        GridBagConstraints gbc = g.getConstraints(glPanel);
        gbc.gridx = nSubWindow;
        g.setConstraints(glPanel, gbc);
    }


    public void showPopupMenu(StsMouse mouse)
    {
        tp = new JPopupMenu();
        StsMenuItem removePanel = new StsMenuItem();
        removePanel.setMenuActionListener("Remove Panel", this, "removePanel", null);
        tp.add(removePanel);
        tp.show(glPanel, mouse.getX(), mouse.getY());
    }
	public void cancelPopupMenu(StsMouse mouse)
	{
		   if (tp != null) tp.setVisible(false);
    }
    /** view2d */
    public void doInitialize()
    {

    }

    public void setAxisRanges()
    {

    }

    public StsView copy(StsGLPanel3d in)
    {
        return null;
    }

    public void setDefaultView()
    {

    }

    public void computeProjectionMatrix(GL gl, GLU glu)
    {
        if(axisRanges == null) return;

        //gl.glViewport(Seismic.SeismicView.insets.left, 0, getGLWidth() - Seismic.SeismicView.insets.left, getGLHeight());

        gl.glMatrixMode(GL.GL_PROJECTION);
        gl.glLoadIdentity();
        axisRanges[1][0] = (float)wellViewModel.getWindowMdepthMax();
        axisRanges[1][1] = (float)wellViewModel.getWindowMdepthMin(this);
        glu.gluOrtho2D(axisRanges[0][0], axisRanges[0][1], axisRanges[1][0], axisRanges[1][1]);
        gl.glGetDoublev(GL.GL_PROJECTION_MATRIX, glPanel.projectionMatrix, 0);
        gl.glMatrixMode(GL.GL_MODELVIEW);
    }

	public void display(GLAutoDrawable component)
	{
        if(glPanel.panelViewChanged)
        {
            initializeView();
            glPanel.panelViewChanged = false;
        }
        if(Main.isGLDebug) System.out.println("StsWellSeismicCurtainView.display() called.");
        displayTexture();
	}

    protected void displayTexture()
    {
        if(textureTiles == null)
        {
            if(!initializeTextureTiles()) return;
            textureChanged = true;

        }

        //computeProjectionMatrix(gl, glu); // should call only if view changed
        gl.glDisable(GL.GL_DEPTH_TEST);
        gl.glDisable(GL.GL_LIGHTING);
        gl.glEnable(GL.GL_TEXTURE_2D);
        gl.glEnable(GL.GL_BLEND);
        gl.glTexEnvf(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
        gl.glShadeModel(GL.GL_FLAT);

        if(isPixelMode != getIsPixelMode())
        {
            textureChanged = true;
            isPixelMode = !isPixelMode;
        }

        setGLColorList(gl);

        if(textureChanged)
        {
            textureTiles.deleteTextures(gl);
            textureTiles.displayTiles(this, gl, isPixelMode, getData(), nullByte);
            textureChanged = false;
        }
        else
            textureTiles.displayTiles(this, gl, isPixelMode, (byte[])null, nullByte);

        gl.glDisable(GL.GL_TEXTURE_2D);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL.GL_LIGHTING);

        if (textureTiles.shader != StsJOGLShader.NONE) StsJOGLShader.disableARBShader(gl);
    }

    protected void displayValues(GL gl)
    {
        StsMousePoint mousePoint = glPanel.getMousePoint();
        if(mousePoint == null) return;
        int glMouseX = mousePoint.x;
        int glMouseY = wellViewModel.cursorY;
        int width = getWidth();
        int height = getHeight();
		int locationX = StsMath.minMax(glMouseX, 0, Math.max(0, width - 70));
        int locationY = StsMath.minMax(glMouseY, 2, Math.max(0, height - 10));

        double depth;
        String depthLabel;
        if(wellViewModel.cursorPicked)
        {
            depth = wellViewModel.getCursorMdepth();
            glMouseY = wellViewModel.getCursorY();
            depthLabel = wellViewModel.getMdepthStringFromGLCursorY();
        }
        else
        {
            depth = wellViewModel.getMdepthFromGLMouseY(glMouseY);
            depthLabel = wellViewModel.convertMDepthToZString(depth);
        }

        float fraction = glMouseX / width;
        float xCoordinate = getValueFromPanelXFraction(fraction);
        double yCoordinate = depth;
        String xLabel = getValueLabel(xCoordinate);
        StsColor.BLACK.setGLColor(gl);
        if(glMouseX > width - 150)
        {
            float f = Math.max(0, (width - 150) / width);
            xCoordinate = StsMath.interpolateValue(axisRanges[0][0], axisRanges[0][1], f);
        }
        if(glMouseY > height - 40)
        {
            float f = Math.max(0, (height - 40) / width);
            yCoordinate = StsMath.interpolateValue(axisRanges[1][0], axisRanges[1][1], f);
        }
        //int glMouseY = getHeight() - glMouseY;
        // int glMouseY = glMouseY;
        StsGLDraw.fontHelvetica12(gl, 2, locationY, depthLabel + " " + xLabel);
        float value = compute2dValue(gl, xCoordinate, yCoordinate);
        if(value == StsParameters.nullValue) return;
        String tempString = Float.toString(value);
        int fontHeight = wellViewModel.fontHeight;
        locationY -= 1.5 * fontHeight;
        StsGLDraw.fontHelvetica12(gl, locationX, locationY, tempString);
    }

    public int getValueFromPanelXFraction(double fraction)
    {
        int trace = (int)((nTraces * fraction));
        if (trace < 0) trace = 0;
        if (trace >= nTraces) trace=nTraces-1;
        return trace;
    }

    private boolean initializeTextureTiles()
    {

        float[][] textureRanges;
        textureRanges = totalAxisRanges;
        // if(textureTiles != null) deleteTextureAndSurface(gl);
        textureTiles = StsTextureTiles.constructor(model, this, nTraces, nSamples, true, textureRanges);
        if(textureTiles == null) return false;
        textureChanged = true;
        return true;
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
         int numBytes = getWidth() * getHeight() * 4;
         if (backingBuffer == null || backingBuffer.capacity() < numBytes)
         {
             try {
                 backingBuffer = ByteBuffer.allocateDirect(numBytes);
             } catch (Exception e)
             {
                 return false;
             } catch (OutOfMemoryError e2)
             {
                 return false;
             }
         }
         backingBuffer.rewind();
         gl.glViewport(0,0,getWidth(),getHeight());
         gl.glPixelStorei(GL.GL_PACK_ALIGNMENT, 1);
         gl.glReadPixels(0,0,getWidth(),getHeight(),GL.GL_BGRA,GL.GL_UNSIGNED_BYTE,backingBuffer);
         return true;
     }

    public boolean restoreBacking()
     {
         int numBytes = getWidth() * getHeight() * 4;
         if (backingBuffer == null || backingBuffer.capacity() < numBytes)
             return false;

         backingBuffer.rewind();
         gl.glViewport(0,0,getWidth(),getHeight());
         gl.glRasterPos2i(0,0);
         gl.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1);
         gl.glDrawPixels(getWidth(),getHeight(),GL.GL_BGRA,GL.GL_UNSIGNED_BYTE,backingBuffer);
         return true;
     }

    public void deleteBacking()
     {
         if (backingBuffer != null)
             backingBuffer = null;
     }

    public void setInsets(boolean axisOn)
	{
		if(axisOn)
		{
			int leftInset = halfWidth + majorTickLength + 2 * verticalFontDimension.width + 2 * fontLineSpace;
			;
			int bottomInset = 0;
			int topInset = halfWidth + majorTickLength + 4 * horizontalFontDimension.height + 3 * fontLineSpace;
			int rightInset = 0;
			insets = new Insets(topInset, leftInset, bottomInset, rightInset);
		}
		else
		{
			insets = new Insets(0, 0, 0, 0);
		}
	}

    public void computeProjectionMatrix()
	{

	}

    public void reshape(GLAutoDrawable component, int x, int y, int width, int height)
	{

	}

    protected StsColor getGridColor()
	{
		return StsColor.BLACK;
	}
    /*
    public void setDisplayPanelSize(Dimension size)
	{
		glPanel.setPanelSize(size);
	}
    */
    public void viewChangedRepaint()
	{
		glPanel.repaint();
	}

    public String getViewClassname()
	{
		return StsTimeLogCurveWellView.viewCurtain;
	}

    public void resetToOrigin() {}

    public void actionPerformed(ActionEvent actionEvent)
    {
        textureChanged = true;
        repaint();
    }

    public void drawTextureTileSurface(StsTextureTile tile, GL gl, boolean is3d)
    {
        tile.drawQuadStripSurface2d(gl, xCoordinates, yCoordinates);
    }

    public void addTextureToDeleteList(StsGLPanel glPanel)
    {}

    public void deleteTextureAndSurface(GL gl)
    {}


    public int getWindowX(float xCoordinate)
    {
        float xMin = axisRanges[0][0];
        float xMax = axisRanges[0][1];
        return (int)(getWidth() * (xCoordinate - xMin) / (xMax - xMin));
    }

    public float getXCoordinate(int mouseX)
    {
        float f = (float)mouseX / getWidth();
        return StsMath.interpolateValue(axisRanges[0][0], axisRanges[0][1], f);
    }

    public boolean textureChanged()
    {
        textureChanged = true;
        return true;
    }

    public boolean dataChanged()
    {
        return textureChanged();
    }

    public void geometryChanged()
    {
    }

    /** Called to actually delete the displayables on the delete list. */
    public void deleteTexturesAndDisplayLists(GL gl)
    {
        if (textureTiles == null) return;
        if (debug) StsException.systemDebug(this, "deleteTextureTileSurface");
        textureTiles.deleteTextures(gl);
        textureTiles.deleteDisplayLists(gl);
        textureChanged = true;
    }
}
