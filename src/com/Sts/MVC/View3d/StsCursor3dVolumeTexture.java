package com.Sts.MVC.View3d;

import com.Sts.Utilities.StsException;
import com.Sts.Utilities.StsParameters;
import com.Sts.Utilities.StsToolkit;

import javax.media.opengl.GL;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: 2/28/11
 * Time: 1:26 PM
 * To change this template use File | Settings | File Templates.
 */

/**
 * Subclassed by IJK-orthgonal volumes displayed as textures on cursor sections.
 * StsSeismicCursorSection (seismic volumes), StsFractureCursorSection (stimulated volumes).
 */
abstract public class StsCursor3dVolumeTexture extends StsCursor3dTexture
{
    transient boolean isPixelMode = true;
    /** indicates whether cursorSection is being drawn in time or depth. */
    transient byte zDomain = StsParameters.TD_NONE;

    abstract protected void checkTextureAndGeometryChanges(StsGLPanel3d glPanel3d, boolean is3d);

    abstract protected void displayTiles3d(StsGLPanel3d glPanel3d, GL gl);

    abstract protected void displayTiles2d(StsGLPanel3d glPanel3d, GL gl);

    abstract protected boolean initializeTextureTiles(StsGLPanel3d glPanel3d);

    abstract protected void computeTextureData();

    /** This is a standard pattern used by all IJK-orthogonal volume cursor texture displayable classes.
     *  Currently subclassed by StsSeismicCursorSection and StsFractureCursorSection.
     *
     * @param glPanel3d
     * @param is3d
     * @param cursorSection
     */
    public void displayTexture(StsGLPanel3d glPanel3d, boolean is3d, StsCursorSection cursorSection)
    {
        GL gl = glPanel3d.getGL();

        try
        {
            // check for any visibility limits
            if(!isVisible()) return;

            if(debug) StsException.systemDebug(this, "displayTexture", "seismicCursorSection " + getName());

            // check if texture or geometry has changed
            checkTextureAndGeometryChanges(glPanel3d, is3d);

            // geometry has changed because: domain changed (time <-> depth).
            // Rebuild displayLists.
            if(geometryChanged)
            {
                deleteDisplayLists(gl);
                geometryChanged = false;
            }
            // texture has changed because 1) cursor has moved; 2) texture on cursor has changed; 3) subVolume changed; pixelMode changed.
            // Run same displayList with new texture.
            if(textureChanged)
            {
                deleteTextures(gl);
            }
            // construct textureTiles initially or rebuild them if the size has changed (sets textureChanged = true)
            if(!initializeTextureTiles(glPanel3d)) return;

            // plane for this direction and coordinate can't be drawn because it is cropped out
            if(textureTiles.isDirCoordinateCropped()) return;

            // texture has changed so recompute the data
            if(textureChanged)
            {
                computeTextureData();
            }

            // build displayLists if they are required and don't exist
            textureTiles.checkBuildDisplayLists(gl, is3d);

            // initialize gl parameters including colorList for texture
            if(!enableGLState(glPanel3d, gl, is3d)) return;

            // display the tiles (displayLists or immediate mode)
            displayTiles(glPanel3d, gl, is3d);
        }
        catch(Exception e)
        {
            StsException.outputException("StsSeismicCursorSection.displayTexture() failed.", e, StsException.WARNING);
        }
        finally
        {
            // restore any GL state parameters to their S2S defaults
            disableGLState(glPanel3d, gl, is3d);
        }
    }

    protected boolean enableGLState(StsGLPanel3d glPanel3d, GL gl, boolean is3d)
    {
        gl.glDisable(GL.GL_LIGHTING);
        gl.glEnable(GL.GL_TEXTURE_2D);
        gl.glEnable(GL.GL_BLEND);
        gl.glTexEnvf(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
        gl.glShadeModel(GL.GL_FLAT);
        return true;
    }

    protected void disableGLState(StsGLPanel3d glPanel3d, GL gl, boolean is3d)
    {
        gl.glDisable(GL.GL_TEXTURE_2D);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL.GL_LIGHTING);
        gl.glShadeModel(GL.GL_SMOOTH);
        disableShader(gl);
    }

    protected void displayTiles(StsGLPanel3d glPanel3d, GL gl, boolean is3d)
    {
        if(is3d)
        {
            if(debug)
            {
                if(textureChanged)
                    StsException.systemDebug(this, "displayTexture", "Displaying changed texture 3d in dirNo " + dirNo + ".");
                else
                    StsException.systemDebug(this, "displayTexture", "Displaying current texture 3d in dirNo " + dirNo + ".");
            }
            displayTiles3d(glPanel3d, gl);
            textureChanged = false;
        }
        else
        {
            if(debug)
            {
                if(textureChanged)
                    if(debug)
                        StsException.systemDebug(this, "displayTexture", "Displaying changed texture 2d in dirNo " + dirNo + ".");
                    else if(debug)
                        StsException.systemDebug(this, "displayTexture", "Displaying current texture 2d in dirNo " + dirNo + ".");
            }
            displayTiles2d(glPanel3d, gl);
            textureChanged = false;
        }
    }
}
