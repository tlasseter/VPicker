package com.Sts.Interfaces;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version beta 1.0
 */

import com.Sts.MVC.View3d.*;

import javax.media.opengl.*;

public interface StsTextureSurfaceFace
{
    /** draw the texture on this tiled surface. */
    public void drawTextureTileSurface(StsTextureTile tile, GL gl, boolean is3d);
    /** data has not changed, but texture colors or range has.  Update textures. */
    public boolean textureChanged();
    /** data has changed. Update data and textures */
    public boolean dataChanged();
    /** geometry data and texture are isVisible on has changed. Update the display lists and rerender. */
    public void geometryChanged();
    /** Get the default OpenGL ARB shader */
    public int getDefaultShader();
    /** Indicate whether OpenGL shader is to be used. */
    public boolean getUseShader();
    /** Clean up.  Delete and all textures and display lists. */
	public void deleteTexturesAndDisplayLists(GL gl);
    /** get the Class whose data is being isVisible here. */
    public Class getDisplayableClass();
    /** get name.  Used for debug identification. */
    public String getName();
}
