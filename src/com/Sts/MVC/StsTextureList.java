package com.Sts.MVC;

import com.Sts.Interfaces.*;

import javax.media.opengl.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jan 26, 2009
 * Time: 1:56:35 PM
 *
 * Keeps a list of all textures being used;  when a texture is to be deleted, it is moved to the delete list
 * and all textures are deleted when possible (graphicsContext == GLContext).
 */
public class StsTextureList
{
    /** textures on deleteList to be deleted at next opportunity */
    static ArrayList deleteTextureList = new ArrayList();
    /** textures currently being drawn; all must be deleted if GL is re-inited */
    static ArrayList textureList = new ArrayList();

    static final boolean textureDebug = false;

    static public void addTextureToList(StsTextureSurfaceFace textureTileSurface)
    {
        if(textureList.contains(textureTileSurface)) return;
        if(textureDebug) System.out.println("Added " + textureTileSurface + " to textureList");
        textureList.add(textureTileSurface);
    }

    static public void addTextureToDeleteList(StsTextureSurfaceFace deleteTextureTileSurface)
    {
        if(deleteTextureList.contains(deleteTextureTileSurface)) return;
        if(textureDebug) System.out.println("Added " + deleteTextureTileSurface + " to deleteTextureList");
        textureList.remove(deleteTextureTileSurface);
        deleteTextureList.add(deleteTextureTileSurface);
    }

    static public void deleteAllTextures(GL gl)
    {
        StsTextureSurfaceFace textureTileSurface;
        Iterator iterator = textureList.iterator();
        while(iterator.hasNext())
        {
            textureTileSurface = (StsTextureSurfaceFace)iterator.next();
            if(textureDebug) System.out.println("Deleting " + textureTileSurface + " from deleteTextureList");
            textureTileSurface.deleteTexturesAndDisplayLists(gl);
        }
        textureList.clear();

        iterator = deleteTextureList.iterator();
        while(iterator.hasNext())
        {
            textureTileSurface = (StsTextureSurfaceFace)iterator.next();
            if(textureDebug) System.out.println("Deleting " + textureTileSurface + " from deleteTextureList");
            textureTileSurface.deleteTexturesAndDisplayLists(gl);
        }
        deleteTextureList.clear();
    }

    static public void checkDeleteTextures(GL gl)
    {
        StsTextureSurfaceFace textureSurface;
        Iterator iterator = deleteTextureList.iterator();
        while(iterator.hasNext())
        {
            textureSurface = (StsTextureSurfaceFace)iterator.next();
            textureSurface.deleteTexturesAndDisplayLists(gl);
        }
        deleteTextureList.clear();
    }

    static public boolean deleteTexture(StsTextureSurfaceFace textureSurface, GL gl)
    {
        if(textureSurface == null) return false;
        if(textureDebug) System.out.println("Delete textureSurface "+ textureSurface.toString());
        textureSurface.deleteTexturesAndDisplayLists(gl);
//        textureList.remove(textureTiles);
        return true;
    }

    static public void clearAllTextures()
    {
        deleteTextureList = textureList;
        textureList = new ArrayList();
    }

    static public void deleteViewClassTextures(Class viewClass)
    {
        StsTextureSurfaceFace textureSurface;
        Iterator iterator = deleteTextureList.iterator();
        while(iterator.hasNext())
        {
            textureSurface = (StsTextureSurfaceFace)iterator.next();
            if(textureSurface.getDisplayableClass() == viewClass)
                textureSurface.textureChanged();
        }
    }
}
