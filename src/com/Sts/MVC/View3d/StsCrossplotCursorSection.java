package com.Sts.MVC.View3d;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.*;

public class StsCrossplotCursorSection extends StsCursor3dTexture
{
    public StsCrossplot crossplot = null;
    transient int nRows, nCols;
    transient boolean isPixelMode;
    transient public int nPlane = -1;
    transient byte[] data = null;
    transient boolean textureChanged = true;

    public StsCrossplotCursorSection()
    {
    }
    
    public StsCrossplotCursorSection(StsModel model, StsCrossplot crossplot, StsCursor3d cursor3d, int dir)
    {
        this.crossplot = crossplot;
        initialize(model, cursor3d, dir);
    }

    public boolean initialize(StsModel model, StsCursor3d cursor3d, int dir)
    {
        super.initialize(model, cursor3d, dir);
        StsCrossplotClass crossplotClass = (StsCrossplotClass)StsModel.getCurrentModel().getStsClass(StsCrossplot.class);
        if(crossplotClass == null) return false;
        isPixelMode = crossplotClass.getIsPixelMode();
        if(crossplot == null) return false;
        crossplot = (StsCrossplot)crossplotClass.getFirst();
        // initializeTextureTiles(glPanel3d);
        return true;
    }


    public boolean setDirCoordinate(float dirCoordinate)
    {
        if(this.dirCoordinate == dirCoordinate) return false;
		if(crossplot == null) return false;
        this.dirCoordinate = dirCoordinate;
        textureChanged();
        if(textureTiles != null) textureTiles.setTilesDirCoordinate(dirCoordinate);
        data = crossplot.readBytePlaneData(dirNo, dirCoordinate);
        return true;
    }

    public boolean isDisplayableObject(Object object) { return (object instanceof StsCrossplot); }

    public boolean isDisplayingObject(Object object) 
    {
        if(crossplot == object) return true;
        if(this.textureTiles != null && object == textureTiles.cropVolume) return true;
        return false;
    }

    public boolean setObject(Object object)
	{
       if(object == crossplot) return false;
       if(isDisplayableObject(object))
       {
           crossplot = (StsCrossplot)object;
           isVisible = true;
           textureChanged();
           return true;
       }
       else if(object == null)
       {
           isVisible = false;
           crossplot = null;
           textureChanged();
           return true;
       }
       return false;
    }

    public Object getObject() { return crossplot; }

    public Class getDisplayableClass() { return StsCrossplot.class; }

    public boolean canDisplayClass(Class c) { return StsCrossplot.class.isAssignableFrom(c); }

    public void displayTexture(StsGLPanel3d glPanel3d, boolean is3d, StsCursorSection cursorSection)
    {
        if(crossplot == null) return;
        if(!isVisible) return;
        if(!crossplot.getIsVisibleOnCursor()) return;

        GL gl = glPanel3d.getGL();

        checkTextureAndGeometryChanges();

	    if(!initializeTextureTiles(glPanel3d)) return;

        if(textureChanged)
        {
//			System.out.println("Deleting texture.");
            deleteTexturesAndDisplayLists(gl);
            deleteDisplayLists(gl);
        }

        if(is3d)
        {
            gl.glDepthFunc(GL.GL_LEQUAL);
            glPanel3d.setViewShift(gl, 1.0f);
        }

        gl.glDisable(GL.GL_LIGHTING);
        gl.glEnable(GL.GL_TEXTURE_2D);
        gl.glEnable(GL.GL_BLEND);
        gl.glTexEnvf( GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE );
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
        gl.glShadeModel(GL.GL_FLAT);
        createPolygonColorList(crossplot, gl);

        if(textureChanged)
        {
//			System.out.println("Texture changed, reading new texture.");
            data = crossplot.readBytePlaneData(dirNo, dirCoordinate);
//            textureChanged = false;
            byte[] subVolumePlane = cursorSection.subVolumePlane;
            if(subVolumePlane != null) data = applySubVolume(data, subVolumePlane);
        }

        if(is3d)
		{
			if (textureChanged)
			{
//				System.out.println("Displaying changed texture 3d.");
				textureTiles.displayTiles(this, gl, isPixelMode, data, nullByte);
				textureChanged = false;
			}
			else
				textureTiles.displayTiles(this, gl, isPixelMode, (byte[])null, nullByte);
		}
        else
        {
            StsViewCursor viewCursor = (StsViewCursor)glPanel3d.getView();
			if(textureChanged)
			{
//				System.out.println("Displaying changed texture 2d.");
				textureTiles.displayTiles2d(this, gl, viewCursor.axesFlipped, isPixelMode, data, nullByte);
				textureChanged = false;
			}
			else
			{
//				System.out.println("Displaying old texture 2d.");
				textureTiles.displayTiles2d(this, gl, viewCursor.axesFlipped, isPixelMode, (byte[])null, nullByte);
			}
		}

        gl.glDisable(GL.GL_TEXTURE_2D);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL.GL_LIGHTING);

//        StsGLDraw.disableTransparentOverlay(gl);
        if(is3d)
        {
            gl.glDepthFunc(GL.GL_LESS);
            glPanel3d.resetViewShift(gl);
        }
    }

    private void checkTextureAndGeometryChanges()
    {

        if(isPixelMode != crossplot.getIsPixelMode())
        {
            isPixelMode = !isPixelMode;
            textureChanged = true;
        }

        if(crossplot.getPolygons().getSize() == 0)
        {
            textureChanged = true;
            return;
        }
    }
    
    private boolean initializeTextureTiles(StsGLPanel3d glPanel3d)
    {
        if(crossplot == null) return false;
		// if(!glPanel3d.initialized) return false;
        StsSeismicVolume[] seismicVolumes = crossplot.getVolumes();
        if(seismicVolumes == null || seismicVolumes[0] == null) return false;
        StsSeismicVolume seismicVolume = seismicVolumes[0];

        nRows = seismicVolume.getNCursorRows(dirNo);
        nCols = seismicVolume.getNCursorCols(dirNo);
		nPlane = seismicVolume.getCursorPlaneIndex(dirNo, dirCoordinate);
		if(nPlane == -1) return false;
		data = crossplot.readBytePlaneData(dirNo, dirCoordinate);

        if(textureTiles == null)
        {
            StsCropVolume subVolume = model.getProject().getCropVolume();
            textureTiles = StsTextureTiles.constructor(model, this, dirNo, seismicVolume, isPixelMode, subVolume);
            if(textureTiles == null) return false;
            nTextureRows = textureTiles.nTotalRows;
            nTextureCols = textureTiles.nTotalCols;
            textureChanged = true;
        }
        else if(!textureTiles.isSameSize(seismicVolume))
        {
            textureTiles.constructTiles(seismicVolume);
            nTextureRows = textureTiles.nTotalRows;
            nTextureCols = textureTiles.nTotalCols;
            textureChanged = true;
        }
        textureTiles.setTilesDirCoordinate(dirCoordinate);
		return true;
     }

    protected void deleteDisplayLists(GL gl)
	{
		/*
		if(!usingDisplayLists) return;
		usingDisplayLists = false;
		if(textureTiles != null)
			textureTiles.deleteDisplayList(gl);
		*/
	}

    public int getDefaultShader() { return StsJOGLShader.NONE; }
    public boolean getUseShader() { return false; }

    private byte[] applySubVolume(byte[] planeData, byte[] subVolumePlane)
    {
        if(planeData == null) return null;
        if(subVolumePlane == null) return planeData;
        if(!crossplot.getCrossplotClass().getDisplayOnSubVolumes()) return planeData;

        if(planeData.length != subVolumePlane.length)
        {
            StsException.systemError("StsSeismicCursorSection.applySubVolume() failed.\n" +
                " planeData.length = " + planeData.length + " subVolumePlane.length = " + subVolumePlane.length);
            return null;
        }
        byte[] subVolumeData = new byte[planeData.length];
        for(int n = 0; n < planeData.length; n++)
        {
            if(subVolumePlane[n] == 0) subVolumeData[n] = 0;
            else                       subVolumeData[n] = planeData[n];
        }
        return subVolumeData;
    }

    public void drawTextureTileSurface(StsTextureTile tile, GL gl, boolean is3d)
    {
		if(is3d)
			tile.drawQuadSurface3d(gl, StsProject.TD_TIME);
		else
			tile.drawQuadSurface2d(gl);
    }
/*
	public void drawTextureTileSurface2d(StsTextureTile tile, GL gl)
	{
		tile.drawQuadSurface2d(gl);
	}
*/
    // this could be put in a displayList, but is so small in terms of the number of colors
    // that it is easy to recompute on the fly
    private void createPolygonColorList(StsCrossplot crossplot, GL gl)
    {
        StsObjectRefList polygons = crossplot.getPolygons();
        int nPolygons = polygons.getSize();
        int nColors = nPolygons + 1;
        int nAllocatedColors = StsMath.nextBaseTwoInt(nColors);
        float[][] arrayRGBA = new float[4][nAllocatedColors];

        float[] rgba = new float[4];
//        System.out.println("Number of polygons=" + nPolygons);
        for(int n = 0; n < nPolygons; n++)
        {
            StsXPolygon polygon = (StsXPolygon)polygons.getElement(n);
            Color color = polygon.getColor();
            if(color == null) continue;
            rgba = color.getRGBComponents(rgba);
            for(int c = 0; c < 4; c++) arrayRGBA[c][n+1] = rgba[c];
        }
        gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_R, nAllocatedColors, arrayRGBA[0], 0);
        gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_G, nAllocatedColors, arrayRGBA[1], 0);
        gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_B, nAllocatedColors, arrayRGBA[2], 0);
        gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_A, nAllocatedColors, arrayRGBA[3], 0);
        gl.glPixelTransferf(GL.GL_MAP_COLOR, 1);
    }

    public void display(StsGLPanel3d glPanel3d, boolean is3d)
    {
    }

    public boolean hasData(int nPlane) { return nPlane == this.nPlane; }

	public boolean textureChanged()
	{
		nPlane = -1;
		textureChanged = true;
        return true;
    }

    public boolean dataChanged()
    {
        data = null;
        textureChanged();
        return true;
    }
/*
    public void addTextureToDeleteList()
    {
        nPlane = -1;
        data = null;
        textureChanged();
    }
*/
    public void cropChanged()
    {
        textureTiles = null;
        textureChanged();
    }

    public void subVolumeChanged()
    {
    }

    public String propertyReadout(StsPoint point)
    {
        StsRotatedGridBoundingBox boundingBox = model.getProject().getRotatedBoundingBox();
        if(boundingBox == null) return "";

        if(dirNo == StsCursor3d.XDIR)
        {
            int row = Math.round(boundingBox.getRowCoor(point.v[1]));
            int col = Math.round(boundingBox.getSliceCoor(point.v[2]));
            return getTypeString(row, col);
        }
        else if(dirNo == StsCursor3d.YDIR)
        {
            int row = Math.round(boundingBox.getColCoor(point.v[0]));
            int col = Math.round(boundingBox.getSliceCoor(point.v[2]));
            return getTypeString(row, col);
        }
        else
        {
            int row = Math.round(boundingBox.getRowCoor(point.v[1]));
            int col = Math.round(boundingBox.getColCoor(point.v[0]));
            return getTypeString(row, col);
        }
    }

    private String getTypeString(int row, int col)
    {
        if(data == null || !isInRange(row, col)) return " ";
        int nType = data[row*nCols + col];
        String typeName = crossplot.getTypeLibrary().getTypeName(nType-1);
        return new String(crossplot.getName() + " type: " + typeName + " ");
    }

    public String getName()
    {
        return "Cursor view[" + dirNo + "] of: " + crossplot.getName();
    }
}