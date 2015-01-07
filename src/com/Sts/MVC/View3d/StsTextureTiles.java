package com.Sts.MVC.View3d;

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.nio.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 *
 * @author TJLasseter
 * @version beta 1.0
 */

public class StsTextureTiles
{
    StsModel model;
    /** bounding box which holds this plane of tiles */
    public StsRotatedGridBoundingBox boundingBox = null;
    /** bounding box which trims this plane of tiles */
    public StsCropVolume cropVolume = null;
    /** number of rows covering all tiles */
    public int nTotalRows;
    /** number of cols covering all tiles */
    public int nTotalCols;
    /** inner class tiles */
    transient private StsTextureTile[] tiles;
    /** coordinate direction of plane */
    public int dir = -1;
    /** coordinate value defining plane position */
    public float dirCoordinate = nullValue;
    /** number of background rows (power of 2 > nTextureRows) */
    private int nBackgroundRows;
    /** number of background cols ( power of 2 > nTextureCols) */
    private int nBackgroundCols;
    /** display in pixel mode or not */
    public boolean isPixelMode;
	/* shader contouring */
	public int shader = StsJOGLShader.NONE;
    /** use shader if available and required */
    public boolean useShader = false;
    private boolean cropChanged;
    /** current ranges being used for this set of tiles */
    public float[][] axisRanges = null;
    /** textured surface these textures are on */
    StsTextureSurfaceFace textureSurface;
    /** Display lists should be used (controlled by View:Display Options) */
    transient boolean useDisplayLists;
    /** Display lists are currently being used for surface geometry */
    transient boolean usingDisplayLists = false;
    /** textures have been deleted or don't exist yet */
//    transient public boolean textureDeleted = true;
    /** indicates texture has been changed and needs to be rebuilt */
//    transient boolean textureChanged = true;
    /** flag to indicate textures should be deleted (if they exist) and created */
    //    private boolean deleteTextures = true;
    /** flag to indicate surface displayList (geometry) should be deleted */
//    private boolean deleteDisplayList = true;

    static private boolean runTimer = false;
    static private StsTimer timer = null;
    static public final boolean debug = StsGLPanel.textureDebug;

    static final float nullValue = StsParameters.nullValue;

    private StsTextureTiles(StsModel model, StsTextureSurfaceFace textureSurface, int dirNo, StsRotatedGridBoundingBox boundingBox, boolean isPixelMode, StsCropVolume cropVolume)
    {
        this.model = model;
        this.textureSurface = textureSurface;
        this.dir = dirNo;
        this.isPixelMode = isPixelMode;
        this.cropVolume = cropVolume;
        if (runTimer && timer == null) timer = new StsTimer();
        constructTiles(boundingBox);
        initShader();
        StsTextureList.addTextureToList(textureSurface);
        if (debug) System.out.println("StsTextureTiles.constructor() for " + textureSurface.toString() + " called on thread " + Thread.currentThread().getName() + " dirNo " + dir);
    }

    private StsTextureTiles(StsModel model, StsTextureSurfaceFace textureSurface, int nRows, int nCols, boolean isPixelMode)
    {
        this.model = model;
        this.textureSurface = textureSurface;
        this.isPixelMode = isPixelMode;
        if (runTimer && timer == null) timer = new StsTimer();
        constructTiles(nRows, nCols, StsGLPanel.maxTextureSize);
        initShader();
        StsTextureList.addTextureToList(textureSurface);
        if (debug) System.out.println("StsTextureTiles.constructor(nTextureRows, nTextureCols, ...) for " + textureSurface.toString() + " called on thread " + Thread.currentThread().getName());
    }

    private StsTextureTiles(StsModel model, StsTextureSurfaceFace textureSurface, int nRows, int nCols, boolean isPixelMode, float[][] axisRanges)
    {
        this.model = model;
        this.textureSurface = textureSurface;
        this.isPixelMode = isPixelMode;
        dir = 1;
        if (runTimer && timer == null) timer = new StsTimer();
        constructTiles(nRows, nCols, StsGLPanel.maxTextureSize, axisRanges);
        initShader();
//        cropChanged();
        if (debug) System.out.println("StsTextureTiles.constructor(nTextureRows, nTextureCols, ...) for " + textureSurface.toString() + " called on thread " + Thread.currentThread().getName());
    }

    public void constructTiles(StsRotatedGridBoundingBox boundingBox)
    {
        this.boundingBox = boundingBox;
        setGridSize(dir);
        constructTiles(StsGLPanel.maxTextureSize);
    }

    public void constructTiles(int nRows, int nCols, int maxTextureSize)
    {
        nTotalRows = nRows;
        nTotalCols = nCols;
	    shader = textureSurface.getDefaultShader();
        createTiles(maxTextureSize);
        cropChanged();
    }

    public void constructTiles(int maxTextureSize)
    {
	    shader = textureSurface.getDefaultShader();
        createTiles(maxTextureSize);
        cropChanged();
    }

    public void constructTiles(int maxTextureSize, float[][] axisRanges)
    {
        constructTiles(nTotalRows, nTotalCols, maxTextureSize, axisRanges);
        cropChanged();
    }


    private void constructTiles(int nTotalRows, int nTotalCols, int maxTextureSize, float[][] axisRanges)
    {
        createTiles(nTotalRows, nTotalCols, maxTextureSize, axisRanges);
        cropChanged();
    }

    static public StsTextureTiles constructor(StsModel model, StsTextureSurfaceFace textureSurface, int dirNo, StsRotatedGridBoundingBox boundingBox, boolean isPixelMode, StsCropVolume cropVolume)
    {
        // if (!glPanel.initialized) return null;  // glPanel not initialized yet
        try
        {
            return new StsTextureTiles(model, textureSurface, dirNo, boundingBox, isPixelMode, cropVolume);
        }
        catch (Exception e)
        {
            StsException.systemError("StsTextureTiles.constructor() failed.");
            return null;
        }
    }

    static public StsTextureTiles constructor(StsModel model, StsTextureSurfaceFace textureSurface, int nRows, int nCols, boolean isPixelMode)
    {
        // if (!glPanel.initialized) return null;  // glPanel not initialized yet
        try
        {
            return new StsTextureTiles(model, textureSurface, nRows, nCols, isPixelMode);
        }
        catch (Exception e)
        {
            StsException.systemError("StsTextureTiles.constructor() failed.");
            return null;
        }
    }

    static public StsTextureTiles constructor(StsModel model, StsTextureSurfaceFace textureSurface, int nRows, int nCols, boolean isPixelMode, float[][] axisRanges)
    {
        // if (!glPanel.initialized) return null;  // glPanel not initialized yet
        try
        {
            return new StsTextureTiles(model, textureSurface, nRows, nCols, isPixelMode, axisRanges);
        }
        catch (Exception e)
        {
            StsException.systemError("StsTextureTiles.constructor() failed.");
            return null;
        }
    }

    private void setGridSize(int dir)
    {
        switch (dir)
        {
            case StsCursor3d.XDIR:
                nTotalRows = boundingBox.nRows;
                nTotalCols = boundingBox.nSlices;
                break;
            case StsCursor3d.YDIR:
                nTotalRows = boundingBox.nCols;
                nTotalCols = boundingBox.nSlices;
                break;
            case StsCursor3d.ZDIR:
                nTotalRows = boundingBox.nRows;
                nTotalCols = boundingBox.nCols;
                break;
        }
    }

   private void initShader()
	{
		boolean useShader = textureSurface.getUseShader() && StsJOGLShader.canUseShader;
        setShader(useShader);
    }

    private void setShader(boolean useShader)
    {
        if(useShader)
			shader = textureSurface.getDefaultShader();
		else
			shader = StsJOGLShader.NONE;
//        textureSurface.textureChanged();
    }

	public boolean shaderChanged()
	{
		boolean usingShader = shader != StsJOGLShader.NONE;
		boolean useShader = textureSurface.getUseShader() && StsJOGLShader.canUseShader;
		if(useShader == usingShader)return false;
        setShader(useShader);
        return true;
    }

    public boolean usingShader() { return shader != StsJOGLShader.NONE; }

    private void createTiles(int maxTextureSize)
    {
        int nTileRows; // number of rows of tiles
        int nTileCols; // number of columns of tiles
        int nRowsRemaining; // number of rows in last row of tiles
        int nColsRemaining; // number of columns in last column of tiles
        int nBackgroundRowsRemaining; // number of base rows in last row of tiles ( power of 2 > nRowsRemaining
        int nBackgroundColsRemaining; // number of base cols in last col of tiles ( power of 2 > nColsRemaining

        if (maxTextureSize <= 0) return;
        if (runTimer) timer.start();
        maxTextureSize /= 2;  // divide by 2 to be on the safe side...
        if (maxTextureSize <= 0)
        {
            StsException.systemError("StsTextureTiles.createTiles() failed. Max texture size returned as: " + maxTextureSize);
            return;
        }

        nBackgroundRows = StsMath.nextBaseTwoInt(nTotalRows);
        nBackgroundCols = StsMath.nextBaseTwoInt(nTotalCols);

        if (nBackgroundRows <= maxTextureSize)
        {
            nTileRows = 1;
            nRowsRemaining = nTotalRows;
            nBackgroundRowsRemaining = maxTextureSize;
        }
        else
        {
            nTileRows = StsMath.ceiling(((float) (nTotalRows - 1)) / (maxTextureSize - 1));
            nRowsRemaining = (nTotalRows - 1) % (maxTextureSize - 1) + 1;
            nBackgroundRowsRemaining = StsMath.nextBaseTwoInt(nRowsRemaining);
        }
        if (nBackgroundCols <= maxTextureSize)
        {
            nTileCols = 1;
            nColsRemaining = nTotalCols;
            nBackgroundColsRemaining = maxTextureSize;
        }
        else
        {
            nTileCols = StsMath.ceiling(((float) (nTotalCols - 1)) / (maxTextureSize - 1));
            nColsRemaining = (nTotalCols - 1) % (maxTextureSize - 1) + 1;
            nBackgroundColsRemaining = StsMath.nextBaseTwoInt(nColsRemaining);
        }

        tiles = new StsTextureTile[nTileRows * nTileCols];

        int nTile = 0;
        int row = 0;
        for (int nTileRow = 0; nTileRow < nTileRows - 1; nTileRow++)
        {
            int col = 0;
            for (int nTileCol = 0; nTileCol < nTileCols - 1; nTileCol++, nTile++)
            {
                tiles[nTile] = new StsTextureTile(nTotalRows, nTotalCols, row, col, this, maxTextureSize, nTile);
                col += (maxTextureSize - 1);
            }
            tiles[nTile++] = new StsTextureTile(nTotalRows, nTotalCols, row, col, this, maxTextureSize, nTile);
            row += (maxTextureSize - 1);
        }
        int col = 0;
        for (int nTileCol = 0; nTileCol < nTileCols - 1; nTileCol++, nTile++)
        {
            tiles[nTile] = new StsTextureTile(nTotalRows, nTotalCols, row, col, this, maxTextureSize, nTile);
            col += (maxTextureSize - 1);
        }
        tiles[nTile++] = new StsTextureTile(nTotalRows, nTotalCols, row, col, this, maxTextureSize, nTile);

        if (runTimer) timer.stopPrint("Construct tiles for direction " + dir + ".");
        //       deleteTextures = true;
    }

    private void createTiles(int nRows, int nCols, int maxTextureSize, float[][] axisRanges)
    {
        int nTileRows; // number of rows of tiles
        int nTileCols; // number of columns of tiles
        int nRowsRemaining; // number of rows in last row of tiles
        int nColsRemaining; // number of columns in last column of tiles
        int nBackgroundRowsRemaining; // number of base rows in last row of tiles ( power of 2 > nRowsRemaining
        int nBackgroundColsRemaining; // number of base cols in last col of tiles ( power of 2 > nColsRemaining

        this.nTotalRows = nRows;
        this.nTotalCols = nCols;
        this.axisRanges = axisRanges;

        maxTextureSize /= 2;  // divide by 2 to be on the safe side...
        if (maxTextureSize <= 0) maxTextureSize = 512;
        if (runTimer) timer.start();
        if (maxTextureSize <= 0)
        {
            StsException.systemError("StsTextureTiles.createTiles() failed. Max texture size returned as: " + maxTextureSize);
            return;
        }

        nBackgroundRows = StsMath.nextBaseTwoInt(nTotalRows);
        nBackgroundCols = StsMath.nextBaseTwoInt(nTotalCols);

        if (nBackgroundRows <= maxTextureSize)
        {
            nTileRows = 1;
            nRowsRemaining = nTotalRows;
            nBackgroundRowsRemaining = maxTextureSize;
        }
        else
        {
            nTileRows = StsMath.ceiling(((float) (nTotalRows - 1)) / (maxTextureSize - 1));
            nRowsRemaining = (nTotalRows - 1) % (maxTextureSize - 1) + 1;
            nBackgroundRowsRemaining = StsMath.nextBaseTwoInt(nRowsRemaining);
        }
        if (nBackgroundCols <= maxTextureSize)
        {
            nTileCols = 1;
            nColsRemaining = nTotalCols;
            nBackgroundColsRemaining = maxTextureSize;
        }
        else
        {
            nTileCols = StsMath.ceiling(((float) (nTotalCols - 1)) / (maxTextureSize - 1));
            nColsRemaining = (nTotalCols - 1) % (maxTextureSize - 1) + 1;
            nBackgroundColsRemaining = StsMath.nextBaseTwoInt(nColsRemaining);
        }

        tiles = new StsTextureTile[nTileRows * nTileCols];

        int nTile = 0;
        int row = 0;
        for (int nTileRow = 0; nTileRow < nTileRows - 1; nTileRow++)
        {
            int col = 0;
            for (int nTileCol = 0; nTileCol < nTileCols - 1; nTileCol++, nTile++)
            {
                tiles[nTile] = new StsTextureTile(nTotalRows, nTotalCols, row, col, this, maxTextureSize, axisRanges, nTile);
                col += (maxTextureSize - 1);
            }
            tiles[nTile++] = new StsTextureTile(nTotalRows, nTotalCols, row, col, this, maxTextureSize, axisRanges, nTile);
            row += (maxTextureSize - 1);
        }
        int col = 0;
        for (int nTileCol = 0; nTileCol < nTileCols - 1; nTileCol++)
        {
            tiles[nTile++] = new StsTextureTile(nTotalRows, nTotalCols, row, col, this, maxTextureSize, axisRanges, nTile);
            col += (maxTextureSize - 1);
        }
        tiles[nTile++] = new StsTextureTile(nTotalRows, nTotalCols, row, col, this, maxTextureSize, axisRanges, nTile);
        if (runTimer) timer.stopPrint("Construct tiles for direction " + dir + ".");
        //       deleteTextures = true;
    }

    public void displayTiles(StsTextureSurfaceFace textureSurface, GL gl, boolean isPixelMode, byte[] data, byte nullByte)
    {
        if (tiles == null) return;
        checkPixelMode(isPixelMode);
        if (runTimer) timer.start();
        for (int n = 0; n < tiles.length; n++)
            tiles[n].display(textureSurface, data, isPixelMode, gl, n, nullByte);
        if (runTimer) timer.stopPrint("    Display tiles.");
    }

    public void displayTiles(StsTextureSurfaceFace textureSurface, GL gl, boolean isPixelMode, float[] data, byte nullByte)
    {
        displayTiles(textureSurface, gl, isPixelMode, data, false, nullByte);
    }

    public void displayTiles(StsTextureSurfaceFace textureSurface, GL gl, boolean isPixelMode, float[] data, boolean reverse, byte nullByte)
    {
        if (tiles == null) return;
        checkPixelMode(isPixelMode);
        if (runTimer) timer.start();
        for (int n = 0; n < tiles.length; n++)
            tiles[n].display(textureSurface, data, isPixelMode, gl, reverse, n, nullByte);
        if (runTimer) timer.stopPrint("    Display tiles.");
    }

    public void displayTiles(StsTextureSurfaceFace textureSurface, GL gl, boolean isPixelMode, ByteBuffer data, byte nullByte)
    {
        if (tiles == null) return;
        checkPixelMode(isPixelMode);
        if (runTimer) timer.start();
        for (int n = 0; n < tiles.length; n++)
            tiles[n].display(textureSurface, data, isPixelMode, gl, nullByte);
        if (runTimer) timer.stopPrint("    Display tiles.");
    }

    public void addData(GL gl, byte[] data)
    {
        for (int n = 0; n < tiles.length; n++)
            tiles[n].addData(gl, data);
    }

    public void setAxesFlipped(boolean axesFlipped)
    {
        for (int n = 0; n < tiles.length; n++)
            tiles[n].setAxesFlipped(axesFlipped);

    }

    public void displayTiles2d(StsTextureSurfaceFace textureSurface, GL gl, boolean axesFlipped, boolean isPixelMode, byte[] planeData, byte nullByte)
    {
        if (tiles == null) return;
        checkPixelMode(isPixelMode);
        if (runTimer) timer.start();
        for (int n = 0; n < tiles.length; n++)
            tiles[n].display2d(textureSurface, isPixelMode, planeData, gl, axesFlipped, n, nullByte);
        if (runTimer) timer.stopPrint("    Display 2d tiles.");
    }

    public void displayTiles2d(StsTextureSurfaceFace textureSurface, GL gl, boolean axesFlipped, boolean isPixelMode, ByteBuffer planeData, byte nullByte)
    {
        if (tiles == null) return;
        checkPixelMode(isPixelMode);
        if (runTimer) timer.start();
        for (int n = 0; n < tiles.length; n++)
            tiles[n].display2d(textureSurface, isPixelMode, planeData, gl, axesFlipped, n, nullByte);
        if (runTimer) timer.stopPrint("    Display 2d tiles.");
    }

    private void checkPixelMode(boolean isPixelMode)
    {
        if (isPixelMode == this.isPixelMode) return;
        this.isPixelMode = isPixelMode;
        //         deleteTextures = true;
    }

    public boolean isDirCoordinateCropped()
    {
        if (cropVolume == null) return false;
        if (!cropVolume.getApplyCrop()) return false;
        if (dirCoordinate == nullValue) return false;

        switch (dir)
        {
            case StsCursor3d.XDIR:
                return dirCoordinate < cropVolume.xMin ||
                    dirCoordinate > cropVolume.xMax;
            case StsCursor3d.YDIR:
                return dirCoordinate < cropVolume.yMin ||
                    dirCoordinate > cropVolume.yMax;
            case StsCursor3d.ZDIR:
                return dirCoordinate < cropVolume.zMin ||
                    dirCoordinate > cropVolume.zMax;
            default:
                return false;
        }
    }

    public boolean deleteTextures(GL gl)
    {
//        if(textureDeleted) return;

        if (tiles == null || tiles[0] == null || tiles[0].texture == 0) return false;
        for (int n = 0; n < tiles.length; n++)
            tiles[n].deleteTexture(gl, n);

        if(debug) System.out.println(textureSurface.getName() + "deleting tiles 0 thru " + ( tiles.length-1 ) + " cursorSection[ " + dir + "] on thread " + Thread.currentThread().getName());

//        textureDeleted = true;
        textureSurface.textureChanged();
        return true;
    }

    public void deleteDisplayLists(GL gl)
    {
        if(!usingDisplayLists) return;
        usingDisplayLists = false;
        boolean displayListsDeleted = false;
        if(tiles == null) return;
        for (int n = 0; n < tiles.length; n++)
            if(tiles[n].deleteDisplayList(gl, n)) displayListsDeleted = true;

 //       if(debug && displayListsDeleted)
 //           System.out.println("Tile displayLists deleted on thread " + Thread.currentThread().getName() + " dirNo " + dir);
    }

    public boolean checkBuildDisplayLists(GL gl, boolean is3d)
    {
        useDisplayLists = model.useDisplayLists;
        if(useDisplayLists == usingDisplayLists) return false;
        if (!useDisplayLists && usingDisplayLists)
            deleteDisplayLists(gl);
        else if(useDisplayLists && !usingDisplayLists)
            constructSurface(textureSurface, gl, useDisplayLists, is3d);
        usingDisplayLists = useDisplayLists;
        return true;
    }

    public void setTilesDirCoordinate(float dirCoordinate)
    {
        if (this.dirCoordinate == dirCoordinate) return;
        //       deleteTextures = true;
        this.dirCoordinate = dirCoordinate;
        if (tiles == null) return;
        for (int n = 0; n < tiles.length; n++)
            if (tiles[n] != null) tiles[n].setDirCoordinate(dirCoordinate);
    }

    public void cropChanged()
    {
        int cropRowMin, cropRowMax, cropColMin, cropColMax;

        if (cropVolume == null) return;
        boolean applyCrop = cropVolume.getApplyCrop();
        if (applyCrop)
        {
            float cropRowMinCoor = cropVolume.getCursorRowMinCoor(dir);
            cropRowMin = boundingBox.getCursorRow(dir, cropRowMinCoor);

            float cropRowMaxCoor = cropVolume.getCursorRowMaxCoor(dir);
            cropRowMax = boundingBox.getCursorRow(dir, cropRowMaxCoor);

            float cropColMinCoor = cropVolume.getCursorColMinCoor(dir);
            cropColMin = boundingBox.getCursorCol(dir, cropColMinCoor);

            float cropColMaxCoor = cropVolume.getCursorColMaxCoor(dir);
            cropColMax = boundingBox.getCursorCol(dir, cropColMaxCoor);

            //           deleteTextures = true;
        }
        else
        {
            if (boundingBox == null) return;
            cropRowMin = boundingBox.getCursorRowMin(dir);
            cropRowMax = boundingBox.getCursorRowMax(dir);
            cropColMin = boundingBox.getCursorColMin(dir);
            cropColMax = boundingBox.getCursorColMax(dir);
        }
        if (tiles == null) return;
        for (int n = 0; n < tiles.length; n++)
        {
            if (tiles[n] != null)
                tiles[n].adjust(cropRowMin, cropRowMax, cropColMin, cropColMax);
        }
    }

    public boolean isSameSize(StsRotatedGridBoundingBox boundingBox)
    {
        switch (dir)
        {
            case StsCursor3d.XDIR:
                if (nTotalRows != boundingBox.nRows) return false;
                if (nTotalCols != boundingBox.nSlices) return false;
                if (this.boundingBox.yMin != boundingBox.yMin) return false;
                if (this.boundingBox.yMax != boundingBox.yMax) return false;
                if (this.boundingBox.zMin != boundingBox.zMin) return false;
                if (this.boundingBox.zMax != boundingBox.zMax) return false;
                break;
            case StsCursor3d.YDIR:
                if (nTotalRows != boundingBox.nCols) return false;
                if (nTotalCols != boundingBox.nSlices) return false;
                if (this.boundingBox.xMin != boundingBox.xMin) return false;
                if (this.boundingBox.xMax != boundingBox.xMax) return false;
                if (this.boundingBox.zMin != boundingBox.zMin) return false;
                if (this.boundingBox.zMax != boundingBox.zMax) return false;
                break;
            case StsCursor3d.ZDIR:
                if (nTotalRows != boundingBox.nRows) return false;
                if (nTotalCols != boundingBox.nCols) return false;
                if (this.boundingBox.xMin != boundingBox.xMin) return false;
                if (this.boundingBox.xMax != boundingBox.xMax) return false;
                if (this.boundingBox.yMin != boundingBox.yMin) return false;
                if (this.boundingBox.yMax != boundingBox.yMax) return false;
                break;
        }
        return true;
    }

    public boolean isSameSize(int nRows, int nCols)
    {
        return nTotalRows == nRows && nTotalCols == nCols;
    }

    public void constructSurface(StsTextureSurfaceFace surface, GL gl, boolean useDisplayLists, boolean is3d)
    {
        if(debug) StsException.systemDebug(this, "constructSurface", "for " + surface.getName());
        if (tiles == null) return;
        this.useDisplayLists = useDisplayLists;
        usingDisplayLists = useDisplayLists;
        for (int n = 0; n < tiles.length; n++)
            if (tiles[n] != null) tiles[n].constructSurface(surface, gl, useDisplayLists, is3d, n);
    }

//    public void setShader(int shader) { this.shader = shader; }

//    public int getEnableShader() { return this.shader; }
/*
	protected void finalize() throws Throwable
	{
		try
		{
		  if(tiles == null) return;
	      for(int n = 0; n < tiles.length; n++)
		    if(tiles[n] != null) System.out.println("tile not null "+n);

		}
		finally
		{
			super.finalize();
		}
    }
*/
}
