
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.util.*;

/** An StsZoneBlock is the intersection of a zone and a block.  The zoneBlock is bound by a topGrid and a botGrid.
 *  The top and bottom grids are horizons interpreted by the geophysicist generally from seismic data.
 *  Between these grids, the zoneBlock is divided into a set of subZones.  These are defined by the geologist
 *  generally from well data based on log information.  These subZones have a layering style: proportional, onlap, offlap, or top-bot-truncated.
 *  The geologist can use subhorizons for defining the subzone gridding as well.  Subhorizons are optional boundary surfaces between subzones.
 *  Each subzone can then be divided into proportional layers.  The number of layers is defined by the simulation engineer.
 */
public class StsZoneBlock extends StsGridBoundingBox implements StsSelectable
{
    transient protected StsZone zone;
    transient private StsBlock block;
    transient private StsColor stsColor;
    transient private StsList zoneBlockSides;
    transient private ArrayList<StsPropertyVolume> zoneProperties;
    transient protected StsBlockGrid topGrid;
    transient protected StsBlockGrid botGrid;
    /** convenience copy; value defined by StsZone.nSubZones */
    transient protected int nSubZones;
    transient int cursorRow = -1, cursorCol = -1, cursorSlice = -1;
    transient float currentCursorXCoor = StsParameters.nullValue;
    transient float currentCursorYCoor = StsParameters.nullValue;
    transient float currentCursorZCoor = StsParameters.nullValue;
    transient StsList cursorRowGridLines, cursorColGridLines, cursorSliceGridLines;
    transient StsList cursorRowEdgeLoops, cursorColEdgeLoops, cursorSliceEdgeLoops;

    transient public int simulationRowMin;

    transient float[][][] subZoneGridZ;
    transient float[][][][] layerGridZs;

    transient float[][][] subZoneVerticalGridZs;

    /** cell types are either full, empty, or truncated. */
    transient byte[][][] cellTypes;
    transient int currentNSlice = -1;
    transient float[][][] sliceNormals;

	transient private int sliceSurfDisplayListNum = 0; // display list number (> 0)
	transient private int sliceGridDisplayListNum = 0; // display list number (> 0)

    static public int nextColorIndex = 0;

    static final int ROW = StsParameters.ROW;
    static final int COL = StsParameters.COL;
    static final int PLUS = StsParameters.PLUS;
    static final int MINUS = StsParameters.MINUS;

    static final float nullValue = StsParameters.nullValue;

    static public final float[] rowSectionNormal = new float[] { 0.0f, -1.0f, 0.0f };
    static public final float[] colSectionNormal = new float[] { -1.0f, 0.0f, 0.0f };

    public StsZoneBlock()
    {
    }

    public StsZoneBlock(StsZone zone, StsBlock block, StsBlockGrid botGrid, StsBlockGrid topGrid)
    {
        super(false);
        this.zone = zone;
        this.block = block;
        this.botGrid = botGrid;
        this.topGrid = topGrid;
        initialize(block);
        nSubZones = zone.getNSubZones();
        constructZoneBlockSides();
    }

    /** Required since this is an StsObject; currently StsZoneBlocks are not written to the database and
     *  are rebuilt when the grid model is constructed.
     */
    public boolean initialize(StsModel model)
    {
        return true;
    }

    public void constructGrids()
    {
        constructGrid();
    }

    // called when a subZone is changed
    public void reconstructSubZoneGrid(StsSubZone subZone)
    {
//        constructZoneGrid();
        if(zoneBlockSides != null) zoneBlockSides.forEach("reconstructSubZoneGrid", subZone);
    }

    public StsZone getZone() { return zone; }
    public StsBlockGrid getBottomGrid() { return this.botGrid; }
    public StsBlockGrid getTopGrid() { return this.topGrid; }
    public Object[] getSides() { return zoneBlockSides.getList(); }

    public StsColor getStsColor()
    {
        if(stsColor != null) return stsColor;
        stsColor = StsColor.colors32[nextColorIndex++%32];
        return stsColor;
    }

    public float getZoneGridZ(int nSubZone, float rowF, float colF)
    {
        return getSubZoneGridZ(nSubZone, 0, rowF, colF);
    }

    public float getSubZoneGridZ(int nZoneLayer, int row, int col)
    {
        int[] nSubZoneAndLayer = zone.getSubZoneAndLayer(nZoneLayer);
        if(nSubZoneAndLayer == null) return nullValue;
        int nSubZone = nSubZoneAndLayer[0];
        int nLayer = nSubZoneAndLayer[1];
        return getSubZoneGridZ(nSubZone, nLayer, row, col);
    }

    public float getSubZoneGridZ(int nSubZone, int nSlice, int row, int col)
    {
        float[][][] subZoneGridZ = getSubZoneGridZs(nSubZone);
        if(subZoneGridZ != null)
        {
            if(!isInsideRowCol(row, col))
            {
                StsException.systemError(this, "getSubZoneGridZ", " row " + row + " col " + col + " is outside range " + getLabel());
                return nullValue;
            }
            return subZoneGridZ[nSlice][row - rowMin][col - colMin];
        }
        else
        {
            if(this.subZoneGridZ == null) return nullValue;

            float topZ = getZoneGridZ(nSubZone, row, col);
            if(topZ == nullValue) return nullValue;
            if(nSlice == 0) return topZ;

            float botZ = getZoneGridZ(nSubZone+1, row, col);
            if(botZ == nullValue) return nullValue;

            StsObjectRefList subZones = zone.getSubZones();
            StsSubZone subZone = (StsSubZone)subZones.getElement(nSubZone);
            if(subZone == null) return nullValue;

            int subZoneType = subZone.getSubZoneType();
            int nLayers = subZone.getNLayers();

            float f, dZ, thickness;
            switch(subZoneType)
            {
                case StsZone.SUBZONE_UNIFORM:
                    f = ((float)nSlice)/nLayers;
                    return topZ + f*(botZ - topZ);
                case StsZone.SUBZONE_OFFLAP:
                    thickness = subZone.getLayerThickness();
                    dZ = (nLayers - nSlice)*thickness;
                    return Math.max(botZ - dZ, topZ);
                case StsZone.SUBZONE_ONLAP:
                    thickness = subZone.getLayerThickness();
                    dZ = nSlice*thickness;
                    return Math.min(topZ + dZ, botZ);
                default:
                    return nullValue;
            }
        }
    }

    // get z value midway between in subZoneLayer
    public float getSubZoneGridSliceZ(int nSubZone, int nSlice, int row, int col)
    {
        float[][][] subZoneGridZ = getSubZoneGridZs(nSubZone);
        if(subZoneGridZ != null)
        {
            float zTop = subZoneGridZ[nSlice][row - rowMin][col - colMin];
            float zBot = subZoneGridZ[nSlice+1][row - rowMin][col - colMin];
            return 0.5f*(zTop + zBot);
        }
        else
        {
            if(this.subZoneGridZ == null) return nullValue;

            float topZ = getZoneGridZ(nSubZone, row, col);
            if(topZ == nullValue) return nullValue;
            if(nSlice == 0) return topZ;

            float botZ = getZoneGridZ(nSubZone+1, row, col);
            if(botZ == nullValue) return nullValue;

            StsObjectRefList subZones = zone.getSubZones();
            StsSubZone subZone = (StsSubZone)subZones.getElement(nSubZone);
            if(subZone == null) return nullValue;

            int subZoneType = subZone.getSubZoneType();
            int nLayers = subZone.getNLayers();

            float f, dZ, thickness;
            switch(subZoneType)
            {
                case StsZone.SUBZONE_UNIFORM:
                    f = (nSlice+0.5f)/nLayers;
                    return topZ + f*(botZ - topZ);
                case StsZone.SUBZONE_OFFLAP:
                    thickness = subZone.getLayerThickness();
                    dZ = (nLayers - nSlice - 0.5f)*thickness;
                    return Math.max(botZ - dZ, topZ);
                case StsZone.SUBZONE_ONLAP:
                    thickness = subZone.getLayerThickness();
                    dZ = (nSlice + 0.5f)*thickness;
                    return Math.min(topZ + dZ, botZ);
                default:
                    return nullValue;
            }
        }
    }

    private float[][][] getSubZoneGridZs(int nSubZone)
    {
        if(layerGridZs == null) return null;
        if(nSubZone >= layerGridZs.length) return null;
        return layerGridZs[nSubZone];
    }
/*
        else
        {

    }

    public float getSubZoneGridZ(int nSubZone, int nSlice, int subZoneType, int nLayers, int row, int col)
    {
        if(subZoneGridZs != null && subZoneGridZs[nSubZone] != null)
            return subZoneGridZs[nSubZone][nSlice][row - boundingBox.rowMin][col - boundingBox.colMin];
        else
        {
            if(gridZ == null) return nullValue;

            float topZ = gridZ[nSubZone][row][col];
            float botZ = gridZ[nSubZone+1][row][col];
            float f, dZ;

            switch(subZoneType)
            {
                case StsZone.SUBZONE_UNIFORM:
                    f = ((float)nSlice)/nSubZones;
                    return botZ + f*(topZ - botZ);
                case StsZone.SUBZONE_OFFLAP:
                    dZ = nSlice*zone.getSubZoneThickness();
                    dZ = Math.min(dZ, botZ - topZ);
                    return botZ - dZ;
                case StsZone.SUBZONE_ONLAP:
                    dZ = (nLayers - nSlice)*zone.getSubZoneThickness();
                    dZ = Math.min(dZ, botZ - topZ);
                    return topZ + dZ;
                default:
                    return nullValue;
            }
        }
    }

    public float getAbsGridZ(int nLayer, int row, int col)
    {
        float z = StsParameters.nullValue;
        try { z = gridZ[nSlice][row][col]; }
        catch(Exception e) { e.printStackTrace(); }
        return z;
    }

    public float getRelGridZ(int nSlice, int row, int col)
    {
        return gridZ[nSlice][row - boundingBox.rowMin][col - boundingBox.colMin];
    }
*/
    // ij is normal row column ordering; k is ordered from the bottom up
    public boolean constructGrid()
    {
        int i, j, k, n;

        try
        {
            zone.checkLoadSubHorizonGrids();

            int nRows = getNRows();
            int nCols = getNCols();
            nSubZones = zone.getNSubZones();

            if( nSubZones == 0 )
            {
                StsException.systemError(this, "constructGrid", "No subzones in ZoneBlock.constructSubZoneGrid()");
                return false;
            }
            // we will refer to subHorizonSurfaces as subHorizons
            StsObjectRefList subHorizons = zone.getSubHorizonSurfaces();
            int nSubHorizons = subHorizons == null ? 0 : subHorizons.getSize();
            /*
            if(nSubHorizons != nSubZones+1)
                StsException.systemError(this, "constructZoneGrid",
                    "Number of subzones (" + nSubZones + ") != number of zone subHorizons + 1 (" + nSubHorizons + "+1)");
            */
            subZoneGridZ = new float[nSubZones+1][][];

            if( botGrid == null ||  topGrid == null )
                StsException.systemError(this, "constructGrid",  "Bottom or top grid is null in ZoneBlock.constructSubZoneGrid()");

            float[][] topPoints = topGrid.getPointsZ();
            float[][] botPoints = botGrid.getPointsZ();

            if( botPoints == null ||  topPoints == null )
                StsException.systemError(this, "constructGrid", "Bottom or top points is null in ZoneBlock.constructSubZoneGrid()");

            subZoneGridZ[0] = topPoints;
            subZoneGridZ[nSubZones] = botPoints;

            for(k = 1; k < nSubZones; k++)
                subZoneGridZ[k] = new float[nRows][nCols];

            if(nSubHorizons <= 2)
            {
                for(i = 0; i < nRows; i++)
                    for(j = 0; j < nCols; j++)
                        constructZoneGrid(i, j, nSubZones, subZoneGridZ);
                return true;
            }

            float[] zIJs = new float[nSubZones+1];  // scratch vector of Zs at IJ location
            float[][][] subHorizonZs = new float[nSubHorizons][][];
            for(n = 0; n < nSubHorizons; n++)
            {
                StsModelSurface surface = (StsModelSurface)subHorizons.getElement(n);
                subHorizonZs[n] = surface.getPointsZ();
            }
            
            for(i = 0; i < nRows; i++)
                for(j = 0; j < nCols; j++)
                    if(!constructZoneGrid(i, j, subHorizonZs, nSubHorizons, zIJs, subZoneGridZ))
                        constructZoneGrid(i, j, nSubZones, subZoneGridZ);
            return true;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "constructSubZoneGrid", e);
            return false;
        }
    }

    public float getZoneGridZ(int nSlice, int row, int col)
    {
        if(subZoneGridZ == null) return nullValue;
        if(subZoneGridZ.length < nSlice+1) return nullValue;
        return subZoneGridZ[nSlice][row - rowMin][col - colMin];
    }

    public boolean constructLayerGrids()
    {
        if(subZoneGridZ == null && !constructGrid()) return false;

        try
        {
            StsObjectRefList subZones = zone.getSubZones();
            int nSubZones = subZones.getSize();
            layerGridZs = new float[nSubZones][][][];

            for(int sz = 0; sz < nSubZones; sz++)
            {
                StsSubZone subZone = (StsSubZone)subZones.getElement(sz);
                int nLayers = subZone.getNLayers();
                float[][][] subZoneGridZ = new float[nLayers+1][][];
                layerGridZs[sz] = subZoneGridZ;
                subZoneGridZ[0] = this.subZoneGridZ[sz];
                subZoneGridZ[nLayers] = this.subZoneGridZ[sz+1];
                if(nLayers > 1)
                {
                    int subZoneType = subZone.getSubZoneType();
                    constructLayerGrid(nLayers, subZoneGridZ, subZoneType);
                }
            }
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsZoneBlock.constructLayerGrid() failed.", e, StsException.WARNING);
            return false;
        }
    }

    private void constructLayerGrid(int nLayers, float[][][] subZoneGridZ, int subZoneType)
    {
        int nRows = getNRows();
        int nCols = getNCols();
        for(int i = 0; i < nRows; i++)
        {
            for(int j = 0; j < nCols; j++)
            {
                float topZ = subZoneGridZ[0][i][j];
                float botZ = subZoneGridZ[nLayers][i][j];

                if(botZ == nullValue || topZ == nullValue)
                {
                    for(int k = 0; k <= nLayers; k++)
                        subZoneGridZ[k][i][j] = nullValue;
                    return;
                }

                float dZ, z;

                switch(subZoneType)
                {
                    case StsZone.SUBZONE_UNIFORM:
                        dZ = (botZ - topZ)/nLayers;
                        z = topZ + dZ;
                        for(int k = 1; k < nLayers; k++)
                        {
                            subZoneGridZ[k][i][j] = z;
                            z += dZ;
                        }
                        break;
                    case StsZone.SUBZONE_OFFLAP:
                        dZ = zone.getSubZoneThickness();
                        z = botZ - dZ;
                        for(int k = nLayers-1; k >= 1 ; k--)
                        {
                            z = Math.max(z, topZ);
                            subZoneGridZ[k][i][j] = z;
                            z -= dZ;
                        }
                        break;
                    case StsZone.SUBZONE_ONLAP:
                        dZ = zone.getSubZoneThickness();
                        z = topZ + dZ;
                        for(int k = 1; k < nLayers; k++)
                        {
                            z = Math.min(z, botZ);
                            subZoneGridZ[k][i][j] = z;
                            z += dZ;
                        }
                        break;
                }
            }
        }
    }

    private void constructZoneGrid(int i, int j, int nSubZones, float[][][] gridZ)
    {
        float topZ = gridZ[0][i][j];
        float botZ = gridZ[nSubZones][i][j];

        if(botZ == nullValue || topZ == nullValue)
        {
            for(int k = 1; k < nSubZones; k++)
                gridZ[k][i][j] = nullValue;
            return;
        }

        float dZ, z;
        switch(zone.getSubZoneType())
        {
            case StsZone.SUBZONE_UNIFORM:
                dZ = (botZ - topZ)/(nSubZones);
                z = topZ + dZ;
                for(int k = 1; k < nSubZones; k++)
                {
                    gridZ[k][i][j] = z;
                    z += dZ;
                }
                break;
            case StsZone.SUBZONE_OFFLAP:
                dZ = zone.getSubZoneThickness();
                z = botZ - dZ;
                for(int k = nSubZones; k >= 1 ; k--)
                {
                    z = Math.max(z, topZ);
                    gridZ[k][i][j] = z;
                    z -= dZ;
                }
                break;
            case StsZone.SUBZONE_ONLAP:
                z = topZ;
                dZ = zone.getSubZoneThickness();
                for(int k = 1; k < nSubZones; k++)
                {
                    gridZ[k][i][j] = z;
                    z += dZ;
                    z = Math.min(z, botZ);
                }
                break;
        }
    }

    // This column of Z values is defined by the subHorizons, so there is no changes for
    // subZoneType as this is ignored.
    private boolean constructZoneGrid(int i, int j, float[][][] subHorizonZs, int nSubHorizons,
                                      float[] zIJs, float[][][] zoneGridZ)
    {
        int k;

        if(subHorizonZs == null) return false;

        try
        {
            float topZ = zoneGridZ[0][i][j];
            float botZ = zoneGridZ[nSubHorizons-1][i][j];
            if(botZ == nullValue || topZ == nullValue)
            {
                for(k = 0; k < nSubHorizons; k++)
                    zoneGridZ[k][i][j] = nullValue;
                return true;
            }

            // get Z values from voxelGrids
            for(k = 0; k < nSubHorizons; k++)
                zIJs[k] = subHorizonZs[k][i][j];

            // adjust so no dZ is negative
            for(k = 1; k < nSubHorizons; k++)
                if(zIJs[k] < zIJs[k-1]) zIJs[k] = zIJs[k-1];

            float topHorizonZ = zIJs[0];
            float botHorizonZ = zIJs[nSubHorizons-1];
            if(botHorizonZ <= topHorizonZ) return false;

            // now rescale so that they match top and bot blockGrid z values
            float offset = topZ - zIJs[0];
            float scale = (botZ - topZ)/(botHorizonZ - topHorizonZ);
            if(StsMath.sameAs(offset, 0.0f) && StsMath.sameAs(scale, 1.0f))
            {
                for(k = 1; k < nSubHorizons-1; k++)
                    zoneGridZ[k][i][j] = zIJs[k];
            }
            else
            {
                zoneGridZ[0][i][j] = topZ;

                for(k = 1; k < nSubHorizons-1; k++)
                    zoneGridZ[k][i][j] = topZ + (zIJs[k] - topHorizonZ)*scale;

                zoneGridZ[nSubHorizons-1][i][j] = botZ;
            }
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsZoneBlock.constructZoneGrid() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    public StsPropertyVolume getPropertyVolume(StsPropertyType propertyType)
    {
        String propertyName = propertyType.name;
        int nPropertyVolumes = zoneProperties.size();
        for(int n = 0; n < nPropertyVolumes; n++)
        {
            StsPropertyVolume propertyVolume = zoneProperties.get(n);
            if(propertyVolume.getName().equals(propertyName))
                return propertyVolume;
        }
        return null;
    }

    public StsPropertyVolume getPropertyVolume(String propertyName)
    {
        if(zoneProperties == null) return null;
        int nPropertyVolumes = zoneProperties.size();
        for(int n = 0; n < nPropertyVolumes; n++)
        {
            StsPropertyVolume propertyVolume = zoneProperties.get(n);
            if(propertyVolume.getName().equals(propertyName))
                return propertyVolume;
        }
        return null;
    }

    public StsPropertyVolume getEclipsePropertyVolume(String eclipsePropertyName)
    {
        if(zoneProperties == null)
            return zone.getEclipsePropertyVolume(eclipsePropertyName);
        int nPropertyVolumes = zoneProperties.size();
        for(int n = 0; n < nPropertyVolumes; n++)
        {
            StsPropertyVolume propertyVolume = zoneProperties.get(n);
            if(propertyVolume.propertyType.eclipseName.equals(eclipsePropertyName))
                return propertyVolume;
        }
        return zone.getEclipsePropertyVolume(eclipsePropertyName);
    }

    public void displayGridSlice(StsGLPanel3d glPanel3d, int nSlice, boolean displayGaps)
    {
/*
        int i, ii, j, jj;
        float[][] sliceZs;

        try
        {
            // borrow quadStrips from bot blockGrid (has the same ij info)
            StsList quadStrips = botGrid.getQuadStrips();
            if(quadStrips == null) return;

            GL gl = glPanel3d.getGL();

            if(nSlice != currentNSlice && sliceSurfDisplayListNum > 0)
            {
                gl.glDeleteLists(sliceSurfDisplayListNum, 1 );
                sliceSurfDisplayListNum = 0;
            }

            currentNSlice = nSlice;

            boolean useDisplayLists = currentModel.getBooleanProperty("Use Display Lists");

            if(sliceSurfDisplayListNum > 0 && useDisplayLists)
            {
                gl.glCallList( sliceSurfDisplayListNum );
                //timer.stop("display list surface draw: ");
                return;
            }

            sliceZs = getSliceZs(nSlice);
            if(sliceZs == null) return;

            sliceNormals = StsGridUtilities.computeGridNormals(sliceZs, boundingBox, sliceNormals);

            if (sliceSurfDisplayListNum == 0 && useDisplayLists)  // build display list
            {
                sliceSurfDisplayListNum = gl.glGenLists(1);
                if(sliceSurfDisplayListNum == 0)
                {
                    currentModel.model.statusBarMessage("System Error in StsZoneBlock.displayGridSlice(): " +
                                           "Failed to allocate a display list");
                    return;
                }

                gl.glNewList(sliceSurfDisplayListNum, GL.GL_COMPILE_AND_EXECUTE);
                drawGridSliceSurface(gl, nSlice, sliceZs, quadStrips, displayGaps);
                gl.glEndList();
                //timer.stop("display list surface setup: ");
                sliceZs = null;
                sliceNormals = null;
            }
            else // immediate mode draw
            {
                if(sliceSurfDisplayListNum > 0)
                {
                    gl.glDeleteLists(sliceSurfDisplayListNum, 1 );
                    sliceSurfDisplayListNum = 0;
                }
                drawGridSliceSurface(gl, nSlice, sliceZs, quadStrips, displayGaps);
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsZoneBlock.displayGridSlice() failed.",
                e, StsException.WARNING);
        }
*/
    }

    public float[][] getSliceZs(int nSlice)
    {
        int i, ii, j, jj;

        try
        {
            int[] subZoneAndLayer = zone.getSubZoneAndLayer(nSlice);
            if(subZoneAndLayer == null) return null;

            int nSubZone = subZoneAndLayer[0];
            int nSubZoneSlice = subZoneAndLayer[1];
            int nRows = getNRows();
            int nCols = getNCols();
            float[][] sliceZs = new float[nRows][nCols];
            for(i = rowMin, ii = 0; i <= rowMax; i++, ii++)
                for(j = colMin, jj = 0; j <= colMax; j++, jj++)
                    sliceZs[ii][jj] = getSubZoneGridSliceZ(nSubZone, nSubZoneSlice, i, j);

            return sliceZs;
        }
        catch(Exception e)
        {
            StsException.outputException("StsZoneBlock.getSliceZs() failed.",
                e, StsException.WARNING);
            return null;
        }
    }
/*
    private void drawGridSliceSurface(GL gl, int nSlice, float[][] sliceZs, StsList quadStrips, boolean displayGaps)
    {

        StsPropertyVolumeOld pv = zone.getPropertyVolume();
        if(pv != null && pv.scaledColorsAreValid())
        {
            StsColor[][] colors = pv.getSliceColors(boundingBox, nSlice);
            if(colors == null) return;
            StsGLDraw.drawQuadStripProps(gl, boundingBox, quadStrips, sliceZs,
                sliceNormals, displayGaps, nSlice, colors);
        }
        else
        {
            zone.getStsColor().setGLColor(gl);
            StsGLDraw.drawQuadStrips(gl, boundingBox, quadStrips, sliceZs, sliceNormals, displayGaps);
        }
    }
*/
    public void deleteDisplayLists()
    {
        GL gl = currentModel.getWin3dGL();
        if(gl == null) return;

        deleteSliceSurfDisplayList(gl);
//        if(zoneBlockSides != null) zoneBlockSides.forEach("deleteDisplayLists", gl);
        if(zoneBlockSides == null) return;
        int nSides = zoneBlockSides.getSize();
        for(int n = 0; n < nSides; n++)
        {
            StsZoneSide zoneSide = (StsZoneSide)zoneBlockSides.getElement(n);
            zoneSide.deleteDisplayLists(gl);
        }
    }

    public void deleteSliceSurfDisplayList(GL gl)
    {
        if(sliceSurfDisplayListNum <= 0) return;

        gl.glDeleteLists(sliceSurfDisplayListNum, 1 );
        sliceSurfDisplayListNum = 0;
    }

    private void constructZoneBlockSides()
    {
        StsSurfaceEdge topEdge, botEdge;
        StsModelSurface topSurface, botSurface;
        StsZoneSide zoneSide;

        try
        {
            StsObjectRefList blockSides = block.getBlockSides();
            if(blockSides == null) return;
            int nBlockSides = blockSides.getSize();
            if(nBlockSides == 0) return;

            if(zoneBlockSides == null) zoneBlockSides = new StsObjectList(nBlockSides, 1);

            topSurface = zone.getTopModelSurface();
            botSurface = zone.getBaseModelSurface();

            for(int n = 0; n < nBlockSides; n++)
            {
                StsBlockSide blockSide = (StsBlockSide)blockSides.getElement(n);
                topEdge = topGrid.getSurfaceEdge(blockSide.section);
                if(topEdge == null) continue;
                botEdge = botGrid.getSurfaceEdge(blockSide.section);
                if(botEdge == null) continue;
                zoneSide = new StsZoneSide(zone, this, topEdge, botEdge);
                zoneBlockSides.add(zoneSide);
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsZoneBlock.constructZoneBlockSides() failed.",
                e, StsException.WARNING);
        }
    }

    // nSlice is index into zone of all subZone layers. Get subZone and subZone slice.
    private float getSubZoneGridZ(int nSlice, float rowF, float colF)
    {
        int[] subZoneAndSlice = zone.getSubZoneAndLayer(nSlice);
        if(subZoneAndSlice == null) return nullValue;

        int nSubZone = subZoneAndSlice[0];
        nSlice = subZoneAndSlice[1];
        return getSubZoneGridZ(nSubZone, nSlice, rowF, colF);
    }

    public float getSubZoneF(float[] xyz)
    {
        float rowF = topGrid.getRowCoor(xyz);
        float colF = topGrid.getColCoor(xyz);
        return getZoneF(rowF, colF, xyz[2]);
    }
/*
    public float getSubZoneF(float rowF, float colF, float z)
    {
        int i, j;
        float dx, dy;
        float zTop, zBot;
        float w, weight = 0.0f;
        float zTopW = 0.0f, zBotW = 0.0f;
        boolean allNulls = true;

        try
        {
            if(boundingBox == null) boundingBox = block.getBoundingBox();
//            if(!isInsideRowCol(rowF, colF)) return nullValue;

            i = (int)rowF;
            if(i >= rowMax) i = rowMax - 1;
            else if(i < rowMin) i = rowMin;
            dy = rowF - i;

            j = (int)colF;
            if(j >= colMax) j = colMax - 1;
            else if(j < colMin) j = colMin;
            dx = colF - j;

            zTop = topGrid.getComputeExtrapolatedZ(i, j);
            zBot = botGrid.getComputeExtrapolatedZ(i, j);
            if(zTop != nullValue && zBot != nullValue)
            {
                w = (1.0f-dy) * (1.0f-dx);
                weight += w;
                zTopW += w*zTop;
                zBotW += w*zBot;
                allNulls = false;
            }

            zTop = topGrid.getComputeExtrapolatedZ(i+1, j);
            zBot = botGrid.getComputeExtrapolatedZ(i+1, j);
            if(zTop != nullValue && zBot != nullValue)
            {
                w = dy * (1.0f-dx);
                weight += w;
                zTopW += w*zTop;
                zBotW += w*zBot;
                allNulls = false;
            }

            zTop = topGrid.getComputeExtrapolatedZ(i, j+1);
            zBot = botGrid.getComputeExtrapolatedZ(i, j+1);
            if(zTop != nullValue && zBot != nullValue)
            {
                w = (1.0f-dy) * dx;
                weight += w;
                zTopW += w*zTop;
                zBotW += w*zBot;
                allNulls = false;
            }

            zTop = topGrid.getComputeExtrapolatedZ(i+1, j+1);
            zBot = botGrid.getComputeExtrapolatedZ(i+1, j+1);
            {
                w = dy*dx;
                weight += w;
                zTopW += w*zTop;
                zBotW += w*zBot;
                allNulls = false;
            }

            if (allNulls) return nullValue;

            float f;
            if(weight == 1.0f)
                f = (z - zTopW)/(zBotW - zTopW);
            else
                f = (z*weight - zTopW)/(zBotW - zTopW);

            if(f < 0.001f) f = 0.0f;
            else if(f > 0.999f) f = 1.0f;
            return f;
        }
        catch(Exception e)
        {
            StsException.outputException("StsZoneBlock.getSubZoneF() failed.",
                e, StsException.WARNING);
            return nullValue;
        }
    }
*/
    

    public float getZoneF(float rowF, float colF, float z)
    {
        int i, j;
        float dx, dy;
        float zTop, zBot;
        float w, weight = 0.0f;
        float zTopW = 0.0f, zBotW = 0.0f;
        boolean allNulls = true;

        try
        {
            i = (int)rowF;
            if(i >= rowMax) i = rowMax - 1;
            else if(i < rowMin) i = rowMin;
            dy = rowF - i;

            j = (int)colF;
            if(j >= colMax) j = colMax - 1;
            else if(j < colMin) j = colMin;
            dx = colF - j;

            // check if point is on corners
            if(dx == 0.0f)
            {
               if(dy == 0.0f)
                {
                    zTop = topGrid.getComputeExtrapolatedZ(i, j);
                    zBot = botGrid.getComputeExtrapolatedZ(i, j);
                    return checkRoundoffZoneF(z, zTop, zBot);
                }
                else if(dy == 1.0f)
                {
                    zTop = topGrid.getComputeExtrapolatedZ(i+1, j);
                    zBot = botGrid.getComputeExtrapolatedZ(i+1, j);
                    return checkRoundoffZoneF(z, zTop, zBot);
                }
            }
            else if(dx == 1.0f)
            {
               if(dy == 0.0f)
                {
                    zTop = topGrid.getComputeExtrapolatedZ(i, j+1);
                    zBot = botGrid.getComputeExtrapolatedZ(i, j+1);
                    return checkRoundoffZoneF(z, zTop, zBot);
                }
                else if(dy == 1.0f)
                {
                    zTop = topGrid.getComputeExtrapolatedZ(i+1, j+1);
                    zBot = botGrid.getComputeExtrapolatedZ(i+1, j+1);
                    return checkRoundoffZoneF(z, zTop, zBot);
                }
            }

            if(dy != 1.0f && dx != 1.0f)
            {
                zTop = topGrid.getComputeExtrapolatedZ(i, j);
                zBot = botGrid.getComputeExtrapolatedZ(i, j);
                if(zTop != nullValue && zBot != nullValue)
                {
                    w = (1.0f-dy) * (1.0f-dx);
                    weight += w;
                    zTopW += w*zTop;
                    zBotW += w*zBot;
                    allNulls = false;
                }
            }

            if(dy != 0.0f && dx != 1.0f)
            {
                zTop = topGrid.getComputeExtrapolatedZ(i+1, j);
                zBot = botGrid.getComputeExtrapolatedZ(i+1, j);
                if(zTop != nullValue && zBot != nullValue)
                {
                    w = dy * (1.0f-dx);
                    weight += w;
                    zTopW += w*zTop;
                    zBotW += w*zBot;
                    allNulls = false;
                }
            }

            if(dy != 1.0f && dx != 0.0f)
            {
                zTop = topGrid.getComputeExtrapolatedZ(i, j+1);
                zBot = botGrid.getComputeExtrapolatedZ(i, j+1);
                if(zTop != nullValue && zBot != nullValue)
                {
                    w = (1.0f-dy) * dx;
                    weight += w;
                    zTopW += w*zTop;
                    zBotW += w*zBot;
                    allNulls = false;
                }
            }

            if(dy != 0.0f && dx != 0.0f)
            {
                zTop = topGrid.getComputeExtrapolatedZ(i+1, j+1);
                zBot = botGrid.getComputeExtrapolatedZ(i+1, j+1);
                {
                    w = dy*dx;
                    weight += w;
                    zTopW += w*zTop;
                    zBotW += w*zBot;
                    allNulls = false;
                }
            }

            if (allNulls) return nullValue;

            float f;
                return checkRoundoffZoneF(z*weight, zTopW, zBotW);
        }
        catch(Exception e)
        {
            StsException.outputException("StsZoneBlock.getSubZoneF() failed.",
                e, StsException.WARNING);
            return nullValue;
        }
    }

    private float checkRoundoffZoneF(float f)
    {
        if(f < 0.001f) f = 0.0f;
        else if(f > 0.999f) f = 1.0f;
        return f;
    }

    private float checkRoundoffZoneF(float z, float zTop, float zBot)
    {
        float f;
        if(zBot > zTop)
        {
            f = (z - zTop)/(zBot - zTop);
            if(f < 0.001f) f = 0.0f;
            else if(f > 0.999f) f = 1.0f;
            return f;
        }
        else
            return 0.0f;
    }

    /** For now we are assuming that subZones are equally spaced in zone and each subZone has only one layer */
    public float getLayerF(float rowF, float colF, float z)
    {
        float f = getZoneF(rowF, colF, z);
        if(f == nullValue)
            return 0.0f;
        return zone.getLayerF(f);
    }
    // nSlice is the index in the subZone from 0 to nLayers+1
    private float getSubZoneGridZ(int nSubZone, int nSlice, float rowF, float colF)
    {
        int i, j;
        float dx, dy;
        float z, w, weight = 0, zWeighted = 0;
        boolean allNulls = true;

        try
        {
            if(!isInsideRowCol(rowF, colF)) return nullValue;

            i = (int)rowF;
            if(i == rowMax) i--;
            dy = rowF - i;

            j = (int)colF;
            if(j == colMax) j--;
            dx = colF - j;

            if( (z = getSubZoneGridZ(nSubZone, nSlice, i, j)) != nullValue)
    //        if((w = (1.0f-dy) * (1.0f-dx)) > StsParameters.roundOff &&
    //           ((z = getGridZ(nSlice, i, j)) != nullValue))
            {
                w = (1.0f-dy) * (1.0f-dx);
                weight += w;
                zWeighted += w*z;
                allNulls = false;
            }
            if( (z = getSubZoneGridZ(nSubZone, nSlice, i+1, j)) != nullValue)
    //        if((w = dy * (1.0f-dx)) > StsParameters.roundOff &&
    //           ((z = getGridZ(nSlice, i+1, j)) != nullValue))
            {
                w = dy * (1.0f-dx);
                weight += w;
                zWeighted += w*z;
                allNulls = false;
            }
            if( (z = getSubZoneGridZ(nSubZone, nSlice, i, j+1)) != nullValue)
    //        if((w = (1.0f-dy) * dx) > StsParameters.roundOff &&
    //           ((z = getGridZ(nSlice, i, j+1)) != nullValue))
            {
                w = (1.0f-dy) * dx;
                weight += w;
                zWeighted += w*z;
                allNulls = false;
            }
            if( (z = getSubZoneGridZ(nSubZone, nSlice, i+1, j+1)) != nullValue)
    //        if((w = dy*dx) > StsParameters.roundOff &&
    //           ((z = getGridZ(nSlice, i+1, j+1)) != nullValue))
            {
                w = dy*dx;
                weight += w;
                zWeighted += w*z;
                allNulls = false;
            }

            if (allNulls) return nullValue;

            if(weight == 1.0f)
                return zWeighted;
            else
                return zWeighted / weight;
        }
        catch(Exception e)
        {
            StsException.outputException("StsZoneBlock.getInterpolatedZ() failed.",
                e, StsException.WARNING);
            return nullValue;
        }
    }

    public boolean fillZoneBlockArray(StsZoneBlock[][] zoneBlockArray,
                        int rowMin, int rowMax, int colMin, int colMax)
    {
        if(topGrid == null) return false;
        return topGrid.fillZoneBlockArray(zoneBlockArray, rowMin, rowMax, colMin, colMax, this);
    }

    public void display(StsModel model, StsGLPanel3d glPanel3d, String displayMode, StsPropertyVolume propertyVolume)
    {
        if(zoneBlockSides == null) return;
        if(!block.getIsVisible()) return;

        StsColor stsColor = null;
        boolean displayLayers = false;
        boolean displayProperties = false;
        if(displayMode == StsBuiltModelClass.displayBlocksString)
            stsColor = getStsColor();
        else if(displayMode == StsBuiltModelClass.displayZonesString)
            stsColor = zone.getStsColor();
        else if(displayMode == StsBuiltModelClass.displayLayersString)
        {
            displayLayers = true;
        }
        else if(displayMode == StsBuiltModelClass.staticGetDisplayModeString())
        {
            displayLayers = true;
            displayProperties = true;
        }

        int nZoneBlockSides = zoneBlockSides.getSize();
        for(int n = 0; n < nZoneBlockSides; n++)
        {
            StsZoneSide zoneSide = (StsZoneSide)zoneBlockSides.getElement(n);
            zoneSide.display(model, glPanel3d, displayMode, displayLayers, displayProperties, stsColor, propertyVolume);
        }
    }

    private StsColor getBlockColor()
    {
        return getStsColor();
    }

    public void drawCursorSection(StsModel model, StsGLPanel3d glPanel3d, int dir, float dirCoordinate, StsPoint point0, StsPoint point1,
                                  boolean isDragging, String displayMode)
    {
        GL gl = glPanel3d.getGL();
        if(gl == null) return;

        if(!block.isVisible) return;

        boolean displayFill = zone.getDisplayFill();
        boolean displayGrid = zone.getDisplayGrid();

        if(dir == StsCursor3d.XDIR)
	    {
            StsGridPoint p0 = new StsGridPoint(point0, botGrid);
            int col = p0.col;

            if(currentCursorXCoor != dirCoordinate)
            {
                currentCursorXCoor = dirCoordinate;
                cursorCol = col;
                cursorColGridLines = constructCursorGridLines(COL, col);
                cursorColEdgeLoops = null;
            }

            if(isDragging)
                drawCursorGridLines(cursorColGridLines, gl);
            else
            {
                if(cursorColEdgeLoops == null)
                    cursorColEdgeLoops = constructCursorEdgeLoops(COL, col, cursorColGridLines, colSectionNormal);
                StsPropertyVolume propertyVolume = block.getCurrentPropertyVolume();
                drawCursorEdgeLoops(model, cursorColEdgeLoops, glPanel3d, displayFill, displayGrid, displayMode, propertyVolume);
            }
        }
        else if(dir == StsCursor3d.YDIR)
        {
            StsGridPoint p0 = new StsGridPoint(point0, botGrid);
            int row = p0.row;

            if(currentCursorYCoor != dirCoordinate)
            {
                currentCursorYCoor = dirCoordinate;
                cursorRow = row;
                cursorRowGridLines = constructCursorGridLines(ROW, row);
                cursorRowEdgeLoops = null;
            }

            if(isDragging)
                drawCursorGridLines(cursorRowGridLines, gl);
            else
            {
                if(cursorRowEdgeLoops == null)
                    cursorRowEdgeLoops = constructCursorEdgeLoops(ROW, row, cursorRowGridLines, rowSectionNormal);
                StsPropertyVolume propertyVolume = block.getCurrentPropertyVolume();
                drawCursorEdgeLoops(model, cursorRowEdgeLoops, glPanel3d, displayFill, displayGrid, displayMode, propertyVolume);
            }
        }
        else if(dir == StsCursor3d.ZDIR)
        {
            /*
            float z = point0.getZorT();
            int slice = StsBuiltModel.builtModel.modelBoundingBox.getNearestSliceCoor(z);

            if(currentCursorZCoor != dirCoordinate)
            {
                currentCursorZCoor = dirCoordinate;
                cursorSlice = slice;
                cursorSliceGridLines = constructCursorSliceLines(ROW, slice);
                cursorRowEdgeLoops = null;
            }

            if(isDragging)
                drawCursorGridLines(cursorRowGridLines, gl);
            else
            {
                if(cursorRowEdgeLoops == null)
                    cursorRowEdgeLoops = constructCursorEdgeLoops(ROW, slice, cursorRowGridLines, null);
                drawCursorEdgeLoops(cursorSliceEdgeLoops, glPanel3d, displayFill, displayGrid, displayBlockColors, displayZoneColors);
            }
            */
        }
    }

    private void drawCursorEdgeLoops(StsModel model, StsList cursorEdgeLoops, StsGLPanel3d glPanel3d,
                                     boolean displayFill, boolean displayGrid, String displayMode, StsPropertyVolume propertyVolume)
    {
        StsEdgeLoop cursorEdgeLoop;

        if(cursorEdgeLoops == null) return;
        int nEdgeLoops = cursorEdgeLoops.getSize();
        if(nEdgeLoops == 0) return;

        boolean displaySubPolygons = false;
        boolean displayProperties = false;
        StsColor stsColor = null;
        if(displayMode == StsBuiltModelClass.displayBlocksString)
            stsColor = getStsColor();
        else if(displayMode == StsBuiltModelClass.displayZonesString)
            stsColor = zone.getStsColor();
        else if(displayMode == StsBuiltModelClass.displayLayersString)
        {
            displaySubPolygons = true;
        }
        else if(displayMode == StsBuiltModelClass.staticGetDisplayModeString())
        {
            displaySubPolygons = true;
            displayProperties = true;
        }

 //       boolean useAdjustedPoints = currentModel.project.getUseAdjustedPoints();
        for(int n = 0; n < nEdgeLoops; n++)
        {
            cursorEdgeLoop = (StsEdgeLoop)cursorEdgeLoops.getElement(n);
            cursorEdgeLoop.checkSetPropertyChanged(displayMode);
            if(displaySubPolygons && !cursorEdgeLoop.subPolygonsConstructed)
                cursorEdgeLoop.constructCursorSubPolygons(this);
            cursorEdgeLoop.display(model, glPanel3d, displayFill, displayGrid, displaySubPolygons, displayProperties, stsColor, propertyVolume, true);
        }
    }

    private StsList constructCursorGridLines(int rowOrCol, int rowCol)
    {
        StsList cursorGridLines;
        Object[] gridLines;
        StsEdgeLinkable gridLine;

        if(topGrid == null || botGrid == null) return null;
        StsEdgeLoop topGridEdgeLoop = topGrid.getEdgeLoop();
        if(topGridEdgeLoop == null) return null;
        StsEdgeLoop botGridEdgeLoop = botGrid.getEdgeLoop();
        if(botGridEdgeLoop == null) return null;
        try
        {
            cursorGridLines = new StsList(4, 4);

            gridLines = topGridEdgeLoop.getGridLines(rowOrCol, rowCol);
            cursorGridLines.add(gridLines);

            gridLines = botGridEdgeLoop.getGridLines(rowOrCol, rowCol);
            cursorGridLines.add(gridLines);

            if(zoneBlockSides != null)
            {
                int nZoneSides = zoneBlockSides.getSize();
                for(int n = 0; n < nZoneSides; n++)
                {
                    StsZoneSide zoneSide = (StsZoneSide)zoneBlockSides.getElement(n);
                    gridLines = zoneSide.getEdgeLoop().getGridLines(rowOrCol, rowCol);
                    cursorGridLines.add(gridLines);
                }
            }
            return cursorGridLines;
        }
        catch(Exception e)
        {
            StsException.outputException("StsZoneBlock.constructCursorGridLines() failed.",
                    e, StsException.WARNING);
            return null;
        }
    }
    /** not finished */
    private StsList constructCursorSliceLines(int rowOrCol, int slice)
    {
        StsList cursorSliceLines;
        Object[] sliceLines;

        try
        {
            cursorSliceLines = new StsList(4, 4);

            int nBlockSides = zoneBlockSides.getSize();
            for(int n = 0; n < nBlockSides; n++)
            {
                StsBlockSide blockSide = (StsBlockSide)zoneBlockSides.getElement(n);
            }
            return cursorSliceLines;
        }
        catch(Exception e)
        {
            StsException.outputException("StsZoneBlock.constructCursorGridLines() failed.",
                    e, StsException.WARNING);
            return null;
        }
    }

    private StsList constructCursorEdgeLoops(int rowOrCol, int rowCol, StsList cursorGridLines, float[] normal)
    {
        StsEdgeLinkable gridLine;
        StsEdgeLoop cursorEdgeLoop;

        try
        {
            StsList cursorEdgeLoops = new StsList(1,2);

            while( (gridLine = getNextGridLine(cursorGridLines)) != null )
            {
                cursorEdgeLoop = constructCursorEdgeLoop(rowOrCol, rowCol, gridLine, cursorGridLines, normal);
                if(cursorEdgeLoop != null) cursorEdgeLoops.add(cursorEdgeLoop);
            }
            return cursorEdgeLoops;
        }
        catch(Exception e)
        {
            StsException.outputException("StsZoneBlock.constructCursorEdgeLoops() failed.",
                    e, StsException.WARNING);
            return null;
        }
    }

    private StsEdgeLinkable getNextGridLine(StsList cursorGridLines)
    {
        if(cursorGridLines == null) return null;
        int nGridLines = cursorGridLines.getSize();
        if(nGridLines <= 0) return null;

        StsEdgeLinkable nextGridLine = (StsEdgeLinkable)cursorGridLines.getFirst();
        cursorGridLines.delete(nextGridLine);
        return nextGridLine;
    }

    private StsEdgeLoop constructCursorEdgeLoop(int rowOrCol, int rowCol, StsEdgeLinkable firstGridLine,
                                StsList cursorGridLines, float[] normal)
    {
        StsEdgeLoop cursorEdgeLoop;
        StsDirectedEdge directedEdge;

        try
        {
            StsCursorSection cursorSection;
            if(rowOrCol == ROW)
                cursorSection = currentModel.getCursor3d().getRowCursorSection();
            else           
                cursorSection = currentModel.getCursor3d().getColCursorSection();
            
            cursorEdgeLoop = new StsEdgeLoop(cursorSection);
            cursorEdgeLoop.setNormal(normal, true);

            directedEdge = new StsDirectedEdge(firstGridLine, MINUS);
            cursorEdgeLoop.addDirectedEdge(directedEdge);

            while( ( directedEdge = getNextDirectedEdge(directedEdge, cursorGridLines) ) != null)
            {
                cursorEdgeLoop.addDirectedEdge(directedEdge);
            }
        /*
            if(!cursorEdgeLoop.isClosed())
            {
//                StsException.systemError("StsZoneBlock.constructCursorEdgeLoop() failed. " +
//                    cursorEdgeLoop.getLabel() + " is not closed.");
                return null;
            }
        */
            if(!cursorEdgeLoop.constructEdgeLinks()) return null;
            if(!cursorEdgeLoop.constructCursorSection(rowOrCol, rowCol, topGrid)) return null;
            return cursorEdgeLoop;
        }
        catch(Exception e)
        {
            StsException.outputException("StsZoneBlock.constructCursorEdgeLoop() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    private StsDirectedEdge getNextDirectedEdge(StsDirectedEdge directedEdge, StsList edges)
    {
        StsEdgeLinkable edge;
        StsPoint prevPoint;

        edge = directedEdge.getEdge();
        int direction = directedEdge.getDirection();
        if(direction == PLUS)
            prevPoint = edge.getLastGridPoint().getPoint();
        else
            prevPoint = edge.getFirstGridPoint().getPoint();

        int nEdges = edges.getSize();
        for(int n = 0; n < nEdges; n++)
        {
            edge = (StsEdgeLinkable)edges.getElement(n);
            if(edge.getFirstGridPoint().getPoint().sameAs(prevPoint))
            {
                edges.delete(edge);
                return new StsDirectedEdge(edge, PLUS);
            }
            if(edge.getLastGridPoint().getPoint() ==  prevPoint)
            {
                edges.delete(edge);
                return new StsDirectedEdge(edge, MINUS);
            }
        }
        return null;
    }

    private void drawCursorGridLines(StsList gridLines, GL gl)
    {
        try
        {
            if(gridLines == null) return;

            zone.getStsColor().setGLColor(gl);

            int nGridLines = gridLines.getSize();
            for(int n = 0; n < nGridLines; n++)
            {
                StsEdgeLinkable gridLine = (StsEdgeLinkable)gridLines.getElement(n);
                float[][] xyzPoints = gridLine.getXYZPoints();
                StsGLDraw.drawSegmentedLine(gl, xyzPoints);
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsZoneBlock.drawCursorGridLines() failed.",
                e, StsException.WARNING);
        }
    }
/*
    private void drawCursorSection(StsGLPanel3d glPanel3d, StsLinkGridLine botGridLine, StsLinkGridLine topGridLine,
                    float startF, float endF, GL gl, float[] sectionNormal)
    {
        int lowerBound = Math.max(botGridLine.getLowerBound(), topGridLine.getLowerBound());
        int upperBound = Math.min(botGridLine.getUpperBound(), topGridLine.getUpperBound());

        float[][] topXYZs = topGrid.getLineXYZs(topGridLine, startF, endF, lowerBound, upperBound);
        if(topXYZs == null) return;
        float[][] botXYZs = botGrid.getLineXYZs(botGridLine, startF, endF, lowerBound, upperBound);
        if(botXYZs == null) return;

        subZoneVerticalGridZs = null;
        float[][] lineRowColFs = null;

        boolean displaySubZoneColors = currentModel.getDisplaySubZoneColors();
        if(zone.getDisplayFill())
        {
            lineRowColFs = topGrid.getLineRowColFs(topGridLine, startF, endF);
            StsPropertyVolumeOld pv = zone.getPropertyVolume();
            if(pv != null && pv.scaledColorsAreValid())
                drawCursorSectionSurface(botXYZs, topXYZs, lineRowColFs, gl, sectionNormal, pv);
            else
                drawCursorSectionSurface(botXYZs, topXYZs, lineRowColFs, gl, sectionNormal, displaySubZoneColors);
        }
        if(zone.getDisplayGrid())
        {
            if(lineRowColFs == null) lineRowColFs = topGrid.getLineRowColFs(topGridLine, startF, endF);
            drawCursorSectionGrid(glPanel3d, botXYZs, topXYZs, lineRowColFs, gl);
        }

    }
*/
    private void drawCursorSectionSurface(float[][] botXYZs, float[][] topXYZs,
                    float[][] rowColFs, GL gl, float[] sectionNormal, boolean displaySubZoneColors)
    {
        try
        {
            if(!displaySubZoneColors)
            {
                zone.getStsColor().setGLColor(gl);
                gl.glBegin(GL.GL_QUAD_STRIP);

                int nPnts = topXYZs.length;
                for(int n = 0; n < nPnts; n++)
                {
                    gl.glNormal3fv(sectionNormal, 0);
                    gl.glVertex3fv(topXYZs[n], 0);
                    gl.glVertex3fv(botXYZs[n], 0);
                }
                gl.glEnd();
            }
            else
            {
                if(subZoneVerticalGridZs == null) computeSubZoneVerticalGridZs(rowColFs);
                StsObjectRefList subZones = zone.getSubZones();
                int nSubZones = subZones.getSize();
                for(int sz = 0; sz < nSubZones; sz++)
                {
                    StsSubZone subZone = (StsSubZone)subZones.getElement(sz);
                    StsColor color = subZone.getStsColor();
		            color.setGLColor(gl);

                    float[][] gridZ = subZoneVerticalGridZs[sz];
                    if(gridZ == null)
                    if(gridZ == null) continue;
                    int nSubZoneLayers = gridZ.length-1;
                    int nPnts = gridZ[0].length;

                    for(int k = 0; k < nSubZoneLayers; k++)
                    {
                        gl.glBegin(GL.GL_QUAD_STRIP);
                        for(int n = 0; n < nPnts; n++)
                        {
                            float[] xyz = botXYZs[n];
                            gl.glNormal3fv(sectionNormal, 0);
                            xyz[2] = gridZ[k][n];
                            gl.glVertex3fv(xyz, 0);
                            xyz[2] = gridZ[k+1][n];
                            gl.glVertex3fv(xyz, 0);
                        }
                        gl.glEnd();
                    }
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsZoneBlock.drawCursorSectionSurface() failed.",
                e, StsException.WARNING);
            gl.glEnd();
        }
    }

    private void drawCursorSectionSurface(float[][] botXYZs, float[][] topXYZs,
                    float[][] rowColFs, GL gl, float[] sectionNormal, StsPropertyVolumeOld pv)
    {
        try
        {
            computeSubZoneVerticalGridZs(rowColFs);
            if(subZoneVerticalGridZs == null) return;

            int nSubZones = zone.getNSubZones();

            float[][] tXYZs = StsMath.copyFloatArray(topXYZs);
            float[][] bXYZs = StsMath.copyFloatArray(botXYZs);

            int nZoneLayer = 0;
            for(int sz = 0; sz < nSubZones; sz++)
            {
                float[][] gridZ = subZoneVerticalGridZs[sz];
                if(gridZ == null) continue;
                int nSubZoneLayers = gridZ.length-1;
                int nPnts = gridZ[0].length;
                float[][] values = pv.getValuesArray(rowColFs, nZoneLayer, nSubZoneLayers);

                for(int k = 0; k < nSubZoneLayers; k++, nZoneLayer++)
                {
                    for(int n = 0; n < nPnts; n++)
                    {
                        tXYZs[n][2] = gridZ[k][n];
                        bXYZs[n][2] = gridZ[k+1][n];
                    }

                    gl.glBegin(GL.GL_QUADS);
                    for(int n = 0; n < nPnts-1; n++)
                    {
                        pv.getScaledColor(values[n][k]).setGLColor(gl);
                        gl.glNormal3fv(sectionNormal, 0);
                        gl.glVertex3fv(tXYZs[n], 0);
                        gl.glVertex3fv(bXYZs[n], 0);
                        gl.glVertex3fv(bXYZs[n+1], 0);
                        gl.glVertex3fv(tXYZs[n+1], 0);
                    }
                    gl.glEnd();
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsBlockSide.displayVerticalZoneGrid(...StsPropertyVolumeOld) failed.",
                                         e, StsException.WARNING);
            gl.glEnd();
        }
    }
/*
    private void drawCursorSection(StsGLPanel3d glPanel3d, StsLinkGridLine botGridLine, StsLinkGridLine topGridLine,
                    float startF, float endF, GL gl, float[] sectionNormal)
    {
        int lowerBound = Math.max(botGridLine.getLowerBound(), topGridLine.getLowerBound());
        int upperBound = Math.min(botGridLine.getUpperBound(), topGridLine.getUpperBound());

        float[][] topXYZs = topGrid.getLineXYZs(topGridLine, startF, endF, lowerBound, upperBound);
        if(topXYZs == null) return;
        float[][] botXYZs = botGrid.getLineXYZs(botGridLine, startF, endF, lowerBound, upperBound);
        if(botXYZs == null) return;

        subZoneVerticalGridZs = null;
        float[][] lineRowColFs = null;

        boolean displaySubZoneColors = currentModel.getDisplaySubZoneColors();
        if(zone.getDisplayFill())
        {
            lineRowColFs = topGrid.getLineRowColFs(topGridLine, startF, endF);
            StsPropertyVolumeOld pv = zone.getPropertyVolume();
            if(pv != null && pv.scaledColorsAreValid())
                drawCursorSectionSurface(botXYZs, topXYZs, lineRowColFs, gl, sectionNormal, pv);
            else
                drawCursorSectionSurface(botXYZs, topXYZs, lineRowColFs, gl, sectionNormal, displaySubZoneColors);
        }
        if(zone.getDisplayGrid())
        {
            if(lineRowColFs == null) lineRowColFs = topGrid.getLineRowColFs(topGridLine, startF, endF);
            drawCursorSectionGrid(glPanel3d, botXYZs, topXYZs, lineRowColFs, gl);
        }

    }
*/

    private void drawCursorSectionGrid(StsGLPanel3d glPanel3d,
                        float[][] botXYZs, float[][] topXYZs, float[][] rowColFs, GL gl)
    {
        try
        {
            computeSubZoneVerticalGridZs(rowColFs);

            gl.glDisable(GL.GL_LIGHTING);
            gl.glLineWidth(StsGraphicParameters.gridLineWidth);
            glPanel3d.setViewShift(gl, StsGraphicParameters.gridShift);
            StsColor drawColor = zone.getDisplayFill() ? StsColor.BLACK : zone.getStsColor();
            drawColor.setGLColor(gl);

            int nPnts = topXYZs.length;

            gl.glBegin(GL.GL_LINES);
            // draw lines down
            for(int n = 0; n < nPnts; n++)
            {
                gl.glVertex3fv(botXYZs[n], 0);
                gl.glVertex3fv(topXYZs[n], 0);
            }
            gl.glEnd();

            float[][] XYZs = StsMath.copyFloatArray(botXYZs);

            for(int sz = 0; sz < nSubZones; sz++)
            {
                float[][] gridZ = subZoneVerticalGridZs[sz];
                int nLines = gridZ.length;
                if(sz < nSubZones-1) nLines--;

                for(int k = 0; k < nLines; k++)
                {
                    gl.glBegin(GL.GL_LINE_STRIP);
                    for(int n = 0; n < nPnts; n++)
                    {
                        float[] xyz = XYZs[n];
                        xyz[2] = gridZ[k][n];
                        gl.glVertex3fv(xyz, 0);
                    }
                    gl.glEnd();
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsZoneBlock.drawCursorSectionGrid() failed.",
                e, StsException.WARNING);
            gl.glEnd();
        }
        finally
        {
            gl.glEnable(GL.GL_LIGHTING);
            glPanel3d.resetViewShift(gl);
        }
    }

    private void computeSubZoneVerticalGridZs(float[][] subZoneGridZ, int gridCol, int nLayers,
                                 int subZoneType, float thickness)
    {
        int n, k;

        try
        {
            float topZ = subZoneGridZ[0][gridCol];
            float botZ = subZoneGridZ[nLayers][gridCol];
            float dZ;
            switch(subZoneType)
            {
                case StsZone.SUBZONE_UNIFORM:
                    float df = 1.0f/nLayers;
                    float f = df;
                    for(k = 1; k < nLayers; k++)
                    {
                        subZoneGridZ[k][gridCol] = topZ + f*(botZ - topZ);
                        f += df;
                    }
                    break;
                case StsZone.SUBZONE_ONLAP:
                    dZ = thickness;
                    for(k = 1; k < nLayers; k++)
                    {
                        subZoneGridZ[k][gridCol] = Math.min(topZ + dZ, botZ);
                        dZ += thickness;
                    }
                    break;
                case StsZone.SUBZONE_OFFLAP:
                    dZ = thickness;
                    for(k = nLayers-1; k >= 1; k--)
                    {
                        subZoneGridZ[k][gridCol] = Math.max(botZ - dZ, topZ);
                        dZ += thickness;
                    }
                    break;
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsZoneBlock.computeSubZoneVerticalGridZs() failed.",
                e, StsException.WARNING);
        }

    }

    private void computeSubZoneVerticalGridZs(float[][] rowColFs)
    {
        try
        {
            if(subZoneVerticalGridZs != null) return;

            StsObjectRefList subZones = zone.getSubZones();
            int nSubZones = subZones.getSize();
            subZoneVerticalGridZs = new float[nSubZones][][];
            for(int sz = 0; sz < nSubZones; sz++)
            {
                StsSubZone subZone = (StsSubZone)subZones.getElement(sz);
                computeSubZoneVerticalGridZ(subZone, rowColFs);
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsZoneBlock.computeSubZoneVerticalGridZs() failed.",
                e, StsException.WARNING);
        }
    }

    private void computeSubZoneVerticalGridZ(StsSubZone subZone, float[][] rowColFs)
    {
        float[] botZs, topZs;

        try
        {
            int subZoneNumber = subZone.getSubZoneNumber();
            int nLayers = subZone.getNLayers();
            int subZoneType = subZone.getSubZoneType();
            float thickness = subZone.getLayerThickness();
            int nPnts = rowColFs.length;
            float[][] subZoneGridZ = new float[nLayers+1][nPnts];
            subZoneVerticalGridZs[subZoneNumber] = subZoneGridZ;
            for(int n = 0; n < nPnts; n++)
            {
                float rowF = rowColFs[n][0];
                float colF = rowColFs[n][1];

                subZoneGridZ[0][n] = getZoneGridZ(subZoneNumber, rowF, colF);
                subZoneGridZ[nLayers][n] = getZoneGridZ(subZoneNumber+1, rowF, colF);

                computeSubZoneVerticalGridZs(subZoneGridZ, n, nLayers, subZoneType, thickness);
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsZoneBlock.computeSubZoneVerticalGridZs() failed.",
                e, StsException.WARNING);
        }
    }

    public boolean delete()
    {
        super.delete();
        if(zoneBlockSides != null)
        {
            zoneBlockSides.deleteAll();
            zoneBlockSides = null;
        }
        return true;
    }

    public String getName() { return new String("ZoneBlock-" + getIndex()); }

    public StsBlock getBlock()
    {
        return block;
    }

    public StsList getZoneBlockSides()
    {
        return zoneBlockSides;
    }

    public boolean hasLayer(int nLayer)
    {
        return zone.hasLayer(nLayer);
    }

    public boolean hasLayerGrid(int nLayerGrid)
    {
        return zone.hasLayerGrid(nLayerGrid);
    }

    public float[][] getLayerGrid(int nLayerGrid)
    {
        int zoneTopLayerNumber = zone.topLayerNumber;
        int nZoneLayerGrid = nLayerGrid - zoneTopLayerNumber;
        int[] subZoneBotLayerGridNumbers = zone.subZoneBotLayerGridNumbers;
        int nSubZones = subZoneBotLayerGridNumbers.length;
        int subZoneBotLayerGridNumber = 0;
        for(int n = 0; n < nSubZones; n++)
        {
            int subZoneTopLayerGridNumber = subZoneBotLayerGridNumber;
            subZoneBotLayerGridNumber = subZoneBotLayerGridNumbers[n];
            if(nZoneLayerGrid >= subZoneTopLayerGridNumber && nZoneLayerGrid <= subZoneBotLayerGridNumber)
            {
                int subZoneLayerGridNumber = nZoneLayerGrid - subZoneTopLayerGridNumber;
                return this.layerGridZs[n][subZoneLayerGridNumber];
            }
        }
        StsException.systemError(this, "getLayerGrid", "Failed to find layer grid " + nLayerGrid + " in " + toDetailString());
        return null;
    }

    public ArrayList<StsPropertyVolume> getZoneProperties()
    {
        if(zoneProperties == null)
            zoneProperties = new ArrayList<StsPropertyVolume>();
        return zoneProperties;
    }

    public void propertyTypeChanged()
    {
        int nZoneBlockSides = zoneBlockSides.getSize();
        for(int n = 0; n < nZoneBlockSides; n++)
        {
            StsZoneSide zoneSide = (StsZoneSide)zoneBlockSides.getElement(n);
            zoneSide.propertyTypeChanged();
        }
        topGrid.propertyTypeChanged();
        botGrid.propertyTypeChanged();
        propertyTypeChanged(cursorRowEdgeLoops);
        propertyTypeChanged(cursorColEdgeLoops);
    }

    private void propertyTypeChanged(StsList cursorEdgeLoops)
    {
        int nLoops = cursorEdgeLoops.getSize();
        for(int n = 0; n < nLoops; n++)
        {
            StsEdgeLoop edgeLoop = (StsEdgeLoop)cursorEdgeLoops.getElement(n);
            edgeLoop.propertyTypeChanged();
        }
    }
}

