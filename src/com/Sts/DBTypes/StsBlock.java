//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Types.Simulation.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;
import com.Sts.Actions.Export.*;

import java.util.*;

public class StsBlock extends StsGridBoundingBox implements StsTreeObjectI, StsXYGridable, StsGrid2dIterable<StsBlock.BlockCellColumn>
{
    protected StsObjectRefList blockSides;
    protected StsObjectRefList connectedBlocks;  // blocks connected across an auxiliary section
    protected byte blockType = NONE;
    protected boolean hasAuxiliarySections = false;
    /** @param hasAuxiliarySections indicates blockSide is an auxiliary section */

    transient protected StsList blockGrids;
    transient protected StsList zoneBlocks;
    transient private StsPropertyVolume currentPropertyVolume;
    // transient public ArrayList<StsPropertyVolume> propertyVolumes;
    transient protected StsPropertyType currentPropertyType;
    /** @param type type may be LOOP or OPEN_LOOP * */


    /**
     * For simulation output, we create a grid of blockCells.  Each block cell is FULL, EMPTY, or TRUNCATED.  If a blockCell is TRUNCATED,
     * it contains a set of polygons which are the faces of the truncated portion of the cell.
     * The number of cell rows and columns is one less in each direction than the grid rows and columns for the block.
     */
    transient private BlockCellColumn[][] blockCellColumns;
    transient public int nCellRows;
    transient public int nCellCols;
    transient public int nTotalLayers;
    transient public float actualMergedCellMinFraction = 1.0f;
    transient public int numberMergedCells = 0;
    transient public int numberUnmergedSmallCells = 0;
    transient public ArrayList<StsNNC> neighborNncList = new ArrayList<StsNNC>();
    transient public TreeSet<BlockCellColumn.GridCell> truncatedCells;
    transient public ArrayList<BlockCellColumn.ParentCell> parentCells;
    transient public int[] indexMap;
    /**
     * simulation rowMin.  Blocks are laid out in a sequence with each block occupying a set of rows with a null row between blocks.
     * So the first block will be from 0 to block[0].nRows-1.  We then skip a row.  The second block is from row block[0].nRows + 1 to block[0].nRows + block[1].nRows.
     */
    transient public int eclipseRowMin;
    transient public double cellXYArea;

    static public StsGridBoundingBox modelBoundingBox;
    static protected boolean hasModelAuxiliarySections = false;

    static int nXPlus, nXMinus, nYPlus, nYMinus;

    static public float minPoreVolume = 0.0f;

    static StsObjectPanel objectPanel = null;
    static public StsFieldBean[] displayFields = null;
    static public StsFieldBean[] propertyFields = null;

    static public final byte NONE = StsParameters.NONE;
    static public final byte LOOP = StsParameters.LOOP;
    static public final byte OPEN_LOOP = StsParameters.OPEN_LOOP;

    static public final byte EMPTY = StsParameters.CELL_EMPTY;
    static public final byte FULL = StsParameters.CELL_FULL;
    static public final byte EDGE = StsParameters.CELL_EDGE;
    static public final String[] cellTypeStrings = new String[]{"EMPTY", "FULL", "EDGE"};

    static public final byte INSIDE_ROW = 3;

    static final int LEFT = StsParameters.LEFT;
    static final int RIGHT = StsParameters.RIGHT;

    static final int ABOVE = -1;
    static final int BETWEEN = 0;
    static final int BELOW = 1;

    static public final byte CELL_NONE = -1;
    static public final byte CELL_PLUS_DX = 0;
    static public final byte CELL_MINUS_DX = 1;
    static public final byte CELL_PLUS_DY = 2;
    static public final byte CELL_MINUS_DY = 3;
    static public final byte CELL_PLUS_DZ = 4;
    static public final byte CELL_MINUS_DZ = 5;

    static final float[] NORMAL_PLUS_DX = new float[]{1.0f, 0.0f, 0.0f};
    static final float[] NORMAL_MINUS_DX = new float[]{-1.0f, 0.0f, 0.0f};
    static final float[] NORMAL_PLUS_DY = new float[]{0.0f, 1.0f, 0.0f};
    static final float[] NORMAL_MINUS_DY = new float[]{0.0f, -1.0f, 0.0f};
    static final float[] NORMAL_PLUS_DZ = new float[]{0.0f, 0.0f, 1.0f};
    static final float[] NORMAL_MINUS_DZ = new float[]{0.0f, 0.0f, -1.0f};
    static float[][] cellOrthogonalNormals = new float[][]{NORMAL_PLUS_DX, NORMAL_MINUS_DX, NORMAL_PLUS_DY, NORMAL_MINUS_DY, NORMAL_PLUS_DZ, NORMAL_MINUS_DZ};
    static final String[] cellDirectionStrings = new String[]{"+X", "-X", "+Y", "-Y", "+Z", "-Z"};
    static final String[] cellParameterStrings = new String[]{"dX+/2", "dY+/2", "dZ+/2", "dX-/2", "dY-/2", "dZ-/2"};

    static final float largeFloat = StsParameters.largeFloat;

    static float blockLineSortMinRowF;
    static float blockLineSortMinColF;

    /** dRow, dCol index changes to neighbors CELL_PLUS_DX, CELL_MINUS_DX, CELL_PLUS_DY, CELL_MINUS_DY */
    static final int[][] neighborDRowCol = new int[][]{{0, 1}, {0, -1}, {1, 0}, {-1, 0}};

    static final boolean debugCell = true;
    static int debugCellRow = 67;
    static int debugCellCol = 0;
    static int debugCellLayer = 3;
    static int debugCellBlockIndex = 1;
    static byte debugDirection = CELL_PLUS_DX;
    // display fields
    static public final StsFieldBean[] fieldBeans =
        {
            new StsBooleanFieldBean(StsBlock.class, "isVisible", "Enable"),
        };

    public StsBlock()
    {
//        setName("Block-" + getIndex());
    }

    static StsBlock constructBlock()
    {
        try
        {
            StsBlock block = new StsBlock();
            block.setName("Block-" + block.getIndex());
            //			block.addToModel();
            //			block.refreshObjectPanel();
            return block;
        }
        catch (Exception e)
        {
            StsException.systemError("StsBlock.constructor() failed.");
            return null;
        }
    }

    public boolean initialize(StsModel model)
    {
        nCellRows = getNRows() - 1;
        nCellCols = getNCols() - 1;
        return true;
    }

    public String getName()
    {
        return super.getName();
    }

    public void setIsVisible(boolean b)
    {
        isVisible = b;
    }

    public StsObjectRefList getBlockSides()
    {
        checkInitializeBlockSides();
        return blockSides;
    }

    public boolean hasAuxiliarySections()
    {
        return hasAuxiliarySections;
    }

    public StsObjectRefList getConnectedBlocks()
    {
        return connectedBlocks;
    }

    static public boolean hasModelAuxiliarySections()
    {
        return hasModelAuxiliarySections;
    }

    static public boolean constructBlocks(StsModel model)
    {
        StsClass sections, blocks;
        int nBlocks;
        boolean constructionOK = true;
        StsLine line;
        boolean isFullyConnected;
        int n;

        try
        {
            // check if blocks already built
            /*
            blocks = model.getCreateStsClass(StsBlock.class);
            nBlocks = blocks.getSize();
            if(nBlocks > 0) return true;
            */
            blocks = model.getCreateStsClass(StsBlock.class);
            blocks.deleteAll();
            model.getCreateStsClass(StsBlockSide.class).deleteAll();
            // model.getCreateStsClass(StsBlockGrid.class).deleteAll();
            sections = model.getCreateStsClass(StsSection.class);
            int nSections = sections.getSize();
            for (n = 0; n < nSections; n++)
            {
                StsSection section = (StsSection) sections.getElement(n);
                line = section.getFirstLine();
                isFullyConnected = line.isFullyConnected();
                if (!isFullyConnected)
                {
                    StsException.systemError("StsBlock.constructBlocks() failed." +
                        section.getLabel() + " does not form a closed loop." +
                        " delete section or complete it. ");
                }
                constructionOK = constructionOK && isFullyConnected;

                line = section.getLastLine();
                isFullyConnected = line.isFullyConnected();
                if (!isFullyConnected)
                {
                    StsException.systemError("StsBlock.constructBlocks() failed." +
                        section.getLabel() + " does not form a closed loop." +
                        " delete section or complete it. ");
                }
                constructionOK = constructionOK && isFullyConnected;
            }
            if (!constructionOK) return false;


            for (n = 0; n < nSections; n++)
            {
                StsSection section = (StsSection) sections.getElement(n);
                section.constructBlockSides();
            }

            StsClass blockSides = model.getCreateStsClass(StsBlockSide.class);
            int nBlockSides = blockSides.getSize();
            for (n = 0; n < nBlockSides; n++)
            {
                StsBlockSide blockSide = (StsBlockSide) blockSides.getElement(n);
                if (blockSide.getBlock() == null)
                    constructionOK = constructionOK && constructBlock(blockSide, model);
            }

            // remove incomplete blocks
            nBlocks = blocks.getSize();
            if (!constructionOK)
            {
                for (n = nBlocks - 1; n >= 0; n--)
                {
                    StsBlock block = (StsBlock) blocks.getElement(n);
                    if (block.getType() == OPEN_LOOP) blocks.delete(block);
                }

                if (!StsYesNoDialog.questionValue(model.win3d,
                    "Do you wish to continue construction with bad blocks removed?"))
                    return false;
            }

            nBlocks = blocks.getSize();
            for (n = 0; n < nBlocks; n++)
            {
                StsBlock block = (StsBlock) blocks.getElement(n);
                block.addConnectedBlocks();
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsBlock.constructBlocks() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    static private boolean constructBlock(StsBlockSide side, StsModel model)
    {
        StsBlockSide firstSide, prevSide, nextSide;

        StsBlock block = StsBlock.constructBlock();

        /** Back-up until we find the same edge, an edge already assigned to a
         *  prefix, or we terminate for some reason (bad edge prefix).
         */

        firstSide = side;
        /*
            prevSide = (StsBlockSide)side.getPrevSide();
            while(prevSide != null && prevSide != side && prevSide.getBlock() == null)
            {
                firstSide = prevSide;
                prevSide = (StsBlockSide)prevSide.getPrevBlockSide();
            }
        */
        block.addBlockSide(firstSide);

        nextSide = (StsBlockSide) firstSide.getNextBlockSide();
        while (nextSide != null && nextSide != firstSide && nextSide.getBlock() == null)
        {
            StsBlockSide otherSide = block.overlapsSideOnOtherSide(nextSide);
            if (otherSide != null)
            {

                new StsMessage(model.win3d, StsMessage.WARNING,
                    "Two sides from same block are on opposite sides of same section:\n" +
                        " First  Side: " + nextSide.getLabel() + ".\n" +
                        " Second Side: " + otherSide.getLabel() + ".\n" +
                        "Insert an auxiliary section to split this block.");
            }

            block.addBlockSide(nextSide);
            nextSide = (StsBlockSide) nextSide.getNextBlockSide();
        }

        if (nextSide != null && nextSide == firstSide)
        {
            block.setType(LOOP);
            block.checkForAuxiliarySections();
            return true;
        }
        else
        {
            nextSide = (StsBlockSide) block.getBlockSides().getLast();
            new StsMessage(model.win3d, StsMessage.WARNING,
                "Sides for this block do not form a closed loop.\n" +
                    "Insert an auxiliary section to join:\n" +
                    " First  Side: " + firstSide.getLabel() + ".\n" +
                    " Last Side: " + nextSide.getLabel() + ".\n");
            block.setType(OPEN_LOOP);
            return false;
        }
    }

    private void addBlockSide(StsBlockSide blockSide)
    {
        checkInitializeBlockSides();
        blockSides.add(blockSide);
        blockSide.setBlock(this);
    }

    private void checkInitializeBlockSides()
    {
        if (blockSides != null) return;
        blockSides = StsObjectRefList.constructor(4, 2, "blockSides", this);

    }

    public StsBlockSide overlapsSideOnOtherSide(StsBlockSide side)
    {
        StsSection section = side.getSection();
        int sectionSide = side.getSide();

        int nBlockSides = blockSides.getSize();

        for (int n = 0; n < nBlockSides; n++)
        {
            StsBlockSide blockSide = (StsBlockSide) blockSides.getElement(n);
            if (blockSide.getSection() == section && sectionSide != blockSide.getSide()) return blockSide;
        }
        return null;
    }

    public void checkForAuxiliarySections()
    {
        int nSides = blockSides.getSize();
        for (int n = 0; n < nSides; n++)
        {
            StsBlockSide side = (StsBlockSide) blockSides.getElement(n);
            if (side.getSection().isAuxiliary())
            {
                hasAuxiliarySections = true;
                hasModelAuxiliarySections = true;
                return;
            }
        }
    }

    // blockGrids are in top down sequence; when we add a new one, set references
    // with blockGrid above
    public void addBlockGrid(StsBlockGrid blockGrid)
    {
        if (blockGrids == null) blockGrids = new StsList(4, 2);
        StsBlockGrid blockGridAbove = (StsBlockGrid) blockGrids.getLast();

        blockGrids.add(blockGrid);
        if (blockGridAbove == null) return;
        blockGrid.setBlockGridAbove(blockGridAbove);
        blockGridAbove.setBlockGridBelow(blockGrid);
    }

    public StsBlockGrid getBlockGrid(StsModelSurface surface)
    {
        if (blockGrids == null) return null;
        for (int n = 0; n < blockGrids.getSize(); n++)
        {
            StsBlockGrid blockGrid = (StsBlockGrid) blockGrids.getElement(n);
            if (blockGrid.getSurface() == surface) return blockGrid;
        }
        return null;
    }

    public StsBlockGrid[] getBlockGrids()
    {
        return (StsBlockGrid[]) blockGrids.getCastList(StsBlockGrid.class);
    }

    static public StsGridBoundingBox computeModelBoundingBox(StsBlock[] blocks)
    {
        try
        {
            modelBoundingBox = new StsGridBoundingBox(false);
            for (StsBlock block : blocks)
            {
                block.computeBoundingBox();
                modelBoundingBox.addBoundingBox(block);
            }
            StsProject project = currentModel.getProject();
            StsRotatedGridBoundingBox projectBoundingBox = project.getRotatedBoundingBox();
            modelBoundingBox.xInc = projectBoundingBox.xInc;
            modelBoundingBox.yInc = projectBoundingBox.yInc;
            modelBoundingBox.zMin = project.getZorTMin(isDepth);
            modelBoundingBox.zMax = project.getZorTMax(isDepth);
            modelBoundingBox.zInc = project.getZorTInc(isDepth);
            modelBoundingBox.computeSliceRange();
            /*
            for(StsBlock block : blocks)
                block.loopBoundingBox.addBorder(5, currentModel.getGridDefinition());
            */
            return modelBoundingBox;
        }
        catch (Exception e)
        {
            StsException.outputException("StsBlock.computeBlockBoundingBoxes() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    public void computeBoundingBox()
    {
        if (blockGrids == null) return;

        int nBlockGrids = blockGrids.getSize();
        for (int n = 0; n < nBlockGrids; n++)
        {
            StsBlockGrid blockGrid = (StsBlockGrid) blockGrids.getElement(n);
            addBoundingBox(blockGrid.getGridBoundingBox());
        }
        nCellRows = getNRows() - 1;
        nCellCols = getNCols() - 1;
        //        boundingBox.addBorder(5, currentModel.getGridDefinition());
    }

    public float getBlockXCoor(int blockCol)
    {
        blockCol = StsMath.minMax(blockCol, 0, getNCols() - 1);
        return blockCol * xInc + xMin;
    }

    public float getBlockYCoor(int blockRow)
    {
        blockRow = StsMath.minMax(blockRow, 0, getNRows() - 1);
        return blockRow * yInc + yMin;
    }

    public int getBlockRowMin()
    {
        return rowMin;
    }

    public int getBlockRowMax()
    {
        return rowMax;
    }

    public int getBlockColMin()
    {
        return colMin;
    }

    public int getBlockColMax()
    {
        return colMax;
    }

    public String getLabel()
    {
        return new String("block-" + getIndex());
    }

    private void addConnectedBlocks()
    {
        if (blockSides == null) return;
        int nBlockSides = blockSides.getSize();
        for (int n = 0; n < nBlockSides; n++)
        {
            StsBlockSide blockSide = (StsBlockSide) blockSides.getElement(n);
            addConnectedBlocks(blockSide);
        }
    }

    private void addConnectedBlocks(StsBlockSide blockSide)
    {
        StsObjectRefList sectionBlockSides;

        StsSection section = blockSide.getSection();
        if (!section.isAuxiliary()) return;
        if (connectedBlocks == null) connectedBlocks = StsObjectRefList.constructor(2, 1, "connectedBlocks", this);
        sectionBlockSides = section.getLeftBlockSides();
        addConnectedBlocks(blockSide, sectionBlockSides);
        sectionBlockSides = section.getRightBlockSides();
        addConnectedBlocks(blockSide, sectionBlockSides);
    }

    private void addConnectedBlocks(StsBlockSide insideBlockSide, StsObjectRefList otherBlockSides)
    {
        int nOtherBlockSides = otherBlockSides.getSize();
        for (int n = 0; n < nOtherBlockSides; n++)
        {
            StsBlockSide otherBlockSide = (StsBlockSide) otherBlockSides.getElement(n);
            if (otherBlockSide != insideBlockSide)
            {
                StsBlock otherBlock = otherBlockSide.getBlock();
                if (!connectedBlocks.contains(otherBlock))
                    connectedBlocks.add(otherBlock);
            }
        }
    }

    public void addBlockZone(StsZoneBlock zoneBlock)
    {
        if (zoneBlocks == null) zoneBlocks = new StsList(4, 2);
        zoneBlocks.add(zoneBlock);
    }

    public StsZoneBlock[] getZoneBlocks()
    {
        if(zoneBlocks == null) return null;
        return (StsZoneBlock[]) zoneBlocks.getCastList(StsZoneBlock.class);
    }

    public StsZoneBlock getZoneBlock(int nLayer)
    {
        int nZoneBlocks = zoneBlocks.getSize();
        for (int n = 0; n < nZoneBlocks; n++)
        {
            StsZoneBlock zoneBlock = (StsZoneBlock) zoneBlocks.getElement(n);
            if (zoneBlock.hasLayer(nLayer))
                return zoneBlock;
        }
        return null;
    }

    public StsZoneBlock getZoneBlock(StsZone zone)
    {
        int nZoneBlocks = zoneBlocks.getSize();
        for (int n = 0; n < nZoneBlocks; n++)
        {
            StsZoneBlock zoneBlock = (StsZoneBlock) zoneBlocks.getElement(n);
            if (zoneBlock.zone == zone)
                return zoneBlock;
        }
        return null;
    }

    /*
        public void clearProperties()
        {
            propertyVolumes = null;
        }
    */
    public void constructBasicProperties(StsZone[] zones, StsEclipseModel eclipseModel)
    {
        constructBasicProperty(eclipseModel.porosity, zones, eclipseModel);
        constructBasicProperty(eclipseModel.permX, zones, eclipseModel);
        constructBasicProperty(eclipseModel.permY, zones, eclipseModel);
        constructBasicProperty(eclipseModel.permZ, zones, eclipseModel);
    }

    private void constructBasicProperty(StsPropertyType propertyType, StsZone[] zones, StsEclipseModel eclipseModel)
    {
        StsPropertyVolume blockPropertyVolume = eclipseModel.getBlockPropertyVolume(this, propertyType);
        if (blockPropertyVolume != null) return;

        byte distributionType = getZonesPropertyType(propertyType, zones);
        blockPropertyVolume = eclipseModel.createBlockPropertyVolume(this, propertyType, distributionType);
        int topLayerNumber = 0;
        StsPropertyVolume zonePropertyVolume;
        for (StsZone zone : zones)
        {
            StsZoneBlock zoneBlock = getZoneBlock(zone);
            zonePropertyVolume = zoneBlock.getPropertyVolume(propertyType.name);
            if (zonePropertyVolume == null)
                zonePropertyVolume = zone.getPropertyVolume(propertyType.name);
            if (zonePropertyVolume == null)
            {
                StsException.systemError(this, "constructBasicProperty", "Couldn't find zonePropertyVolume for this block " + toDetailString() + " and zone " + zone.toDetailString());
            }
            int nLayers = zone.getNLayers();
            blockPropertyVolume.setValues(zonePropertyVolume, topLayerNumber, nLayers);
            topLayerNumber += nLayers;
        }
    }

    private byte getZonesPropertyType(StsPropertyType propertyType, StsZone[] zones)
    {
        int nZones = zones.length;
        byte distributionType = zones[0].getPropertyVolume(propertyType.name).getDistributionType();
        for (int n = 1; n < nZones; n++)
        {
            StsZone zone = zones[n];
            byte zoneType = zone.getPropertyVolume(propertyType.name).getDistributionType();
            distributionType = StsPropertyVolume.getDistributionType(distributionType, zoneType);
        }
        return distributionType;
    }

    public void constructCellColumns(int nTotalLayers)
    {
        int nRows = getNRows();
        int nCols = getNCols();
        nCellRows = nRows - 1;
        nCellCols = nCols - 1;
        int nZoneBlocks = zoneBlocks.getSize();
        for (int n = 0; n < nZoneBlocks; n++)
        {
            StsZoneBlock zoneBlock = (StsZoneBlock) zoneBlocks.getElement(n);
            zoneBlock.constructGrids();
        }
        blockCellColumns = new BlockCellColumn[nCellRows][nCellCols];
    }

    public void clearTransientArrays()
    {
        blockCellColumns = null;
        neighborNncList = null;
        truncatedCells = null;
        parentCells = null;
        indexMap = null;
        actualMergedCellMinFraction = 1.0f;
        numberMergedCells = 0;
        numberUnmergedSmallCells = 0;
    }

    /**
     * These polygons are on blockZoneSides and define the intersections with grid row and column vertical planes.
     * These polygons are subdivided into layerPolygons.  Determine the layer range and add this range to the blockCell
     * at this row and column.
     */
    /*
    public void addGridPolygons(byte geometryType, Iterator gridPolygonIterator, StsZoneBlock zoneBlock)
    {
        while(gridPolygonIterator.hasNext())
        {
            StsPolygon polygon = (StsPolygon)gridPolygonIterator.next();
            addCellColumnGridPolygon(polygon, zoneBlock, geometryType);
        }
    }
    */
    static public boolean addCellColumnPolygons(StsModel model, StsBlock[] blocks)
    {
        for (StsBlock block : blocks)
        {
            block.addCellColumnQuads();
            block.addCellColumnsBlockGridPolygons();
        }
        StsObject[] sections = model.getObjectList(StsSection.class);
        int nSections = sections.length;
        for (int s = 0; s < nSections; s++)
        {
            StsSection section = (StsSection) sections[s];
            byte geometryType = section.getGeometryType();
            if (geometryType == StsSection.GEOM_UNALIGNED)
                addCellColumnZoneSidePolygons(section);
        }
        for (StsBlock block : blocks)
            block.addEmptyCellColumns();
        return true;
    }

    static private void addCellColumnZoneSidePolygons(StsSection section)
    {
        addCellColumnZoneSidePolygons(section.getLeftZoneSides(), LEFT);
        addCellColumnZoneSidePolygons(section.getRightZoneSides(), RIGHT);
    }

    static private void addCellColumnZoneSidePolygons(StsList zoneSides, int side)
    {
        for (int n = 0; n < zoneSides.getSize(); n++)
        {
            StsZoneSide zoneSide = (StsZoneSide) zoneSides.getElement(n);
            addCellColumnZoneSidePolygons(zoneSide, side);
        }
    }

    /**
     * These polygons are on blockZoneSides and define the intersections with grid row and column vertical planes.
     * These polygons are subdivided into layerPolygons.  Determine the layer range and add this range to the blockCell
     * at this row and column.
     *
     * @param zoneSide zoneBlockSide containing intersections
     */
    static private void addCellColumnZoneSidePolygons(StsZoneSide zoneSide, int side)
    {
        Iterator gridPolygonIterator = getGridPolygonIterator(zoneSide);
        if (gridPolygonIterator == null) return;
        StsZoneBlock zoneBlock = zoneSide.getZoneBlock();
        StsBlock block = zoneBlock.getBlock();
        int blockNumber = block.getIndex();
        while (gridPolygonIterator.hasNext())
        {
            StsPolygon polygon = (StsPolygon) gridPolygonIterator.next();
            block.addCellColumnZoneSidePolygon(polygon, blockNumber);
        }
    }

    static private Iterator getGridPolygonIterator(StsZoneSide zoneSide)
    {
        if (zoneSide == null) return null;
        StsEdgeLoop edgeLoop = zoneSide.getEdgeLoop();
        if (edgeLoop == null)
        {
            StsException.systemError(StsBlock.class, "getGridPolygonIterator", "edgeLoop is null for this zoneSide " + zoneSide.toString());
            return null;
        }
        return edgeLoop.getZoneSideGridPolygonIterator();
    }

    public void addCellColumnZoneSidePolygon(StsPolygon polygon, int blockNumber)
    {
        int row = polygon.row;
        int col = polygon.col;
        int blockRow = row - rowMin;
        int blockCol = col - colMin;
        polygon.nBlock = blockNumber;
        if (!isInsideRowCol(row, col))
        {
            StsException.systemError(this, "addCellColumnZoneSidePolygon", " row " + " col " + col + " outside boundingBox range " + toDetailString());
            return;
        }
        BlockCellColumn blockCellColumn = getBlockZoneSideCellColumn(blockRow, blockCol);
        if (blockCellColumn == null) return;
        blockCellColumn.addGridPolygon(polygon);
        if (debugCell && debugBlockIJB(blockRow, blockCol))
        {
            blockCellColumn.blockColumnDebug("adding zoneSide polygon: " + polygon.toString());
        }
    }

    public void addCellColumnBlockGridPolygon(StsPolygon polygon)
    {
        int row = polygon.row;
        int col = polygon.col;
        int blockRow = row - rowMin;
        int blockCol = col - colMin;
        BlockCellColumn blockCellColumn = getBlockGridEdgeCellColumn(blockRow, blockCol);
        if (blockCellColumn == null) return;
        blockCellColumn.addGridPolygon(polygon);
        if (debugCell && debugBlockIJB(blockRow, blockCol))
            blockCellColumn.blockColumnDebug("adding block grid polygon: " + polygon.toString());
    }

    static public boolean constructBlockCellColumns(StsModel model, StsBlock[] blocks)
    {
        for (StsBlock block : blocks)
            block.constructBlockCellColumns();
        return true;
    }

    static public boolean constructLayerColumnGrids(StsModel model, StsBlock[] blocks)
    {
        for (StsBlock block : blocks)
            block.constructLayerColumnGrids();
        return true;
    }

    public boolean constructLayerColumnGrids()
    {
        for (int row = 0; row < nCellRows; row++)
            for (int col = 0; col < nCellCols; col++)
                blockCellColumns[row][col].constructLayerGrids(sliceMax);
        return true;
    }

    static public boolean constructEdgeGridCells(StsModel model, StsBlock[] blocks)
    {
        for (StsBlock block : blocks)
            block.constructEdgeGridCells();
        return true;
    }

    public boolean constructEdgeGridCells()
    {
        for (int row = 0; row < nCellRows; row++)
            for (int col = 0; col < nCellCols; col++)
                blockCellColumns[row][col].constructGridCells();
        return true;
    }

    static public void computeCellProperties(StsZone[] zones, StsBlock[] blocks, StsEclipseModel eclipseModel)
    {
        int nLayer;

        int nZones = zones.length;
        int nBlocks = blocks.length;
        for (int z = 0; z < nZones; z++)
        {
            StsZone zone = zones[z];
            StsZoneBlock[] zoneBlocks = zone.getZoneBlocks();
            int zoneTopLayer = zone.getTopLayerNumber();
            int zoneBotLayer = zone.getBottomLayerNumber();
            for (nLayer = zoneTopLayer; nLayer <= zoneBotLayer; nLayer++)
            {
                for (int n = 0; n < nBlocks; n++)
                {
                    StsBlock block = blocks[n];
                    StsPropertyVolume[] blockPermeabilityVolumes = block.getPermeabilityVolumes(eclipseModel);
                    if (blockPermeabilityVolumes == null)
                        blockPermeabilityVolumes = block.getZoneBlockPermeabilityVolumes(zoneBlocks[n]);
                    block.computeCellPermeabilities(blockPermeabilityVolumes);
                    StsPropertyVolume blockPorosityVolume = eclipseModel.getBlockPropertyVolume(block, eclipseModel.porosity);
                    if (blockPorosityVolume == null)
                        blockPorosityVolume = zoneBlocks[n].getPropertyVolume(eclipseModel.porosity);
                    block.computeCellPorosity(blockPorosityVolume);
                }
            }
        }
    }

    private StsPropertyVolume[] getPermeabilityVolumes(StsEclipseModel eclipseModel)
    {
        StsPropertyVolume permXVolume = eclipseModel.getBlockPropertyVolume(this, eclipseModel.permX);
        if (permXVolume == null) return null;
        StsPropertyVolume permYVolume = eclipseModel.getBlockPropertyVolume(this, eclipseModel.permY);
        if (permYVolume == null) return null;
        StsPropertyVolume permZVolume = eclipseModel.getBlockPropertyVolume(this, eclipseModel.permZ);
        if (permZVolume == null) return null;
        return new StsPropertyVolume[]{permXVolume, permYVolume, permZVolume};
    }

    private StsPropertyVolume[] getZoneBlockPermeabilityVolumes(StsZoneBlock zoneBlock)
    {
        StsPropertyVolume permXVolume = zoneBlock.getEclipsePropertyVolume("PERMX");
        StsPropertyVolume permYVolume = zoneBlock.getEclipsePropertyVolume("PERMY");
        StsPropertyVolume permZVolume = zoneBlock.getEclipsePropertyVolume("PERMZ");
        return new StsPropertyVolume[]{permXVolume, permYVolume, permZVolume};
    }

    protected void computeCellPermeabilities(StsPropertyVolume[] blockPermeabilityVolumes)
    {
        for (int row = 0; row < nCellRows; row++)
            for (int col = 0; col < nCellCols; col++)
            {
                BlockCellColumn cellColumn = blockCellColumns[row][col];
                cellColumn.computeCellPermeabilities(blockPermeabilityVolumes);
            }
    }

    protected void computeCellPorosity(StsPropertyVolume blockPorosityVolume)
    {
        for (int row = 0; row < nCellRows; row++)
            for (int col = 0; col < nCellCols; col++)
            {
                BlockCellColumn cellColumn = blockCellColumns[row][col];
                cellColumn.computeCellPorosity(blockPorosityVolume);
            }
    }

    /*
        public void addBlockGridCells(StsCellTypeGrid cellTypeGrid)
        {
            int rowMin = cellTypeGrid.getRowMin();
            int rowMax = cellTypeGrid.getRowMax();
            int colMin = cellTypeGrid.getColMin();
            int colMax = cellTypeGrid.getColMax();
            byte[][] cellTypes = cellTypeGrid.getCellTypes();
            for(int row = rowMin, r = 0; row <= rowMax; row++, r++)
                for(int col = colMin, c = 0; col <= colMax; col++, c++)
                {
                    byte cellType = cellTypes[r][c];
                    if(cellType == EMPTY) continue;
                }
        }
    */
    private BlockCellColumn getBlockZoneSideCellColumn(int blockRow, int blockCol)
    {
        try
        {
            if (debugCell && debugBlockIJB(blockRow, blockCol))
                StsException.systemDebug(this, "getBlockZoneSideCellColumn");
            BlockCellColumn blockCellColumn = getBlockCellColumn(blockRow, blockCol);
            if (blockCellColumn != null)
            {
                blockCellColumn.type = EDGE;
                return blockCellColumn;
            }
            blockCellColumn = constructEdgeBlockCellColumn(blockRow, blockCol);
            blockCellColumns[blockRow][blockCol] = blockCellColumn;
            return blockCellColumn;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "getBlockZoneSideCellColumn", e);
            return null;
        }
    }

    private BlockCellColumn getBlockGridEdgeCellColumn(int blockRow, int blockCol)
    {
        try
        {
            BlockCellColumn blockCellColumn = getBlockCellColumn(blockRow, blockCol);
            if (blockCellColumn != null)
            {
                blockCellColumn.type = EDGE;
                return blockCellColumn;
            }

            blockCellColumn = constructEdgeBlockCellColumn(blockRow, blockCol);
            blockCellColumns[blockRow][blockCol] = blockCellColumn;
            return blockCellColumn;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "getBlockGridEdgeCellColumn", e);
            return null;
        }
    }

    /**
     * slice index is index of vertical grid points, so there is one less layer than slices in block.
     * So if we have slices 0 thru 10, there are 11 slices and 10 layers.
     */
    public void setLayerRange(int sliceMax)
    {
        this.sliceMin = 0;
        this.sliceMax = sliceMax;
        dbFieldChanged("sliceMin", sliceMin);
        dbFieldChanged("sliceMax", sliceMax);
    }

    public float[][] getLayerGrid(int nLayerGrid)
    {
        int nZoneBlocks = zoneBlocks.getSize();
        for (int n = 0; n < nZoneBlocks; n++)
        {
            StsZoneBlock zoneBlock = (StsZoneBlock) zoneBlocks.getElement(n);
            if (zoneBlock.hasLayerGrid(nLayerGrid))
                return zoneBlock.getLayerGrid(nLayerGrid);
        }
        StsException.systemError(this, "getLayerGrid", "Couldn't find zoneBlock for layer " + nLayerGrid);
        return null;
    }

    public void fillInsideCells()
    {
        for (int row = rowMin + 1, r = 1; row < rowMax; row++, r++)
        {
            int firstInsideCol = -1;
            for (int col = colMin + 1, c = 1; col < colMax; col++, c++)
            {
                byte cellType = getCellType(r, c);
                if (cellType == EMPTY)
                {
                    if (getCellType(r, c - 1) != EMPTY && firstInsideCol == -1)
                        firstInsideCol = c;
                    if (getCellType(r, c + 1) != EMPTY && firstInsideCol != -1)
                    {
                        int lastInsideCol = c;
                        for (int ic = firstInsideCol; ic <= lastInsideCol; ic++)
                            setCellType(r, ic, INSIDE_ROW);
                        firstInsideCol = -1;
                    }
                }
            }
        }
        for (int col = colMin + 1, c = 1; col < colMax; col++, c++)
        {
            int firstInsideRow = -1;
            for (int row = rowMin + 1, r = 1; row < rowMax; row++, r++)
            {
                byte cellType = getCellType(r, c);
                if (cellType == INSIDE_ROW)
                {
                    byte prevRowCellType = getCellType(r - 1, c);
                    if (prevRowCellType != INSIDE_ROW && prevRowCellType != EMPTY && firstInsideRow == -1)
                        firstInsideRow = r;
                    byte nextRowCellType = getCellType(r + 1, c);
                    if (nextRowCellType != INSIDE_ROW && nextRowCellType != EMPTY && firstInsideRow != -1)
                    {
                        int lastInsideRow = r;
                        for (int ir = firstInsideRow; ir <= lastInsideRow; ir++)
                            setCellType(ir, c, FULL);
                        firstInsideRow = -1;
                    }
                }
            }
        }
    }

    public byte getCellType(int row, int col)
    {
        try
        {
            return blockCellColumns[row][col].type;
        }
        catch (Exception e)
        {
            StsException.systemError(this, "getCellType", "row " + row + " col " + col);
            return EMPTY;
        }
    }

    public void setCellType(int row, int col, byte type)
    {
        try
        {
            if (blockCellColumns[row][col] == null)
                blockCellColumns[row][col] = new BlockCellColumn(type, row, col);
        }
        catch (Exception e)
        {
            StsException.systemError(this, "setCellType", "row " + row + " col " + col);
        }
    }

    public void addCellColumnQuads()
    {
        int nBlockGrids = blockGrids.getSize();
        for (int n = 0; n < nBlockGrids; n++)
        {
            StsBlockGrid blockGrid = (StsBlockGrid) blockGrids.getElement(n);
            StsEdgeLoop edgeLoop = blockGrid.getEdgeLoop();
            StsList quadStripList = edgeLoop.getQuadStrips();
            int nQuadStrips = quadStripList.getSize();
            for (int i = 0; i < nQuadStrips; i++)
            {
                QuadStrip quadStrip = (QuadStrip) quadStripList.getElement(i);
                int row = quadStrip.rowNumber;
                int firstCol = quadStrip.firstCol;
                int lastCol = quadStrip.lastCol;
                for (int col = firstCol; col < lastCol; col++)
                    addCellColumnQuad(blockGrid, row, col);
            }
        }
    }

    public void addEmptyCellColumns()
    {
        for (int blockRow = 0; blockRow < nCellRows; blockRow++)
            for (int blockCol = 0; blockCol < nCellCols; blockCol++)
                if (blockCellColumns[blockRow][blockCol] == null)
                    blockCellColumns[blockRow][blockCol] = constructEmptyBlockCellColumn(blockRow, blockCol);
    }

    public void addCellColumnsBlockGridPolygons()
    {
        int nBlockGrids = blockGrids.getSize();
        for (int n = 0; n < nBlockGrids; n++)
        {
            StsBlockGrid blockGrid = (StsBlockGrid) blockGrids.getElement(n);
            StsEdgeLoop edgeLoop = blockGrid.getEdgeLoop();
            Iterator<StsPolygon> polygonIterator = edgeLoop.getPolygonIterator();
            while (polygonIterator.hasNext())
            {
                StsPolygon polygon = polygonIterator.next();
                addCellColumnBlockGridPolygon(polygon);
            }
        }
    }

    private void addCellColumnQuad(StsBlockGrid blockGrid, int row, int col)
    {
        int blockRow = row - rowMin;
        int blockCol = col - colMin;
        BlockCellColumn blockColumn = getBlockCellColumn(blockRow, blockCol);
        if (blockColumn == null)
            blockColumn = blockCellColumns[blockRow][blockCol] = constructFullBlockCellColumn(blockRow, blockCol);
        if (debugCell && debugBlockIJB(blockRow, blockCol))
            blockColumn.blockColumnDebug("adding quad polygon");
        StsPolygon polygon = blockGrid.getQuadPolygon(row, col);
        blockColumn.addGridPolygon(polygon);
    }

    private boolean debugBlockIJB(int blockRow, int blockCol)
    {
        return blockRow + rowMin == debugCellRow && blockCol + colMin == debugCellCol && getIndex() == debugCellBlockIndex;
    }

    private boolean debugIJB(int row, int col)
    {
        return row == debugCellRow && col == debugCellCol && getIndex() == debugCellBlockIndex;
    }

    static public boolean debugIJK(int row, int col, int nLayer)
    {
        boolean debugOn = row == debugCellRow && col == debugCellCol && nLayer == debugCellLayer;
        if (debugOn)
            System.out.println("debug is on.");
        return debugOn;
    }


    public boolean debugBlockIJK(int blockRow, int blockCol, int nLayer)
    {
        boolean debugOn = blockRow + rowMin == debugCellRow && blockCol + colMin == debugCellCol && nLayer == debugCellLayer;
        if (debugOn)
            System.out.println("debug is on.");
        return debugOn;
    }

    static public boolean debugIJK(BlockCellColumn.GridCell gridCell)
    {
        return debugIJK(gridCell.getRow(), gridCell.getCol(), gridCell.nLayer);
    }

    private void constructBlockCellColumns()
    {
        cellXYArea = xInc * yInc;

        StsZoneBlock[] zoneBlocks = (StsZoneBlock[]) this.zoneBlocks.getCastList(StsZoneBlock.class);
        for (int blockRow = 0; blockRow < nCellRows; blockRow++)
            for (int blockCol = 0; blockCol < nCellCols; blockCol++)
                blockCellColumns[blockRow][blockCol].constructPrisms(zoneBlocks);
    }

    private BlockCellColumn constructEdgeBlockCellColumn(int blockRow, int blockCol)
    {
        return new BlockCellColumn(EDGE, blockRow, blockCol);
    }

    private BlockCellColumn constructEmptyBlockCellColumn(int blockRow, int blockCol)
    {
        return new BlockCellColumn(EMPTY, blockRow, blockCol);
    }

    private BlockCellColumn constructFullBlockCellColumn(int blockRow, int blockCol)
    {
        return new BlockCellColumn(FULL, blockRow, blockCol);
    }

    public boolean isCellActive(int blockRow, int blockCol, int nLayer)
    {
        BlockCellColumn cellColumn = getBlockCellColumn(blockRow, blockCol);
        if (cellColumn == null) return false;
        return cellColumn.isCellActive(nLayer);
    }

    public float computeCellDepth(int blockRow, int blockCol, int nLayer, float[][] topLayerGrid, float[][] botLayerGrid, float fillZ)
    {
        BlockCellColumn cellColumn = getBlockCellColumn(blockRow, blockCol);
        if (cellColumn == null || cellColumn.type == EMPTY) return fillZ;
        if (cellColumn.type == FULL)
            return computeFullCellDepth(blockRow, blockCol, topLayerGrid, botLayerGrid, fillZ);
        else // EDGE cellColumn
            return cellColumn.computeEdgeCellDepth(nLayer, fillZ);

    }

    private float computeFullCellDepth(int blockRow, int blockCol, float[][] topLayerGrid, float[][] botLayerGrid, float fillZ)
    {
        float topZ = getLayerAvgZ(blockRow, blockCol, topLayerGrid);
        if (topZ == nullValue) return fillZ;
        float botZ = getLayerAvgZ(blockRow, blockCol, botLayerGrid);
        if (botZ == nullValue) return fillZ;
        return (topZ + botZ) / 2;
    }

    private float getLayerAvgZ(int i, int j, float[][] layerGrid)
    {
        float zSum = 0.0f;
        float nValues = 0;
        float z;
        z = layerGrid[i][j];
        if (z != nullValue)
        {
            zSum += z;
            nValues++;
        }
        z = layerGrid[i + 1][j];
        if (z != nullValue)
        {
            zSum += z;
            nValues++;
        }
        z = layerGrid[i][j + 1];
        if (z != nullValue)
        {
            zSum += z;
            nValues++;
        }
        z = layerGrid[i + 1][j + 1];
        if (z != nullValue)
        {
            zSum += z;
            nValues++;
        }
        if (nValues == 0) return nullValue;
        return zSum / nValues;
    }
/*
    public float[] getCellCenter(int row, int col, int nLayer, float[][] topLayerGrid, float[][] botLayerGrid)
    {
        int blockRow = row - rowMin;
        int blockCol = col = colMin;
        BlockCellColumn cellColumn = getBlockCellColumn(blockRow, blockCol);
        if (cellColumn == null || cellColumn.type == EMPTY) return null;
        if (cellColumn.type == FULL)
            return getFullCellCenter(blockRow, blockCol, topLayerGrid, botLayerGrid);
        else // EDGE cellColumn
            return cellColumn.getEdgeCellCenter(nLayer);
    }
*/
/*
    public float getCellAreaXPlus(int blockRow, int blockCol, int nLayer, float[][] topLayerGrid, float[][] botLayerGrid)
    {
        BlockCellColumn cellColumn = getBlockCellColumn(blockRow, blockCol);
        if (cellColumn == null || cellColumn.type == EMPTY) return 0.0f;
        BlockCellColumn.GridCell gridCell = getBlockGridCell(blockRow, blockCol, nLayer);
        return gridCell.getFaceArea(CELL_PLUS_DX);
    }

    public float getCellAreaXMinus(int blockRow, int blockCol, int nLayer, float[][] topLayerGrid, float[][] botLayerGrid)
    {
        return getCellAreaXPlus(blockRow, blockCol - 1, nLayer, topLayerGrid, botLayerGrid);
    }

    public float getCellAreaYPlus(int blockRow, int blockCol, int nLayer, float[][] topLayerGrid, float[][] botLayerGrid)
    {
        BlockCellColumn cellColumn = getBlockCellColumn(blockRow, blockCol);
        if (cellColumn == null || cellColumn.type == EMPTY) return 0.0f;
        BlockCellColumn.GridCell gridCell = getBlockGridCell(blockRow, blockCol, nLayer);
        return gridCell.getFaceArea(CELL_PLUS_DY);
    }

    public float getCellAreaYMinus(int blockRow, int blockCol, int nLayer, float[][] topLayerGrid, float[][] botLayerGrid)
    {
        return getCellAreaYPlus(blockRow - 1, blockCol, nLayer, topLayerGrid, botLayerGrid);
    }
*/

    /*
        public float getCellAreaZPlus(int blockRow, int blockCol, int nLayer)
        {
            BlockCellColumn cellColumn = getBlockCellColumn(blockRow, blockCol);
            if (cellColumn == null || cellColumn.type == EMPTY) return 0.0f;
            return cellColumn.getCellAreaZPlus(nLayer);
        }

        public float computeCellAreaZMinus(int blockRow, int blockCol, int nLayer)
        {
            BlockCellColumn cellColumn = getBlockCellColumn(blockRow, blockCol);
            if (cellColumn == null || cellColumn.type == EMPTY) return 0.0f;
            return cellColumn.getCellAreaZPlus(nLayer - 1);
        }
    */
    static public byte getOtherDirection(byte direction)
    {
        switch (direction)
        {
            case CELL_PLUS_DX:
                return CELL_MINUS_DX;
            case CELL_MINUS_DX:
                return CELL_PLUS_DX;
            case CELL_PLUS_DY:
                return CELL_MINUS_DY;
            case CELL_MINUS_DY:
                return CELL_PLUS_DY;
            case CELL_PLUS_DZ:
                return CELL_MINUS_DZ;
            case CELL_MINUS_DZ:
                return CELL_PLUS_DZ;
            default:
                return -1;
        }
    }

    public float computeTranXPlus(int blockRow, int blockCol, int nLayer)
    {
        BlockCellColumn.GridCell gridCellMinus = getBlockGridCell(blockRow, blockCol, nLayer);
        BlockCellColumn.GridCell gridCellPlus = getBlockGridCell(blockRow, blockCol + 1, nLayer);
        return effectiveTrans(gridCellMinus, gridCellPlus, CELL_MINUS_DX, CELL_PLUS_DX);
    }

    public float computeTranYPlus(int blockRow, int blockCol, int nLayer)
    {
        BlockCellColumn.GridCell gridCellMinus = getBlockGridCell(blockRow, blockCol, nLayer);
        BlockCellColumn.GridCell gridCellPlus = getBlockGridCell(blockRow + 1, blockCol, nLayer);
        return effectiveTrans(gridCellMinus, gridCellPlus, CELL_MINUS_DY, CELL_PLUS_DY);
    }

    public float computeTranZPlus(int blockRow, int blockCol, int nLayer)
    {
        BlockCellColumn.GridCell gridCellMinus = getBlockGridCell(blockRow, blockCol, nLayer);
        BlockCellColumn.GridCell gridCellPlus = getBlockGridCell(blockRow, blockCol, nLayer + 1);
        return effectiveTrans(gridCellMinus, gridCellPlus, CELL_MINUS_DZ, CELL_PLUS_DZ);
    }

    private float effectiveTrans(BlockCellColumn.GridCell gridCellMinus, BlockCellColumn.GridCell gridCellPlus, byte minusDirection, byte plusDirection)
    {
        if (gridCellMinus == null || gridCellPlus == null) return 0.0f;
        float area = (gridCellMinus.faceAreas[plusDirection] + gridCellPlus.faceAreas[minusDirection]) / 2;
        float dxMinus = gridCellMinus.halfSizes[plusDirection];
        float dxPlus = gridCellPlus.halfSizes[minusDirection];
        float kMinus = gridCellMinus.getPermeability(plusDirection);
        float kPlus = gridCellPlus.getPermeability(minusDirection);
        double kAvgDx = kMinus * dxPlus + kPlus * dxMinus;
        if (kAvgDx == 0.0) return 0.0f;
        float tran = (float) (area * kPlus * kMinus / kAvgDx);
        if (tran < 0.0f)
        {
            StsException.systemError(this, "effectiveTrans", "Tran is negative: " + tran + " minusDx: " + dxMinus + " plusDx: " + dxPlus + " kMinus: " + kMinus + " kPlus: " + kPlus + " area: " + area);
            tran = 0.0f;
        }
        if(gridCellMinus.debugBlockIJK())
            StsException.systemDebug(this, "effectiveTrans", "Tran " + tran + " grid cell " + gridCellMinus.toString() + " to " + gridCellPlus.toString());
        return tran;
    }

    public void setGridCellPoreVolume(int blockRow, int blockCol, int nLayer, float poreVolume)
    {
        BlockCellColumn.GridCell gridCell = getBlockGridCell(blockRow, blockCol, nLayer);
        if (gridCell == null) return;
        gridCell.setGridCellPoreVolume(poreVolume);
    }

    public void setGridCellTran(int blockRow, int blockCol, int nLayer, byte direction, float trans, float area)
    {
        BlockCellColumn.GridCell gridCell = getBlockGridCell(blockRow, blockCol, nLayer);
        if (gridCell == null) return;
        gridCell.setGridCellTran(direction, trans);
    }

    public void setPropertyVolumeBlockValue(StsPropertyType propertyType, int blockRow, int blockCol, int nLayer, float value, StsEclipseModel eclipseModel)
    {
        StsPropertyVolume propertyVolume = eclipseModel.getBlockPropertyVolume(this, propertyType);
        propertyVolume.setBlockValue(blockRow, blockCol, nLayer, value);
    }

    public void addPropertyVolumeBlockValue(StsPropertyType propertyType, int blockRow, int blockCol, int nLayer, float value, StsEclipseModel eclipseModel)
    {
        StsPropertyVolume propertyVolume = eclipseModel.getBlockPropertyVolume(this, propertyType);
        propertyVolume.addBlockValue(blockRow, blockCol, nLayer, value);
    }

    private float effectiveTrans(double minusDx, double plusDx, double kMinus, double kPlus, double area)
    {
        double kAvgDx = kMinus * plusDx + kPlus * minusDx;
        if (kAvgDx == 0.0) return 0.0f;
        float tran = (float) (area * kPlus * kMinus / kAvgDx);
        if (tran < 0.0f)
        {
            StsException.systemError(this, "effectiveTrans", "Tran is negative: " + tran + " minusDx: " + minusDx + " plusDx: " + plusDx + " kMinus: " + kMinus + " kPlus: " + kPlus + " area: " + area);
            tran = 0.0f;
        }
        return tran;
    }

    /** compute the thickknesses at a center point on the centerLayerGrid up to the topLayerGrid and down to the botLayerGrid */
    public float[] getAvgThicknesses(int blockRow, int blockCol, float[][] topLayerGrid, float[][] centerLayerGrid, float[][] botLayerGrid)
    {
        float zTop = getLayerAvgZ(blockRow, blockCol, topLayerGrid);
        if (zTop == nullValue) return null;
        float zCen = getLayerAvgZ(blockRow, blockCol, centerLayerGrid);
        if (zCen == nullValue) return null;
        float zBot = getLayerAvgZ(blockRow, blockCol, botLayerGrid);
        if (zBot == nullValue) return null;
        return new float[]{zCen - zTop, zBot - zCen};
    }

    public float getAvgGridCellThickness(int blockRow, int blockCol, float[][] topLayerGrid, float[][] botLayerGrid)
    {
        float zTop = getLayerAvgZ(blockRow, blockCol, topLayerGrid);
        if (zTop == nullValue) return 0.0f;
        float zBot = getLayerAvgZ(blockRow, blockCol, botLayerGrid);
        if (zBot == nullValue) return 0.0f;
        return zBot - zTop;
    }

    private float getFullCellAreaXPlus(int blockRow, int blockCol, float[][] topLayerGrid, float[][] botLayerGrid)
    {
        try
        {
            float dz1 = botLayerGrid[blockRow][blockCol + 1] - topLayerGrid[blockRow][blockCol + 1];
            float dz2 = botLayerGrid[blockRow + 1][blockCol + 1] - topLayerGrid[blockRow + 1][blockCol + 1];
            return yInc * (dz1 + dz2) / 2;
        }
        catch (Exception e)
        {
            StsException.systemError(this, "getFullCellAreaXPlus");
            return 0.0f;
        }
    }

    private float getFullCellAreaXMinus(int blockRow, int blockCol, float[][] topLayerGrid, float[][] botLayerGrid)
    {
        try
        {
            float dz1 = botLayerGrid[blockRow][blockCol] - topLayerGrid[blockRow][blockCol];
            float dz2 = botLayerGrid[blockRow + 1][blockCol] - topLayerGrid[blockRow + 1][blockCol];
            return yInc * (dz1 + dz2) / 2;
        }
        catch (Exception e)
        {
            StsException.systemError(this, "getFullCellAreaXPlus");
            return 0.0f;
        }
    }

    private float getFullCellAreaYPlus(int blockRow, int blockCol, float[][] topLayerGrid, float[][] botLayerGrid)
    {
        float dz1 = botLayerGrid[blockRow + 1][blockCol] - topLayerGrid[blockRow + 1][blockCol];
        float dz2 = botLayerGrid[blockRow + 1][blockCol + 1] - topLayerGrid[blockRow + 1][blockCol + 1];
        return xInc * (dz1 + dz2) / 2;
    }

    private float getFullCellAreaYMinus(int blockRow, int blockCol, float[][] topLayerGrid, float[][] botLayerGrid)
    {
        float dz1 = botLayerGrid[blockRow][blockCol] - topLayerGrid[blockRow][blockCol];
        float dz2 = botLayerGrid[blockRow][blockCol + 1] - topLayerGrid[blockRow][blockCol + 1];
        return xInc * (dz1 + dz2) / 2;
    }

    private float[] getFullCellCenter(int blockRow, int blockCol, float[][] topLayerGrid, float[][] botLayerGrid)
    {
        float topZ = getLayerAvgZ(blockRow, blockCol, topLayerGrid);
        if (topZ == nullValue) return null;
        float botZ = getLayerAvgZ(blockRow, blockCol, botLayerGrid);
        if (botZ == nullValue) return null;
        float avgZ = (topZ + botZ) / 2;
        float xCenter = getXCoor(blockCol + 0.5f);
        float yCenter = getYCoor(blockRow + 0.5f);
        return new float[]{xCenter, yCenter, avgZ};
    }

    public BlockCellColumn.GridCell getGridCell(int[] ijkb)
    {
        int blockRow = ijkb[0] - rowMin;
        int blockCol = ijkb[1] - colMin;
        int nLayer = ijkb[2];
        return getBlockGridCell(blockRow, blockCol, nLayer);
    }

    public BlockCellColumn.GridCell getGridCell(int row, int col, int nLayer)
    {
        int blockRow = row - rowMin;
        int blockCol = col - colMin;
        BlockCellColumn blockCellColumn = getBlockCellColumn(blockRow, blockCol);
        if (blockCellColumn == null) return null;
        return blockCellColumn.getGridCell(nLayer);
    }

    public BlockCellColumn.GridCell getParentCellOrGridCell(int row, int col, int nLayer)
    {
        int blockRow = row - rowMin;
        int blockCol = col - colMin;
        BlockCellColumn blockCellColumn = getBlockCellColumn(blockRow, blockCol);
        if (blockCellColumn == null) return null;
        return blockCellColumn.getParentCellOrGridCell(nLayer);
    }

    public BlockCellColumn.GridCell getBlockGridCell(int blockRow, int blockCol, int nLayer)
    {
        BlockCellColumn blockCellColumn = getBlockCellColumn(blockRow, blockCol);
        if (blockCellColumn == null) return null;
        return blockCellColumn.getGridCell(nLayer);
    }

    public BlockCellColumn.GridCell getBlockGridCell(int[] ijk)
    {
        BlockCellColumn blockCellColumn = getBlockCellColumn(ijk[0], ijk[1]);
        if (blockCellColumn == null) return null;
        return blockCellColumn.getGridCell(ijk[2]);
    }

    public void deleteBlockGridCell(int blockRow, int blockCol, int nLayer)
    {
        BlockCellColumn blockCellColumn = getBlockCellColumn(blockRow, blockCol);
        BlockCellColumn.GridCell gridCell = blockCellColumn.gridCells[nLayer];
        gridCell.poreVolume = 0.0f;
        gridCell.volume = 0.0f;
        // blockCellColumn.gridCells[nLayer] = null;
    }

    public boolean getBlockGridCellOk(int blockRow, int blockCol, int nLayer)
    {
        BlockCellColumn blockCellColumn = getBlockCellColumn(blockRow, blockCol);
        if (blockCellColumn == null) return false;
        if (blockCellColumn.type == EMPTY) return false;
        if (blockCellColumn.type == FULL) return true;
        else return blockCellColumn.getGridCell(nLayer) != null;
    }

    public double getBlockGridCellParameter(int blockRow, int blockCol, int nLayer, byte direction, float[][] topLayerGrid, float[][] botLayerGrid)
    {
        BlockCellColumn blockCellColumn = getBlockCellColumn(blockRow, blockCol);
        if (blockCellColumn == null)
        {
            StsException.systemError(this, "getBlockGridCellParameter", "Tried to get gridCell direction from null column.\n" +
                "    block: " + getIndex() + " row: " + (blockRow + rowMin) + " col " + (blockCol + colMin) + " nLayer " + nLayer + " direction " + cellParameterStrings[direction]);
            return 0.0;
        }
        byte columnType = blockCellColumn.type;
        if (columnType == EMPTY)
        {
            StsException.systemError(this, "getBlockGridCellParameter", "Tried to get gridCell direction from empty cell.\n" +
                "    block: " + getIndex() + " row: " + (blockRow + rowMin) + " col " + (blockCol + colMin) + " nLayer " + nLayer + " direction " + cellParameterStrings[direction]);
            return 0.0;
        }
        else if (columnType == FULL)
        {
            switch (direction)
            {
                case CELL_PLUS_DX:
                case CELL_MINUS_DX:
                    return xInc / 2;
                case CELL_PLUS_DY:
                case CELL_MINUS_DY:
                    return yInc / 2;
                case CELL_PLUS_DZ:
                case CELL_MINUS_DZ:
                    return getAvgGridCellThickness(blockRow, blockCol, topLayerGrid, botLayerGrid) / 2;
                default:
                    StsException.systemError(this, "getBlockGridCellParameter", "Undefined Cell direction " + direction);
                    return 0.0;
            }
        }
        else // columnType == EDGE
            return blockCellColumn.getGridCellParameter(nLayer, direction);
    }

    public BlockCellColumn getCellColumn(int row, int col)
    {
        if (!this.isInsideRowCol(row, col)) return null;
        return blockCellColumns[row - rowMin][col - colMin];
    }

    public BlockCellColumn getBlockCellColumn(int blockRow, int blockCol)
    {
        if (!this.isInsideBlockCellRowCol(blockRow, blockCol)) return null;
        try
        {
            return blockCellColumns[blockRow][blockCol];
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "getBlockCellColumn", e);
            return null;
        }
    }

    /** There are colMax cellColumns in this block.  Return false if this blockCol number exceeds this. */
    public boolean canConnectXPlus(int blockRow, int blockCol, int nLayer)
    {
        if (blockCol + 1 >= nCellCols) return false;
        return getBlockGridCellOk(blockRow, blockCol, nLayer) && getBlockGridCellOk(blockRow, blockCol + 1, nLayer);
    }

    /** There are rowMax cellRows in this block.  Return false if this blockRow number exceeds this. */
    public boolean canConnectYPlus(int blockRow, int blockCol, int nLayer)
    {
        if (blockRow + 1 >= nCellRows) return false;
        return getBlockGridCellOk(blockRow, blockCol, nLayer) && getBlockGridCellOk(blockRow + 1, blockCol, nLayer);
    }

    /** There are nTotalLayers layers in this block.  Return false if this layer+1 exceeds this. */
    public boolean canConnectZPlus(int blockRow, int blockCol, int nLayer)
    {
        if (nLayer + 1 >= sliceMax) return false;
        return getBlockGridCellOk(blockRow, blockCol, nLayer) && getBlockGridCellOk(blockRow, blockCol, nLayer + 1);
    }

    /*
        public void createPropertyVolume(StsPropertyType propertyType, StsSpectrum propertySpectrum)
        {
            addPropertyVolume(new StsPropertyVolume(propertyType, this, propertySpectrum, false));
        }

        public void createPropertyVolume(StsPropertyType propertyType)
        {
            addPropertyVolume(new StsPropertyVolume(propertyType, this, null, false));
        }

        public void addPropertyVolume(StsPropertyVolume propertyVolume)
        {
            if(propertyVolumes == null)
                propertyVolumes = new ArrayList<StsPropertyVolume>();
            propertyVolumes.add(propertyVolume);
        }

        public StsPropertyVolume getPropertyVolume(StsPropertyType propertyType)
        {
            if(propertyType == null) return null;
            if(propertyVolumes == null) return null;
            String propertyName = propertyType.name;
            int nPropertyVolumes = propertyVolumes.size();
            for(int n = 0; n < nPropertyVolumes; n++)
            {
                StsPropertyVolume propertyVolume = propertyVolumes.get(n);
                if(propertyVolume.getName().equals(propertyName))
                    return propertyVolume;
            }
            return null;
        }
    */
    public float getPropertyVolumeValue(StsPropertyVolume propertyVolume, int blockRow, int blockCol, int nLayer)
    {
        return propertyVolume.getBlockValue(blockRow, blockCol, nLayer);
    }

    public int[] getEclipseRowCol(int row, int col)
    {
        if (!isInsideCellRowCol(row, col)) return null;
        int blockRow = row - rowMin;
        int blockCol = col - colMin;
        BlockCellColumn blockCellColumn = blockCellColumns[blockRow][blockCol];
        byte type = blockCellColumn.type;
        if (type == EMPTY) return null;
        int eclipseRow = eclipseRowMin + rowMax - row;
        int eclipseCol = blockCol + 1;
        return new int[]{eclipseRow, eclipseCol};
    }

    public int[] getEclipseRowColF(float rowF, float colF)
    {
        int row = Math.round(rowF);
        int col = Math.round(colF);
        return getEclipseRowCol(row, col);
    }

    public int[] getEclipseRowColFromBlockRowCol(int blockRow, int blockCol)
    {
        return getEclipseRowCol(blockRow + rowMin, blockCol + colMin);
    }

    public BlockCellColumn getGridObject(int row, int col)
    {
        return getCellColumn(row, col);
    }

    public void createTruncatedCellListComputeTruncatedFraction(float mergeCellMinFraction)
    {
        actualMergedCellMinFraction = 1.0f;
        truncatedCells = new TreeSet<BlockCellColumn.GridCell>(getTruncatedFractionComparator());

        // set time constants for all cells
        Iterator<BlockCellColumn.GridCell> gridCellIterator = getGridCellIterator();

        while (gridCellIterator.hasNext())
        {
            StsBlock.BlockCellColumn.GridCell gridCell = gridCellIterator.next();
            if (gridCell == null) continue;
            if (gridCell.debugBlockIJK() || gridCell.debugBlockIJKNeighbor(debugDirection))
                StsException.systemDebug(this, "createTruncatedCellListComputeTruncatedFraction");
            // if(!gridCell.isEdgeCell() || gridCell.volume <= 0.0f) continue;
            if (gridCell.volume <= 0.0f) continue;
            float truncatedFraction = gridCell.truncatedFraction;
            if (truncatedFraction > mergeCellMinFraction) continue;
            actualMergedCellMinFraction = Math.min(actualMergedCellMinFraction, truncatedFraction);
            truncatedCells.add(gridCell);
        }
    }

    public int getNTruncatedCells()
    {
        if (truncatedCells == null) return 0;
        return truncatedCells.size();
    }

    public void checkClearParentCells()
    {
        parentCells = null;
    }

    public void recomputeParentCellTruncatedCellFractions(float mergeCellMinFraction)
    {
        actualMergedCellMinFraction = 1.0f;
        numberUnmergedSmallCells = 0;
        numberMergedCells = 0;
        for (BlockCellColumn.ParentCell parentCell : parentCells)
        {
            parentCell.computeTruncatedFraction();
            if (parentCell.truncatedFraction < mergeCellMinFraction)
            {
                numberUnmergedSmallCells += parentCell.children.size();
                StsException.systemDebug(this, "recomputeParentCellTruncatedCellFractions", "ParentCell truncated fraction too small " + parentCell.toString());
            }
            else
                numberMergedCells += parentCell.children.size();
            actualMergedCellMinFraction = Math.min(actualMergedCellMinFraction, parentCell.truncatedFraction);
        }

        for (BlockCellColumn.GridCell truncatedCell : truncatedCells)
        {
            if (truncatedCell.parentCell != null) continue;
            if (truncatedCell.truncatedFraction < mergeCellMinFraction)
            {
                numberUnmergedSmallCells++;
                StsException.systemDebug(this, "recomputeParentCellTruncatedCellFractions", "ParentCell truncated fraction too small " + truncatedCell.toString());
            }
            actualMergedCellMinFraction = Math.min(actualMergedCellMinFraction, truncatedCell.truncatedFraction);
        }
    }

    public void recomputeAdjustedTruncatedCellFractions(float mergeCellMinFraction)
    {
        actualMergedCellMinFraction = 1.0f;
        numberUnmergedSmallCells = 0;
        numberMergedCells = 0;
        for (BlockCellColumn.GridCell truncatedCell : truncatedCells)
        {
            if (truncatedCell.parentCell != null) continue;
            if (truncatedCell.truncatedFraction < mergeCellMinFraction)
            {
                numberUnmergedSmallCells++;
                StsException.systemDebug(this, "recomputeAdjustedTruncatedCellFractions", "tCell truncated fraction too small " + truncatedCell.toString());
            }
            actualMergedCellMinFraction = Math.min(actualMergedCellMinFraction, truncatedCell.truncatedFraction);
        }
    }

    public void checkMergeCells(StsEclipseOutput eclipseOutput)
    {
        parentCells = new ArrayList<BlockCellColumn.ParentCell>();
        // clear parentCells set on previous merge operation
        for (BlockCellColumn.GridCell gridCell : truncatedCells)
            gridCell.parentCell = null;

        // create list of parentCells whose truncated children are merged into it
        for (BlockCellColumn.GridCell gridCell : truncatedCells)
        {
            if (gridCell.truncatedFraction >= eclipseOutput.mergeCellMinFraction) break;
            gridCell.checkSetParentCell();
        }
        /*
            for(BlockCellColumn.ParentCell parentCell : parentCells)
            {
                parentCell.computeCellGeometry();
            }
            StsEclipseOutput.EclipseWriter outputWriter = eclipseOutput.outputWriter;
            for(BlockCellColumn.ParentCell parentCell : parentCells)
                parentCell.connectParentCell(outputWriter);
        */
    }

    public void merge(TreeSet<StsBlock.BlockCellColumn.GridCell> debugCellList, StsEclipseModel eclipseModel, StsEclipseOutput eclipseOutput)
    {
        if (parentCells == null) return;

        for (BlockCellColumn.ParentCell parentCell : parentCells)
        {
            parentCell.mergeOrAdjust(eclipseModel, eclipseOutput);
            debugCellList.add(parentCell);
        }
        for (BlockCellColumn.GridCell truncatedCell : truncatedCells)
            if (truncatedCell.parentCell == null) debugCellList.add(truncatedCell);
    }

    public boolean delete()
    {
        if (blockSides != null)
        {
            blockSides.deleteAll();
            blockSides = null;
        }
        if (blockGrids != null)
            blockGrids = null;
        if (connectedBlocks != null)
        {
            connectedBlocks.delete();
            connectedBlocks = null;
        }
        super.delete();
        return true;
    }

    public void deleteTransientArrays()
    {
        blockGrids = null;
        zoneBlocks = null;
        deleteLineZones();
    }

    private void deleteLineZones()
    {
        StsObject[] lineObjects;
        lineObjects = currentModel.getObjectList(StsLine.class);
        deleteLineZones(lineObjects);
        lineObjects = currentModel.getObjectList(StsFaultLine.class);
        deleteLineZones(lineObjects);      
    }

    private void deleteLineZones(StsObject[] lineObjects)
    {
        int nObjects = lineObjects.length;
        for(int n = 0; n < nObjects; n++)
        {
            StsLine line = (StsLine)lineObjects[n];
            line.deleteTransients();
        }
    }

    public void setCurrentPropertyType(StsPropertyType propertyType)
    {
        if (currentPropertyType == propertyType) return;
        currentPropertyType = propertyType;
        propertyTypeChanged();
    }

    public void propertyTypeChanged()
    {
        int nZoneBlocks = zoneBlocks.getSize();
        for (int n = 0; n < nZoneBlocks; n++)
        {
            StsZoneBlock zoneBlock = (StsZoneBlock) zoneBlocks.getElement(n);
            zoneBlock.propertyTypeChanged();
        }
    }

    public StsPropertyVolume getCurrentPropertyVolume(){ return currentPropertyVolume; }

    public boolean anyDependencies()
    {
        return false;
    }

    public StsFieldBean[] getDisplayFields()
    {
        if (displayFields == null)
        {
            displayFields = new StsFieldBean[]
                {
                    new StsBooleanFieldBean(StsBlock.class, "isVisible", "Enable")
                };
        }
        return displayFields;
    }

    public StsFieldBean[] getPropertyFields()
    {
        return propertyFields;
    }

    public Object[] getChildren()
    {
        return new Object[0];
    }

    public StsObjectPanel getObjectPanel()
    {
        if (objectPanel == null) objectPanel = StsObjectPanel.constructor(this, true);
        return objectPanel;
    }

    public void treeObjectSelected()
    {
        currentModel.getCreateStsClass("com.Sts.DBTypes.StsBlock").selected(this);
    }

    public ColumnEdgeGridCellIterator getEdgeGridCellIterator()
    {
        return new ColumnEdgeGridCellIterator(this);
    }

    public ColumnGridCellIterator getGridCellIterator()
    {
        return new ColumnGridCellIterator(this);
    }

    public BlockCellColumn[][] getBlockCellColumns()
    {
        return blockCellColumns;
    }

    public void setBlockCellColumns(BlockCellColumn[][] blockCellColumns)
    {
        this.blockCellColumns = blockCellColumns;
    }

    public void setBlockCellColumn(int row, int col, BlockCellColumn blockCellColumn)
    {
        blockCellColumns[row][col] = blockCellColumn;
    }

    public void setCurrentPropertyVolume(StsPropertyVolume currentPropertyVolume)
    {
        if (this.currentPropertyVolume == currentPropertyVolume) return;
        this.currentPropertyVolume = currentPropertyVolume;
        propertyTypeChanged();
    }
/*
    public void createIndexMap(StsPropertyVolume indexVolume, int nEclipseRows)
    {
        int blockRow, blockCol, slice;
        int index;
        float[] indexMap = indexVolume.getValues();
        Arrays.fill(indexMap, -1.0f);
        int nCellRows = getNCellRows();
        int nCellCols = getNCellCols();
        int nCellSlices = getNCellSlices();
        try
        {
            for(blockRow = 0; blockRow < nCellRows; blockRow++)
            {
                for(blockCol = 0; blockCol < nCellCols; blockCol++)
                {
                    BlockCellColumn column = this.getBlockCellColumn(blockRow, blockCol);
                    if(column != null && column.type == EMPTY) continue;
                    for(slice = 0; slice < nCellSlices; slice++)
                    {
                        BlockCellColumn.GridCell gridCell = column.getGridCell(slice);
                        if(gridCell == null || gridCell.poreVolume <= 0.0f) continue;
                        int eclipseIndex = gridCell.getEclipseIndex(blockRow, blockCol, slice, nEclipseRows);
                        index = indexVolume.getBlockIndex3d(blockRow, blockCol, slice);
                        indexMap[index] = eclipseIndex;
                    }
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "createIndexMap", e);
        }
    }
 */

    class ColumnGridCellIterator implements Iterator<BlockCellColumn.GridCell>
    {
        StsGridBoundingBox.Cell2dIterator<BlockCellColumn> blockColumnIterator;
        Iterator<BlockCellColumn.GridCell> gridCellIterator = null;
        BlockCellColumn.GridCell nextCell;

        public ColumnGridCellIterator(StsBlock block)
        {
            blockColumnIterator = new StsGridBoundingBox.Cell2dIterator<BlockCellColumn>(block);
        }

        public boolean hasNext()
        {
            if (gridCellIterator != null && gridCellIterator.hasNext())
                return true;
            else
            {
                while (blockColumnIterator.hasNext())
                {
                    BlockCellColumn blockColumn = blockColumnIterator.next();
                    if (blockColumn.type == EMPTY) continue;
                    gridCellIterator = blockColumn.getGridCellIterator();
                    if (gridCellIterator.hasNext())
                        return true;
                }
            }
            return false;
        }

        public BlockCellColumn.GridCell next(){ return gridCellIterator.next(); }

        public void remove(){ }
    }

    class ColumnEdgeGridCellIterator extends ColumnGridCellIterator implements Iterator<BlockCellColumn.GridCell>
    {
        public ColumnEdgeGridCellIterator(StsBlock block)
        {
            super(block);
        }

        public boolean hasNext()
        {
            if (gridCellIterator != null && gridCellIterator.hasNext())
                return true;
            else
            {
                while (blockColumnIterator.hasNext())
                {
                    BlockCellColumn blockColumn = blockColumnIterator.next();
                    if (blockColumn.type != EDGE) continue;
                    gridCellIterator = blockColumn.getGridCellIterator();
                    return true;
                }
            }
            return false;
        }
    }

    static public Comparator<BlockCellColumn.GridCell> getTruncatedFractionComparator()
    {
        return new TruncatedFractionComparator();
    }

    static class TruncatedFractionComparator implements Comparator<BlockCellColumn.GridCell>
    {
        public int compare(BlockCellColumn.GridCell cellA, BlockCellColumn.GridCell cellB)
        {
            return Float.compare(cellA.truncatedFraction, cellB.truncatedFraction);
        }
    }

    public BlockCellColumn constructor(byte type, int blockRow, int blockCol)
    {
        return new BlockCellColumn(type, blockRow, blockCol);
    }

    /**
     * A BlockCellColumn is a vertical prism of blockLines which trace down the edges of a column where a single grid cell intersects this block.
     * Most edges are at the corners or on the sides of the bounding grid cell.
     */
    public class BlockCellColumn
    {
        /** relative row location in block */
        public int blockRow;
        /** relative col location in block */
        public int blockCol;
        /** column type is either FULL, EMPTY, or an EDGE type.  If EDGE, then it contains some truncating grid polygons */
        byte type;
        /** gridPolygons which intersect this column. Two types: SURFACE which are from horizon intersections and GRID which are zoneSide intersections */
        ArrayList<StsPolygon> gridPolygons = new ArrayList<StsPolygon>();
        /** hash map of polygonLinkPoints in column keyed by gridRowCol */
        HashMap<StsGridRowCol, PolygonLinkPoint> linkHashMap = new HashMap<StsGridRowCol, PolygonLinkPoint>();
        /** A set of lines which trace down the edges of the column prism. */
        ArrayList<BlockLine> blockLines = new ArrayList<BlockLine>();
        /** Points from column polygons which don't fall on a row and or column grid. */
        ArrayList<PolygonLinkPoint> missedPoints = new ArrayList<PolygonLinkPoint>();
        /** Polygons for each layer at intersection with layer grids (nLayers+1 polygons). */
        StsPolygon[] layerPolygons;
        /** Grid cell geometry information */
        GridCell[] gridCells;

        public BlockCellColumn(byte type, int blockRow, int blockCol)
        {
            this.type = type;
            this.blockRow = blockRow;
            this.blockCol = blockCol;
            if (debugCell && debugBlockIJB(blockRow, blockCol))
                blockColumnDebug("creating column of type: " + cellTypeStrings[type]);
        }

        void addGridPolygon(StsPolygon gridPolygon)
        {
            if (!debugCheck(gridPolygon)) return;
            gridPolygons.add(gridPolygon);
        }

        private boolean debugCheck(StsPolygon gridPolygon)
        {
            if (debugCell && debugBlockIJB(blockRow, blockCol))
                blockColumnDebug("adding polygon: " + gridPolygon.toString());
            if (gridPolygon.row - rowMin == blockRow && gridPolygon.col - colMin == blockCol)
                return true;
            StsException.systemError(this, "debugPolygonCheck", "Attempting to add polygon row " + gridPolygon.row + " col " + gridPolygon.col +
                " to column at row " + blockRow + " col " + blockCol);
            return false;
        }

        void constructPrisms(StsZoneBlock[] zoneBlocks)
        {
            if (type != EDGE) return;
            if (debugCell && debugBlockIJB(blockRow, blockCol))
                blockColumnDebug("constructing prism. ");

            Collections.sort(gridPolygons);
            for (StsPolygon polygon : gridPolygons)
            {
                ArrayList<PolygonLinkPoint> links = createLinkPoints(polygon);
                // addRowAndColDuplicatePoints(links);
                for (PolygonLinkPoint link : links)
                    checkAddToBlockLine(link);
            }
            for (BlockLine blockLine : blockLines)
                blockLine.sortPoints();

            sortBlockLines();
            addMissingPoints();
            //TODO deal with any polygon points that weren't added: listed as missedPoints
            computeBlockLineGridPoints(zoneBlocks);
        }

        private void blockColumnDebug(String message)
        {
            StsException.systemDebug("BLOCK COLUMN DEBUG: block-" + getIndex() + " row " + debugCellRow + " col " + debugCellCol + " type " + cellTypeStrings[type] + " " + message);
        }

        private void addMissingPoints()
        {

        }

        ArrayList<PolygonLinkPoint> createLinkPoints(StsPolygon polygon)
        {
            int nPoints = polygon.nPoints;
            ArrayList<PolygonLinkPoint> linkPoints = new ArrayList<PolygonLinkPoint>(nPoints);
            for (int n = 0; n < nPoints; n++)
            {
                StsGridRowCol gridRowCol = polygon.gridRowCols[n];
                PolygonLinkPoint existingLink = linkHashMap.get(gridRowCol);
                if (existingLink != null)
                    linkPoints.add(existingLink);
                else
                {
                    PolygonLinkPoint link = new PolygonLinkPoint(polygon, gridRowCol, polygon.pntNrmls[n]);
                    linkHashMap.put(gridRowCol, link);
                    linkPoints.add(link);
                    /*
                    if(gridRowCol.isRowAndCol())
                    {
                        PolygonLinkPoint copyLink = new PolygonLinkPoint(link);
                        linkPoints.add(copyLink);
                    }
                    */
                }
            }
            connectLinks(polygon, linkPoints);
            return linkPoints;
        }

        private void addRowAndColDuplicatePoints(ArrayList<PolygonLinkPoint> links)
        {
            PolygonLinkPoint nextLink, link;

            PolygonLinkPoint firstLink = links.get(0);
            int nPoints = links.size();
            nextLink = firstLink;
            for (int n = 1; n < nPoints; n++)
            {
                link = nextLink;
                nextLink = links.get(n);
                addRowAndColDuplicatePoints(link, nextLink, links);
            }
            link = nextLink;
            nextLink = firstLink;
            addRowAndColDuplicatePoints(link, nextLink, links);
        }

        void addRowAndColDuplicatePoints(PolygonLinkPoint link, PolygonLinkPoint nextLink, ArrayList<PolygonLinkPoint> links)
        {
            PolygonLinkPoint newLink;

            if (!link.isRowAndCol()) return;
            int linkDirection = getUpOrDownLinkDirection(link, nextLink);
            // direction is down:
            // insert newLink after link and make newLink.downLink = nextLink
            // set link.downLink to null to effectively terminate this blockLine at this rowCol

            if (linkDirection == 1)
            {
                newLink = new PolygonLinkPoint(link);
                int index = links.indexOf(link);
                links.add(index + 1, newLink);
                link.setCwLink(newLink);
                newLink.setCcwLink(link);
                newLink.setDownLink(nextLink);
                nextLink.setUpLink(newLink);
                link.setDownLink(null);
            }
            // direction is up:
            // insert newLink after link and make newLink.upLink = nextLink
            // set link.upLink to null to effectively terminate this blockLine at this rowCol
            else if (linkDirection == -1) // direction is up
            {
                newLink = new PolygonLinkPoint(link);
                int index = links.indexOf(link);
                links.add(index + 1, newLink);
                link.setCwLink(newLink);
                newLink.setCcwLink(link);
                newLink.setUpLink(nextLink);
                nextLink.setDownLink(newLink);
                link.setUpLink(null);
            }
        }

        private void connectLinks(StsPolygon polygon, ArrayList<PolygonLinkPoint> links)
        {
            PolygonLinkPoint link, nextLink;
            PolygonLinkPoint firstLink = links.get(0);
            int nPoints = links.size();
            if (polygon.getType() == StsPolygon.SURFACE)
            {
                nextLink = firstLink;
                for (int n = 1; n < nPoints; n++)
                {
                    link = nextLink;
                    nextLink = links.get(n);
                    link.setCwLink(nextLink);
                    nextLink.setCcwLink(link);
                }
                link = nextLink;
                nextLink = firstLink;
                link.setCwLink(nextLink);
                nextLink.setCcwLink(link);
            }
            else // polygon.type == GRID
            {
                nextLink = firstLink;
                for (int n = 1; n < nPoints; n++)
                {
                    link = nextLink;
                    nextLink = links.get(n);
                    connectLinks(link, nextLink, links);
                }
                link = nextLink;
                nextLink = firstLink;
                connectLinks(link, nextLink, links);
            }
        }

        void connectLinksNew(PolygonLinkPoint link, PolygonLinkPoint nextLink, ArrayList<PolygonLinkPoint> links)
        {
            if (isCrossLink(link, nextLink))
            {
                if (link.getCwLink() != null)
                {
                    StsException.systemError(this, "connectLinks", link.toString() + " already has cwLink. Can't connect to: " + nextLink.toString());
                    return;
                }
                link.setCwLink(nextLink);
                nextLink.setCcwLink(link);
                return;
            }
            int linkDirection = getUpOrDownLinkDirection(link, nextLink);
            if (linkDirection == 1) // direction is down
            {
                nextLink.setUpLink(link);
                link.setDownLink(nextLink);
            }
            else if (linkDirection == -1)
            {
                nextLink.setDownLink(link);
                link.setUpLink(nextLink);
            }
        }

        void connectLinks(PolygonLinkPoint link, PolygonLinkPoint nextLink, ArrayList<PolygonLinkPoint> links)
        {
            if (isCrossLink(link, nextLink))
            {
                if (link.getCwLink() != null)
                {
                    StsException.systemError(this, "connectLinks", link.toString() + " already has cwLink. Can't connect to: " + nextLink.toString());
                    return;
                }
                link.setCwLink(nextLink);
                nextLink.setCcwLink(link);
                // return;
            }
            int linkDirection = getUpOrDownLinkDirection(link, nextLink);
            if (linkDirection == 1) // direction is down
            {
                if (link.getDownLink() != null)
                {
                    PolygonLinkPoint newLink = new PolygonLinkPoint(link);
                    int index = links.indexOf(link);
                    links.add(index + 1, newLink);
                    link.setCwLink(newLink);
                    newLink.setCcwLink(link);
                    newLink.setDownLink(nextLink);
                    nextLink.setUpLink(newLink);
                }
                else if (nextLink.getUpLink() != null)
                {
                    PolygonLinkPoint newLink = new PolygonLinkPoint(nextLink);
                    int index = links.indexOf(link);
                    links.add(index, newLink);
                    nextLink.setCcwLink(newLink);
                    newLink.setCwLink(nextLink);
                    newLink.setUpLink(link);
                    link.setDownLink(newLink);
                }
                else
                {
                    nextLink.setUpLink(link);
                    link.setDownLink(nextLink);
                }
            }
            else if (linkDirection == -1)
            {
                if (link.getUpLink() != null)
                {
                    PolygonLinkPoint newLink = new PolygonLinkPoint(link);
                    int index = links.indexOf(link);
                    links.add(index + 1, newLink);
                    link.setCwLink(newLink);
                    newLink.setCcwLink(link);
                    newLink.setUpLink(nextLink);
                    nextLink.setDownLink(newLink);
                }
                else if (nextLink.getDownLink() != null)
                {
                    PolygonLinkPoint newLink = new PolygonLinkPoint(nextLink);
                    int index = links.indexOf(link);
                    links.add(index, newLink);
                    link.setCcwLink(newLink);
                    newLink.setCwLink(link);
                    newLink.setDownLink(link);
                    link.setUpLink(newLink);
                }
                else
                {
                    nextLink.setDownLink(link);
                    link.setUpLink(nextLink);
                }
            }
        }

        boolean isCrossLink(PolygonLinkPoint link, PolygonLinkPoint nextLink)
        {
            return link.isRowAndCol() && nextLink.isRowAndCol();
        }

        int getUpOrDownLinkDirection(PolygonLinkPoint link, PolygonLinkPoint nextLink)
        {
            if (!link.isSameRowOrCol(nextLink)) return 0;
            float z = link.xyz[2];
            float nextZ = nextLink.xyz[2];
            return Float.compare(nextZ, z);
        }

        /**
         * At this point, all blockLines are on either a row or column or both.
         * There will be at most four groups corresponding to four sides of cells.
         * Sort them into the order: rowMin-increasingCol, colMax-increasingRow, rowMax-decreasingCol, colMin-decreasing row.
         */
        void sortBlockLines()
        {
            if (blockLines.size() == 0) return;
            // first find row and column minimum
            blockLineSortMinRowF = largeFloat;
            blockLineSortMinColF = largeFloat;
            for (BlockLine blockLine : blockLines)
            {
                blockLineSortMinRowF = Math.min(blockLineSortMinRowF, blockLine.rowF);
                blockLineSortMinColF = Math.min(blockLineSortMinColF, blockLine.colF);
            }
            Collections.sort(blockLines);
        }

        void checkAddToBlockLine(PolygonLinkPoint link)
        {
            StsGridRowCol gridRowCol;

            BlockLine blockLine;

            if (link.isRowAndCol())
            {
                blockLine = getBlockLine(link.gridRowCol);
                PolygonLinkPoint newLink = new PolygonLinkPoint(link);
                blockLine.addLinkPoint(newLink);
            }

            PolygonLinkPoint upLink = link.getUpLink();
            PolygonLinkPoint downLink = link.getDownLink();

            if (upLink == null && downLink == null)
            {
                gridRowCol = link.gridRowCol;
                if (!gridRowCol.isRowOrCol())
                {
                    missedPoints.add(link);
                    return;
                }
                blockLine = getBlockLine(gridRowCol);
                // if (blockLine.hasRowOrCol(gridRowCol)) return;
                blockLine.addLinkPoint(link);
            }
            else if (upLink != null)
            {
                gridRowCol = StsGridRowCol.getCommonGridRowCol(link.gridRowCol, upLink.gridRowCol);
                blockLine = upLink.getBlockLine();
                if (blockLine == null || !blockLine.matches(gridRowCol))
                {
                    blockLine = getBlockLine(gridRowCol);
                    blockLine.addLinkPoint(upLink);
                }
                blockLine.addLinkPoint(link);
            }
            else if (downLink != null)
            {
                gridRowCol = StsGridRowCol.getCommonGridRowCol(link.gridRowCol, downLink.gridRowCol);
                blockLine = downLink.getBlockLine();  // shouldn't exist debug check it
                if (blockLine == null || !blockLine.matches(gridRowCol))
                {
                    blockLine = getBlockLine(gridRowCol);
                    blockLine.addLinkPoint(link);
                    blockLine.addLinkPoint(downLink);
                }
                else
                {
                    blockLine.addLinkPoint(downLink);
                }
            }
        }

        /*
        void checkAddToBlockLine(StsPolygon polygon, StsGridRowCol gridRowCol, float[] xyz)
        {
            PolygonLinkPoint linkPoint;


            if (!gridRowCol.isRowOrCol())
            {
                linkPoint = new PolygonLinkPoint(polygon, gridRowCol, xyz);
                missedPoints.add(linkPoint);
                return;
            }
            BlockLine blockLine = getBlockLine(gridRowCol);
            if (!blockLine.hasRowOrCol(gridRowCol))
            {
                linkPoint = new PolygonLinkPoint(polygon, gridRowCol, xyz);
                blockLine.addLinkPoint(linkPoint);
            }
        }
        */
        BlockLine getBlockLine(PolygonLinkPoint link)
        {
            StsGridRowCol gridRowCol = link.gridRowCol;
            return getBlockLine(gridRowCol);
        }

        BlockLine getBlockLine(StsGridRowCol gridRowCol)
        {
            for (BlockLine blockLine : blockLines)
            {
                if (blockLine.matches(gridRowCol)) return blockLine;
            }
            BlockLine blockLine = new BlockLine(gridRowCol);
            blockLines.add(blockLine);
            return blockLine;
        }

        void computeBlockLineGridPoints(StsZoneBlock[] zoneBlocks)
        {
            for (BlockLine blockLine : blockLines)
            {
                blockLine.computeGridCrossingPoints();
                blockLine.computeLayerFs(zoneBlocks);
            }
        }

        void constructLayerGrids(int nTotalLayers)
        {
            try
            {
                if (type != EDGE) return;
                if (debugCell && debugBlockIJB(blockRow, blockCol))
                    StsException.systemDebug(this, "constructLayerGrids");
                for (BlockLine blockLine : blockLines)
                    blockLine.initializeLayerPoints();
                layerPolygons = new StsPolygon[nTotalLayers + 1];
                int nBlockLines = blockLines.size();
                BlockLine[] blockLineArray = new BlockLine[nBlockLines];
                blockLines.toArray(blockLineArray);
                float[] minBlockLineFs = new float[nBlockLines];
                float[] maxBlockLineFs = new float[nBlockLines];
                float minBlockLineF;
                float maxBlockLineF;
                for (int i = 0; i < nBlockLines; i++)
                {
                    minBlockLineFs[i] = blockLineArray[i].getMinBlockLineF();
                    maxBlockLineFs[i] = blockLineArray[i].getMaxBlockLineF();
                }
                minBlockLineF = StsMath.min(minBlockLineFs);
                maxBlockLineF = StsMath.max(maxBlockLineFs);
                if (minBlockLineF == nullValue || maxBlockLineF == nullValue) return;

                // possibly add a partial layer at the top
                int minPolygon = (int) minBlockLineF;
                int maxPolygon = (int) maxBlockLineF;
                if (minPolygon < minBlockLineF)
                {
                    layerPolygons[minPolygon] = constructLayerPolygon(minBlockLineF, blockLineArray);
                    minPolygon++;
                }
                for (int n = minPolygon; n <= maxPolygon; n++)
                    layerPolygons[n] = constructLayerPolygon(n, blockLineArray);
                if (maxPolygon < maxBlockLineF)
                    layerPolygons[maxPolygon + 1] = constructLayerPolygon(maxBlockLineF, blockLineArray);
            }
            catch (Exception e)
            {
                StsException.outputWarningException(this, "constructLayerGrids", e);
            }
        }

        StsPolygon constructLayerPolygon(float layerF, BlockLine[] blockLines)
        {
            int nBlockLines = blockLines.length;
            int nNonNullPolygonPoints = 0;
            float[][] polygonPoints = new float[nBlockLines][];
            for (int i = 0; i < nBlockLines; i++)
            {
                float[] xyz = blockLines[i].getLayerPoint(layerF);
                if (xyz != null) nNonNullPolygonPoints++;
                polygonPoints[i] = xyz;
            }
            if (nNonNullPolygonPoints == 0) return null;
            if (nNonNullPolygonPoints < nBlockLines)
                interpolateNullPoints(polygonPoints, nNonNullPolygonPoints);
            StsPolygon layerPolygon = new StsPolygon();
            layerPolygon.setPntNrmls(polygonPoints);
            return layerPolygon;
        }

        void constructLayerGridsOld(int nTotalLayers)
        {
            if (type != EDGE) return;
            for (BlockLine blockLine : blockLines)
                blockLine.initializeLayerPoints();
            layerPolygons = new StsPolygon[nTotalLayers + 1];
            int nBlockLines = blockLines.size();
            BlockLine[] blockLineArray = new BlockLine[nBlockLines];
            blockLines.toArray(blockLineArray);

            for (int n = 0; n <= nTotalLayers; n++)
            {
                float[][] polygonPoints = new float[nBlockLines][];
                int nNonNullPolygonPoints = 0;
                for (int i = 0; i < nBlockLines; i++)
                {
                    float[] xyz = blockLineArray[i].getLayerPoint(n);
                    if (xyz != null) nNonNullPolygonPoints++;
                    polygonPoints[i] = xyz;
                }
                if (nNonNullPolygonPoints == 0) continue;
                if (nNonNullPolygonPoints < nBlockLines)
                    interpolateNullPoints(polygonPoints, nNonNullPolygonPoints);
                StsPolygon layerPolygon = new StsPolygon();
                layerPolygon.pntNrmls = polygonPoints;
                layerPolygons[n] = layerPolygon;
            }
        }

        void constructGridCells()
        {
            if (type == EDGE)
            {
                if (debugCell && debugBlockIJB(blockRow, blockCol))
                    StsException.systemDebug(this, "constructEdgeGridCells", "called.");
                gridCells = new GridCell[sliceMax];
                StsPolygon botLayerPolygon = layerPolygons[0];
                for (int n = 0; n < sliceMax; n++)
                {
                    StsPolygon topLayerPolygon = botLayerPolygon;
                    botLayerPolygon = layerPolygons[n + 1];
                    if (topLayerPolygon == null || botLayerPolygon == null) continue;
                    gridCells[n] = new GridCell(n, topLayerPolygon, botLayerPolygon);
                }
            }
            else if (type == FULL)
            {
                gridCells = new GridCell[sliceMax];
                for (int n = 0; n < sliceMax; n++)
                    gridCells[n] = new GridCell(n);
            }
        }

        GridCell getGridCell(int nLayer)
        {
            if (type == EMPTY) return null;
            if (nLayer < 0 || nLayer >= gridCells.length)
            {
                StsException.systemError(this, "getGridCell", "nLayer outside range 0 to " + (gridCells.length - 1));
                return null;
            }
            return gridCells[nLayer];
        }

        GridCell getParentCellOrGridCell(int nLayer)
        {
            if (type == EMPTY) return null;
            GridCell gridCell = gridCells[nLayer];
            if (gridCell == null) return null;
            if (gridCell.parentCell != null)
                return gridCell.parentCell;
            else
                return gridCell;
        }

        double getGridCellParameter(int nLayer, byte direction)
        {
            if (gridCells[nLayer] == null) return 0.0;
            return gridCells[nLayer].getHalfLength(direction);
        }

        private void interpolateNullPoints(float[][] points, int nNonNullPolygonPoints)
        {
            int nPoints = points.length;
            IndexPoint[] indexPoints = new IndexPoint[nNonNullPolygonPoints];
            int nNonNull = 0;
            for (int n = 0; n < nPoints; n++)
                if (points[n] != null) indexPoints[nNonNull++] = new IndexPoint(n, points[n]);
            if (nNonNull == 1)
            {
                float[] nonNullPoint = indexPoints[0].point;
                for (int n = 0; n < nPoints; n++)
                    points[n] = nonNullPoint;
            }
            for (int n = 0; n < nPoints; n++)
            {
                if (points[n] == null)
                {
                    IndexPoint[] interpolatingPoints = getInterpolatingIndexPoints(n, indexPoints, nPoints);
                    if (interpolatingPoints == null) continue;
                    float di = interpolatingPoints[1].index - interpolatingPoints[0].index;
                    if (di < 0) di += nPoints;
                    float dj = n - interpolatingPoints[0].index;
                    if (dj < 0) dj += nPoints;
                    float f = dj / di;
                    points[n] = StsMath.interpolate(interpolatingPoints[0].point, interpolatingPoints[1].point, f);
                }
            }
        }

        IndexPoint[] getInterpolatingIndexPoints(int index, IndexPoint[] indexPoints, int nPoints)
        {
            IndexPoint prevPoint, nextPoint;
            int nIndexPoints = indexPoints.length;
            nextPoint = indexPoints[nIndexPoints - 1];
            for (int n = 0; n < nIndexPoints; n++)
            {
                prevPoint = nextPoint;
                nextPoint = indexPoints[n];
                int prevIndex = prevPoint.index;
                int nextIndex = nextPoint.index;
                int adjustedIndex = index;
                if (nextIndex < prevIndex)
                    nextIndex += nPoints;
                if (adjustedIndex < prevIndex)
                    adjustedIndex += nPoints;
                if (prevIndex < adjustedIndex && nextIndex > adjustedIndex)
                    return new IndexPoint[]{prevPoint, nextPoint};
            }
            StsException.systemError(this, "getInterpolatingIndexPoints", "failed to find interpolating points for index " + index);
            return null;
        }

        class IndexPoint
        {
            float[] point;
            int index;

            IndexPoint(int index, float[] point)
            {
                this.index = index;
                this.point = point;
            }
        }

        private float[] getPrevNonNullPoint(int index, float[][] points)
        {
            for (int n = index - 1; n >= 0; n--)
                if (points[n] != null) return points[n];
            int nPoints = points.length;
            for (int n = nPoints - 1; n > index; n--)
                if (points[n] != null) return points[n];
            return null;
        }

        private float[] getNextNonNullPoint(int index, float[][] points)
        {
            int nPoints = points.length;
            for (int n = index + 1; n < nPoints; n++)
                if (points[n] != null) return points[n];
            for (int n = 0; n < index; n++)
                if (points[n] != null) return points[n];
            return null;
        }

        private boolean isCellActive(int nLayer)
        {
            if (type == FULL) return true;
            else if (type == EMPTY) return false;
            else return layerPolygons[nLayer] != null && layerPolygons[nLayer + 1] != null;
        }

        public float computeEdgeCellDepth(int nLayer, float fillZ)
        {
            float depth = (float) computeFrustumCenterZ(layerPolygons[nLayer], layerPolygons[nLayer + 1]);
            if (depth == 0.0f) depth = fillZ;
            return depth;
        }

        double computeFrustumCenterZ(StsPolygon topPolygon, StsPolygon botPolygon)
        {
            if (topPolygon == null || botPolygon == null) return 0.0;
            double topZ = topPolygon.computeAverageZ();
            double botZ = botPolygon.computeAverageZ();
            return (topZ + botZ) / 2;
        }

        public double computeEdgeCellLayerPolygonArea(int nLayer)
        {
            return layerPolygons[nLayer].computePolygonArea();
        }

        public float getEdgeCellAreaXPlus(int nLayer)
        {
            GridCell gridCell = getGridCell(nLayer);
            if (gridCell == null) return nullValue;
            return gridCell.faceAreas[CELL_PLUS_DX];
        }

        public float getEdgeCellAreaYPlus(int nLayer)
        {
            GridCell gridCell = getGridCell(nLayer);
            if (gridCell == null) return nullValue;
            return gridCell.faceAreas[CELL_PLUS_DY];
        }

        public float getCellAreaZ(int nLayer)
        {
            if (type == FULL)
                return xInc * yInc;
            else // EDGE cellColumn
                return (float) computeEdgeCellLayerPolygonArea(nLayer);
        }

        public float getCellAreaZPlus(int nLayer)
        {
            if (type == FULL) return xInc * yInc;
            GridCell gridCell = getGridCell(nLayer);
            if (gridCell == null) return nullValue;
            float area = gridCell.faceAreas[CELL_PLUS_DZ];
            if (area == nullValue)
                StsException.systemError(this, "getComputeCellAreaZPlus", "value is null.");
            return area;
        }

        public float[] getEdgeCellCenter(int nLayer)
        {
            float[] topAvgCenter = layerPolygons[nLayer].computeCenter();
            float[] botAvgCenter = layerPolygons[nLayer + 1].computeCenter();
            return StsMath.interpolate(topAvgCenter, botAvgCenter, 0.5f);
        }

        public void computeCellPermeabilities(StsPropertyVolume[] blockPermeabilityVolumes)
        {
            if (type == EMPTY) return;
            StsPropertyVolume permXVolume = blockPermeabilityVolumes[0];
            StsPropertyVolume permYVolume = blockPermeabilityVolumes[1];
            StsPropertyVolume permZVolume = blockPermeabilityVolumes[2];

            float kx, ky, kz;
            for (int nLayer = 0; nLayer < sliceMax; nLayer++)
            {
                if (gridCells == null) continue;
                GridCell gridCell = gridCells[nLayer];
                if (gridCell == null) continue;
                kx = getPropertyVolumeValue(permXVolume, blockRow, blockCol, nLayer);
                ky = getPropertyVolumeValue(permYVolume, blockRow, blockCol, nLayer);
                kz = getPropertyVolumeValue(permZVolume, blockRow, blockCol, nLayer);
                gridCell.setKxyz(kx, ky, kz);
            }
        }

        public void computeCellPorosity(StsPropertyVolume blockPorosityVolume)
        {
            if (type == EMPTY) return;
            for (int nLayer = 0; nLayer < sliceMax; nLayer++)
            {
                GridCell gridCell = gridCells[nLayer];
                if (gridCell == null) continue;
                float porosity = getPropertyVolumeValue(blockPorosityVolume, blockRow, blockCol, nLayer);
                gridCell.setPorosity(porosity);
            }
        }

        /*
            public float computeEdgeCellVolume(int nLayer)
            {
                StsPolygon topPolygon = layerPolygons[nLayer];
                if (topPolygon == null) return 0.0f;
                StsPolygon botPolygon = layerPolygons[nLayer + 1];
                if (botPolygon == null) return 0.0f;
                int nPoints = topPoints.length;
                float[][] topPoints = topPolygon.pntNrmls;
                float[][] botPoints = botPolygon.pntNrmls;
                int nGoodPoints = 0;
                float[][] topGoodPoints = new float[nPoints][];
                float[][] botGoodPoints = new float[nPoints][];
                for(int n = 0; n < nPoints; n++)
                {
                    if(topPoints[n] != null || botPoints[n] != null)
                    {
                        topGoodPoints[nGoodPoints] = topPoints[n];
                        botGoodPoints[nGoodPoints] = botPoints[n];
                        nGoodPoints++;
                    }
                }
                float[] topPoint = topGoodPoints[nGoodPoints-1];
                float[] nextTopPoint = topGoodPoints[0];
                float[] botPoint = botGoodPoints[nGoodPoints-1];
                float[] nextBotPoint = botGoodPoints[0];
                for(int n = 0; n < nGoodPoints; n++)
                {
                    float[] prevTopPoint = topPoint;
                    topPoint = nextTopPoint;
                    nextTopPoint = topGoodPoints[(n+1)%nGoodPoints];
                    float[] prevBotPoint = botPoint;
                    botPoint = nextBotPoint;
                    nextBotPoint = botGoodPoints[(n+1)%nGoodPoints];
                    if(topGoodPoints[n] == null)

                }
                StsMath.copy(topPolygon.pntNrmls);
                float[][] botPoints = StsMath.copy(botPolygon.pntNrmls);

                return StsPolygon.computePolygonFrustumVolume(topPolygon, botPolygon);
            }
        */
        public float getCellAreaX(int nLayer)
        {
            StsPolygon topPolygon = layerPolygons[nLayer];
            if (topPolygon == null) return 0.0f;
            StsPolygon botPolygon = layerPolygons[nLayer + 1];
            if (botPolygon == null) return 0.0f;
            return 0.0f;
        }

        public String toString()
        {
            return StsParameters.cellTypeString(type);
        }

        class BlockLine implements Comparable<BlockLine>
        {
            int rowOrCol;
            float rowF;
            float colF;
            float rowMinF; // used in sorting BlockLines in clockwise order around square
            float colMinF; // used in sorting BlockLines in clockwise order around square
            StsPoint pointAbove; // layer number at indexAboveLast point
            StsPoint pointBelow; // layer number at indexAboveLast+1 point
            int indexAbove;
            /**
             * points along a BlockLine which define gridCrossing intersections. Used in layerGrid intersection calculations.
             * A gridCrossing.point coordinates are x,y,z, and f where f is the layer number
             */
            StsGridCrossingPoint[] gridCrossingPoints;
            int nGridCrossingPoints;

            ArrayList<PolygonLinkPoint> linkPoints = new ArrayList<PolygonLinkPoint>();

            BlockLine(StsGridRowCol rowCol)
            {
                this.rowOrCol = rowCol.rowOrCol;
                this.rowF = rowCol.getRowF();
                this.colF = rowCol.getColF();
            }

            boolean hasLinkPoint(PolygonLinkPoint link)
            {
                for (PolygonLinkPoint linkPoint : linkPoints)
                    if (linkPoint.gridRowCol == link.gridRowCol) return true;
                return false;
            }

            boolean hasRowOrCol(StsGridRowCol gridRowCol)
            {
                for (PolygonLinkPoint linkPoint : linkPoints)
                    if (linkPoint.gridRowCol == gridRowCol) return true;
                return false;
            }

            void addLinkPoint(PolygonLinkPoint linkPoint)
            {
                if (hasLinkPoint(linkPoint)) return;
                linkPoints.add(linkPoint);
                linkPoint.setBlockLine(this);
            }

            boolean matches(StsGridRowCol gridRowCol)
            {
                if (gridRowCol.rowOrCol != rowOrCol) return false;
                if (rowOrCol == StsGridRowCol.ROW)
                    return gridRowCol.rowF == rowF;
                else if (rowOrCol == StsGridRowCol.COL)
                    return gridRowCol.colF == colF;
                else if (rowOrCol == StsGridRowCol.ROWCOL)
                    return gridRowCol.rowF == rowF && gridRowCol.colF == colF;
                else
                    return false;
            }

            void sortPoints()
            {
                Collections.sort(linkPoints);
            }

            public int compareTo(BlockCellColumn.BlockLine other)
            {
                float perimeter = getPerimeter(this, blockRow, blockCol);
                float otherPerimeter = getPerimeter(other, blockRow, blockCol);
                int compare = Float.compare(perimeter, otherPerimeter);
                if (compare == 0)
                    StsException.systemError(this, "compareTo", "Failed comparing blockLines " + toString() + " and " + other.toString());
                return compare;
            }

            /** At this point, all blockLines being sorted are on either a row or column */
            int getGroup(BlockLine blockLine)
            {
                if (StsRowCol.isCol(blockLine.rowOrCol))
                {
                    if (blockLine.colF == blockLineSortMinColF) return 0;
                    else return 2;
                }
                else if (StsRowCol.isRow(blockLine.rowOrCol))
                {
                    if (blockLine.rowF == blockLineSortMinRowF) return 3;
                    else return 1;
                }
                else
                    return -1;
            }

            float getPerimeter(BlockCellColumn.BlockLine blockLine, int row, int col)
            {
                float rowFactor = blockLine.rowF - (row + rowMin);
                float colFactor = blockLine.colF - (col + colMin);
                if (rowFactor < 0.0f || rowFactor > 1.0f || colFactor < 0.0f || colFactor > 1.0f)
                {
                    StsException.systemError(this, "getPerimeter", "rowFactor: " + rowFactor + " colFactor: " + colFactor);
                    return 0.0f;
                }
                if (colFactor == 0.0f)
                    return rowFactor;
                else if (colFactor == 1.0f)
                    return 3.0f - rowFactor;
                else if (rowFactor == 0.0f)
                    return 4.0f - colFactor;
                else if (rowFactor == 1.0f)
                    return 1.0f + colFactor;
                else
                {
                    StsException.systemError(this, "getPerimeter", "Failed for blockLine: " + blockLine.toString());
                    return 0.0f;
                }
            }

            void computeGridCrossingPoints()
            {
                int nPoints = linkPoints.size();
                StsGridCrossingPoint[] linePoints = new StsGridCrossingPoint[nPoints];
                for (int n = 0; n < nPoints; n++)
                {
                    PolygonLinkPoint linkPoint = linkPoints.get(n);
                    StsGridRowCol gridPoint = linkPoint.gridRowCol;
                    linePoints[n] = new StsGridCrossingPoint(gridPoint.rowF, gridPoint.colF, linkPoint.xyz, StsBlock.this);
                }

                ArrayList<StsGridCrossingPoint> gridCrossingPointsArray = StsGridCrossings.computeGridCrossingPoints(linePoints, StsBlock.this);
                nGridCrossingPoints = gridCrossingPointsArray.size();
                gridCrossingPoints = new StsGridCrossingPoint[nGridCrossingPoints];
                gridCrossingPointsArray.toArray(gridCrossingPoints);
            }


            public int getPointPosition(float topZ, float botZ, float z)
            {
                if (z < topZ && !StsMath.sameAs(z, topZ))
                    return ABOVE;
                else if (z > botZ && !StsMath.sameAs(z, botZ))
                    return BELOW;
                else
                    return BETWEEN;
            }

            void computeLayerFs(StsZoneBlock[] zoneBlocks)
            {
                for (StsGridCrossingPoint gridCrossingPoint : gridCrossingPoints)
                {
                    float rowF = gridCrossingPoint.iF;
                    float colF = gridCrossingPoint.jF;
                    float z = gridCrossingPoint.point.getZ();
                    if (debugCell && debugBlockIJB(blockRow, blockCol))
                        StsException.systemDebug(this, "computeLayerFs");
                    gridCrossingPoint.point.setF(nullValue);
                    int nZoneBlocks = zoneBlocks.length;
                    for (int n = 0; n < nZoneBlocks; n++)
                    {
                        StsZoneBlock zoneBlock = zoneBlocks[n];
                        float topZ = zoneBlock.topGrid.getPointZ(rowF, colF);
                        if (topZ == nullValue) continue;
                        float botZ = zoneBlock.botGrid.getPointZ(rowF, colF);
                        if (botZ == nullValue) continue;
                        if (z <= botZ || n == nZoneBlocks - 1)
                        {
                            float layerF = zoneBlock.getLayerF(rowF, colF, z);
                            gridCrossingPoint.point.setF(layerF);
                            break;
                        }
                    }
                }
            }

            /*
                        void computeLayerFs(StsZoneBlock[] zoneBlocks)
                        {
                            int nZoneBlocks = zoneBlocks.length;
                            int nZoneBlock = 0;
                            StsZoneBlock zoneBlock;

                            for (StsGridCrossingPoint gridCrossingPoint : gridCrossingPoints)
                            {
                                float rowF = gridCrossingPoint.iF;
                                float colF = gridCrossingPoint.jF;
                                float z = gridCrossingPoint.point.getZ();
                                if (debugCell && debugBlockIJB(blockRow, blockCol))
                                    StsException.systemDebug(this, "computeLayerFs");
                                nZoneBlock = 0;
                                while(true)
                                {
                                    zoneBlock = zoneBlocks[nZoneBlock];
                                    float topZ = zoneBlock.topGrid.getComputePointZ(rowF, colF);
                                    float botZ = zoneBlock.botGrid.getComputePointZ(rowF, colF);
                                    int pointPosition = getPointPosition(topZ, botZ, z);

                                    if(pointPosition == BETWEEN)
                                        break;
                                    else if(pointPosition == ABOVE)
                                    {
                                        if(nZoneBlock == 0)
                                        {
                                            StsException.systemError(this, "computeLayerFs", " z " + z + " <  topZ " + topZ);
                                            z = topZ;
                                            break;
                                        }
                                        else nZoneBlock--;
                                    }
                                    else if(pointPosition == BELOW)
                                    {
                                        if(nZoneBlock == nZoneBlocks-1)
                                        {
                                            StsException.systemError(this, "computeLayerFs", " z " + z + " > botZ " + botZ);
                                            z = botZ;
                                            break;
                                        }
                                        else nZoneBlock++;
                                    }
                                }

                                float layerF = zoneBlock.getLayerF(rowF, colF, z);
                                gridCrossingPoint.point.setF(layerF);
                            }
                        }
            */
            void initializeLayerPoints()
            {
                indexAbove = 0;
                pointAbove = gridCrossingPoints[indexAbove].point;
                if (nGridCrossingPoints < 2)
                    pointBelow = pointAbove;
                else
                    pointBelow = gridCrossingPoints[indexAbove + 1].point;
            }

            float[] getLayerPoint(float layerF)
            {
                float fAbove = pointAbove.getF();
                if (fAbove > layerF) return null;
                float fBelow = pointBelow.getF();
                while (layerF > fBelow)
                {
                    indexAbove++;
                    if (indexAbove >= nGridCrossingPoints) return null;
                    pointAbove = pointBelow;
                    fAbove = fBelow;
                    pointBelow = gridCrossingPoints[indexAbove].point;
                    fBelow = pointBelow.getF();
                }
                if (layerF >= fBelow)
                    return pointBelow.v;
                else
                {
                    float f = (layerF - fAbove) / (fBelow - fAbove);
                    return StsMath.interpolate(pointAbove, pointBelow, f);
                }
            }

            float[] getLayerPointOld(int nLayer)
            {
                float fAbove = pointAbove.getF();
                if (fAbove > nLayer) return null;
                if (fAbove == nLayer) return pointAbove.v;
                float fBelow = pointBelow.getF();
                while (fBelow < nLayer)
                {
                    indexAbove++;
                    if (indexAbove >= nGridCrossingPoints) return null;
                    pointAbove = pointBelow;
                    fAbove = fBelow;
                    pointBelow = gridCrossingPoints[indexAbove].point;
                    fBelow = pointBelow.getF();
                }
                if (fBelow <= fAbove)
                {
                    StsException.systemError(this, "getLayerPoint", "fBelow " + fBelow + " <= fAbove " + fAbove);
                    return null;
                }
                float f = (nLayer - fAbove) / (fBelow - fAbove);
                return StsMath.interpolate(pointAbove, pointBelow, f);
            }

            float getMinBlockLineF()
            {
                for (int n = 0; n < gridCrossingPoints.length; n++)
                {
                    float pointF = gridCrossingPoints[n].point.getF();
                    if (pointF != nullValue) return pointF;
                }
                return nullValue;
            }

            float getMaxBlockLineF()
            {
                for (int n = gridCrossingPoints.length - 1; n >= 0; n--)
                {
                    float pointF = gridCrossingPoints[n].point.getF();
                    if (pointF != nullValue) return pointF;
                }
                return nullValue;
            }

            public String toString()
            {
                return new String("blockLine at: rowF " + rowF + " colF: " + colF + " on: " + StsParameters.rowCol(rowOrCol));
            }
        }

        /** used for testing only */
        public void constructGridCells(int nCells)
        {
            gridCells = new GridCell[nCells];
            for (int n = 0; n < nCells; n++)
                gridCells[n] = new GridCell(n);
        }

        public Iterator<GridCell> getGridCellIterator()
        {
            ArrayList<GridCell> gridCellList = new ArrayList<GridCell>(Arrays.asList(gridCells));
            return gridCellList.iterator();
        }

        public class GridCell
        {
            public byte type = FULL;
            public int nLayer;
            public float[] cellCenter;
            public float topCenterZ;
            public float botCenterZ;
            public float volume;
            public float poreVolume;
            public float[] faceAreas;
            public float[] halfSizes;
            public float[] tranXYZ;
            public float truncatedFraction = 1.0f;
            private float[] kxyz;
            private float porosity;
            public ParentCell parentCell;
            public ArrayList<StsNNC> nncList;

            static final float largeFloat = StsParameters.largeFloat;
            static final float nullValue = StsParameters.nullValue;

            public GridCell()
            {
            }

            /** constructed used for full column grid cells. */
            public GridCell(int nLayer)
            {
                this.nLayer = nLayer;
                initializeTransArray();
                initializeGeometryArrays();
                computeFullCellGeometry();
            }

            /** constructor used for edge column grid cells */
            public GridCell(int nLayer, StsPolygon topLayerPolygon, StsPolygon botLayerPolygon)
            {
                try
                {
                    type = EDGE;
                    this.nLayer = nLayer;
                    initializeTransArray();
                    initializeGeometryArrays();
                    computeEdgeCellGeometry();
                }
                catch (Exception e)
                {
                    StsException.outputWarningException(this, "constructor", e);
                }
            }

            protected void initializeTransArray()
            {
                tranXYZ = new float[6];
                // Arrays.fill(tranXYZ, nullValue);
            }

            protected void initializeGeometryArrays()
            {
                halfSizes = new float[6];
                // Arrays.fill(halfSizes, nullValue);
                faceAreas = new float[6];
                // Arrays.fill(faceAreas, nullValue);
            }

            public StsBlock getBlock(){ return StsBlock.this; }

            void computeCellGeometry()
            {
                if (type == FULL)
                    computeFullCellGeometry();
                else
                    computeEdgeCellGeometry();
            }

            void computeEdgeCellGeometry()
            {
                if (debugBlockIJK())
                    StsException.systemDebug(this, "computeEdgeCellGeometry");
                if (layerPolygons == null) return;
                StsPolygon topPolygon = layerPolygons[nLayer];
                if (topPolygon == null) return;
                StsPolygon botPolygon = layerPolygons[nLayer + 1];
                if (botPolygon == null) return;
                PolygonGeometry topPolygonGeometry = new PolygonGeometry(topPolygon);
                topCenterZ = topPolygonGeometry.center[2];
                setFaceArea(CELL_MINUS_DZ, topPolygonGeometry.area);
                PolygonGeometry botPolygonGeometry = new PolygonGeometry(botPolygon);
                botCenterZ = botPolygonGeometry.center[2];
                setFaceArea(CELL_PLUS_DZ, botPolygonGeometry.area);
                cellCenter = StsMath.interpolate(topPolygonGeometry.center, botPolygonGeometry.center, 0.5f);
                if (Float.isNaN(cellCenter[0]))
                    StsException.systemError(this, "computeCellGeometry", "value is a NaN");
                computeEdgeCellXFaces();
                computeEdgeCellYFaces();

                double topArea = topPolygonGeometry.area;
                double botArea = botPolygonGeometry.area;
                double topZ = topPolygon.computeAverageZ();
                double botZ = botPolygon.computeAverageZ();
                volume = (float) ((1.0 / 3.0) * (topArea + botArea + Math.sqrt(topArea * botArea)) * (botZ - topZ));
                if (volume < 0.0)
                {
                    StsException.systemError(this, "computeEdgeCellGeometry", "volume is negative");
                    return;
                }

                adjustHalfSizes();

                truncatedFraction = (float) ((topArea + botArea) / cellXYArea / 2);
                if (truncatedFraction > 0.99f) truncatedFraction = 1.0f;
            }

            protected void computeFullCellGeometry()
            {
                if (debugBlockIJK())
                    StsException.systemDebug(this, "computeFullCellGeometry");
                float[][] topLayerGrid = StsBlock.this.getLayerGrid(nLayer);
                float[][] botLayerGrid = StsBlock.this.getLayerGrid(nLayer + 1);
                faceAreas = new float[6];
                halfSizes = new float[6];
                setFaceArea(0, StsBlock.this.getFullCellAreaXPlus(blockRow, blockCol, topLayerGrid, botLayerGrid));
                setFaceArea(1, StsBlock.this.getFullCellAreaXMinus(blockRow, blockCol, topLayerGrid, botLayerGrid));
                setFaceArea(2, StsBlock.this.getFullCellAreaYPlus(blockRow, blockCol, topLayerGrid, botLayerGrid));
                setFaceArea(3, StsBlock.this.getFullCellAreaYMinus(blockRow, blockCol, topLayerGrid, botLayerGrid));
                setFaceArea(4, xInc * yInc);
                setFaceArea(5, xInc * yInc);
                setHalfSize(0, xInc/2);
                setHalfSize(1, xInc/2);
                setHalfSize(2, yInc/2);
                setHalfSize(3, yInc/2);
                topCenterZ = getLayerAvgZ(blockRow, blockCol, topLayerGrid);
                botCenterZ = getLayerAvgZ(blockRow, blockCol, botLayerGrid);
                setHalfSize(4, (botCenterZ - topCenterZ) / 2);
                setHalfSize(5, (botCenterZ - topCenterZ) / 2);
                volume = xInc * yInc * (botCenterZ - topCenterZ);
            }

            void addVolume(float dVolume)
            {
                if(debugBlockIJK())
                     StsException.systemDebug(this, "addVolume", "adjusted from " + volume + " to " + (volume + dVolume));
                volume += dVolume;
                poreVolume = porosity*volume;
                truncatedFraction = (float)(volume/((botCenterZ - topCenterZ)*cellXYArea));
                if(dVolume < 0.0f && truncatedFraction < 0.1f)
                    StsException.systemDebug(this, "addVolume", "truncated Fraction small " + truncatedFraction + " for cell " + toString());
            }

            void setFaceArea(int direction, float area)
            {
                if(debugBlockIJK())
                     StsException.systemDebug(this, "setFaceArea", "area: " + area);
                if (area < 0.0f)
                    StsException.systemError(this, "setFaceArea", "faceArea < 0");
                faceAreas[direction] = area;
            }

            public void addFaceArea(byte direction, float area)
            {
                if(debugBlockIJK())
                    StsException.systemDebug(this, "addFaceArea", "adjusted from " + faceAreas[direction] + " to " + (faceAreas[direction] + area));
                if (area < 0.0f)
                    StsException.systemError(this, "setFaceArea", "faceArea < 0");
                faceAreas[direction] += area;
            }

            public float getFaceArea(byte direction)
            {
                if (faceAreas != null) return faceAreas[direction];
                StsException.systemError(this, "getFaceArea", "face area undefined");
                return 0.0f;
            }

            void adjustHalfSizes()
            {
                adjustHalfSize(CELL_PLUS_DX, CELL_MINUS_DX);
                adjustHalfSize(CELL_PLUS_DY, CELL_MINUS_DY);
                adjustHalfSize(CELL_PLUS_DZ, CELL_MINUS_DZ);
            }

            void adjustHalfSize(byte plusDirection, byte minusDirection)
            {
                float avgArea = (faceAreas[plusDirection] + faceAreas[minusDirection]) / 2;
                if (avgArea <= 0.0f && volume > 0.0f)
                    StsException.systemError(this, "adjustHalfSize", "avgArea <= 0.0f, but volume > 0.0f");
                if (avgArea == 0.0f) return;
                float size = volume / avgArea;
                halfSizes[plusDirection] = halfSizes[minusDirection] = size / 2;
            }

            void setHalfSize(int direction, float halfSize)
            {
                if (halfSize < 0.0f)
                    halfSize = 0.0f;
                halfSizes[direction] = halfSize;
            }

            void computeEdgeCellXFaces()
            {
                StsPolygon topLayerPolygon = layerPolygons[nLayer];
                if (topLayerPolygon == null) return;
                StsPolygon botLayerPolygon = layerPolygons[nLayer + 1];
                if (botLayerPolygon == null) return;

                float[][] topColSegment = topLayerPolygon.getColSegment(blockCol + 1, StsBlock.this);
                if (topColSegment == null) return;
                float[][] botColSegment = botLayerPolygon.getColSegment(blockCol + 1, StsBlock.this);
                if (botColSegment == null) return;
                float topYLength = Math.abs(topColSegment[0][1] - topColSegment[1][1]);
                float botYLength = Math.abs(botColSegment[0][1] - botColSegment[1][1]);
                float dz0 = botColSegment[0][2] - topColSegment[0][2];
                float dz1 = botColSegment[1][2] - topColSegment[1][2];
                float area = (topYLength + botYLength) * (dz0 + dz1) / 4;
                setFaceArea(CELL_PLUS_DX, area);

                topColSegment = topLayerPolygon.getColSegment(blockCol, StsBlock.this);
                if (topColSegment == null) return;
                botLayerPolygon = layerPolygons[nLayer + 1];
                if (botLayerPolygon == null) return;
                botColSegment = botLayerPolygon.getColSegment(blockCol, StsBlock.this);
                if (botColSegment == null) return;
                topYLength = Math.abs(topColSegment[0][1] - topColSegment[1][1]);
                botYLength = Math.abs(botColSegment[0][1] - botColSegment[1][1]);
                dz0 = botColSegment[0][2] - topColSegment[0][2];
                dz1 = botColSegment[1][2] - topColSegment[1][2];
                area = (topYLength + botYLength) * (dz0 + dz1) / 4;
                setFaceArea(CELL_MINUS_DX, area);
            }

            public void computeEdgeCellYFaces()
            {
                StsPolygon topLayerPolygon = layerPolygons[nLayer];
                if (topLayerPolygon == null) return;
                StsPolygon botLayerPolygon = layerPolygons[nLayer + 1];
                if (botLayerPolygon == null) return;

                float[][] topRowSegment = topLayerPolygon.getRowSegment(blockRow + 1, StsBlock.this);
                if (topRowSegment == null) return;
                float[][] botRowSegment = botLayerPolygon.getRowSegment(blockRow + 1, StsBlock.this);
                if (botRowSegment == null) return;
                float topXLength = Math.abs(topRowSegment[0][0] - topRowSegment[1][0]);
                float botXLength = Math.abs(botRowSegment[0][0] - botRowSegment[1][0]);
                float dz0 = botRowSegment[0][2] - topRowSegment[0][2];
                float dz1 = botRowSegment[1][2] - topRowSegment[1][2];
                float area = (topXLength + botXLength) * (dz0 + dz1) / 4;
                setFaceArea(CELL_PLUS_DY, area);

                topRowSegment = topLayerPolygon.getRowSegment(blockRow, StsBlock.this);
                if (topRowSegment == null) return;
                botRowSegment = botLayerPolygon.getRowSegment(blockRow, StsBlock.this);
                if (botRowSegment == null) return;
                topXLength = Math.abs(topRowSegment[0][0] - topRowSegment[1][0]);
                botXLength = Math.abs(botRowSegment[0][0] - botRowSegment[1][0]);
                dz0 = botRowSegment[0][2] - topRowSegment[0][2];
                dz1 = botRowSegment[1][2] - topRowSegment[1][2];
                area = (topXLength + botXLength) * (dz0 + dz1) / 4;
                setFaceArea(CELL_MINUS_DY, area);
            }

            public ParentCell getCreateParentCell()
            {
                if (parentCell != null) return parentCell;
                return createParentCell();
            }

            public ParentCell createParentCell()
            {
                parentCell = new ParentCell(this);
                parentCells.add(parentCell);
                parentCell.addChild(this);
                return parentCell;
            }

            public boolean checkSetParentCell()
            {
                if (debugBlockIJK() || debugBlockIJKNeighbor(debugDirection))
                    StsException.systemDebug(this, "checkSetParentCell");
                if (isParent()) return false;
                GridCell bestNeighborCell = getBestNeighborCell();
                if (bestNeighborCell == null)
                {
                    StsException.systemError(this, "checkSetParentCell", "failed to find cell to merge with " + toString());
                    return false;
                }
                ParentCell parentCell = bestNeighborCell.getCreateParentCell();
                parentCell.addChild(this);
                return true;
            }

            private void mergeNeighborParentCell(GridCell bestNeighborCell)
            {
                ParentCell neighborParentCell = bestNeighborCell.parentCell;
                parentCells.remove(neighborParentCell);
                ArrayList<GridCell> neighborChildren = neighborParentCell.children;
                for (GridCell neighborChild : neighborChildren)
                    parentCell.addChild(neighborChild);
            }

            protected boolean isParent(){ return false; }

            protected boolean debugBlockIJK()
            {
                return debugBlockIJK(blockRow, blockCol, nLayer);
            }

            /** Debug neighbor cell in defined direction from current debug row, col,layer. */
            protected boolean debugBlockIJKNeighbor(byte direction)
            {
                switch (direction)
                {
                    case CELL_NONE:
                        return false;
                    case CELL_PLUS_DX:
                        return debugBlockIJK(blockRow, blockCol - 1, nLayer);
                    case CELL_MINUS_DX:
                        return debugBlockIJK(blockRow, blockCol + 1, nLayer);
                    case CELL_PLUS_DY:
                        return debugBlockIJK(blockRow - 1, blockCol, nLayer);
                    case CELL_MINUS_DY:
                        return debugBlockIJK(blockRow + 1, blockCol, nLayer);
                    default:
                        return debugBlockIJK(blockRow, blockCol, nLayer);
                }
            }

            protected int[] getNeighborBlockIJK(byte direction)
            {
                switch (direction)
                {
                    case CELL_PLUS_DX:
                        return new int[]{blockRow, blockCol + 1, nLayer};
                    case CELL_MINUS_DX:
                        return new int[]{blockRow, blockCol - 1, nLayer};
                    case CELL_PLUS_DY:
                        return new int[]{blockRow + 1, blockCol, nLayer};
                    case CELL_MINUS_DY:
                        return new int[]{blockRow - 1, blockCol, nLayer};
                    case CELL_PLUS_DZ:
                        return new int[]{blockRow, blockCol, nLayer + 1};
                    case CELL_MINUS_DZ:
                        return new int[]{blockRow, blockCol, nLayer - 1};
                    default:
                        return new int[]{blockRow, blockCol, nLayer};
                }
            }

            private boolean debugBlockIJK(int blockRow, int blockCol, int nLayer)
            {
                return StsBlock.debugIJK(blockRow + rowMin, blockCol + colMin, nLayer);
            }

            /**
             * First try to merge with a neighbor cell which is not already a parent: keeps children to a max of two;
             * if this fails, then find the best parent to merge with (parent with fewest children)
             */
            private GridCell getBestNeighborCell()
            {
                GridCell bestNeighborCell = getBestNeighborCell(false);
                if (bestNeighborCell != null) return bestNeighborCell;
                return getBestNeighborCell(true);
            }

            private GridCell getBestNeighborCell(boolean parentsOnly)
            {
                float bestTruncatedFraction = 0.0f;
                GridCell bestNeighborCell = null;
                for (byte direction = 0; direction < 4; direction++)
                {
                    GridCell gridCell = getRowColNeighborCell(direction);
                    if (gridCell == null) continue;
                    if (parentsOnly != gridCell.isParent()) continue;
                    // float trans = getTrans(direction);
                    // if(trans == 0.0f) continue;
                    float truncatedFraction = gridCell.truncatedFraction;
                    if (truncatedFraction > bestTruncatedFraction)
                    {
                        bestNeighborCell = gridCell;
                        bestTruncatedFraction = truncatedFraction;
                    }
                }
                return bestNeighborCell;
            }

            private GridCell getBestNeighborChildCell()
            {
                float bestTruncatedFraction = 0.0f;
                GridCell bestNeighborCell = null;
                for (byte direction = 0; direction < 4; direction++)
                {
                    GridCell gridCell = getRowColNeighborCell(direction);
                    if (gridCell == null) continue;
                    if (gridCell.isParent()) continue;
                    // float trans = getTrans(direction);
                    // if(trans == 0.0f) continue;
                    float truncatedFraction = gridCell.truncatedFraction;
                    if (truncatedFraction > bestTruncatedFraction)
                    {
                        bestNeighborCell = gridCell;
                        bestTruncatedFraction = truncatedFraction;
                    }
                }
                return bestNeighborCell;
            }

            private GridCell getBestNeighborParentCell()
            {
                int bestMinChildren = 1000;
                GridCell bestNeighborCell = null;
                for (byte direction = 0; direction < 4; direction++)
                {
                    GridCell gridCell = getRowColNeighborCell(direction);
                    if (gridCell == null) continue;
                    if (!gridCell.isParent()) continue;
                    // float trans = getTrans(direction);
                    // if(trans == 0.0f) continue;
                    int nChildren = ((ParentCell) gridCell).children.size();
                    if (nChildren < bestMinChildren)
                    {
                        bestNeighborCell = gridCell;
                        bestMinChildren = nChildren;
                    }
                }
                return bestNeighborCell;
            }

            private GridCell getRowColNeighborCell(byte direction)
            {
                int row = blockRow + neighborDRowCol[direction][0];
                int col = blockCol + neighborDRowCol[direction][1];
                return getBlockGridCell(row, col, nLayer);
            }

            public boolean isNeighbor(GridCell neighborCell)
            {
                int dBlockRow = neighborCell.getBlockRow() - getBlockRow();
                int dBlockCol = neighborCell.getBlockCol() - getBlockCol();
                int dBlockLayer = neighborCell.nLayer - nLayer;

                if (dBlockRow == 1 || dBlockRow == -1)
                    return dBlockCol == 0 && dBlockLayer == 0;
                else if (dBlockCol == 1 || dBlockCol == -1)
                    return dBlockRow == 0 && dBlockLayer == 0;
                else if (dBlockLayer == 1 || dBlockLayer == -1)
                    return dBlockRow == 0 && dBlockCol == 0;
                else
                    return false;
            }

            /**
             * Merge this child to the parent: 1) add volume and poreVolume; 2) set child inactive; 3) get trans from parent to child;
             * 3) harmonic avg this trans with child NNCs, and connect these NNCs with parent.
             */
            /*
            public void mergeChild(GridCell childCell)
            {
                childCell.nncList = null;
                childCell.setCellInactive();
            }
            */
            public void connectNncToParent(StsNNC nnc, GridCell parentCell, float trans)
            {
                if (nnc.leftGridCell == this)
                    nnc.leftGridCell = parentCell;
                else if (nnc.rightGridCell == this)
                    nnc.rightGridCell = parentCell;
                nnc.trans = nnc.trans * trans / (nnc.trans + trans);
                if (nnc.trans < 0.0f)
                {
                    StsException.systemError(this, "connectNncToParent", "nnc.trans < 0.0f");
                    nnc.trans = 0.0f;
                }
                if (StsNNC.isDebug(nnc))
                {
                    StsException.systemDebug(this, "connectNncToParent", nnc.toString());
                }
            }

            protected void setCellInactive(StsEclipseModel eclipseModel)
            {
                StsPropertyVolume actnumVolume = eclipseModel.getBlockPropertyVolume(StsBlock.this, eclipseModel.actnum);
                actnumVolume.setBlockValue(blockRow, blockCol, nLayer, 0.0f);
                StsPropertyVolume poreVolume = eclipseModel.getBlockPropertyVolume(StsBlock.this, eclipseModel.poreVolume);
                poreVolume.setBlockValue(blockRow, blockCol, nLayer, 0.0f);
            }

            protected void setCellActive(StsEclipseModel eclipseModel)
            {
                StsPropertyVolume actnumVolume = eclipseModel.getBlockPropertyVolume(StsBlock.this, eclipseModel.actnum);
                actnumVolume.setBlockValue(blockRow, blockCol, nLayer, 1.0f);
                StsPropertyVolume poreVolume = eclipseModel.getBlockPropertyVolume(StsBlock.this, eclipseModel.poreVolume);
                poreVolume.setBlockValue(blockRow, blockCol, nLayer, this.poreVolume);
            }

            protected void setCellDepth(StsEclipseModel eclipseModel)
            {
                StsPropertyVolume depthVolume = eclipseModel.getBlockPropertyVolume(StsBlock.this, eclipseModel.depth);
                depthVolume.setBlockValue(blockRow, blockCol, nLayer, cellCenter[2]);
            }

            protected float getTrans(byte direction, StsEclipseModel eclipseModel)
            {
                StsPropertyVolume transVolume;

                if (debugBlockIJK(blockRow, blockCol, nLayer) || debugBlockIJK(blockRow, blockCol - 1, nLayer))
                {
                    StsException.systemDebug(this, "getTrans");
                }
                switch (direction)
                {
                    case CELL_PLUS_DX:
                        transVolume = eclipseModel.getBlockPropertyVolume(StsBlock.this, eclipseModel.tranX);
                        return transVolume.getBlockValue(blockRow, blockCol, nLayer);
                    case CELL_MINUS_DX:
                        transVolume = eclipseModel.getBlockPropertyVolume(StsBlock.this, eclipseModel.tranX);
                        return transVolume.getBlockValue(blockRow, blockCol - 1, nLayer);
                    case CELL_PLUS_DY:
                        transVolume = eclipseModel.getBlockPropertyVolume(StsBlock.this, eclipseModel.tranY);
                        return transVolume.getBlockValue(blockRow, blockCol, nLayer);
                    case CELL_MINUS_DY:
                        transVolume = eclipseModel.getBlockPropertyVolume(StsBlock.this, eclipseModel.tranY);
                        return transVolume.getBlockValue(blockRow - 1, blockCol, nLayer);
                    case CELL_PLUS_DZ:
                        transVolume = eclipseModel.getBlockPropertyVolume(StsBlock.this, eclipseModel.tranZ);
                        return transVolume.getBlockValue(blockRow, blockCol, nLayer);
                    case CELL_MINUS_DZ:
                        transVolume = eclipseModel.getBlockPropertyVolume(StsBlock.this, eclipseModel.tranZ);
                        return transVolume.getBlockValue(blockRow, blockCol, nLayer - 1);
                    default:
                        return 0.0f;
                }
            }

            protected void setTrans(byte direction, float value, StsEclipseModel eclipseModel)
            {
                StsPropertyVolume transVolume;

                switch (direction)
                {
                    case CELL_PLUS_DX:
                        transVolume = eclipseModel.getBlockPropertyVolume(StsBlock.this, eclipseModel.tranX);
                        transVolume.setBlockValue(blockRow, blockCol, nLayer, value);
                        break;
                    case CELL_MINUS_DX:
                        transVolume = eclipseModel.getBlockPropertyVolume(StsBlock.this, eclipseModel.tranX);
                        transVolume.setBlockValue(blockRow, blockCol - 1, nLayer, value);
                        break;
                    case CELL_PLUS_DY:
                        transVolume = eclipseModel.getBlockPropertyVolume(StsBlock.this, eclipseModel.tranY);
                        transVolume.setBlockValue(blockRow, blockCol, nLayer, value);
                        break;
                    case CELL_MINUS_DY:
                        transVolume = eclipseModel.getBlockPropertyVolume(StsBlock.this, eclipseModel.tranY);
                        transVolume.setBlockValue(blockRow - 1, blockCol, nLayer, value);
                        break;
                    case CELL_PLUS_DZ:
                        transVolume = eclipseModel.getBlockPropertyVolume(StsBlock.this, eclipseModel.tranZ);
                        transVolume.setBlockValue(blockRow, blockCol, nLayer, value);
                        break;
                    case CELL_MINUS_DZ:
                        transVolume = eclipseModel.getBlockPropertyVolume(StsBlock.this, eclipseModel.tranZ);
                        transVolume.setBlockValue(blockRow, blockCol, nLayer - 1, value);
                        break;
                    default:
                        StsException.systemError(this, "setTrans", "undefined direction " + direction);
                }
            }

            public void computeCellPermeability()
            {
                GridCell gridCell = gridCells[nLayer];
                if (gridCell == null) return;

                StsPropertyVolume permXVolume, permYVolume, permZVolume;
                StsZoneBlock zoneBlock = StsBlock.this.getZoneBlock(nLayer);
                permXVolume = zoneBlock.getEclipsePropertyVolume("PERMX");
                permYVolume = zoneBlock.getEclipsePropertyVolume("PERMY");
                permZVolume = zoneBlock.getEclipsePropertyVolume("PERMZ");
                float kx = getPropertyVolumeValue(permXVolume, blockRow, blockCol, nLayer);
                float ky = getPropertyVolumeValue(permYVolume, blockRow, blockCol, nLayer);
                float kz = getPropertyVolumeValue(permZVolume, blockRow, blockCol, nLayer);
                gridCell.setKxyz(kx, ky, kz);
            }

            private boolean isEdgeCell()
            {
                return nncList != null;
            }

            protected byte getChildDirection(GridCell childCell)
            {
                if (childCell.debugBlockIJK())
                {
                    StsException.systemDebug(this, "getChildDirection");
                }
                if (blockRow == childCell.getBlockRow())
                {
                    if (blockCol < childCell.getBlockCol())
                        return CELL_PLUS_DX;
                    else
                        return CELL_MINUS_DX;
                }
                else if (blockCol == childCell.getBlockCol())
                {
                    if (blockRow < childCell.getBlockRow())
                        return CELL_PLUS_DY;
                    else
                        return CELL_MINUS_DY;
                }
                return 0;
            }

            public int getBlockRow(){ return blockRow; }

            public int getBlockCol(){ return blockCol; }

            public int getRow(){ return blockRow + rowMin; }

            public int getCol(){ return blockCol + colMin; }

            public int[] getIJKB()
            {
                return new int[]{blockRow + rowMin, blockCol + colMin, nLayer, StsBlock.this.getIndex()};
            }

            public int[] getBlockIJK()
            {
                return new int[]{blockRow, blockCol, nLayer};
            }

            public int getEclipseIndex(int blockRow, int blockCol, int slice, int nEclipseRows)
            {
                if (parentCell != null)
                    return parentCell.getEclipseIndex(blockRow, blockCol, slice, nEclipseRows);
                else
                {
                    int[] eclipseRowCol = getEclipseRowColFromBlockRowCol(blockRow, blockCol);
                    int eclipseRow = eclipseRowCol[0] - 1;
                    int eclipseCol = eclipseRowCol[1] - 1;
                    return slice * nCellCols * nEclipseRows + nCellCols * eclipseRow + eclipseCol;
                }
            }

            class PolygonGeometry
            {
                float area;
                float[] center;

                PolygonGeometry(StsPolygon polygon)
                {
                    computePolygonGeometry(polygon);
                }

                void computePolygonGeometry(StsPolygon polygon)
                {
                    float[][] points = polygon.pntNrmls;
                    float xc = 0, yc = 0, zc = 0;
                    int nPoints = points.length;
                    float[] point = points[nPoints - 1];
                    float x = 0.0f;
                    float y = 0.0f;
                    float x0 = point[0];
                    float y0 = point[1];
                    float z = point[2];
                    float zSum = 0.0f;
                    area = 0.0f;
                    for (int n = 0; n < nPoints; n++)
                    {
                        float xp = x;
                        float yp = y;
                        point = points[n];
                        x = point[0] - x0;
                        y = point[1] - y0;
                        z = point[2];
                        float ya = (y + yp);
                        float xa = (x + xp);
                        float xy = x * yp - xp * y;
                        xc += xa * xy;
                        yc += ya * xy;
                        area += xy;
                        zSum += z;
                    }
                    if (area > 0.0f)
                    {
                        area /= 2;
                        xc /= (6 * area);
                        xc += x0;
                        yc /= (6 * area);
                        yc += y0;
                    }
                    else
                    {
                        area = 0.0f;
                        center = points[0];
                        return;
                    }
                    zc = zSum / nPoints;
                    center = new float[]{xc, yc, zc};
                }
            }

            private boolean isBetween(float f)
            {
                return f >= 0.0f && f < 1.0f;
            }

            public void setKxyz(float kx, float ky, float kz)
            {
                if (debugBlockIJK())
                    StsException.systemDebug(this, "setPermeabilities");
                kxyz = new float[]{kx, ky, kz};
            }

            public float getPermeability(byte direction)
            {
                switch (direction)
                {
                    case CELL_PLUS_DX:
                    case CELL_MINUS_DX:
                        return kxyz[0];
                    case CELL_PLUS_DY:
                    case CELL_MINUS_DY:
                        return kxyz[1];
                    case CELL_PLUS_DZ:
                    case CELL_MINUS_DZ:
                        return kxyz[2];
                    default:
                        return kxyz[0];
                }
            }

            public void setPorosity(float porosity)
            {
                this.porosity = porosity;
                poreVolume = porosity * volume;
            }

            public float getPorosity(){ return porosity; }

            public void setKxyz(float[] kxyz)
            {
                this.kxyz = kxyz;
            }

            public float[] getKxyz()
            {
                if (kxyz != null) return kxyz;
                computeCellPermeability();
                return kxyz;
            }

            public float[] getCellCenter()
            {
                if (cellCenter == null)
                {
                    StsBlock block = StsBlock.this;
                    float[][] topLayerGrid = block.getLayerGrid(nLayer);
                    float[][] botLayerGrid = block.getLayerGrid(nLayer + 1);
                    topCenterZ = getLayerAvgZ(blockRow, blockCol, topLayerGrid);
                    if (topCenterZ == nullValue) return null;
                    botCenterZ = getLayerAvgZ(blockRow, blockCol, botLayerGrid);
                    if (botCenterZ == nullValue) return null;
                    float avgZ = (topCenterZ + botCenterZ) / 2;
                    float xCenter = getXCoor(blockCol + 0.5f);
                    float yCenter = getYCoor(blockRow + 0.5f);
                    cellCenter = new float[]{xCenter, yCenter, avgZ};
                }
                return cellCenter;
            }

            public float getCellDepth()
            {
                float[] cellCenter = getCellCenter();
                return cellCenter[2];
            }

            public float getCellSize(byte direction)
            {
                if (halfSizes == null)
                {
                    halfSizes = new float[6];
                    halfSizes[0] = halfSizes[1] = xInc;
                    halfSizes[2] = halfSizes[3] = yInc;
                }
                if (direction == CELL_MINUS_DX || direction == CELL_PLUS_DX)
                    return halfSizes[0] + halfSizes[1];
                else if (direction == CELL_MINUS_DY || direction == CELL_PLUS_DY)
                    return halfSizes[2] + halfSizes[3];
                else
                    return halfSizes[4] + halfSizes[5];
            }

            public float[] getCellCenterVector(GridCell neighborCell)
            {
                float[] neighborCenter = neighborCell.getCellCenter();
                return StsMath.subtract(neighborCenter, getCellCenter());
            }

            public double getHalfLength(byte direction)
            {
                switch (direction)
                {
                    case CELL_PLUS_DX:
                    case CELL_MINUS_DX:
                    case CELL_PLUS_DY:
                    case CELL_MINUS_DY:
                        return halfSizes[direction];
                    case CELL_PLUS_DZ:
                        return botCenterZ - cellCenter[2];
                    case CELL_MINUS_DZ:
                        return cellCenter[2] - topCenterZ;
                    default:
                        StsException.systemError(this, "getParameter", "Undefined Cell direction " + direction);
                        return 0.0;
                }
            }

            public double getHalfLength(byte parameter, GridCell childCell)
            {
                int blockCol = childCell.getBlockCol();
                int blockRow = childCell.getBlockRow();
                switch (parameter)
                {
                    case CELL_PLUS_DX:
                        return getBlockXCoor(blockCol + 1) - cellCenter[0];
                    case CELL_MINUS_DX:
                        return cellCenter[0] - getBlockXCoor(blockCol);
                    case CELL_PLUS_DY:
                        return getBlockYCoor(blockRow + 1) - cellCenter[1];
                    case CELL_MINUS_DY:
                        return cellCenter[1] - getBlockYCoor(blockRow);
                    case CELL_PLUS_DZ:
                        return botCenterZ - cellCenter[2];
                    case CELL_MINUS_DZ:
                        return cellCenter[2] - topCenterZ;
                    default:
                        StsException.systemError(this, "getParameter", "Undefined Cell parameter " + parameter);
                        return 0.0;
                }
            }

            public float computeVolume()
            {
                StsPolygon topPolygon = layerPolygons[nLayer];
                if (topPolygon == null) return 0.0f;
                StsPolygon botPolygon = layerPolygons[nLayer + 1];
                if (botPolygon == null) return 0.0f;
                double topArea = topPolygon.computePolygonArea();
                double botArea = botPolygon.computePolygonArea();
                double topZ = topPolygon.computeAverageZ();
                double botZ = botPolygon.computeAverageZ();
                volume = (float) ((1.0 / 3.0) * (topArea + botArea + Math.sqrt(topArea * botArea)) * (botZ - topZ));
                if (volume < 0.0)
                {
                    StsException.systemError(StsPolygon.class, "computePolygonFrustumVolume", "volume is negative");
                    return 0.0f;
                }
                truncatedFraction = (float) ((topArea + botArea) / cellXYArea / 2);
                if (truncatedFraction > 0.99f) truncatedFraction = 1.0f;
                return volume;
            }

            public void addGridCellTran(byte direction, float trans, StsEclipseModel eclipseModel)
            {
                switch (direction)
                {
                    case CELL_PLUS_DX:
                        addPropertyVolumeBlockValue(eclipseModel.tranX, blockRow, blockCol, nLayer, trans, eclipseModel);
                        break;
                    case CELL_MINUS_DX:
                        addPropertyVolumeBlockValue(eclipseModel.tranX, blockRow, blockCol - 1, nLayer, trans, eclipseModel);
                        break;
                    case CELL_PLUS_DY:
                        addPropertyVolumeBlockValue(eclipseModel.tranY, blockRow, blockCol, nLayer, trans, eclipseModel);
                        break;
                    case CELL_MINUS_DY:
                        addPropertyVolumeBlockValue(eclipseModel.tranY, blockRow - 1, blockCol, nLayer, trans, eclipseModel);
                        break;
                    case CELL_PLUS_DZ:
                        addPropertyVolumeBlockValue(eclipseModel.tranY, blockRow, blockCol, nLayer, trans, eclipseModel);
                        break;
                    case CELL_MINUS_DZ:
                        addPropertyVolumeBlockValue(eclipseModel.tranY, blockRow - 1, blockCol, nLayer, trans, eclipseModel);
                        break;
                }
            }

            public void setGridCellPoreVolume(float poreVolume)
            {
                this.poreVolume = poreVolume;
            }

            public void setGridCellTran(byte direction, float trans)
            {
                if (debugBlockIJK() || debugBlockIJKNeighbor(debugDirection))
                    StsException.systemDebug(this, "setGridCellTran");
                // if(tranXYZ[direction] != 0.0f)
                //    StsException.systemError(this, "setGridCellTran", "Tran value already set.");
                tranXYZ[direction] = trans;
            }

            public void addNNC(StsNNC nnc)
            {
                if (nncList == null) nncList = new ArrayList<StsNNC>();
                nncList.add(nnc);
            }

            public String toString()
            {
                return "row " + (blockRow + rowMin) + " col " + (blockCol + colMin) + " layer " + nLayer + " block " + StsBlock.this.getIndex() + " volume " + volume + " truncated Fraction " + truncatedFraction;
            }

            public String toIJKBString()
            {
                return new String("" + (blockRow + rowMin) + "," + (blockCol + colMin) + "," + nLayer + "," + (StsBlock.this.getIndex()));
            }
        }

        public class ParentCell extends GridCell
        {
            public ArrayList<GridCell> children = new ArrayList<GridCell>();
            ArrayList<StsNNC> neighborNncList;
            CellAdjustment[] cellAdjustments;

            ParentCell(GridCell initialCell)
            {
                blockRow = initialCell.getBlockRow();
                blockCol = initialCell.getBlockCol();
                this.nLayer = initialCell.nLayer;
                setKxyz(initialCell.getKxyz());
                cellCenter = initialCell.cellCenter;
                topCenterZ = initialCell.topCenterZ;
                botCenterZ = initialCell.botCenterZ;
            }

            void addChild(GridCell child)
            {
                if (children == null) children = new ArrayList<GridCell>();
                children.add(child);
                child.parentCell = this;
            }

            protected boolean isParent(){ return true; }

            public void mergeOrAdjust(StsEclipseModel eclipseModel, StsEclipseOutput eclipseOutput)
            {
                if (eclipseOutput.adjustCells)
                    adjust(eclipseModel, eclipseOutput);
                else
                    merge(eclipseModel, eclipseOutput);
            }

            public void adjust(StsEclipseModel eclipseModel, StsEclipseOutput eclipseOutput)
            {
                if (debugBlockIJK() || debugBlockIJKNeighbor(debugDirection))
                    StsException.systemDebug(this, "merge");

                try
                {
                    // for (GridCell child : children)
                    //    child.computeCellGeometry();
                    createCellAdjustments();

                    if (poreVolume < minPoreVolume)
                    {
                        StsMessageFiles.logMessage("Merged cell pore volume still < min (" + minPoreVolume + "): " + toString() + " making inactive.");
                        setCellInactive(eclipseModel);
                        return;
                    }
                }
                catch (Exception e)
                {
                    StsException.outputWarningException(this, "merge", e);
                }
            }

            private boolean createCellAdjustments()
            {
                int nChildren = children.size();
                float totalVolume = 0.0f;
                for (GridCell child : children)
                    totalVolume += child.volume;
                float avgVolume = totalVolume / nChildren;

                int nAdjustments = children.size() - 1;
                cellAdjustments = new CellAdjustment[nAdjustments];
                GridCell parentCell = children.get(0);
                float parentVolume = parentCell.volume;
                // parentVolume can't be less than avgVolume
                if (parentVolume <= avgVolume) return false;
                for (int n = 0; n < nAdjustments; n++)
                {
                    GridCell childCell = children.get(n + 1);
                    float dVolume = avgVolume - childCell.volume;
                    float fraction = 1.0f - dVolume / parentVolume;
                    cellAdjustments[n] = new CellAdjustment(parentCell, childCell, fraction);
                }
                return true;
            }

            /**
             * Interface between parent and child cell has been moved in one of four directions: + or - X or Y.
             * The parent cell has been shrunk in that direction and the child cell has grown.
             * The change is expressed as a fraction of the dimension of the parent size, f.
             * Interface change distance: delta = f*(parentCell.minusHalfSize + parentCell.plusHalfSize)
             * Reduce each of the parent halfSizes by delta/2.
             * Shift the parent center by delta/2.
             * Reduce parent volume and poreVolume by f.
             * Increase the child halfSizes by delta/2.
             * Shift the child center by delta/2.
             * Increase the child volume and poreVolume by amount parent is reduced.
             * Decrease parent top and bottom areas by f.
             * Increase child top and bottom areas by f.
             * Recompute parent and child trans and NNCs.
             */
            class CellAdjustment
            {
                byte direction;
                float fraction;

                CellAdjustment(GridCell parentCell, GridCell childCell, float fraction)
                {
                    direction = parentCell.getChildDirection(childCell);
                    this.fraction = fraction;
                    float dVolume = parentCell.volume * (1.0f - fraction);
                    parentCell.addVolume(-dVolume);
                    childCell.addVolume(dVolume);
                    float directionFraction = adjustFaceAreas(parentCell, childCell, direction, fraction);
                    float delta = parentCell.getCellSize(direction) * directionFraction;
                    adjustLateralFaces(parentCell, childCell, direction, directionFraction);
                    shiftCenter(parentCell, childCell, direction, delta / 2);
                    parentCell.adjustHalfSizes();
                    childCell.adjustHalfSizes();
                }

                float adjustFaceAreas(GridCell parentCell, GridCell childCell, byte direction, float f)
                {
                    if (direction == CELL_PLUS_DX)
                        return adjustFaceArea(parentCell, childCell, CELL_PLUS_DX, CELL_MINUS_DX, f);
                    else if (direction == CELL_MINUS_DX)
                        return adjustFaceArea(parentCell, childCell, CELL_MINUS_DX, CELL_PLUS_DX, f);
                    else if (direction == CELL_PLUS_DY)
                        return adjustFaceArea(parentCell, childCell, CELL_PLUS_DY, CELL_MINUS_DY, f);
                    else if (direction == CELL_MINUS_DY)
                        return adjustFaceArea(parentCell, childCell, CELL_MINUS_DY, CELL_PLUS_DY, f);
                    else
                        return 0.0f;
                }

                float adjustFaceArea(GridCell parentCell, GridCell childCell, byte plusDirection, byte minusDirection, double f)
                {
                    float a0 = parentCell.faceAreas[plusDirection];
                    float a1 = parentCell.faceAreas[minusDirection];
                    if(StsMath.sameAs(a0, a1)) return (float)f;
                    float area = (float) Math.sqrt(a1 * a1 * (1 - f) + a0 * a0 * f);
                    parentCell.setFaceArea(plusDirection, area);
                    childCell.setFaceArea(minusDirection, area);
                    return (area - a0) / (a1 - a0);
                }


                public void shiftCenter(GridCell parentCell, GridCell childCell, byte direction, double adjustment)
                {
                    if (direction == CELL_PLUS_DX)
                    {
                        parentCell.cellCenter[0] -= adjustment;
                        childCell.cellCenter[0] -= adjustment;
                    }
                    else if (direction == CELL_MINUS_DX)
                    {
                        parentCell.cellCenter[0] += adjustment;
                        childCell.cellCenter[0] += adjustment;
                    }
                    else if (direction == CELL_PLUS_DY)
                    {
                        parentCell.cellCenter[1] -= adjustment;
                        childCell.cellCenter[1] -= adjustment;
                    }
                    else if (direction == CELL_MINUS_DY)
                    {
                        parentCell.cellCenter[1] += adjustment;
                        childCell.cellCenter[1] += adjustment;
                    }
                }

                void adjustLateralFaces(GridCell parentCell, GridCell childCell, byte direction, float fraction)
                {
                    adjustLateralFace(parentCell, childCell, CELL_MINUS_DZ, fraction);
                    adjustLateralFace(parentCell, childCell, CELL_PLUS_DZ, fraction);
                    if (direction == CELL_PLUS_DX || direction == CELL_MINUS_DX)
                    {
                        adjustLateralFace(parentCell, childCell, CELL_PLUS_DY, fraction);
                        adjustLateralFace(parentCell, childCell, CELL_MINUS_DY, fraction);
                    }
                    else if (direction == CELL_PLUS_DY || direction == CELL_MINUS_DY)
                    {
                        adjustLateralFace(parentCell, childCell, CELL_PLUS_DX, fraction);
                        adjustLateralFace(parentCell, childCell, CELL_MINUS_DX, fraction);
                    }
                }

                void adjustLateralFace(GridCell parentCell, GridCell childCell, byte direction, float fraction)
                {
                    float area = parentCell.faceAreas[direction];
                    float dArea = area * (1.0f - fraction);
                    parentCell.addFaceArea(direction, -dArea);
                    childCell.addFaceArea(direction, dArea);
                }
            }

            public void merge(StsEclipseModel eclipseModel, StsEclipseOutput eclipseOutput)
            {
                if (debugBlockIJK() || debugBlockIJKNeighbor(debugDirection))
                    StsException.systemDebug(this, "merge");

                try
                {
                    computeCellGeometry();
                    connectMergeParentCell(eclipseOutput.outputWriter);

                    // merge children to parent (combine volumes and adjust child NNCs)
                    nncList = new ArrayList<StsNNC>();
                    for (GridCell child : children)
                    {
                        if (child.nncList != null)
                        {
                            nncList.addAll(child.nncList);
                            child.nncList = null;
                        }
                        child.setCellInactive(eclipseModel);
                    }

                    if (poreVolume < minPoreVolume)
                    {
                        StsMessageFiles.logMessage("Merged cell pore volume still < min (" + minPoreVolume + "): " + toString() + " making inactive.");
                        setCellInactive(eclipseModel);
                        return;
                    }
                    setCellActive(eclipseModel);
                    setCellDepth(eclipseModel);

                    for (byte direction = 0; direction < 6; direction++)
                        if (tranXYZ[direction] != nullValue)
                            setTrans(direction, tranXYZ[direction], eclipseModel);

                    // add the child neighbors which are neighbor NNCs relative to parent
                    if (neighborNncList != null)
                    {
                        if (StsBlock.this.neighborNncList == null)
                            StsBlock.this.neighborNncList = new ArrayList<StsNNC>();
                        StsBlock.this.neighborNncList.addAll(neighborNncList);
                        neighborNncList = null;
                    }
                }
                catch (Exception e)
                {
                    StsException.outputWarningException(this, "merge", e);
                }
            }

            void computeCellGeometry()
            {
                if (debugBlockIJK() || debugBlockIJKNeighbor(debugDirection))
                    StsException.systemDebug(this, "computeCellGeometry");
                try
                {
                    // weight cell center and top and bottom face centers by each child poreVolume
                    cellCenter = new float[3];
                    double[] center = new double[3];
                    double topCenter = 0.0;
                    double botCenter = 0.0;
                    for (GridCell child : children)
                    {
                        poreVolume += child.poreVolume;
                        volume += child.volume;
                        float[] childCellCenter = child.getCellCenter();
                        for (int i = 0; i < 3; i++)
                            center[i] += child.volume * childCellCenter[i];
                        topCenter += child.volume * child.topCenterZ;
                        botCenter += child.volume * child.botCenterZ;
                    }
                    for (int i = 0; i < 3; i++)
                        cellCenter[i] = (float) (center[i] / volume);
                    topCenterZ = (float) (topCenter / volume);
                    botCenterZ = (float) (botCenter / volume);
                }
                catch (Exception e)
                {
                    StsException.outputWarningException(this, "computeCellGeometry", e);
                }
            }

            public void connectMergeParentCell(StsEclipseOutput.EclipseWriter outputWriter)
            {
                if (debugBlockIJK() || debugBlockIJKNeighbor(debugDirection))
                    StsException.systemDebug(this, "connectParentCell");

                for (GridCell child : children)
                {
                    for (byte n = 0; n < 6; n++)
                        addChildTran(n, child, outputWriter);
                }
            }

            private void addChildTran(byte childNeighborDirection, GridCell child, StsEclipseOutput.EclipseWriter outputWriter)
            {
                float childNeighborTran = child.tranXYZ[childNeighborDirection];
                if (childNeighborTran == 0.0f) return;
                GridCell childNeighborCell = getChildNeighbor(child, childNeighborDirection);
                if (childNeighborCell == null) return;
                if (hasChildCell(childNeighborCell)) return;

                float[] childNeighborVector = child.getCellCenterVector(childNeighborCell);
                float[] parentChildVector = getCellCenterVector(child);
                float[] parentNeighborVector = StsMath.getAddVectors(parentChildVector, childNeighborVector);
                double parentNeighborLength = StsMath.normalize(parentNeighborVector);
                double childNeighborLength = StsMath.normalize(childNeighborVector);
                float lengthCorrection = (float) (childNeighborLength / parentNeighborLength);

                float[] tranNormal = getCellOrthogonalNormal(childNeighborDirection);
                float areaCorrection = Math.abs(StsMath.dot(parentNeighborVector, tranNormal));
                float correctionFactor = lengthCorrection * areaCorrection;
                float trans = childNeighborTran * correctionFactor;

                double childHalfLength = child.getHalfLength(childNeighborDirection);
                double parentHalfLength = getHalfLength(childNeighborDirection, child);

                float oldCorrectionFactor;
                if (parentHalfLength <= 0.0f)
                {
                    StsException.systemError(this, "addChildTran", "parentHalfLength<=0.0f for parent cell: " + toString());
                    oldCorrectionFactor = 1.0f;
                }
                else
                {
                    if (childHalfLength < 0.0)
                    {
                        StsException.systemError(this, "addChildTran", "childHalfLength<0.0f for parent cell: " + toString());
                        childHalfLength = 0.0;
                    }
                    oldCorrectionFactor = (float) (childHalfLength / parentHalfLength);
                }
                outputWriter.println("Parent Cell " + toString());
                if (isRegularConnection(childNeighborCell))
                {
                    if (tranXYZ[childNeighborDirection] != nullValue)
                    {
                        StsException.systemError(this, "addChild", "trans already set for cell " + toString() + " direction " + childNeighborDirection);
                        outputWriter.println("    trans already set: " + tranXYZ[childNeighborDirection]);
                    }
                    outputWriter.println("    add reg trans " + cellDirectionStrings[childNeighborDirection] + " " + trans + " child trans " + childNeighborTran + " length " + parentNeighborLength + " length factor: " + lengthCorrection + " area factor: " + areaCorrection);
                    tranXYZ[childNeighborDirection] = trans;
                }
                else
                {
                    GridCell neighborNncCell = getNeighborConnection(childNeighborCell);
                    if (debugBlockIJK() || childNeighborCell.debugBlockIJK())
                        StsException.systemDebug(this, "addChildTran: create Interblock NNC");
                    StsNNC nnc = new StsNNC(this, neighborNncCell, trans);
                    if (!hasNeighborNNC(neighborNncCell))
                    {
                        addNeighborNNC(nnc);
                        outputWriter.println("    add neighbor NNC " + nnc.toString());
                        outputWriter.println("   " + cellDirectionStrings[childNeighborDirection] + " child trans " + childNeighborTran + " length " + parentNeighborLength + " length factor: " + lengthCorrection + " area factor: " + areaCorrection);
                    }
                    else
                    {
                        outputWriter.println("    couldn't add neighbor NNC; already exists: " + nnc);
                    }
                }
            }

            private float[] getCellOrthogonalNormal(int direction)
            {
                return cellOrthogonalNormals[direction];
            }

            private GridCell getChildNeighbor(GridCell child, byte direction)
            {
                int[] neighborBlockIJK = child.getNeighborBlockIJK(direction);
                return StsBlock.this.getBlockGridCell(neighborBlockIJK);
            }

            private boolean hasChildCell(GridCell childCell)
            {
                for (GridCell child : children)
                    if (child == childCell) return true;
                return false;
            }

            private boolean isRegularConnection(GridCell neighborCell)
            {
                if (neighborCell.parentCell != null)
                    return isRegularConnection(neighborCell.parentCell);
                return isNeighbor(neighborCell);
            }

            private GridCell getNeighborConnection(GridCell neighborCell)
            {
                if (neighborCell.parentCell == null)
                    return neighborCell;
                else
                    return neighborCell.parentCell;
            }

            private boolean hasNeighborNNC(GridCell otherCell)
            {
                if (neighborNncList != null)
                {
                    for (StsNNC nnc : neighborNncList)
                        if (nnc.rightGridCell == otherCell || nnc.leftGridCell == otherCell) return true;
                }
                if (otherCell instanceof ParentCell)
                    return ((ParentCell) otherCell).hasParentInterblockNNC(this);

                return false;
            }

            private boolean hasParentInterblockNNC(GridCell otherCell)
            {
                if (neighborNncList != null)
                {
                    for (StsNNC nnc : neighborNncList)
                        if (nnc.rightGridCell == otherCell || nnc.leftGridCell == otherCell) return true;
                }
                return false;
            }

            private void addNeighborNNC(StsNNC nnc)
            {
                if (nnc.leftGridCell.debugBlockIJK() || nnc.rightGridCell.debugBlockIJK())
                    StsException.systemDebug(this, "addInterblockNNC");

                if (neighborNncList == null) neighborNncList = new ArrayList<StsNNC>();
                neighborNncList.add(nnc);
            }

            public float computeTruncatedFraction()
            {
                truncatedFraction = 0.0f;
                for (GridCell child : children)
                    truncatedFraction += child.truncatedFraction;
                return truncatedFraction;
            }

            public String toString()
            {
                if (children.size() == 2)
                {
                    GridCell first = children.get(0);
                    int firstRow = first.getRow();
                    int firstCol = first.getCol();
                    GridCell second = children.get(1);
                    int secondRow = second.getRow();
                    int secondCol = second.getCol();
                    if (firstRow == secondRow)
                        return " row " + firstRow + " cols " + firstCol + " - " + secondCol + " layer " + nLayer + " block " + StsBlock.this.getIndex() + " volume " + volume + " truncated fraction " + truncatedFraction;
                    else
                        return " col " + firstCol + " rows " + firstRow + " - " + secondRow + " layer " + nLayer + " block " + StsBlock.this.getIndex() + " volume " + volume + " truncated fraction " + truncatedFraction;
                }
                else
                {
                    GridCell first = children.get(0);
                    int rowMin = first.getBlockRow();
                    int rowMax = rowMin;
                    int colMin = first.getBlockCol();
                    int colMax = colMin;
                    int nChildren = children.size();
                    for (int n = 1; n < nChildren; n++)
                    {
                        GridCell child = children.get(n);
                        int row = child.getRow();
                        int col = child.getCol();
                        if (row < rowMin) rowMin = row;
                        else if (row > rowMax) rowMax = row;
                        if (col < colMin) colMin = col;
                        else if (col > colMax) colMax = col;
                    }
                    return nChildren + " children: rows " + rowMin + "-" + rowMax + " cols " + colMin + "-" + colMax + " layer " + nLayer + " block " + StsBlock.this.getIndex() + " volume " + volume + " truncated fraction " + truncatedFraction;
                }
            }

            public String toIJKBString()
            {
                if (children.size() == 2)
                {
                    GridCell first = children.get(0);
                    int firstRow = first.getRow();
                    int firstCol = first.getCol();
                    GridCell second = children.get(1);
                    int secondRow = second.getRow();
                    int secondCol = second.getCol();
                    if (firstRow == secondRow)
                        return firstRow + "," + firstCol + "-" + secondCol + "," + nLayer + "," + StsBlock.this.getIndex();
                    else
                        return firstRow + "-" + secondRow + "," + firstCol + "," + nLayer + "," + StsBlock.this.getIndex();
                }
                else
                {
                    GridCell first = children.get(0);
                    int firstRow = first.getRow();
                    int lastRow = firstRow;
                    int firstCol = first.getCol();
                    int lastCol = firstCol;
                    int nChildren = children.size();
                    for (int n = 1; n < nChildren; n++)
                    {
                        GridCell child = children.get(n);
                        int row = child.getRow();
                        int col = child.getCol();
                        if (row < firstRow) firstRow = row;
                        else if (row > lastRow) lastRow = row;
                        if (col < firstCol) firstCol = col;
                        else if (col > lastCol) lastCol = col;
                    }
                    return firstRow + "-" + lastRow + "," + firstCol + "-" + lastCol + "," + nLayer + "," + StsBlock.this.getIndex() + "(" + nChildren + ")";
                }
            }
        }
    }

    class PolygonLinkPoint implements Comparable<PolygonLinkPoint>
    {
        StsGridRowCol gridRowCol;
        StsPolygon polygon;
        float[] xyz;
        private PolygonLinkPoint downLink;
        private PolygonLinkPoint upLink;
        private PolygonLinkPoint cwLink;
        private PolygonLinkPoint ccwLink;
        private BlockCellColumn.BlockLine blockLine;

        PolygonLinkPoint(StsPolygon polygon, StsGridRowCol gridRowCol, float[] xyz)
        {
            this.polygon = polygon;
            this.gridRowCol = gridRowCol;
            this.xyz = xyz;
        }

        PolygonLinkPoint(PolygonLinkPoint linkPoint)
        {
            gridRowCol = linkPoint.gridRowCol;
            polygon = linkPoint.polygon;
            xyz = linkPoint.xyz;
        }

        public String toString()
        {
            if (blockLine != null)
                return polygon.toString() + " " + gridRowCol.toString() + " " + blockLine.toString() + " z: " + xyz[2];
            else
                return polygon.toString() + " " + gridRowCol.toString() + " z: " + xyz[2];
        }

        public boolean isRowAndCol()
        {
            return gridRowCol.isRowAndCol();
        }

        public boolean isSameRowOrCol(PolygonLinkPoint nextLink)
        {
            return gridRowCol.isSameRowOrCol(nextLink.gridRowCol);
        }

        public int compareTo(PolygonLinkPoint other)
        {
            float z = xyz[2];
            float zOther = other.xyz[2];
            return Float.compare(z, zOther);
        }

        public PolygonLinkPoint getDownLink()
        {
            return downLink;
        }

        public void setDownLink(PolygonLinkPoint downLink)
        {
            if (this.downLink != null)
            {
                StsException.systemError(this, "setDownLink", "not null");
                return;
            }
            this.downLink = downLink;
        }

        public PolygonLinkPoint getUpLink()
        {
            return upLink;
        }

        public void setUpLink(PolygonLinkPoint upLink)
        {
            if (this.upLink != null)
            {
                StsException.systemError(this, "setUpLink", "not null");
                return;
            }
            this.upLink = upLink;
        }

        public PolygonLinkPoint getCwLink()
        {
            return cwLink;
        }

        public void setCwLink(PolygonLinkPoint cwLink)
        {
            if (this.cwLink != null)
            {
                StsException.systemError(this, "setCwLink", "not null");
                return;
            }
            this.cwLink = cwLink;
        }

        public PolygonLinkPoint getCcwLink()
        {
            return ccwLink;
        }

        public void setCcwLink(PolygonLinkPoint ccwLink)
        {
            if (this.ccwLink != null)
            {
                StsException.systemError(this, "setCcwLink", "not null");
                return;
            }
            this.ccwLink = ccwLink;
        }

        public BlockCellColumn.BlockLine getBlockLine()
        {
            return blockLine;
        }

        public void setBlockLine(BlockCellColumn.BlockLine blockLine)
        {
            /*
            if (this.blockLine != null && this.blockLine != blockLine)
            {
                StsException.systemError(this, "setBlockLine", "Block line already set for this link: " + toString());
            }
            */
            this.blockLine = blockLine;
        }
    }
}
