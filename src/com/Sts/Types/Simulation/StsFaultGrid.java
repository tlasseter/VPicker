package com.Sts.Types.Simulation;

import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;
import com.Sts.Actions.Export.*;
import com.Sts.UI.*;
import com.Sts.Types.*;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Feb 11, 2010
 * Time: 7:25:19 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsFaultGrid
{
    StsSection section;
    public int nRows;
    public int nCols;
    public SectionCell[][] sectionCells;
    public boolean isFaulted;
    public float thickness;
    public float faultPerm;
    public ArrayList<StsNNC> nncList = new ArrayList<StsNNC>();

    static double tranMultiplier = StsEclipseOutput.tranMultiplier;

    static final int LEFT = StsParameters.LEFT;
    static final int RIGHT = StsParameters.RIGHT;

    static final boolean debugNoThrow = false;

    public StsFaultGrid(StsSection section)
    {
        this.section = section;
        isFaulted = section.isFault();
        nRows = section.getNRows();
        nCols = section.getNCols();
        if (!checkConstructLayerPolygons()) return;
        sectionCells = new SectionCell[nRows - 1][nCols - 1];
        addFaultCellPolygons();
    }

    private boolean checkConstructLayerPolygons()
    {
        if (!checkConstructLayerPolygons(section.getLeftZoneSides(), LEFT)) return false;
        if (!checkConstructLayerPolygons(section.getRightZoneSides(), RIGHT)) return false;
        return true;
    }

    private boolean checkConstructLayerPolygons(StsList zoneSides, int side)
    {
        for (int n = 0; n < zoneSides.getSize(); n++)
        {
            StsZoneSide zoneSide = (StsZoneSide) zoneSides.getElement(n);
            if(zoneSide == null) continue;
            if (!zoneSide.checkConstructLayerPolygons()) return false;
        }
        return true;
    }

    private Iterator getGridPolygonIterator(StsZoneSide zoneSide)
    {
        StsEdgeLoop edgeLoop = zoneSide.getEdgeLoop();
        return edgeLoop.getZoneSideGridPolygonIterator();
    }

    /**
     * The fault grid is a 2D grid of quad cells with horizontal edge rows and monotonically in z increasing rib columns.
     * The edgeloops for the zoneSides on each side of the fault section contain quadstrips with a set of quad cells from
     * this section, and trimmed polygons on the edges of the loop.  The edgeloop quad-cells and polygons have a row and column
     * number corresponding to the quad cell indexes of the section (type SECTION).  These quad-cells and polygons are further subdivided into
     * grid polygons each with the grid row and column (type grid). These grid polygons are further divided layer subpolygons.
     * For each grid polygon, there are two sets of layer polygons for the left and right sides of the section (types LAYER_LEFT and
     * LAYER_RIGHT).  For each of these layer subpolygons, the row and column number correspond to the grid cell it belongs to and the layer
     * number is the global layer it belongs to.  These are used in simulation grid construction, where we wish to compute NNCs
     * (non-neighbor connections) from the overlap areas of layer subpolygons on each side of a grid polygon.
     */
    private void addFaultCellPolygons()
    {
        addFaultCellPolygons(section.getLeftZoneSides(), LEFT);
        addFaultCellPolygons(section.getRightZoneSides(), RIGHT);
    }

    private void addFaultCellPolygons(StsList zoneSides, int side)
    {
        for (int n = 0; n < zoneSides.getSize(); n++)
        {
            StsZoneSide zoneSide = (StsZoneSide) zoneSides.getElement(n);
            if(zoneSide == null || zoneSide.getSection() == null) continue;
            byte geometryType = zoneSide.getGeometryType();
            if (geometryType == StsSection.GEOM_UNALIGNED)
                addFaultCellPolygons(zoneSide, side);
        }
    }

    private void addFaultCellPolygons(StsZoneSide zoneSide, int side)
    {
        StsEdgeLoop edgeLoop = zoneSide.getEdgeLoop();
        StsList quadStrips = edgeLoop.getQuadStrips();
        int nQuadStrips = quadStrips.getSize();
        int blockNumber = zoneSide.getBlock().getIndex();
        for (int n = 0; n < nQuadStrips; n++)
        {
            QuadStrip quadStrip = (QuadStrip) quadStrips.getElement(n);
            int row = quadStrip.rowNumber;
            int firstCol = quadStrip.firstCol;
            int lastCol = quadStrip.lastCol;
            for (int col = firstCol; col < lastCol; col++)
            {
                SectionCell sectionCell = getSectionCell(row, col);
                Object[] quadGridPolygons = quadStrip.polygons[col - firstCol];
                sectionCell.addGridPolygons(quadGridPolygons, side, blockNumber);
            }
        }
        StsList edgePolygons = edgeLoop.getPolygons();
        int nEdgePolygons = edgePolygons.getSize();
        for (int n = 0; n < nEdgePolygons; n++)
        {
            StsPolygon edgePolygon = (StsPolygon) edgePolygons.getElement(n);
            int row = edgePolygon.row;
            int col = edgePolygon.col;
            SectionCell sectionCell = getSectionCell(row, col);
            StsList gridPolygonsList = edgePolygon.getSubPolygons();
            sectionCell.addGridPolygons(gridPolygonsList, side, blockNumber);
        }
    }

    private SectionCell getSectionCell(int row, int col)
    {
        if (row < 0 || row >= nRows - 1 || col < 0 || col >= nCols - 1)
        {
            StsException.systemError(this, "getSectionCell", "Requested grid cell row " + row + " or col " + col + "out of range");
            return null;
        }
        SectionCell sectionCell = sectionCells[row][col];
        if (sectionCell != null) return sectionCell;
        sectionCell = new SectionCell();
        sectionCells[row][col] = sectionCell;
        return sectionCell;
    }

    public boolean constructNNCs(StsBlock[] blocks, StsEclipseOutput eclipseOutput, StsStatusUI status, float dProgress)
    {
        int row = -1, col = -1;
        dProgress /= nRows;
        try
        {
            for (row = 0; row < nRows - 1; row++)
            {
                eclipseOutput.progress += dProgress;
                status.setProgress(eclipseOutput.progress);
                for (col = 0; col < nCols - 1; col++)
                {
                    SectionCell sectionCell = sectionCells[row][col];
                    if (sectionCell == null) continue;
                    ArrayList<FaultGridCell> faultGridCells = sectionCell.faultGridCells;
                    for (FaultGridCell faultGridCell : faultGridCells)
                        if (!faultGridCell.constructNNCs(blocks)) return false;
                }
            }
            Collections.sort(nncList);
            // checkCombineNNCs();
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "outputNNCs", "row: " + row + " col: " + col, e);
            return false;
        }
    }

    public void checkCombineNNCs()
    {
        Iterator<StsNNC> iterator = nncList.iterator();
        StsNNC combineNNC = null;
        while (iterator.hasNext())
        {
            StsNNC nnc = iterator.next();
            // if NNCs have the same connection, sum the transmissiblity
            // and remove this redundant connection
            if (combineNNC != null)
            {
                if (nnc.sameNNC(combineNNC))
                {
                    combineNNC.trans += nnc.trans;
                    combineNNC.area += nnc.area;
                    iterator.remove();
                }
                else
                    combineNNC = nnc;
            }
            else
                combineNNC = nnc;
        }
    }
/*
    public void checkMergeCells(float minTimeConstant)
    {
        for(StsNNC nnc : nncList)
        {
            checkMergeCell(nnc.leftGridCell, minTimeConstant);
            checkMergeCell(nnc.rightGridCell, minTimeConstant);
        }
        for(StsNNC nnc : nncList)
        {
            nnc.leftGridCell.merge(nnc);
            nnc.rightGridCell.merge(nnc);
        }
    }

    private void checkMergeCell(StsBlock.BlockCellColumn.GridCell gridCell, float minTimeConstant)
    {
        if(gridCell.timeConstant >= minTimeConstant) return;
        gridCell.checkSetParentCell();
    }
*/
    class TransmissibilityComparator implements Comparator<StsNNC>
    {
        public TransmissibilityComparator(){ }

        public int compare(StsNNC nnc1, StsNNC nnc2)
        {
            if (nnc1.trans < nnc2.trans) return -1;
            else if (nnc1.trans > nnc2.trans) return 1;
            else return 0;
        }
    }

    class SectionCell
    {
        public ArrayList<FaultGridCell> faultGridCells;

        SectionCell()
        {
            faultGridCells = new ArrayList<FaultGridCell>();
        }

        void addGridCell(FaultGridCell faultGridCell)
        {

        }

        void addGridPolygons(Object[] gridPolygons, int side, int blockNumber)
        {
            for (int n = 0; n < gridPolygons.length; n++)
            {
                StsPolygon gridPolygon = (StsPolygon) gridPolygons[n];
                int row = gridPolygon.row;
                int col = gridPolygon.col;
                FaultGridCell faultGridCell = getGridCell(row, col);
                faultGridCell.addGridPolygon(gridPolygon, side, blockNumber);
            }
        }

        void addGridPolygons(StsList gridPolygonsList, int side, int blockNumber)
        {
            addGridPolygons(gridPolygonsList.getTrimmedList(), side, blockNumber);
        }

        private FaultGridCell getGridCell(int row, int col)
        {
            for (FaultGridCell faultGridCell : faultGridCells)
            {
                if (faultGridCell.row == row && faultGridCell.col == col)
                    return faultGridCell;
            }
            FaultGridCell faultGridCell = new FaultGridCell(row, col);
            faultGridCells.add(faultGridCell);
            return faultGridCell;
        }
    }

    class FaultGridCell
    {
        int row, col;

        ArrayList<StsPolygon> leftPolygons;
        ArrayList<StsPolygon> rightPolygons;

        static final int LEFT = StsParameters.LEFT;
        static final int RIGHT = StsParameters.RIGHT;
        static final byte LAYER_RIGHT = StsPolygon.LAYER_RIGHT;
        static final byte LAYER_LEFT = StsPolygon.LAYER_LEFT;

        FaultGridCell(int row, int col)
        {
            this.row = row;
            this.col = col;
            leftPolygons = new ArrayList<StsPolygon>();
            rightPolygons = new ArrayList<StsPolygon>();
        }

        void addGridPolygon(StsPolygon gridPolygon, int side, int blockNumber)
        {
            byte type = gridPolygon.type;
            if (type == StsPolygon.GRID)
            {
                StsList subPolygonsList = gridPolygon.getSubPolygons();
                addLayerPolygons(subPolygonsList, side);
                return;
            }
            else if (side == LEFT && type == LAYER_LEFT)
                addLayerPolygon(gridPolygon, LEFT, blockNumber);
            else if (side == RIGHT && type == LAYER_RIGHT)
                addLayerPolygon(gridPolygon, RIGHT, blockNumber);
            else
                StsException.systemError(this, "addGridPolygon", "attempted to add to layer polygon list polygon " + gridPolygon);
        }

        void addLayerPolygons(StsList layerPolygonsList, int side)
        {
            Object[] layerPolygons = layerPolygonsList.getList();
            if (side == LEFT)
                addLayerPolygons(leftPolygons, layerPolygons);
            else
                addLayerPolygons(rightPolygons, layerPolygons);
        }

        void addLayerPolygons(ArrayList<StsPolygon> polygonList, Object[] polygons)
        {
            for (int n = 0; n < polygons.length; n++)
                polygonList.add((StsPolygon) polygons[n]);
        }

        void addLayerPolygon(StsPolygon layerPolygon, int side, int blockNumber)
        {
            layerPolygon.nBlock = blockNumber;
            if (side == LEFT)
                leftPolygons.add(layerPolygon);
            else
                rightPolygons.add(layerPolygon);
        }

        public boolean constructNNCs(StsBlock[] blocks)
        {
            for (StsPolygon rightPolygon : rightPolygons)
            {
                double[] normal = rightPolygon.getNormal();
                StsMath.normalize(normal);
                for (StsPolygon leftPolygon : leftPolygons)
                {
                    double area = StsPolygon.computeIntersectionAreaNew(rightPolygon, leftPolygon, normal);
                    if (area <= 0.0) continue;

                    int[] ijkbRight = getCellIJKB(rightPolygon);
                    int[] ijkbLeft = getCellIJKB(leftPolygon);
                    StsBlock.BlockCellColumn.GridCell rightCell = getGridCell(ijkbRight, blocks);
                    if (rightCell == null)
                    {
                        StsException.systemError(this, "computeNNCs", "Failed to find cell for " + StsNNC.getStringIJKB(ijkbRight) + " area: " + area);
                        continue;
                    }
                    if(rightCell.poreVolume <= 0.0f)
                        StsException.systemDebug(this, "computeNNCs", "NNC rightCell has zero pore volume " + rightCell.toString());

                    StsBlock.BlockCellColumn.GridCell leftCell = getGridCell(ijkbLeft, blocks);
                    if (leftCell == null)
                    {
                        StsException.systemError(this, "computeNNCs", "Failed to find cell for " + StsNNC.getStringIJKB(ijkbLeft) + " area: " + area);
                        continue;
                    }
                    if(leftCell.poreVolume <= 0.0f)
                        StsException.systemDebug(this, "computeNNCs", "NNC leftCell has zero pore volume " + leftCell.toString());

                    StsNNC nnc = StsNNC.constructor(rightPolygon.computeCenter(), normal, rightCell, leftCell, area);
                    if(debugNoThrow && !StsNNC.ijkCompare(ijkbLeft, ijkbRight))
                    {
                        StsException.systemDebug(this, "computeNNCs", "ijkLeft != ijkRight: " + StsNNC.getStringIJKB(ijkbLeft) + " " + StsNNC.getStringIJKB(ijkbRight) + " area: " + area + " nnc.area: " + nnc.area);
                    }
                    if (nnc != null)
                    {
                        if (Double.isNaN(nnc.trans) || Double.isInfinite(nnc.trans))
                            StsException.systemError(this, "computeNNCs", "Transmissiblility is bad " + nnc.trans);
                        else
                            nncList.add(nnc);
                    }
                }
            }
            return true;
        }

        int[] getCellIJKB(StsPolygon polygon)
        {
            return polygon.getIJKB();
        }

        StsBlock.BlockCellColumn.GridCell getGridCell(int[] ijkb, StsBlock[] blocks)
        {
            int blockIndex = ijkb[3];
            StsBlock block = getBlock(blockIndex, blocks);
            if (block == null) return null;
            return block.getGridCell(ijkb);
        }

        StsBlock getBlock(int blockIndex, StsBlock[] blocks)
        {
            for (int n = 0; n < blocks.length; n++)
                if (blocks[n].getIndex() == blockIndex)
                    return blocks[n];
            StsException.systemError(this, "getBlock", "Failed to find block with index: " + blockIndex);
            return null;
        }
    }
}
