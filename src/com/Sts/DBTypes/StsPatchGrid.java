package com.Sts.DBTypes;

import com.Sts.Actions.Wizards.SurfaceCurvature.*;
import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.Utilities.Seismic.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.*;
import java.util.*;

/**
 * A patchGrid has a XY boundingBox containing a row-col array of transient patchPoints. It is constructed by 2D and 3D auto-pickers.
 * For the volume auto-picker, there may be a group of patchGrids which are connected together and overlapping.  The patchGrid with
 * the smallest ID serves as the parent and the others are children in a linked list with the first being the childGrid.
 */
public class StsPatchGrid extends StsXYGridBoundingBox implements Comparable<StsPatchGrid>, StsSerializable, StsXYSurfaceLinkGridable
{
	private static final long serialVersionUID = 1L;
	/** this is the ID of this  patchGrid; it will be reset to a rowSorted sequence index when patch is completed */
	public int id;
	/** parent grid of this grid; null if this is the root grid */
	public StsPatchGrid parentGrid = null;
	/** childGrid is the next grid in the linked-list. */
	public StsPatchGrid childGrid = null;
	/** the original id of this patch in gridList; saved when patch is re-indexed and is used for debugging purposes */
	int originalID = -1;
	/** type of this grid: Min, Max +Zero-crossing, or -Zero-crossing */
	byte patchType;
	boolean isVisible = true;
	StsPatchVolume patchVolume;
	int nPatchPoints;
	float dataMin;
	float dataMax;
	float zMin = StsParameters.largeFloat;
	float zMax = -StsParameters.largeFloat;

	float[][] pointsZ;
	float[][] rowCorrels;
	float[][] colCorrels;

	transient PatchPoint[][] patchPoints = null;
	/** flag indicating this patchGrid has been added to current patchVolume.rowGrids hashmap list; avoids overhead of trying to re-add */
	transient public boolean rowGridAdded = false;

	transient RowColGrid rowColGrid = new RowColGrid(rowMin, rowMax, colMin, colMax);

	transient float[][] values;

	transient int[][] colorIndices;

	transient StsDiamondStrips diamondStrips;

	transient public int nValuePatchPoints = 0;
	transient public double sum = 0;
	/** points on this grid are in a hashMap with a hash key whose value is the grid index: col + row*patchVolume.nCols */
	// transient HashMap<Integer, StsPatchPoint> patchPointsHashMap;
	/** max size of all patches generated; printed out for diagnostics */
	static int maxGridSize = 0;
	/**
	 * index to be assigned to next patchGrid created; this is the index during construction for a multilayered patch;
	 * incremented for each one; reset on initialization
	 */
	static int nextPatchID = 0;
	/** final index in gridList reset on initialization */
	static int nextFinalPatchID = 0;

	static final float nullValue = StsParameters.nullValue;

	// Multiply # pts in SVD to get ChiSqr limit
	static final private double chiSqrMultiplyer = 2;
	//static final private double stdDevFactor = 1.5;
	static final byte FILTER_NONE = 0;
	static final byte FILTER_ON_CHI_SQ = 1;
	static final byte FILTER_ON_STD_DEV = 2;

	static public byte filterType = FILTER_ON_CHI_SQ;

	static public final float badCurvature = StsQuadraticCurvature.badCurvature;
	static public final float curvatureTest = StsQuadraticCurvature.curvatureTest;

	static boolean sortRowFirst = true;

	static final int largeInt = Integer.MAX_VALUE;

	/*---------------------------- SYSTEM DEBUG FLAGS (DON'T EDIT) ---------------------------------------------------*/
	/** a static final: if false, all blocks bounded by this boolean will not be compiled or checked */
	static public final boolean debug = StsPatchVolume.debug;
	/** standard no-debug flag */
	static public final int NO_DEBUG = -1;
	/** patchPoint patch debug. This will be set to true if debug==true and debugPatchInitialID is set. */
	static public boolean debugPatchGrid;
	/** debugPoint is true if we have debugRow && debugCol set and either debugPatchGrid is set or debugSlice is set */
	static boolean debugPoint;
	/*---------------------------- USER DEBUG FLAGS (SET IN PatchPickPanel GUI (DON't EDIT) --------------------------*/
	/** various debugPatchGrid print of patches in rowGrid and prevRowGrid arrays; check if this is not -1 and prints if id matches this value.  Make this -1 if you want no debugPatchGrid prints. */
	static public int debugPatchInitialID = NO_DEBUG;
	/** debugPatchID may change if original patch is merged into a new one; when this occurs, set debugCurrentPatchID to new one and track it */
	static int debugPatchID = NO_DEBUG;
	/** patchPoint row debug; used when row & ool and either slice or patchID are set. Set if you want a specific row & col at either a slice or patch debugged. */
	static int debugPointRow = NO_DEBUG;
	/** patchPoint col debug; used when row & ool and either slice or patchID are set. Set if you want a specific row & col at either a slice or patch debugged. */
	static int debugPointCol = NO_DEBUG;
	/** patchPoint slice debug. Set if you want a specific row & col & slice debugged. */
	static int debugPointSlice = NO_DEBUG;
	/*----------------------------------------------------------------------------------------------------------------*/


	public StsPatchGrid()
	{
	}

	private StsPatchGrid(StsPatchVolume patchVolume, byte patchType)
	{
		this.patchVolume = patchVolume;
		this.patchType = patchType;
	}

	static public StsPatchGrid construct(StsPatchVolume patchVolume, byte patchType)
	{
		StsPatchGrid patchGrid = new StsPatchGrid(patchVolume, patchType);
		patchGrid.setup();
		return patchGrid;
	}

	static public void staticInitialize()
	{
		nextPatchID = 0;
		nextFinalPatchID = 0;
		debugPatchID = debugPatchInitialID;
	}

	private void setup()
	{
		id = nextPatchID++;
		if (debugPatchID != NO_DEBUG && id == debugPatchID)
			StsException.systemDebug(this, "setup", "patch " + id + " initialized");
	}

	static public void initializeDebug(StsPatchPickPanel pickPanel)
	{
		debugPatchInitialID = pickPanel.patchId;
		debugPatchID = debugPatchInitialID;
		debugPatchGrid = debug && debugPatchInitialID != NO_DEBUG;
		debugPointRow = pickPanel.pointRow;
		debugPointCol = pickPanel.pointCol;
		debugPointSlice = pickPanel.pointSlice;
		debugPoint = debug && debugPointRow != NO_DEBUG && debugPointCol != NO_DEBUG && (debugPatchGrid || debugPointSlice != NO_DEBUG);
	}

	public boolean debug()
	{
		return debugPatchID != -1 && id == debugPatchID;
	}

	/**
	 * give two different size grids with this intersection, check if we can merge the two;
	 * possible only if they don't have common occupancy of any row-col location
	 * @param patchGrid1 first of patchGrids to merge
	 * @param patchGrid2 second of patchGrids to merge
	 */
	static public boolean mergePatchPointsOK(StsPatchGrid patchGrid1, StsPatchGrid patchGrid2)
	{
		// check for any overlap between this grid and patchPointGrid
		RowColGridFilled intersectGrid = patchGrid1.rowColGridIntersect(patchGrid2);
		int r1 = intersectGrid.rowMin - patchGrid1.rowMin;
		int r2 = intersectGrid.rowMin - patchGrid2.rowMin;
		int c1 = intersectGrid.colMin - patchGrid1.colMin;
		int c2 = intersectGrid.colMin - patchGrid2.colMin;
		for (int row = 0; row < intersectGrid.nRows; row++, r1++, r2++)
		{
			int cc1 = c1;
			int cc2 = c2;
			for (int col = 0; col < intersectGrid.nCols; col++, cc1++, cc2++)
			{
				if (patchGrid1.patchPoints[r1][cc1] != null && patchGrid2.patchPoints[r2][cc2] != null)
					return false;
			}
		}
		return true;
	}

	/**
	 * We wish to merge the points from removedGrid into this one.  Copy both sets of points to a new grid which is union of two.
	 * reset the removedGrid patchPoints.id to this id
	 * @param removedGrid newPatchGrid to be merged to this (otherPatchGrid).
	 * @return true if merged successfully
	 */
	boolean mergePatchPoints(StsPatchGrid removedGrid)
	{
		RowColGrid union = rowColGridUnion(removedGrid); //union of this and newPatchGrid
		PatchPoint[][] newMergedPatchPoints = new PatchPoint[union.nRows][union.nCols]; // create empty mergedPoints grid
		if (!copyPatchPointsTo(newMergedPatchPoints, union)) return false; // copy this patchPoints to mergedPoints
		if (!removedGrid.copyPatchPointsTo(newMergedPatchPoints, union)) return false;
		removedGrid.resetPatchPointsGrid(this);
		resetPatchPoints(union, newMergedPatchPoints);
		nPatchPoints += removedGrid.nPatchPoints;
		return true;
	}

	/**
	 * Create new patchPoint array for partially merged patches which is the union of the two;
	 * as the newConnection is between existing points on the respective patches, it will be included.
	 * @param fromGrid points are to be moved from the fromGrid to "this" toGrid
	 * @param newConnection the new connection between the patches
	 * @param orderReversed false if connection is from the toGrid to the prevPoint on the fromGrid; otherwise reversed and true
	 */
	protected StsPatchGrid moveNonOverlappingPointsFrom(StsPatchGrid fromGrid, Connection newConnection, boolean orderReversed)
	{
		try
		{
			RowColGrid unionGrid = rowColGridUnion(fromGrid); //union of this and otherGrid
			PatchPoint[][] newMergedPatchPoints = new PatchPoint[unionGrid.nRows][unionGrid.nCols]; // create empty mergedPoints grid
			if (!copyPatchPointsTo(newMergedPatchPoints, unionGrid))
				return null; // copy this patchPoints to mergedPoints
			return fromGrid.moveNonOverlappingPointsTo(this, newMergedPatchPoints, unionGrid, newConnection, orderReversed);
			// resetPatchPointsGrid(this);
			//fromGrid.adjustGridSizeDown();
			//resetPatchPoints(unionGrid, newMergedPatchPoints);
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "moveNonOverlappingPointsFrom", "Failed moving points from patch " + fromGrid.toString() +
					" to this grid " + toString(), e);
			return null;
		}
	}

	static public void debugCheckOverlappedGrids(StsPatchGrid fromGrid, StsPatchGrid toGrid)
	{
		String directionString = "FROM grid: " + fromGrid.id + " TO grid: " + toGrid.id;
		fromGrid.debugCheckOverlappedGrid(directionString);
		toGrid.debugCheckOverlappedGrid(directionString);
	}

	static public void debugCheckOverlappedGridPoint(PatchPoint fromPoint, PatchPoint toPoint, byte initialState, byte state)
	{
		String stateString = stateStrings[state];
		if (initialState != state)
			stateString = stateStrings[initialState] + "->" + stateString;

		String fromPointGridID = (fromPoint != null) ? Integer.toString(fromPoint.patchGrid.id) : "NULL";
		String toPointGridID = (toPoint != null) ? Integer.toString(toPoint.patchGrid.id) : "NULL";

		String directionString = "FROM grid: " + fromPointGridID + " TO grid: " + toPointGridID + " state: " + stateString;
		if (fromPoint != null)
			fromPoint.patchGrid.debugCheckOverlappedGridPoint(directionString, fromPoint);
		if (toPoint != null)
			toPoint.patchGrid.debugCheckOverlappedGridPoint(directionString, toPoint);
	}

	private void debugCheckOverlappedGrid(String directionString)
	{
		for (int row = 0; row < nRows; row++)
			for (int col = 0; col < nCols; col++)
			{
				PatchPoint patchPoint = getGridPatchPoint(row, col);
				if (patchPoint == null) continue;
				debugCheckOverlappedGridPoint(directionString, patchPoint);
			}
	}

	private void debugCheckOverlappedGridPoint(String directionString, PatchPoint patchPoint)
	{
		if (patchPoint == null) return;
		int row = getGridRow(patchPoint);
		int col = getGridCol(patchPoint);
		Connection connection = patchPoint.getRowConnection();
		if (connection != null)
		{
			PatchPoint prevRowPoint = getGridPatchPoint(row, col - 1);
			if (prevRowPoint == null)
				StsException.systemError(this, "debugCheckOverlappingGrid", directionString + " grid: " + id +
						" point row: " + (row + rowMin) + " col: " + (col + colMin) + " HAS ROW CONNECTION BUT NO CONNECTED POINT. rowConnection: " + connection.toString());
		}
		connection = patchPoint.getColConnection();
		if (connection != null)
		{
			PatchPoint prevColPoint = getGridPatchPoint(row - 1, col);
			if (prevColPoint == null)
				StsException.systemError(this, "debugCheckOverlappingGrid", directionString + " grid: " + id +
						" point row: " + row + " col: " + col + " HAS COL CONNECTION BUT NO CONNECTED POINT. colConnection: " + connection.toString());
		}
	}

	/**
	 * this patchGrid is being merged to another patchGrid with this ID. Reset all the merged patchPoints to this id.
	 * @param newPatchGrid patch which this patch is being merged into
	 */
	private void resetPatchPointsGrid(StsPatchGrid newPatchGrid)
	{
		for (int row = 0; row < nRows; row++)
			for (int col = 0; col < nCols; col++)
				if (patchPoints[row][col] != null) patchPoints[row][col].setPatchGrid(newPatchGrid);
	}

	RowColGridFilled rowColGridIntersect(StsPatchGrid otherGrid)
	{
		int rowMin, rowMax, colMin, colMax;
		rowMin = Math.max(this.rowMin, otherGrid.rowMin);
		rowMax = Math.min(this.rowMax, otherGrid.rowMax);
		colMin = Math.max(this.colMin, otherGrid.colMin);
		colMax = Math.min(this.colMax, otherGrid.colMax);
		return new RowColGridFilled(rowMin, rowMax, colMin, colMax);
	}

	RowColGrid rowColGridUnion(StsPatchGrid otherGrid)
	{
		int rowMin, rowMax, colMin, colMax;
		rowMin = Math.min(this.rowMin, otherGrid.rowMin);
		rowMax = Math.max(this.rowMax, otherGrid.rowMax);
		colMin = Math.min(this.colMin, otherGrid.colMin);
		colMax = Math.max(this.colMax, otherGrid.colMax);
		return new RowColGrid(rowMin, rowMax, colMin, colMax);
	}


	static final int maxChildCount = 10000;

	/**
	 * Combine link-list from two grids together into a new link-list attached to one of the two grids.
	 * If a grid has a non-null parentGrid, then it belongs to the childGrid list for that parentGrid; it is a childGrid.
	 * if a grid nas a null parent and a non-null child, then it has a list of childGrids; it is a parentGrid.
	 * If a grid has null parentGrid and null childGrid, then it is brand new.
	 * So: if one grid is a parent and the other a child, add child to this parent.
	 * if both grids are parents, merge the two lists together with the one with the lowest ID as the new parent.
	 * if one grid is a child and the other new, then add the new one to the child.parent list.
	 * if both grids are children, use the parent of the two with the lowest ID as the parent for the other.
	 * Given the first grid as thisGrid and the second as newChildGrid, we have nine possible combinations:
	 * parent-parent, parent-child, parent-new, child-parent, new-parent, child-child, child-new, new-child, and new-new
	 */
	public void combineChildGrids(StsPatchGrid newChildGrid)
	{
		if (debugPatchID != NO_DEBUG && (newChildGrid.id == debugPatchID || id == debugPatchID))
			StsException.systemDebug(this, "combineChildGrids", "DEBUG COMBINE GRID: " + newChildGrid.toFamilyString() + " AND GRID: " + toFamilyString());

		// check for parent-parent, parent-child, parent-new, child-parent, new-parent, child-child, child-new, new-child, and new-new combinations
		if (isParent()) // this grid is a parent
		{
			if (newChildGrid.isParent()) // parent-parent: assume children are mutually exclusive
			{
				checkAddGridsID(newChildGrid);
			}
			else if (newChildGrid.isChild()) // parent-child: use parent-child.parent
			{
				checkAddGridsID(newChildGrid.parentGrid);
			}
			else // parent-new
			{
				addChildGrids(newChildGrid);
			}
		}
		else if (newChildGrid.isParent())
		{
			if (isChild())  // child-parent
			{
				parentGrid.checkAddGridsID(newChildGrid);
			}
			else // new-parent
			{
				newChildGrid.addChildGrids(this);
			}
		}
		else // neither grid is a parent, so find a child with the best parent
		{
			if (isChild())
			{
				if (newChildGrid.isChild()) // child-child, use parent with lowest ID
				{
					newChildGrid.parentGrid.checkAddGridsID(parentGrid);
				}
				else // child-new
				{
					parentGrid.addChildGrids(newChildGrid);
				}
			}
			else if (newChildGrid.isChild())
			{
				newChildGrid.parentGrid.addChildGrids(this); // new-child

			}
			else // new-new; make parent of one with lowest ID
			{
				checkAddGridsID(newChildGrid);
			}
		}
	}

	public final boolean isParent()
	{
		return parentGrid == null && childGrid != null;
	}

	public final boolean isChild()
	{
		return parentGrid != null;
	}

	public final boolean isNew()
	{
		return parentGrid == null && childGrid == null;
	}

	static private String getGridsRelation(StsPatchGrid grid, StsPatchGrid newGrid)
	{
		return grid.getFamilyTypeString() + "-" + newGrid.getFamilyTypeString();
	}

	public String getFamilyTypeString()
	{
		if (isParent()) return " PARENT ";
		else if (isChild()) return " CHILD ";
		else return " NEW ";
	}

	public String getPatchTypeString()
	{
		return " " + StsTraceUtilities.typeStrings[patchType] + " ";
	}

	public String getGridDescription()
	{
		StsPatchGrid parentGrid;
		String gridFamilyType = getFamilyTypeString();
		if (isNew())
			return " gridType: New grids: 1 nPoints: " + nPatchPoints;

		if (this.parentGrid == null)
			parentGrid = this;
		else
			parentGrid = this.parentGrid;

		int nGrids = 1;
		int nPatchPoints = parentGrid.nPatchPoints;
		StsPatchGrid childGrid = parentGrid.childGrid;
		while (childGrid != null)
		{
			nGrids++;
			nPatchPoints += childGrid.nPatchPoints;
			childGrid = childGrid.childGrid;
		}
		return " gridType: " + gridFamilyType + " grids: " + nGrids + " nPoints: " + nPatchPoints;
	}

	private final void checkAddGridsID(StsPatchGrid newChildGrid)
	{
		if (this == newChildGrid) return;
		if (newChildGrid.id < id)
			newChildGrid.addChildGrids(this);
		else
			addChildGrids(newChildGrid);
	}

	/**
	 * Add this newChildGrid and its link-list to this grid link-list and remove the newChildGrid link-list and set its
	 * parentGrid to this grid.
	 * If a link-list exists, then the grid it's associated with must be a parentGrid, otherwise it is a childGrid.
	 * A grid can exist in only one link-list unless it is a rootGrid (no parent, but has children).
	 * @param newChildGrid grid and its link-list to be added to this grid
	 */
	private void addChildGrids(StsPatchGrid newChildGrid)
	{
		if (debugPatchID != NO_DEBUG && (newChildGrid.id == debugPatchID || id == debugPatchID))
			StsException.systemDebug(this, "addChildGrid", "DEBUG CHILD GRID: " + getGridsRelation(this, newChildGrid) + toFamilyString() + " AND " + newChildGrid.toFamilyString());

		// check if parentGrid already contains this childGrid; if not, add it
		// if parent already has child, method returns false; if false, return with no further action
		if (!checkSetChildGrid(this, newChildGrid)) return;
		// set the parentGrid for these added childGrids to this grid
		resetParentForChildren(newChildGrid, this);
	}

	static private boolean checkSetChildGrid(StsPatchGrid parentGrid, StsPatchGrid newChildGrid)
	{
		if (parentGrid == newChildGrid)
			return true;
		int[] gridIDs;
		int count = 0;
		if (debug)
		{
			gridIDs = new int[nextPatchID];
			Arrays.fill(gridIDs, -1);
			gridIDs[parentGrid.id] = count++;
		}
		// check if this newChildGrid is not already a child of this grid
		StsPatchGrid childGrid = parentGrid.childGrid;
		while (childGrid != null)
		{
			if (childGrid == newChildGrid)
			{
				return false; // already in list
			}
			if (debug)
			{
				if (gridIDs[childGrid.id] != -1)
					StsException.systemDebug(StsPatchGrid.class, "addChildGrids", "CHILD GRID REPEATED " + " id: " + childGrid.id);
				else
					gridIDs[childGrid.id] = count;

				count++;
				if (count > maxChildCount)
				{
					StsException.systemDebug(StsPatchGrid.class, "addChildGrids", "MAX CHILD GRID COUNT " + maxChildCount + " exceeded.");
					return false;
				}
			}
			parentGrid = childGrid;
			childGrid = parentGrid.childGrid;
		}
		// we have found the last non-null child in list (parentGrid); add the newChildGrid to it
		parentGrid.childGrid = newChildGrid;
		return true;
	}

	static public void resetParentForChildren(StsPatchGrid firstChildGrid, StsPatchGrid newParentGrid)
	{
		int count = 0;
		int[] gridIDs;
		if (debug)
		{
			gridIDs = new int[nextPatchID];
			Arrays.fill(gridIDs, -1);
			gridIDs[newParentGrid.id] = count++;
		}

		StsPatchGrid childGrid = firstChildGrid;
		while (childGrid != null)
		{
			if (debugPatchID != NO_DEBUG && childGrid.id == debugPatchID)
				StsException.systemDebug(StsPatchGrid.class, "resetParentForChildren", "GRID PARENT SET. GRID: " + childGrid.toFamilyString() + " NEW PARENT " + newParentGrid.id);

			childGrid.parentGrid = newParentGrid;
			if (debug)
			{
				if (gridIDs[childGrid.id] != -1)
					StsException.systemDebug(StsPatchGrid.class, "addChildGrids", "RESET ID CHILD GRID REPEATED " + " id: " + childGrid.id);
				else
					gridIDs[childGrid.id] = count;

				count++;
				if (count > maxChildCount)
				{
					StsException.systemDebug(StsPatchGrid.class, "addChildGrids", "RESET ID MAX CHILD GRID COUNT " + maxChildCount + " exceeded.");
					return;
				}
			}
			childGrid = childGrid.childGrid;
		}
	}

	/**
	 * This grid is the removedGrid which is being merged with mergedGrid.
	 * removedGrid could be parent, child, or new;
	 * if removedGrid is a parent, then the first child needs to be made the parent and all its children need to have parent reset to first child
	 * if removed grid is a child, it needs to be removed from the parent (deleted from link-list)
	 * if new, we don't have to do anything
	 * @param mergedGrid included for debugging only
	 */
	public void removeChildOrParent(StsPatchGrid mergedGrid)
	{
		if (debugPatchID != NO_DEBUG && (id == debugPatchID || mergedGrid.id == debugPatchID))
			StsException.systemDebug(this, "removeChildOrParent", "DEBUG GRID. REMOVED: " + toString() + " MERGED TO: " + mergedGrid.toString());

		if (isParent())
		{
			if (debugPatchID != NO_DEBUG && (id == debugPatchID || mergedGrid.id == debugPatchID))
				StsException.systemDebug(this, "removeChildOrParent", "GRID PARENT SET TO NULL. GRID: " + childGrid.toFamilyString());
			// parent is removed, and its child is new parent
			childGrid.parentGrid = null;
			StsPatchGrid.resetParentForChildren(childGrid.childGrid, childGrid);

		}
		if (isChild())
		{
			StsPatchGrid parentGrid = this.parentGrid;
			StsPatchGrid childGrid = parentGrid.childGrid;
			while (childGrid != null)
			{
				if (childGrid == this)
				{
					parentGrid.childGrid = childGrid.childGrid;
				}
				parentGrid = childGrid;
				childGrid = childGrid.childGrid;
			}
		}
	}

	private RowColGrid getRowColGrid()
	{
		return new RowColGrid(rowMin, rowMax, colMin, colMax);
	}

	class RowColGrid
	{
		int rowMin = largeInt;
		int rowMax = -largeInt;
		int colMin = largeInt;
		int colMax = -largeInt;
		int nRows = 0, nCols = 0;

		RowColGrid(int rowMin, int rowMax, int colMin, int colMax)
		{
			this.rowMin = rowMin;
			this.rowMax = rowMax;
			this.colMin = colMin;
			this.colMax = colMax;
			nRows = rowMax - rowMin + 1;
			nCols = colMax - colMin + 1;
		}

		public String toString()
		{
			return new String("rowMin: " + rowMin + " rowMax: " + rowMax + " colMin: " + colMin + " colMax: " + colMax);
		}
	}

	class RowColGridFilled extends RowColGrid
	{
		boolean[][] isFilled;
		int nFilled = 0;

		RowColGridFilled(int rowMin, int rowMax, int colMin, int colMax)
		{
			super(rowMin, rowMax, colMin, colMax);
			isFilled = new boolean[nRows][nCols];
		}

		void fill(int row, int col)
		{
			gridFill(row - rowMin, col - colMin);
		}

		boolean gridFill(int row, int col)
		{
			if (isFilled[row][col]) return false;
			isFilled[row][col] = true;
			nFilled++;
			return true;
		}
	}

	/**
	 * return true if patchPoint overlaps this grid; i.e., there already is a point on this grid at the same row and col
	 * @param patchPoint point at whose row and col we want to see if point already exists on this grid
	 * @return true if this point overlaps an existing point on this grid
	 */
	boolean patchPointOverlaps(PatchPoint patchPoint)
	{
		try
		{
			if (getVolPatchPoint(patchPoint) == null) return false;
			if (debug && debugPoint && (doDebugPoint(patchPoint)))
				StsException.systemDebug(this, "patchPointOverlaps", StsPatchVolume.iterLabel + "patchPoint " + patchPoint.toString() +
						" overlaps window " + getVolPatchPoint(patchPoint));
			return true;
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "patchPointOverlaps", "FAILED FOR WINDOW: " + patchPoint.toString() +
					" ON GRID: " + toGridString(), e);
			return false;
		}
	}

	boolean addPatchPoint(PatchPoint patchPoint)
	{
		// if (debugPatchGrid && id == debugPatchID)
		if (debug && debugPoint && (doDebugPoint(patchPoint)))
			StsException.systemDebug(this, "addPatchPoint", StsPatchVolume.iterLabel + "ADD POINT TO GRID: " + id + " POINT: " + patchPoint.toString());

		if (patchPoints == null)
			initializePatchPoints(patchPoint);
		else
			checkAdjustGrid(patchPoint);
		if (!contains(patchPoint))
		{
			StsException.systemError(this, "addPatchPoint", "pointGrid " + rowColGrid.toString() + " doesn't contain window " + patchPoint.toString());
			return false;
		}
		if (patchPoints[patchPoint.getRow() - rowMin][patchPoint.getCol() - colMin] != null)
		{
			StsException.systemError(this, "addPatchPoint", "pointGrid " + rowColGrid.toString() + " ALREADY HAS Point " + patchPoints[patchPoint.getRow() - rowMin][patchPoint.getCol() - colMin].toString() +
					" so can't add " + patchPoint.toString());
			return false;
		}
		patchPoints[patchPoint.getRow() - rowMin][patchPoint.getCol() - colMin] = patchPoint;
		patchPoint.setPatchGrid(this);
		nPatchPoints++;
		return true;
	}

	void setPoint(PatchPoint patchPoint)
	{
		patchPoints[getGridRow(patchPoint)][getGridCol(patchPoint)] = patchPoint;
	}

	final public int getGridRow(PatchPoint patchPoint)
	{
		return patchPoint.getRow() - rowMin;
	}

	final public int getGridCol(PatchPoint patchPoint)
	{
		return patchPoint.getCol() - colMin;
	}

	void checkAdjustGrid(PatchPoint patchPoint)
	{
		int row = patchPoint.getRow();
		int col = patchPoint.getCol();

		int rowMinNew = rowMin, rowMaxNew = rowMax, colMinNew = colMin, colMaxNew = colMax;

		boolean gridChanged = false;
		if (row < this.rowMin)
		{
			rowMinNew = row;
			gridChanged = true;
		}
		if (row > this.rowMax)
		{
			rowMaxNew = row;
			gridChanged = true;
		}
		if (col < this.colMin)
		{
			colMinNew = col;
			gridChanged = true;
		}
		if (col > this.colMax)
		{
			colMaxNew = col;
			gridChanged = true;
		}

		if (!gridChanged) return;

		RowColGrid newRowColGrid = new RowColGrid(rowMinNew, rowMaxNew, colMinNew, colMaxNew);
		copyResetRowColGrid(newRowColGrid);
	}

	void copyResetRowColGrid(RowColGrid newRowColGrid)
	{
		if (debug && debugPatchGrid && id == debugPatchID)
			StsException.systemDebug(this, "copyResetRowColGrid", StsPatchVolume.iterLabel + "grid reset from " + rowColGrid + " to " + newRowColGrid);
		PatchPoint[][] newPatchPoints = copyPatchPoints(newRowColGrid);
		resetPatchPoints(newRowColGrid, newPatchPoints);
	}

	/**
	 * points have been removed from this grid, so resize down.
	 * @return false if grid has zero rows or zer columns.
	 */
	boolean adjustGridSizeDown()
	{
		int newRowMin = rowMin;
		rowMinLoop:
		for (int row = 0; row < nRows; row++, newRowMin++)
			for (int col = 0; col < nCols; col++)
				if (patchPoints[row][col] != null)
					break rowMinLoop;

		int newRowMax = rowMax;
		rowMaxLoop:
		for (int row = nRows - 1; row >= 0; row--, newRowMax--)
			for (int col = 0; col < nCols; col++)
				if (patchPoints[row][col] != null)
					break rowMaxLoop;

		if (newRowMin > newRowMax) return false;

		int newColMin = colMin;
		colMinLoop:
		for (int col = 0; col < nCols; col++, newColMin++)
			for (int row = 0; row < nRows; row++)
				if (patchPoints[row][col] != null)
					break colMinLoop;

		int newColMax = colMax;
		colMaxLoop:
		for (int col = nCols - 1; col >= 0; col--, newColMax--)
			for (int row = 0; row < nRows; row++)
				if (patchPoints[row][col] != null)
					break colMaxLoop;

		if (newColMin > newColMax) return false;

		RowColGrid newRowColGrid = new RowColGrid(newRowMin, newRowMax, newColMin, newColMax);
		PatchPoint[][] newPatchPoints = new PatchPoint[newRowColGrid.nRows][newRowColGrid.nCols];
		copyPatchPointsTo(newPatchPoints, newRowColGrid);
		initializeRowColGrid(newRowColGrid);
		patchPoints = newPatchPoints;
		return true;
	}

	PatchPoint[][] copyPatchPoints(RowColGrid newRowColGrid)
	{
		if (patchPoints == null) return null;
		PatchPoint[][] newPatchPoints = new PatchPoint[newRowColGrid.nRows][newRowColGrid.nCols];
		if (!copyPatchPointsTo(newPatchPoints, newRowColGrid))
			return null;
		else
			return newPatchPoints;
	}

	boolean copyPatchPointsTo(PatchPoint[][] newPatchPoints, RowColGrid newRowColGrid)
	{
		int row = -1, newRow = -1;
		int col = -1, newCol = -1;
		try
		{
			int rowStart = rowMin - newRowColGrid.rowMin;
			int colStart = colMin - newRowColGrid.colMin;
			for (row = 0, newRow = rowStart; row < nRows; row++, newRow++)
				for (col = 0, newCol = colStart; col < nCols; col++, newCol++)
				{
					if (patchPoints[row][col] != null)
					{
						if (newPatchPoints[newRow][newCol] != null)
							StsException.systemError(this, "copyPatchPointsTo", "FAILED COPYING PATCH POINTS from grid: " + id + " " + rowColGrid.toString() +
									" row: " + row + " col: " + col + " to new grid " + newRowColGrid.toString() +
									"row: " + newRow + "col: " + newCol);
						else
							newPatchPoints[newRow][newCol] = patchPoints[row][col];
					}
				}
			return true;
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "copyPatchPointsTo", "Failed copying patch " + rowColGrid.toString() + " row: " + row +
					" to new grid " + newRowColGrid.toString() + "row: " + newRow, e);
			return false;
		}
	}

	/** Overlapping grid flag: point on fromGrid is NULL */
	static final byte GRID_NULL = 0;
	/** FOverlapping grid flag: point exists on toGrid so point can't be moved from fromGrid */
	static final byte GRID_UNCHANGED = 1;
	/** Overlapping grid flag: point doesn't exist on toGrid so point could be moved; connected point(s) on fromGrid can't
	 * so move fromPoint to toGrid and clone and leave on fromGrid with connections that can't be moved */
	static final byte GRID_CLONED = 2;
	/** Overlapping grid flag: point doesn't exist on toGrid and connected ponts are not UNCHANGED so point can be moved to toGrid */
	static final byte GRID_MOVED = 3;
	/** Overlapping grid flag: point on fromGrid is a clone of point on toGrid; if connected points are !UNCHANGED, then the connections
	 * can be moved and the point deleted (state will be set to ALREADY_CLONED_DELETE); otherwise state will be set to UNCHANGED */
	static final byte GRID_ALREADY_CLONED = 4;
	/** Overlapping grid flag: point on fromGrid is a clone of point on toGrid; connections can be moved so this point can be deleted.
	 * see GRID_ALREADY_CLONED */
	static final byte GRID_ALREADY_CLONED_DELETE = 5;
	/** connection to this point is new */
	static final byte GRID_NEW = 6;

	static String[] stateStrings = new String[]{"NULL", "UNCHANGED", "CLONED", "MOVED", "ALREADY_CLONED", "ALREADY_CLONED_DELETE", "NEW"};

	class PointState
	{
		byte state;
		byte rowPlusState = GRID_NULL;
		byte rowMinusState = GRID_NULL;
		byte colPlusState = GRID_NULL;
		byte colMinusState = GRID_NULL;

		PointState(byte state)
		{
			this.state = state;
		}

		void setRowPlusState(byte state) { rowPlusState = state; }
		void setRowMinusState(byte state) { rowMinusState = state; }
		void setColPlusState(byte state) { colPlusState = state; }
		void setColMinusState(byte state) { colMinusState = state; }

		// Connection states can be NULL, UNCHANGED, CLONED, MOVED, ALREADY_CLONED, or NEW
		// If state is UNCHANGED, then we can't move point so leave it alone.
		// If state is MOVED, we can possibly move it:
		//   1) if there is no UNCHANGED connection, then point can be moved
		//   2) otherwise if there is one or more UNCHANGED connections and one or more CLONED or NEW connections,
		//      we need to leave a clone at this location and move it;
		//      so change the state to CLONED;
		// If state is ALREADY_CLONED, this point was cloned from the point on the toGrid
		//   This ALREADY_CLONED point can be deleted if it has no UNCHANGED connection; otherwise this clone will be left as is

		//   1) if there is no UNCHANGED connection, then point can be moved
		//   2) otherwise if there is one or more UNCHANGED connections, we can move it, but we need to leave a clone at this location
		byte adjustState()
		{
			if(state == GRID_NULL || state == GRID_UNCHANGED)  return state;

			if(state == GRID_MOVED)
			{
				if(!hasConnectionState(GRID_UNCHANGED))  return state;
				state = GRID_CLONED;
			}
			else if(state == GRID_ALREADY_CLONED)
			{
				if(!hasConnectionState(GRID_UNCHANGED))
					state = GRID_ALREADY_CLONED_DELETE;
				if(!hasConnectionStateOtherThanNullAnd(GRID_UNCHANGED)) // if only UNCHANGED or NULL connections, can't move any connections
					state = GRID_UNCHANGED;
			}
			else
			{
				StsException.systemError(this, "getAction", "Unhandled state: " + stateStrings[state]);
			}
			return state;
		}

		/** move any connections to unmoved points back to clonedPoint */
		void moveConnections(PatchPoint point, PatchPoint clonedPoint)
		{
			Connection rowConnection = point.getRowConnection();
			if(rowConnection != null && rowMinusState == GRID_UNCHANGED)
			{
				clonedPoint.setRowConnection(rowConnection);
				point.deleteRowConnection();
			}
			Connection colConnection = point.getColConnection();
			if(colConnection != null && colMinusState == GRID_UNCHANGED)
			{
				clonedPoint.setColConnection(colConnection);
				point.deleteColConnection();
			}
		}

		boolean hasConnectionState(byte state)
		{
			if(rowMinusState == state) return true;
			if(colMinusState == state) return true;
			if(rowPlusState == state) return true;
			if(colPlusState == state) return true;
			return false;
		}

		boolean hasConnectionStateOtherThanNullAnd(byte state)
		{
			if(rowMinusState != state && rowMinusState != GRID_NULL) return true;
			if(rowPlusState != state && rowPlusState != GRID_NULL) return true;
			if(colMinusState != state && colMinusState != GRID_NULL) return true;
			if(colPlusState != state && colPlusState != GRID_NULL) return true;
			return false;
		}
		boolean hasMinusConnectionState(byte state)
		{
			if(rowMinusState == state) return true;
			if(colMinusState == state) return true;
			return false;
		}
	}
	/**
	 * We have a newConnection between two different grids and some points on the two grids overlap, so they can't be completely merged,
	 * but we will try to move as many points as possible from the smaller grid to the larger.
	 * We have previously determined that the two connection points don't mutually overlap.
	 * If one point on the newConnection can be moved to the other grid, then the newConnection can be added to that other grid.
	 * If the newConnection cannot be moved, then this partial merge operation will be aborted.
	 * Points are being moved from "from" grid to the "to" grid. The "from" grid has the fewest number of points.
	 * Points on the "from" grid may be deleted or cloned as explained below.  Connections may be moved or retained.
	 * Existing connections on the "from" grid are being used in this process.  We also have a newConnection between the two grids.
	 * Because of this newConnection, we wish to merge these two grids, but overlapping points exist so we will still have two grids.
	 * However, we are trying to move as many points as possible from the smaller grid (the "from" grid) to the larger grid (the "to" grid).
	 * A point on the "from" grid can have four states:
	 * NULL: doesn't exist on "from" grid;
	 * UNCHANGED: point exists at location on "to" grid;
	 * CLONED: point has been moved to "to" grid but connections on "from" grid which aren't moved require a clone there;
	 * MOVED: point is moved to "to" grid and connected points are either moved or cloned so connection(s) have gone with the pair.
	 * There are sixteen combinations between connected pairs, but reversed pairs are the same (U-C is same as C-U, though there might be
	 * an action on only the "C" in this case).  Thus there are 10 unique cases: NN,UU,CC,MM,NU,NC,NM,UC,UM,CM.
	 * Here are the actions for each of these cases:
	 * These four cases won't occur on a new connection since new connection is between the two grids.
	 * NN: both points are on "to" grid, so connection is ok as is. This case won't occur.
	 * UU: both points are on the "from" grid, so connection is ok as is.
	 * CC: both points are on the "to" grid, but clones remain on "from" grid so connection needs to be removed from both clones.
	 * MM: both points are on the "to" grid haveing been movd there, so connection is ok.
	 * These three cases will only occur for a new connection as they involve a connection between the two grids;
	 * NU,UN: try moving the N point on the "to" grid to the "from" grid (either move or clone); if this fails, a new grid is required for the two.
	 * NC,CN: N point is on the "to" grid; C point has been moved to "to" grid, but clone remains to handle other connection which wasn't moved.
	 * Add newConnection to point.
	 * NM,MN: M point has been moved to "to" grid; N point is on "to" grid, so connection is ok as is.
	 * These three cases won't occur on a new connection since new connection is between two grids.
	 * UC, CU: C point has been moved to "to" grid, remove connection from moved point and change connection to cloned point.
	 * UM,MU: Change M point state to C and apply UC logic: point remains on "to" grid, but clone is on "from" grid;
	 * remove connection from moved point and change connection to cloned point.
	 * CM,MC: C point has been moved to "to" grid and cloned on "from" grid; remove connection from C point.
	 * The steps in moving nonoverlapping points:
	 * 1) first set the state of each point in the grid to be moved: NULL, UNCHANGED, CLONED, or MOVED considering
	 * both previous connected points and next connected points (handle latter by reseting prevPoints if connection(s) exist)
	 * 2) add newConnection if possible
	 * @param toGrid grid to which points will be moved; this is the current grid and will have updated points (toPatchPoints) added at end of method
	 * @param toPatchPoints a new array of patchPoints for the toGrid; won't be added to toGrid until completed
	 * @param toRowColGrid describes the geometry of the new toGrid
	 * @param newConnection
	 * @param orderReversed
	 * @return
	 */
	StsPatchGrid moveNonOverlappingPointsTo(StsPatchGrid toGrid, PatchPoint[][] toPatchPoints, RowColGrid toRowColGrid, Connection newConnection, boolean orderReversed)
	{
		int fromRow = -1, toRow = -1;
		int fromCol = -1, toCol = -1;
		PointState[][] pointStates = new PointState[nRows][nCols];
		int toRowStart = rowMin - toRowColGrid.rowMin;
		int toColStart = colMin - toRowColGrid.colMin;
		PatchPoint patchPoint;
		boolean movePoint;
		byte prevColState, prevRowState;
		Connection rowConnection, colConnection;
		PatchPoint newNextPoint = null;
		PatchPoint newPrevPoint = null;
		byte nextPointState, prevPointState;
		PatchPoint fromConnectedPoint, toConnectedPoint;
		PointState fromConnectedState;
		PatchPoint fromPoint, toPoint;
		PointState fromPointState, toPointState;
		// Note that the "from" grid is this grid.
		// increase one grid and decrease the other by number of movedPoints; cloned points don't change count
		int nFromPoints, nToPoints;

		if (StsPatchGrid.debugPatchID != -1 && (id == StsPatchGrid.debugPatchID || toGrid.id == StsPatchGrid.debugPatchID))
			StsException.systemDebug(this, "moveNonOverlappingPointsTo", "MOVE POINTS FROM grid: " + id + " TO grid: " + toGrid.id);

		try
		{
			// save newConnection points as they may be reset by other connections to same point
			newNextPoint = newConnection.window.getCenterPoint();
			newPrevPoint = newConnection.prevWindow.getCenterPoint();

			if (!orderReversed)
			{
				fromConnectedPoint = newNextPoint;
				toConnectedPoint = newPrevPoint;
			}
			else
			{
				fromConnectedPoint = newPrevPoint;
				toConnectedPoint = newNextPoint;
			}
			fromConnectedState = getPointState(pointStates, fromConnectedPoint);

			nFromPoints = nPatchPoints;
			nToPoints = toGrid.nPatchPoints;

			// Set state of newConnection point on fromGrid (which is this grid).
			// fromPoint is newConnection point on fromGrid
			// toPoint is newConnection point on toGrid
			// If fromPoint can be moved, then it doesn't overlap toGrid; set fromState to MOVED;
			// toPoint must overlap the fromGrid (otherwise wouldn't be here).
			// If fromPoint has connection already, we will put a clone on fromGrid, so set the fromState to CLONED.
			// If fromPoint can't be moved, then we will be putting a clone of toPoint on the fromGrid; this will be
			// handled by newConnection section at bottom of routine

			movePoint = toGrid.getVolPatchPoint(fromConnectedPoint) == null;
			if(movePoint)
			{
				if(fromConnectedPoint.hasConnection())
					initializePointState(pointStates, fromConnectedPoint, GRID_CLONED);
				else
					initializePointState(pointStates, fromConnectedPoint, GRID_MOVED);
			}

			// set initial states for fromGrid: NULL, UNCHANGED, MOVED, or ALREADY_CLONED
			// NULL: point doesn't exist
			// UNCHANGED: point exists at this location on toGrid so can't be moved to toGrid
			// MOVED: point doesn't exist at this location on toGrid so point can be moved to toGrid
			// ALREADY_CLONED: point is clone of point on toGrid; point can be moved if connections on fromGrid can also be moved
			// then set the pointState connection flags
			for (fromRow = 0, toRow = toRowStart; fromRow < nRows; fromRow++, toRow++)
			{
				for (fromCol = 0, toCol = toColStart; fromCol < nCols; fromCol++, toCol++)
				{
					patchPoint = patchPoints[fromRow][fromCol];
					if (patchPoint == null) continue; // state flag default is NULL
					byte state;
					movePoint = toPatchPoints[toRow][toCol] == null;
					if (movePoint)
						state = GRID_MOVED;
					else
					{
						if (patchPoint.isCloneFromGrid(toGrid))
							state = GRID_ALREADY_CLONED;
						else
							state = GRID_UNCHANGED;
					}
					PointState pointState = new PointState(state);
					pointStates[fromRow][fromCol] = pointState;

					if (patchPoint.hasRowConnection())
					{
						PointState connectedState = pointStates[fromRow][fromCol - 1];
						connectedState.rowPlusState = state;
						pointState.rowMinusState = connectedState.state;
					}
					if (patchPoint.hasColConnection())
					{
						PointState connectedState = pointStates[fromRow - 1][fromCol];
						connectedState.colPlusState = state;
						pointState.colMinusState = connectedState.state;
					}
				}
			}
			// add the newConnection pointState flags
			int gridRow = fromConnectedPoint.patchGrid.getGridRow(fromConnectedPoint);
			int gridCol = fromConnectedPoint.patchGrid.getGridCol(fromConnectedPoint);
			if(!orderReversed)
			{
				if (newConnection.isRow)
					pointStates[gridRow][gridCol].rowMinusState = GRID_NEW;
				else
					pointStates[gridRow][gridCol].colMinusState = GRID_NEW;
			}
			else
			{
				if (newConnection.isRow)
					pointStates[gridRow][gridCol].rowPlusState = GRID_NEW;
				else
					pointStates[gridRow][gridCol].colPlusState = GRID_NEW;
			}

			// now we can process each of the points based on their pointState flags
			for (fromRow = 0; fromRow < nRows; fromRow++)
			{
				for (fromCol = 0; fromCol < nCols; fromCol++)
				{
					PointState pointState = pointStates[fromRow][fromCol];
					if(pointState == null) continue;
					byte state = pointState.adjustState();
					if(state == GRID_NULL || state == GRID_UNCHANGED)
						continue;
					fromPoint = patchPoints[fromRow][fromCol];
					toPoint = toGrid.getVolPatchPoint(fromPoint);
					if (state == GRID_MOVED)
					{
						moveReplacePoint(fromPoint, null, this, patchPoints, rowColGrid, toGrid, toPatchPoints, toRowColGrid);
						nFromPoints--;
						nToPoints++;
					}
					else if (state == GRID_CLONED)
					{
						PatchPoint clonedPoint = fromPoint.cloneAndClear();
						moveReplacePoint(fromPoint, clonedPoint, this, patchPoints, rowColGrid, toGrid, toPatchPoints, toRowColGrid);
						pointState.moveConnections(fromPoint, clonedPoint);
						nToPoints++;
					}
					else if(state == GRID_ALREADY_CLONED_DELETE)
					{
						pointState.moveConnections(fromPoint, toPoint);
						deletePoint(fromPoint);
						nFromPoints--;
					}
					else if(state == GRID_ALREADY_CLONED)
					{
						pointState.moveConnections(fromPoint, toPoint);
					}
				}
			}

			// MOVED points will change the count, but CLONED points won't
			nPatchPoints = nFromPoints;
			toGrid.nPatchPoints = nToPoints;
			// if points have been moved, this grid will have shrunk, so downsize it
			if (nPatchPoints > 0)
			{
				boolean gridHasPoints = adjustGridSizeDown();
				if (!gridHasPoints)
					patchVolume.removePatchGridFromLists(this);
			}

			StsPatchGrid patchGrid = toGrid;
			// move new toPatchPoints to toGrid as resetConnection() below will use the new points from this toGrid
			toGrid.resetPatchPoints(toRowColGrid, toPatchPoints);
			// if connected points are not on same grid, we need to move the toPoint to the fromGrid if fromState is UNCHANGED
			Connection clonedNewConnection = null;
			if (newNextPoint.patchGrid != newPrevPoint.patchGrid)
			{
				byte fromState = getState(pointStates, fromConnectedPoint);
				byte toState = getState(pointStates, toConnectedPoint);  // here for debugging check
				// fromState is either UNCHANGED, ALREADY_CLONED or MOVED (can't be NULL since point exists on fromGrid)
				// If UNCHANGED: this means the fromPoint overlaps the toGrid which implies that the toPoint doesn't overlap the fromGrid,
				//     so the toPoint will have been moved to fromGrid and may be cloned if any connections already exist to this toPoint on the toGrid.
				// Note fromState might have initially been CLONED; in that case it was reset above to UNCHANGED since CLONED means toPoint exists
				//     and fromPoint must be null since connection points cannot be mutually overlapping;
				//     so toPoint needs to be cloned to fromGrid for newConnection and connection made on fromGrid.
				// If MOVED: then no action is required and newConnection will be added to the new toPoint.

				// toState must be NULL (otherwise connections points would be mutually overlapping which can't be true;
				//     otherwise this routine wouldn't have been called).
				// If fromState is GRID_UNCHANGED: clone the toPoint on toGrid and add to fromGrid (this one).
				// If fromState is GRID_ALREADY_CLONED: we already have a clone of the fromPoint on the toGrid;
				//     use it for the newConnection on the toGrid
				// If fromState is GRID_MOVED or GRID_CLONED, then the fromPoint has been moved to the toGrid so connection is fine.
				if (fromState == GRID_UNCHANGED)
				{
					PatchPoint clonedToPoint = toConnectedPoint.cloneAndClear();
					if (!addPatchPoint(clonedToPoint))
					{
						StsException.systemError(this, "moveNonOverlappingPointsTo", "FAILED TRYING to ADD CLONED POINT TO fromGrid " + toConnectedPoint.toString());
						patchGrid = null;
					}
					else
					{
						clonedNewConnection = newConnection.clone();
						clonedNewConnection.resetConnections(fromConnectedPoint, clonedToPoint);
						patchGrid = fromConnectedPoint.patchGrid;
					}
				}
				else if (fromState == GRID_ALREADY_CLONED)
				{
					clonedNewConnection = newConnection.clone();
					clonedNewConnection.resetConnections(fromConnectedPoint.clonedPoint, toConnectedPoint);
					patchGrid = toConnectedPoint.patchGrid;
				}
				else
				{
					StsException.systemDebug(this, "moveNonOverlappingPointsTo", "NEW CONNECTION state not handled. State: " + stateStrings[fromState] + " CONNECTION: " + newConnection.toString());
					patchGrid = null;
				}
				patchGrid = null;
			}
			nextPointState = getState(pointStates, newNextPoint);
			prevPointState = getState(pointStates, newPrevPoint);
			String directionString = "NEW CONNECTION. FROM grid: " + newNextPoint.getID() + " TO grid: " + newPrevPoint.getID() +
					" FROM STATE: " + stateStrings[nextPointState] + " TO STATE: " + stateStrings[prevPointState];
			newNextPoint.patchGrid.debugCheckOverlappedGridPoint("ORIGINAL CONNECTION nextPoint", newNextPoint);
			newPrevPoint.patchGrid.debugCheckOverlappedGridPoint("ORIGINAL CONNECTION prevPoint", newPrevPoint);
			if (clonedNewConnection != null)
			{
				newNextPoint = clonedNewConnection.window.getCenterPoint();
				newPrevPoint = clonedNewConnection.prevWindow.getCenterPoint();
				newNextPoint.patchGrid.debugCheckOverlappedGridPoint("CLONED CONNECTION nextPoint", newNextPoint);
				newPrevPoint.patchGrid.debugCheckOverlappedGridPoint("CLONED CONNECTION prevPoint", newPrevPoint);
			}
			return patchGrid;
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "moveNonOverlappingPointsTo", "FAILED MOVING POINTS FROM PATCH: " + toString() + " TO PATCH: " + toGrid.toString() +
					" at row: " + (fromRow + rowMin) + " col: " + (fromCol + colMin), e);
			return null;
		}
	}
	/*
	boolean XmoveNonOverlappingPointsTo(StsPatchGrid toGrid, PatchPoint[][] toPatchPoints, RowColGrid toRowColGrid, Connection newConnection, boolean orderReversed)
	{
		int row = -1, newRow = -1;
		int col = -1, newCol = -1;
		byte[][] states = new byte[nRows][nCols];
		int rowStart = rowMin - toRowColGrid.rowMin;
		int colStart = colMin - toRowColGrid.colMin;
		PatchPoint patchPoint;
		boolean movePoint;
		byte prevColState, prevRowState;
		Connection rowConnection, colConnection;
		// increase one grid and decrease the other by number of movedPoints; cloned points don't change count
		int nMovedPoints = 0;
		try
		{
			// first set the state of each point in the grid to be moved: NULL, UNCHANGED, CLONED, or MOVED
			for (row = 0, newRow = rowStart; row < nRows; row++, newRow++)
			{
				for (col = 0, newCol = colStart; col < nCols; col++, newCol++)
				{
					patchPoint = patchPoints[row][col];
					if (patchPoint == null) continue; // state flag default is NULL

					movePoint = toPatchPoints[newRow][newCol] == null;
					prevColState = getPointState(states, row, col - 1); // prevColPoint is in the same row prev col
					prevRowState = getPointState(states, row - 1, col); // prevRowPoint is in the same col prev row

					if (movePoint) // we can move point to unionGrid
					{
						// both connected points are MOVED, CLONED, or NULL: so move this one to unionGrid and delete from this grid
						if (prevColState != GRID_UNCHANGED && prevRowState != GRID_UNCHANGED)
						{
							states[row][col] = GRID_MOVED; // flag as moved: don't actually move until all are checked (in loop below)
							nMovedPoints++;
						}
						// either of two connected points is UNCHANGED, so we need to
						// clone point on this grid and move point to unionGrid
						else
						{
							states[row][col] = GRID_CLONED;
						}
					}
					// can't move point as unionGrid has point at this location
					// if this point has row or col connection and that prevPoint was to be MOVED, then change prevPoint state to CLONED
					else // !movePoint
					{
						states[row][col] = GRID_UNCHANGED;

						if (prevRowState == GRID_MOVED)
						{
							PatchPoint prevRowPoint = patchPoint.getPrevConnectedRowPoint(newConnection);
							if(prevRowPoint != null)
							{
								initializePointState(states, row - 1, col, GRID_CLONED);
								nMovedPoints--;
							}
						}
						if (prevColState == GRID_MOVED)
						{
							PatchPoint prevColPoint = patchPoint.getPrevConnectedColPoint();
							if(prevColPoint != null)
							{
								initializePointState(states, row, col - 1, GRID_CLONED);
								nMovedPoints--;
							}
						}
					}
				}
			}

			PatchPoint nextPoint = newConnection.window.centerPoint;
			PatchPoint prevPoint = newConnection.prevWindow.centerPoint;
			PatchPoint fromPoint, toPoint;
			if(!orderReversed)
			{
				fromPoint = nextPoint;
				toPoint = prevPoint;
			}
			else
			{
				fromPoint = prevPoint;
				toPoint = nextPoint;
			}
			// see if one point of the connection can be moved to the other grid; if not, return false and don't merge nonoverlap points
			byte fromState = getPointState(states, fromPoint);
			byte toState = getPointState(states, toPoint);
			// fromState is either CLONED or MOVED (can't be NULL since point exists), so connection can be added to the toGrid
			if(fromState != GRID_UNCHANGED)
				// this is the "from" grid
				if (!toGrid.patchPointOverlaps(fromPoint))
				{
					if(fromState == GRID_CLONED)
						clonePoint(fromPoint, patchPoints, this.rowColGrid, toGrid, toPatchPoints, toRowColGrid);
					else // fromState == GRID_MOVED
						movePoint(this, fromPoint, patchPoints, this.rowColGrid, toGrid, toPatchPoint, toPatchPoints, toRowColGrid);
				}
				else // !fromGrid.patchPointOverlaps(toPoint);
				{
					if(debug)
					{
						if (!patchPointOverlaps(toPoint))
						{
							StsException.systemError(this, "moveNonOverlappingPointsTo", "POINT SHOULDN'T OVERLAP. point: " + toPoint.toString() +
									"GRID: " + toString());
							return false;
						}
					}
					// Connection can't be moved to toGrid because a point exists where fromPoint is to be moved.
					// But point on toGrid can be moved to fromGrid since state at fromPoint is either CLONED or MOVED.
					// If CLONED: leave a clone at the toPoint and move it to the fromGrid which will move existing connections to the toPoint;
					//            change the state of the fromPoint from NULL to UNCHANGED and check/change states of neighboring points
					// IF MOVED:
					if(fromState == GRID_CLONED)
						clonePoint(fromPoint, patchPoints, this.rowColGrid, toGrid, toPatchPoints, toRowColGrid);
					else // fromState == GRID_MOVED
						movePoint(this, fromPoint, patchPoints, this.rowColGrid, toGrid, toPatchPoint, toPatchPoints, toRowColGrid);

					PatchPoint clonedPoint = smallPatchPoint.clone();
					patchGrid.addPatchPoint(smallPatchPoint);
					smallPatchGrid.replacePoint(clonedPoint);
					connection.resetConnection(smallPatchPoint, orderReversed);
					if (StsPatchVolume.debugOverlapGridGroupOK) patchGrid.combineChildGrids(patchGrid);
					return patchGrid;
				}

			// move nonOverlapping points from this grid to new unionGrid
			// check the newConnection to adjust state of connected points
			setConnectedPointStates(states, newConnection);

			// Having determined the action for each point, MOVE or CLONE the corresponding points
			for (row = 0, newRow = rowStart; row < nRows; row++, newRow++)
			{
				for (col = 0, newCol = colStart; col < nCols; col++, newCol++)
				{
					byte state = states[row][col];
					if (state == GRID_NULL || state == GRID_UNCHANGED) continue;
					patchPoint = patchPoints[row][col];
					if (state == GRID_MOVED)
					{
						movePoint(this, patchPoint, patchPoints, this.rowColGrid, toGrid, toPatchPoint, toPatchPoints, toRowColGrid);
					}
					else if (state == GRID_CLONED)
					{
						clonePoint(patchPoint, patchPoints, this.rowColGrid, toGrid, toPatchPoints, toRowColGrid);
					}
				}
			}
			// move new toPatchPoints to toGrid as resetConnection() below will use the new points from this toGrid
			toGrid.resetPatchPoints(toRowColGrid, toPatchPoints);
			// If this point is UNCHANGED || CLONED, move connections accordingly
			// If a point has been moved, then its connections will have gone with it
			//   (connected prevRow and prevCol points must have been moved as well).
			for (row = 0, newRow = rowStart; row < nRows; row++, newRow++)
			{
				for (col = 0, newCol = colStart; col < nCols; col++, newCol++)
				{
					byte state = states[row][col];
					if (state == GRID_NULL || state == GRID_MOVED) continue;
					PatchPoint fromPoint = patchPoints[row][col];
					PatchPoint toPoint = toPatchPoints[newRow][newCol];
					//patchPoint = patchPoints[row][col];
					//PatchPoint clonedPoint = patchPoint.clonedPoint;
					prevColState = getPointState(states, row, col - 1); // prevColPoint is in the same row and has rowConnection to this point
					prevRowState = getPointState(states, row - 1, col); // prevRowPoint is in the same col and has colConnection to this point
					colConnection = fromPoint.colConnection; // this is the original colConnection
					rowConnection = fromPoint.rowConnection; // this is the original rowConnection

					// point is either UNCHANGED or CLONED

					// if CLONED:  parchPoint has been MOVED to unionGrid;
					//   if prevPoint(s) are UNCHANGED:
					//       connection(s) are from the CLONED fromPoint to the connected point;
					//       remove connections(s) from the toPoint and reset connection.point to CLONED fromPoint
					//   if prevPoint(s) are CLONED:
					//       connection(s) are from toPoint to connected point;
					// 		 remove connection(s) from CLONED fromPoint
					if(state == GRID_CLONED)
					{
						// modify row connection as necessary
						if (rowConnection != null)
						{
							// if prev point in row is UNCHANGED and row connection exists:
							// remove the rowConnection from the moved point (toPoint) and reset connection.point to clonedPoint (fromPoint)
							if (prevColState == GRID_UNCHANGED)
							{
								toPoint.rowConnection = null;
								rowConnection.window.centerPoint = fromPoint;
								// rowConnection.resetConnection(fromPoint);
							}
							// prev point in row is MOVED or CLONED;
							// the connection has been moved: remove the connection with the fromPoint
							else
								fromPoint.rowConnection = null;

						}
						// modify col connection as necessary
						if (colConnection != null)
						{
							// if prev point in col is UNCHANGED and col connection exists:
							// remove the colConnection from the moved point (toPoint) and reset connection.point to clonedPoint (fromPoint)
							if (prevRowState == GRID_UNCHANGED)
							{
								toPoint.colConnection = null;
								colConnection.window.centerPoint = fromPoint;
								//colConnection.resetConnection(fromPoint);
							}
							// prev point in col is MOVED or CLONED;
							// the connection has been moved: remove the connection with the fromPoint
							else
								fromPoint.colConnection = null;
							//colConnection.resetConnection(toPoint);
						}
					}
					// if UNCHANGED:
					//   if prevPoint(s) are UNCHANGED:
					//       connections(s) are ok as is
					//   if prevPoint(s) are CLONED (these point(s) could not have been MOVED (see logic in prev loop):
					//       connection(s) are now to fromPoint (which is UNCHANGED);
					//       so replace corresponding connection.point with fromPoint: FYI: there is no toPoint
					//
					else // point is UNCHANGED
					{
						// if prev point in row is CLONED, replaced connection.otherPoint with this cloned point
						if(prevColState == GRID_CLONED && rowConnection != null)
						{
							PatchPoint prevColPoint = patchPoints[row][col - 1]; // rowConnection connects to point on same row, prev col
							rowConnection.prevWindow.centerPoint = prevColPoint;
							//rowConnection.resetConnection(fromPoint);
						}
						// if prev point in col is CLONED, replaced connection.otherPoint with this cloned point
						if(prevRowState == GRID_CLONED && colConnection != null)
						{
							PatchPoint prevRowPoint = patchPoints[row - 1][col]; // colConnection connects to point on prev row, same col
							colConnection.prevWindow.centerPoint = prevRowPoint;
							//colConnection.resetConnection(fromPoint);
						}
					}
				}
			}
			// MOVED points will change the count, but CLONED points won't
			nPatchPoints -= nMovedPoints;
			toGrid.nPatchPoints += nMovedPoints;
			// if points have been moved, this grid will have shrunk, so downsize it
			if(nPatchPoints > 0) adjustGridSizeDown();
			return true;
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "moveNonOverlappingPointsTo", "FAILED MOVING POINTS FROM PATCH: " + toString()+ " TO PATCH: " + toGrid.toString() +
					" at row: " + (row + rowMin) + " col: " + (col + colMin), e);
			return false;
		}
	}
    */

	/**
	 * clone point fromPatchPoint on fromGrid and replace the point on fromGrid with clonedPoint; move fromPatchPoint to the toGrid
	 * @param fromPatchPoint point on fromGrid to be moved to the toGrid and replaced with a clone of itself
	 * @param fromPatchPoints patchPoints on fromGrid
	 * @param fromRowColGrid grid description of the fromGrid
	 * @param toGrid grid to which fromPatchPoint is to be moved
	 * @param toPatchPoints array of new patchPoints on the toGrid (not yet a member of the toGrid)
	 * @param toRowColGrid grid description of the toGrid
	 */
	static void clonePoint(PatchPoint fromPatchPoint, StsPatchGrid fromGrid, PatchPoint[][] fromPatchPoints, RowColGrid fromRowColGrid, StsPatchGrid toGrid, PatchPoint[][] toPatchPoints, RowColGrid toRowColGrid)
	{
		moveReplacePoint(fromPatchPoint, fromPatchPoint.clone(), fromGrid, fromPatchPoints, fromRowColGrid, toGrid, toPatchPoints, toRowColGrid);
	}

	/**
	 * move point patchPoint from fromGrid to toGrid and null point on fromGrid.
	 * @param fromPatchPoint point on fromGrid to be moved to the toGrid
	 * @param fromPatchPoints patchPoints on fromGrid
	 * @param fromRowColGrid grid description of the fromGrid
	 * @param toGrid grid to which patchPoint is to be moved
	 * @param toPatchPoints array of new patchPoints on the toGrid (not yet a member of the toGrid)
	 * @param toRowColGrid grid description of the toGrid
	 */
	static void movePoint(PatchPoint fromPatchPoint, StsPatchGrid fromGrid, PatchPoint[][] fromPatchPoints, RowColGrid fromRowColGrid, StsPatchGrid toGrid, PatchPoint[][] toPatchPoints, RowColGrid toRowColGrid)
	{
		moveReplacePoint(fromPatchPoint, null, fromGrid, fromPatchPoints, fromRowColGrid, toGrid, toPatchPoints, toRowColGrid);
	}

	/**
	 * move or clone point patchPoint from fromGrid to toGrid and null point on fromGrid if moved. Reset grid id for point on toGrid.
	 * @param fromPatchPoint point on fromGrid to be moved to the toGrid and replaced with replaceFromPatchPoint
	 * @param replaceFromPatchPoint point to replace point on fromGrid: null if moved, clonedPoint if cloned
	 * @param fromPatchPoints patchPoints on fromGrid
	 * @param fromRowColGrid grid description of the fromGrid
	 * @param toGrid grid to which patchPoint is to be moved
	 * @param toPatchPoints array of new patchPoints on the toGrid (not yet a member of the toGrid)
	 * @param toRowColGrid grid description of the toGrid
	 */
	static void moveReplacePoint(PatchPoint fromPatchPoint, PatchPoint replaceFromPatchPoint, StsPatchGrid fromGrid, PatchPoint[][] fromPatchPoints, RowColGrid fromRowColGrid, StsPatchGrid toGrid, PatchPoint[][] toPatchPoints, RowColGrid toRowColGrid)
	{
		fromPatchPoint.patchGrid = toGrid;

		int volRow = fromPatchPoint.getRow();
		int volCol = fromPatchPoint.getCol();

		int fromGridRow = volRow - fromRowColGrid.rowMin;
		int fromGridCol = volCol - fromRowColGrid.colMin;
		fromPatchPoints[fromGridRow][fromGridCol] = replaceFromPatchPoint;
		if (debug && debugPoint && (doDebugPoint(fromPatchPoint)))
			StsException.systemDebug(StsPatchGrid.class, "moveReplacePoint", StsPatchVolume.iterLabel + "patchPoint " + PatchPoint.staticToString(replaceFromPatchPoint));

		int toGridRow = volRow - toRowColGrid.rowMin;
		int toGridCol = volCol - toRowColGrid.colMin;
		toPatchPoints[toGridRow][toGridCol] = fromPatchPoint;
		if (debug && debugPoint && (doDebugPoint(fromPatchPoint)))
			StsException.systemDebug(StsPatchGrid.class, "moveReplacePoint", StsPatchVolume.iterLabel + "patchPoint " + PatchPoint.staticToString(fromPatchPoint));

		if (debug)
			debugCheckMovedClonedPoint(fromGrid, fromPatchPoint, toGrid, toRowColGrid);
	}

	/**
	 * move point patchPoint from toGrid to fromGrid and null point on fromGrid
	 * @param fromPatchPoint point on fromGrid to be moved to the toGrid and replaced with a null
	 * @param fromPatchPoints patchPoints on fromGrid
	 * @param fromRowColGrid grid description of the fromGrid
	 * @param toGrid grid to which patchPoint is to be moved
	 * @param toPatchPoints array of new patchPoints on the toGrid (not yet a member of the toGrid)
	 * @param toRowColGrid grid description of the toGrid
	 */
	static void reverseMovePoint(PatchPoint fromPatchPoint, StsPatchGrid fromGrid, PatchPoint[][] fromPatchPoints, RowColGrid fromRowColGrid, StsPatchGrid toGrid, PatchPoint[][] toPatchPoints, RowColGrid toRowColGrid)
	{
		fromPatchPoint.patchGrid = toGrid;

		int volRow = fromPatchPoint.getRow();
		int volCol = fromPatchPoint.getCol();

		int fromGridRow = volRow - fromRowColGrid.rowMin;
		int fromGridCol = volCol - fromRowColGrid.colMin;
		fromPatchPoints[fromGridRow][fromGridCol] = null;

		int toGridRow = volRow - toRowColGrid.rowMin;
		int toGridCol = volCol - toRowColGrid.colMin;
		toPatchPoints[toGridRow][toGridCol] = fromPatchPoint;

		if (debug)
			debugCheckMovedClonedPoint(fromGrid, fromPatchPoint, toGrid, toRowColGrid);
	}

	static void movePoint(StsPatchGrid fromGrid, PatchPoint fromPatchPoint, PatchPoint[][] fromPatchPoints, RowColGrid fromRowColGrid, StsPatchGrid toGrid, PatchPoint toPatchPoint, PatchPoint[][] toPatchPoints, RowColGrid toRowColGrid)
	{
		fromPatchPoint.patchGrid = toGrid;

		int volRow = fromPatchPoint.getRow();
		int volCol = fromPatchPoint.getCol();

		int fromGridRow = volRow - fromRowColGrid.rowMin;
		int fromGridCol = volCol - fromRowColGrid.colMin;
		fromPatchPoints[fromGridRow][fromGridCol] = fromPatchPoint;

		int toGridRow = volRow - toRowColGrid.rowMin;
		int toGridCol = volCol - toRowColGrid.colMin;
		toPatchPoints[toGridRow][toGridCol] = toPatchPoint;

		if (debug)
			debugCheckMovedClonedPoint(fromGrid, fromPatchPoint, toGrid, toRowColGrid);
	}

	void replacePoint(PatchPoint patchPoint)
	{
		int volRow = patchPoint.getRow();
		int volCol = patchPoint.getCol();

		int fromGridRow = volRow - rowColGrid.rowMin;
		int fromGridCol = volCol - rowColGrid.colMin;
		patchPoints[fromGridRow][fromGridCol] = patchPoint;
	}

	void deletePoint(PatchPoint patchPoint)
	{
		int volRow = patchPoint.getRow();
		int volCol = patchPoint.getCol();

		int fromGridRow = volRow - rowColGrid.rowMin;
		int fromGridCol = volCol - rowColGrid.colMin;
		if (debug && debugPoint)
		{
			PatchPoint deletedPoint = patchPoints[fromGridRow][fromGridCol];
			if (doDebugPoint(deletedPoint))
				StsException.systemDebug(StsPatchGrid.class, "deletePoint", StsPatchVolume.iterLabel + " DELETE POINT: " + PatchPoint.staticToString(deletedPoint));
		}
		patchPoints[fromGridRow][fromGridCol] = null;

	}

	/**
	 * This newConnection is between different grids. For reference this is the grid is the "from" grid: points are being moved from it.
	 * A point on the "from" grid can have four states: NULL (doesn't exist), UNCHANGED (point exists at location on "to" grid),
	 * CLONED (point has been moved to "to" grid but connections on "from" grid which aren't moved require a clone there), and
	 * MOVED (point is moved to "to" grid and connected points are either moved or cloned so connection has gone with the pair).
	 * There are sixteen combinations between connected pairs, but reversed pairs are the same (U-C is same as C-U, though there might be
	 * an action on only the "C" in this case).  Thus there are 10 unique cases: NN,UU,CC,MM,NU,NC,NM,UC,UM,CM.
	 * Here are the actions for each of these cases:
	 * These four cases won't occur on a new connection since new connection is between the two grids.
	 * NN: both points are on "to" grid, so connection is ok as is.
	 * UU: both points are on the "from" grid, so connection is ok as is.
	 * CC: both points are on the "to" grid, but clones remain on "from" grid so connection needs to be removed from both clones.
	 * MM: both points are on the "to" grid haveing been movd there, so connection is ok.
	 * These three cases will only occur for a new connection as they involve a connection between the two grids;
	 * NU,UN: try moving the N point on the "to" grid to the "from" grid (either move or clone); if this fails, a new grid is required for the two.
	 * NC,CN: N point is on the "to" grid; C point has been moved to "to" grid, but clone remains to handle other connection which wasn't moved.
	 * Add newConnection to point.
	 * NM,MN: M point has been moved to "to" grid; N point is on "to" grid, so connection is ok as is.
	 * These three cases won't occur on a new connection since new connection is between two grids.
	 * UC, CU: C point has been moved to "to" grid, remove connection from moved point and change connection to cloned point.
	 * UM,MU: Change M point state to C and apply UC logic: point remains on "to" grid, but clone is on "from" grid;
	 * remove connection from moved point and change connection to cloned point.
	 * CM,MC: C point has been moved to "to" grid and cloned on "from" grid; remove connection from C point.
	 * If a point state is NULL, then it doesn't exist on the "from" grid, but does exist on the "moved" grid.
	 * So if the other point is MOVED or CLONED, then both original points are then on the "to" grid and connection is ok,
	 * i.e., both points are either MOVED/CLONED or NULL and connection is ok as is.
	 * If both points are UNCHANGED, then they are on the "from" grid and connection is ok.
	 * If a point is MOVED and the other is UNCHANGED, then the MOVED point needs to be changed to CLONED so connection stays on the "from" grid;
	 * corresponding connection point needs to be changed to CLONED point.
	 * If both points in connection are MOVED or both are UNCHANGED, connection is fine as is.
	 * If one or the other is NULL, then it doesn't exist on this
	 * If their states are different
	 */

	static private void debugCheckMovedClonedPoint(StsPatchGrid fromGrid, PatchPoint fromPatchPoint, StsPatchGrid toGrid, RowColGrid toRowColGrid)
	{
		if (StsPatchGrid.debugPoint && (StsPatchGrid.doDebugPoint(fromPatchPoint)))
			StsException.systemDebug(fromGrid, "debugCheckMovedClonedPoint", StsPatchVolume.iterLabel +
					" MOVED POINT TO GRID. point: " + fromPatchPoint.toString() + " GRID: " + toGrid.toString());
	}

	/**
	 * given row-col array of states, return state at this row & col
	 * called for point in array or previous row or col point, so currently need only check minimum ranges
	 * @param states 2D byte array of states (NULL, UNChANGED, CLONED, or MOVED.
	 * @param row row of desired state flag
	 * @param col col of desired state flag
	 * @return state at this row and col
	 */
	private final PointState getPointState(PointState[][] states, int row, int col)
	{
		if (col < 0 || col >= states[0].length || row < 0 || row >= states.length) return null;
		return states[row][col];
	}

	private final PointState getPointState(PointState[][] states, PatchPoint point)
	{
		return getPointState(states, getGridRow(point), getGridCol(point));
	}

	private final byte getState(PointState[][] states, int row, int col)
	{
		PointState pointState = getPointState(states, row, col);
		if(pointState == null) return GRID_NULL;
		return pointState.state;
	}

	private final byte getState(PointState[][] states, PatchPoint point)
	{
		return getState(states, getGridRow(point), getGridCol(point));
	}

	/**
	 * given row-col array of states, set state at this row & col
	 * called for point in array or previous row or col point, so currently need only check minimum ranges
	 * @param states 2D byte array of states (NULL, UNChANGED, CLONED, or MOVED.
	 * @param row row of state flag to be set
	 * @param col col of state flag to be set
	 * @param state state flag to be set
	 * @return none
	 */
	private final void initializePointState(PointState[][] states, int row, int col, byte state)
	{
		try
		{
			states[row][col] = new PointState(state);
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "initializePointState", e);
		}
	}

	private final void initializePointState(PointState[][] states, PatchPoint point, byte state)
	{
		initializePointState(states, getGridRow(point), getGridCol(point), state);
	}

	private boolean isMoved(byte state)
	{
		return state == GRID_MOVED;
	}

	private boolean isMovedOrCloned(byte state)
	{
		return state == GRID_MOVED || state == GRID_CLONED;
	}

	void resetPatchPoints(RowColGrid newRowColGrid, PatchPoint[][] newPatchPoints)
	{
		initializeRowColGrid(newRowColGrid);
		patchPoints = newPatchPoints;
	}

	boolean contains(RowColGrid newRowColGrid)
	{
		return newRowColGrid.rowMin >= rowMin && newRowColGrid.rowMax <= rowMax &&
				newRowColGrid.colMin >= colMin && newRowColGrid.colMax <= colMax;
	}

	boolean contains(PatchPoint patchPoint)
	{
		int row = patchPoint.getRow();
		int col = patchPoint.getCol();
		return row >= rowMin && row <= rowMax && col >= colMin && col <= colMax;
	}

	private void initializePatchPoints(PatchPoint patchPoint)
	{
		initializeRowColGrid(patchPoint);
		patchPoints = new PatchPoint[nRows][nCols];
	}

	/**
	 * Called to add a window which doesn't belong to a grid to this patchGrid as long as it doesn't overlap.
	 * If it does overlap, create a new patchGrid and add this window and clone of connectedPatchPoint which
	 * belongs to this grid.
	 */
	public StsPatchGrid checkAddPatchPoint(PatchPoint patchPoint, PatchPoint connectedPatchPoint)
	{
		if (patchPointOverlaps(patchPoint))
		{
			StsPatchGrid newPatchGrid = StsPatchGrid.construct(patchVolume, patchType);
			newPatchGrid.addPatchPoint(patchPoint);
			newPatchGrid.addPatchPoint(connectedPatchPoint.cloneAndClear());
			return newPatchGrid;
		}
		else
		{
			addPatchPoint(patchPoint);
			return this;
		}
	}
/*
	public final void addCorrelation(PatchPoint otherPatchPoint, PatchPoint newPatchPoint, float correl)
	{
		if (otherPatchPoint.getRow() == newPatchPoint.getRow())
			addRowCorrelation(otherPatchPoint, newPatchPoint, correl);
		else
			addColCorrelation(otherPatchPoint, newPatchPoint, correl);
	}

	public final void addRowCorrelation(PatchPoint otherPatchPoint, PatchPoint newPatchPoint, float correl)
	{
		if(debug && debugPoint && (doDebugPoint(otherPatchPoint) || doDebugPoint(newPatchPoint)))
			StsException.systemDebug(this, "addRowCorrelation", StsPatchVolume.iterLabel + "adding row correl to " + otherPatchPoint.toString());
		newPatchPoint.rowConnection =  otherPatchPoint.rowCorrel = correl;
	}

	public final void addColCorrelation(PatchPoint otherPatchPoint, PatchPoint newPatchPoint, float correl)
	{
		if (debug && debugPoint && (doDebugPoint(otherPatchPoint) || doDebugPoint(newPatchPoint)))
			StsException.systemDebug(this, "addColCorrelation", StsPatchVolume.iterLabel + "adding row correl to " + otherPatchPoint.toString());
		otherPatchPoint.colCorrel = correl;
	}
*/

	/** debug this point if point row and col are set and they match debug criteria or point.patchGrid.id matches criteria */
	public static final boolean doDebugPoint(PatchPoint patchPoint)
	{
		if (debugPoint)
		{
			if (patchPoint == null || patchPoint.getRow() != debugPointRow || patchPoint.getCol() != debugPointCol)
				return false;
			if (patchPoint.slice == debugPointSlice) return true;
		}
		return debug && debugPatchGrid && patchPoint.patchGrid != null && patchPoint.patchGrid.id == debugPatchID;
	}

	public static final boolean doDebugPoint(PatchPoint patchPoint, int volRow, int volCol)
	{
		if (debugPoint)
		{
			if (patchPoint == null || volRow != debugPointRow || volCol != debugPointCol)
				return false;
			if (patchPoint.slice == debugPointSlice) return true;
		}
		return debug && debugPatchGrid && patchPoint.patchGrid != null && patchPoint.patchGrid.id == debugPatchID;
	}

	public void clear()
	{
		initializeRowColGrid();
		patchPoints = null;
		if (debugPatchID != NO_DEBUG && id == debugPatchID)
			//if(debugPatchGrid())
			StsException.systemDebug(this, "clear", "clearing patch: " + toString());
	}

	private void initializeRowColGrid()
	{
		initializeRowColGrid(largeInt, -largeInt, largeInt, -largeInt);
	}

	private void initializeRowColGrid(int rowMin, int rowMax, int colMin, int colMax)
	{
		initializeRowColGrid(new RowColGrid(rowMin, rowMax, colMin, colMax));
	}

	private void initializeRowColGrid(PatchPoint patchPoint)
	{
		int row = patchPoint.getRow();
		int col = patchPoint.getCol();
		initializeRowColGrid(row, row, col, col);
	}

	private void initializeRowColGrid(RowColGrid newRowColGrid)
	{
		if (debugPatchID != NO_DEBUG && id == debugPatchID)
			//if(debugPatchGrid())
			StsException.systemDebug(this, "initializeRowColGrid", "FOR GRID: " + id + " RESET rowColGrid FROM " + rowColGrid + " TO " + newRowColGrid);
		rowMin = newRowColGrid.rowMin;
		rowMax = newRowColGrid.rowMax;
		colMin = newRowColGrid.colMin;
		colMax = newRowColGrid.colMax;
		nRows = rowMax - rowMin + 1;
		nCols = colMax - colMin + 1;
		this.rowColGrid = newRowColGrid;
	}

	boolean isDisconnected(int row)
	{
		return rowMax < row;
	}

	boolean isOnePoint()
	{
		return rowMax - rowMin <= 0 && colMax - colMin <= 0;
	}

	boolean isTooSmall(int minNPoints)
	{
		return nPatchPoints < minNPoints;
	}

	public void resetIndex(int index)
	{
		originalID = id;
		id = index;
		if (debugPatchID != -1 && originalID == debugPatchID)
			StsException.systemDebug(this, "resetIndex", "debugPatch id being reset from " + originalID + " to " + id);
	}

	public void finish()
	{
		if (debug && debugPatchGrid && id == debugPatchID)
			StsException.systemDebug(this, "finish", "for patch " + toGridString());

		if (pointsZ != null)
		{
			//	if(debug)
			//		StsException.systemDebug(this, "finish", "PATCH already finished. PATCH: " + toGridString());
			return;
		}
		try
		{
			pointsZ = new float[nRows][nCols];
			rowCorrels = new float[nRows][nCols];
			colCorrels = new float[nRows][nCols];
			for (int row = 0; row < nRows; row++)
			{
				for (int col = 0; col < nCols; col++)
				{
					if (patchPoints[row][col] != null)
					{
						float z = patchPoints[row][col].z;
						pointsZ[row][col] = z;
						zMin = Math.min(zMin, z);
						zMax = Math.max(zMax, z);
						rowCorrels[row][col] = getRowCorrelation(row, col);
						colCorrels[row][col] = getColCorrelation(row, col);
					}
					else
					{
						pointsZ[row][col] = nullValue;
						rowCorrels[row][col] = 0.0f;
						colCorrels[row][col] = 0.0f;
					}
				}
			}
			if (!StsPatchVolume.debug) patchPoints = null;
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "finish", e);
		}
	}

	private float getRowCorrelation(int row, int col)
	{
		if (col >= nCols - 1) return 0.0f;
		PatchPoint otherPoint = patchPoints[row][col + 1];
		if (otherPoint == null) return 0.0f;
		Connection connection = otherPoint.getRowConnection();
		if (connection == null || connection.isMoved) return 0.0f;
		//return connection.stretchCorrelation;
		return connection.window.stretchCorrelation;
	}

	private float getColCorrelation(int row, int col)
	{
		if (row >= nRows - 1) return 0.0f;
		PatchPoint otherPoint = patchPoints[row + 1][col];
		if (otherPoint == null) return 0.0f;
		Connection connection = otherPoint.getColConnection();
		if (connection == null || connection.isMoved) return 0.0f;
		//return connection.stretchCorrelation;
		return connection.window.stretchCorrelation;
	}

	/** get nearest patch whose z at this x,y is just above the z slice plane */
	public float getZDistance(int volumeRow, int volumeCol, float z)
	{
		int patchRowMin = Math.max(0, volumeRow - rowMin - 5);
		int patchRowMax = Math.min(nRows - 1, volumeRow - rowMin + 5);
		int patchColMin = Math.max(0, volumeCol - colMin - 5);
		int patchColMax = Math.min(nCols - 1, volumeCol - colMin + 5);
		float dz = StsParameters.largeFloat;
		for (int patchRow = patchRowMin; patchRow <= patchRowMax; patchRow++)
		{
			for (int patchCol = patchColMin; patchCol <= patchColMax; patchCol++)
			{
				float zPatch = getPatchPointZ(patchRow, patchCol);
				if (zPatch == nullValue) continue;
				float dzPatch = z - zPatch;
				if (dzPatch < 0.0f) continue;
				if (dzPatch < dz) dz = dzPatch;
			}
		}
		return dz;
	}

	/** sort first by rowMin and then by colMin. Return 1 if this rowMin&colMin come after other; 0 if equal; -1 otherwise */
	public int compareTo(StsPatchGrid otherGrid)
	{
		if (sortRowFirst)
		{
			if (rowMin > otherGrid.rowMin) return 1;
			if (rowMin < otherGrid.rowMin) return -1;
			// on the same row
			if (colMin > otherGrid.colMin) return 1;
			if (colMin < otherGrid.colMin) return -1;
			return 0;
		}
		else
		{
			if (colMin > otherGrid.colMin) return 1;
			if (colMin < otherGrid.colMin) return -1;
			// on the same col
			if (rowMin > otherGrid.rowMin) return 1;
			if (rowMin < otherGrid.rowMin) return -1;
			return 0;
		}
	}

	public float[][] getPointsZ()
	{
		return pointsZ;
	}

	public int getGridSize()
	{
		return nRows * nCols;
	}

	public int[] getGridPointsUsed()
	{
		int nUsed = 0;
		int nActualUsed = 0;
		for (int row = 0; row < nRows; row++)
		{
			int colStart = 0, colEnd = nCols - 1;
			for (int col = 0; col < nCols; col++)
			{
				if (pointsZ[row][col] != nullValue)
				{
					colStart = col;
					break;
				}
			}
			for (int col = nCols - 1; col > 0; col--)
			{
				if (pointsZ[row][col] != nullValue)
				{
					colEnd = col;
					break;
				}
			}

			for (int col = colStart; col <= colEnd; col++)
				if (pointsZ[row][col] != nullValue) nActualUsed++;

			nUsed += colEnd - colStart + 1;
		}
		return new int[]{nUsed, nActualUsed};
	}

	public float getDataMin()
	{
		return dataMin;
	}

	public float getDataMax()
	{
		return dataMax;
	}

	public float getzMax()
	{
		return zMax;
	}

	public float getzMin()
	{
		return zMin;
	}

	public boolean computeCurvature(float xInc, float yInc, byte curveType, int filterSize, int minNPoints)
	{
		dataMin = StsPatchVolume.largeFloat;
		dataMax = -StsPatchVolume.largeFloat;

		values = new float[nRows][nCols];
		for (int row = 0; row < nRows; row++)
			Arrays.fill(values[row], nullValue);

		if (nPatchPoints < minNPoints) return false;

		int halfWindow = filterSize / 2;
		// Determine quadratic coefficients for this neighborhood

		float[][] fitPoints = new float[filterSize * filterSize][3];
		for (int volumeRow = rowMin; volumeRow <= rowMax; volumeRow++)
		{
			for (int volumeCol = colMin; volumeCol <= colMax; volumeCol++)
			{
				int nFitPoints = 0;  // number of equations
				int patchPointRow = volumeRow - rowMin;
				int patchPointCol = volumeCol - colMin;
				float zc = getPatchPointZ(patchPointRow, patchPointCol);
				if (zc == StsParameters.nullValue) continue;
				int patchRowMin = Math.max(0, patchPointRow - halfWindow);
				int patchRowMax = Math.min(nRows - 1, patchPointRow + halfWindow);
				int patchColMin = Math.max(0, patchPointCol - halfWindow);
				int patchColMax = Math.min(nCols - 1, patchPointCol + halfWindow);
				float y = (patchRowMin - patchPointRow) * yInc;
				for (int patchRow = patchRowMin; patchRow <= patchRowMax; patchRow++, y += yInc)
				{
					float x = (patchColMin - patchPointCol) * xInc;
					for (int patchCol = patchColMin; patchCol <= patchColMax; patchCol++, x += xInc)
					{
						float z = pointsZ[patchRow][patchCol];
						if (z == StsParameters.nullValue) continue;
						fitPoints[nFitPoints][0] = x;
						fitPoints[nFitPoints][1] = y;
						fitPoints[nFitPoints][2] = z - zc;
						nFitPoints++;
					}
				}
				if (nFitPoints < minNPoints) continue;

				if (!StsQuadraticCurvature.computeSVD(fitPoints, nFitPoints)) continue;

				float val;
				try
				{
					val = StsQuadraticCurvature.getCurvatureComponent(curveType);
				}
				catch (Exception e)
				{
					StsException.systemError(this, "computeCurvature", "getCurvatureComponent failed.");
					continue;
				}

				if (filterType == FILTER_ON_CHI_SQ)
				{
					double chiSqrTest = chiSqrMultiplyer * nFitPoints;
					double chiSqr = StsQuadraticCurvature.computeChiSquared();
					if (chiSqr > chiSqrTest)
					{
						// if(StsPatchVolume.debugPatchGrid) StsException.systemDebug(this, "computeCurvature", "ChiSqr = " + chiSqr + " at volumeRow, volumeCol " + volumeRow + " " + volumeCol);
						//continue;
						if (val > 0) val = badCurvature;
						if (val < 0) val = -badCurvature;
					}
				}

				values[patchPointRow][patchPointCol] = val;

				if (Math.abs(val) > curvatureTest) continue;
				// ChiSqr filtered vals not used for dataMin / dataMax & Statistics
				nValuePatchPoints++;
				sum += val;
				dataMax = Math.max(dataMax, val);
				dataMin = Math.min(dataMin, val);
			}
		}
		return nValuePatchPoints > 0;
	}

	public void drawRow(GL gl, int volumeRow, float y, float xMin, float xInc, StsColorscale colorscale, boolean is3d, boolean displayCurvature, boolean smooth, int boxFilterWidth)
	{
		boolean lineStarted = false;
		float z;
		int col = -1;

		if (volumeRow < rowMin || volumeRow > rowMax) return;
		try
		{
			float x = xMin + colMin * xInc;
			int row = volumeRow - rowMin;
			float[] rowZ;
			if (smooth)
				rowZ = boxFilter2dRow(pointsZ, row, boxFilterWidth, nullValue);
			else
				rowZ = pointsZ[row];

			gl.glDepthFunc(GL.GL_LEQUAL);
			// gl.glLineStipple(1, StsGraphicParameters.dottedLine);
			boolean displayCurvatureColor = (values != null && displayCurvature);
			if (!displayCurvatureColor)
				StsTraceUtilities.getPointTypeColor(patchType).setGLColor(gl);
			for (col = 0; col < nCols; col++, x += xInc)
			{
				z = rowZ[col];
				if (z != nullValue)
				{
					if (displayCurvatureColor)
					{
						// StsTraceUtilities.getPointTypeColor(patchType).setGLColor(gl);
						float v = values[row][col];
						if (v == nullValue)
							StsColor.BLACK.setGLColor(gl);
						else
							colorscale.getStsColor(colorscale.getIndexFromValue(v)).setGLColor(gl);
					}
					if (!lineStarted)
					{
						gl.glBegin(GL.GL_LINE_STRIP);
						lineStarted = true;
					}
					if (is3d)
						gl.glVertex3f(x, y, z);
					else
						gl.glVertex2f(x, z);

					if (rowCorrels[row][col] == 0.0f && lineStarted)
					{
						lineStarted = false;
						gl.glEnd();
					}
				}
				else if (lineStarted)
				{
					lineStarted = false;
					gl.glEnd();
				}
			}
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "drawRow", "Failed for patchGrid " + id + " at row: " + volumeRow + " col: " + col, e);
		}
		finally
		{
			if (lineStarted) gl.glEnd();
			// if(drawingDotted)gl.glDisable(GL.GL_LINE_STIPPLE);
		}
	}

	float[] boxFilter2dRow(float[][] values, int row, int filterWidth, float nullValue)
	{
		FilterVector center = new FilterVector(values[row], nullValue);
		int sideRow;
		for (int i = 0; i < filterWidth; i++)
		{
			// add plus side rows
			sideRow = row + i + 1;
			if (sideRow < nRows)
			{
				FilterVector side = new FilterVector(values[sideRow], nullValue);
				center.add(side);
			}
			sideRow = row - i - 1;
			if (sideRow >= 0)
			{
				FilterVector side = new FilterVector(values[sideRow], nullValue);
				center.add(side);
			}
			center.sum(filterWidth);
		}
		return center.average();
	}

	class FilterVector
	{
		float[] values;
		int[] count;
		int nValues;

		FilterVector(float[] values, float nullValue)
		{
			this.nValues = values.length;
			this.values = new float[nValues];
			count = new int[nValues];
			for (int n = 0; n < nValues; n++)
			{
				if (values[n] != nullValue)
				{
					count[n] = 1;
					this.values[n] = values[n];
				}
			}
		}

		FilterVector(FilterVector vector)
		{
			nValues = vector.nValues;
			values = Arrays.copyOf(vector.values, nValues);
			count = Arrays.copyOf(vector.count, nValues);
		}

		void add(FilterVector side)
		{
			for (int i = 0; i < nValues; i++)
			{
				values[i] += side.values[i];
				count[i] += side.count[i];
			}
		}

		void sum(int filterWidth)
		{
			FilterVector center = new FilterVector(this);
			for (int i = 0; i < filterWidth; i++)
				add(center, i + 1);
		}

		void add(FilterVector shift, int offset)
		{
			//  add values to right of position
			for (int i = offset; i < nValues; i++)
			{
				values[i - offset] += shift.values[i];
				count[i - offset] += shift.count[i];
			}
			// add values to left of position
			for (int i = 0; i < nValues - offset; i++)
			{
				values[i + offset] += shift.values[i];
				count[i + offset] += shift.count[i];
			}
		}

		float[] average()
		{
			for (int n = 0; n < nValues; n++)
				if (count[n] > 0) values[n] /= count[n];
			return values;
		}
	}

	public void drawCol(GL gl, int volumeCol, float x, float yMin, float yInc, StsColorscale colorscale, boolean is3d, boolean displayCurvature)
	{
		boolean lineStarted = false;
		float z;
		int row = -1;

		try
		{
			if (volumeCol < colMin || volumeCol > colMax) return;
			float y = yMin + rowMin * yInc;
			int col = volumeCol - colMin;
			gl.glDepthFunc(GL.GL_LEQUAL);
			boolean displayCurvatureColor = values != null && displayCurvature;
			if (!displayCurvatureColor)
				StsTraceUtilities.getPointTypeColor(patchType).setGLColor(gl);
			for (row = 0; row < nRows; row++, y += yInc)
			{
				z = getPatchPointZ(row, col);

				if (z != StsParameters.nullValue)
				{
					if (displayCurvatureColor)
					{
						float v = values[row][col];
						if (v == nullValue)
							StsColor.BLACK.setGLColor(gl);
						else
							colorscale.getStsColor(colorscale.getIndexFromValue(v)).setGLColor(gl);
					}
					if (!lineStarted)
					{
						gl.glBegin(GL.GL_LINE_STRIP);
						lineStarted = true;
					}
					if (is3d)
						gl.glVertex3f(x, y, z);
					else
						gl.glVertex2f(y, z);

					if (colCorrels[row][col] == 0.0f && lineStarted)
					{
						lineStarted = false;
						gl.glEnd();
					}
				}
				else if (lineStarted)
				{
					lineStarted = false;
					gl.glEnd();
				}
			}
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "drawCol", "Failed for patchGrid " + id + " at row: " + row + " col: " + volumeCol, e);
		}
		finally
		{
			if (lineStarted) gl.glEnd();
		}
	}

	public boolean isPatchGridNearZCursor(float z)
	{
		return z >= zMin && z <= zMax;
	}


	public void drawPatchGrid(GL gl, boolean displayChildPatches, boolean displayCurvature, StsColorscale colorscale)
	{
		draw(gl, displayCurvature, colorscale);
		if (childGrid == null || !displayChildPatches) return;
		childGrid.drawPatchGrid(gl, displayChildPatches, displayCurvature, colorscale);
	}

	public void draw(GL gl, boolean displayCurvature, StsColorscale colorscale)
	{
		if (pointsZ == null) return;
		if (diamondStrips == null) diamondStrips = new StsDiamondStrips(this);
		if (displayCurvature && values != null)
		{
			diamondStrips.setValues(values);
			diamondStrips.drawSurfaceFillWithNulls(gl);
		}
		else
		{
			StsColor patchColor = StsTraceUtilities.getPointTypeColor(patchType);
			patchColor.setGLColor(gl);
			diamondStrips.drawSurfaceFillWithNulls(gl);
		}

		// draw grid lines
		StsColor gridColor = StsColor.BLACK;
		gridColor.setGLColor(gl);
		gl.glDisable(GL.GL_LIGHTING);
		StsGLPanel3d glPanel3d = currentModel.getGlPanel3d();
		gl.glLineWidth(StsGraphicParameters.gridLineWidth);
		glPanel3d.setViewShift(gl, StsGraphicParameters.gridShift);
		diamondStrips.drawGridLines(gl);
		glPanel3d.resetViewShift(gl);
		gl.glEnable(GL.GL_LIGHTING);
		// TriangleStrip.drawSurfaceFillWithNulls(gl, this, tStrips, pointsZ, tStripNormals, true);
	}

	/*
		public void draw(GL gl, float xMin, float xInc, float yMin, float yInc, StsColorscale colorscale)
		{
			float xInc2 = xInc / 2.f;
			float yInc2 = yInc / 2.f;
			float y = yMin + rowMin*yInc;
			for(int row = rowMin; row <= rowMax; row++, y += yInc)
			{
				float[] rowPointsZ1 = pointsZ[row - rowMin];
				float[] rowVals1 = null;
				boolean drawColors = false;
				if (values != null)
				{
					rowVals1 = values[row - rowMin];
					drawColors = true;
				}
				else
					StsColor.GREY.setGLColor(gl);

				gl.glBegin(GL.GL_QUADS);
				float x = xMin + colMin*xInc;
				for (int col = colMin; col <= colMax; col++, x += xInc)
				{
					float z1 = rowPointsZ1[col - colMin];
					if (z1 == StsParameters.nullValue) continue;
					if(drawColors)
					{
						float v = rowVals1[col - colMin];
						StsColor color = colorscale.getStsColor(colorscale.getIndexFromValue(v));
						color.setGLColor(gl);
					}
					gl.glVertex3f(x - xInc2, y - yInc2, z1);
					gl.glVertex3f(x - xInc2, y + yInc2, z1);
					gl.glVertex3f(x + xInc2, y + yInc2, z1);
					gl.glVertex3f(x + xInc2, y - yInc2, z1);
				}
				gl.glEnd();
			}
		}

		public void drawRowColorSurf_busted(GL gl, int row, int planeStartCol, int planeEndCol, float y, float xMin, float xInc, float yInc, float zInc, float z, float curvMin, float curvMax, StsColorscale colorscale)

		{
			int colStart = Math.max(colMin, planeStartCol);
			int colEnd = Math.min(colMax, planeEndCol);

			if (row < rowMin) return;
			if (row > rowMax - 1) return;

			if (zMin > z || zMax < z) return;

			float x = xMin + colStart * xInc;
			float[] rowPointsZ1 = pointsZ[row - rowMin];
			float[] rowPointsZ2 = pointsZ[1 + row - rowMin];
			float[] rowVals1 = null;
			float[] rowVals2 = null;
			if (values != null)
			{
				rowVals1 = values[row - rowMin];
				rowVals2 = values[1 + row - rowMin];
			}
			boolean begin = false;
			int count = 0;
			gl.glShadeModel(GL.GL_SMOOTH);
			gl.glDisable(GL.GL_CULL_FACE);
			gl.glDisable(GL.GL_POLYGON_STIPPLE);

			gl.glPointSize(1.f);
			gl.glColor4f(1.f, 1.f, 1.f, 1.f);
			if ((rowVals1 == null) || (rowVals2 == null))
			{
				return;
			}
			x = xMin + colStart * xInc;


			float v1[] = new float[3];
			float v2[] = new float[3];

			for (int col = colStart; col <= colEnd; col++, x += xInc)
			{

				float z1 = rowPointsZ1[col - colMin];
				float z2 = rowPointsZ2[col - colMin];

				{
					if ((z1 == StsParameters.nullValue) || (z2 == StsParameters.nullValue))
					{
						//System.out.println("end ");
						if (begin)
						{
							if (count == 2)
							{
								gl.glVertex3fv(v2, 0);
								gl.glVertex3fv(v1, 0);
							}

							gl.glEnd();
							begin = false;
							count = 0;
						}
					}
					else
					{
						if (!begin)
						{
							//gl.glBegin(GL.GL_QUAD_STRIP);
							gl.glBegin(GL.GL_LINES);
							begin = true;
							count = 0;
						}
						if (z2 == StsParameters.nullValue)
						{
							v1[0] = v2[0] = x;
							v1[1] = v2[1] = y;
							v1[2] = v2[2] = z1;
							gl.glVertex3fv(v2, 0);
							gl.glVertex3fv(v1, 0);
							count += 2;


						}
						else if (z1 == StsParameters.nullValue)
						{
							v1[0] = v2[0] = x;
							v1[1] = v2[1] = y + yInc;
							v1[2] = v2[2] = z2;
							gl.glVertex3fv(v2, 0);
							gl.glVertex3fv(v1, 0);
							count += 2;


						}
						else
						{

							v2[0] = x;
							v2[1] = y + yInc;
							v2[2] = z2;
							gl.glVertex3fv(v2, 0);

							v1[0] = x;
							v1[1] = y;
							v1[2] = z1;
							gl.glVertex3fv(v1, 0);

							count += 2;
							//System.out.println(" both "+count);
							//System.out.println("    "+v1[0]+" "+v1[1]+" "+v1[2]);
							//System.out.println("    "+v2[0]+" "+v2[1]+" "+v2[2]);
						}
					}

				}
				if (begin)
				{
					if (count == 2)
					{

						gl.glVertex3fv(v2, 0);
						gl.glVertex3fv(v1, 0);
					}

					gl.glEnd();
				}
			}
		}

		transient boolean ppChecked = false;
		transient boolean hasPP = false;

		boolean hasPointParams(GL gl)
		{
			if (ppChecked) return hasPP;

			if (gl.isExtensionAvailable("GL_ARB_point_parameters"))
				hasPP = true;

			if (gl.isExtensionAvailable("GL_EXT_point_parameters"))
				hasPP = true;

			ppChecked = true;

			return hasPP;
		}

		double sqr(double a) { return (a * a); }

		private void setPointParams(GL gl)
		{
			if (!hasPointParams(gl)) return;
			int error = gl.glGetError();
			float maxSize = 50.f;

			gl.glPointParameterf(GL.GL_POINT_SIZE_MIN, 2.0f);
			double[] projectionMatrix = new double[16];
			int[] viewport = new int[4];
			gl.glGetDoublev(GL.GL_PROJECTION_MATRIX, projectionMatrix, 0);
			gl.glGetIntegerv(GL.GL_VIEWPORT, viewport, 0);
			double H = viewport[2];
			double h = 2.0 / projectionMatrix[0];
			double D0 = Math.sqrt(2.0 * H / h);
			double k = 1.0 / (1.0 + 2 * sqr(1. / projectionMatrix[0]));
			float[] atten = new float[3];
			atten[0] = 1.f;
			atten[1] = 0.0f;
			k /= 500.f;
			atten[2] = (float) sqr(1 / D0) * (float) k;
			//System.out.println(" atten "+atten[2]);
			//atten[2] = 0.000001f;

			gl.glPointParameterfv(GL.GL_DISTANCE_ATTENUATION_EXT, atten, 0);
			error = gl.glGetError();
			if (error != 0)
			{
				GLU glu = new GLU();
				System.out.println("pointParams err code " + error + " " + glu.gluErrorString(error));
			}
			System.out.println("atten " + atten[2]);
			float[] ft = new float[1];
			gl.glGetFloatv(GL.GL_POINT_SIZE_MAX, ft, 0);
			System.out.println("window max " + ft[0]);
			gl.glPointSize(ft[0] > maxSize ? maxSize : ft[0]);
			gl.glPointParameterf(GL.GL_POINT_SIZE_MAX, ft[0] > maxSize ? maxSize : ft[0]);
		}
	*/
	public void drawRowVox(GL gl, float yMin, float yInc, float xMin, float xInc, StsColorscale colorscale)
	{
		float y = yMin + rowMin * yInc;
		for (int row = rowMin; row <= rowMax; row++, y += yInc)
		{
			float[] rowPointsZ = pointsZ[row - rowMin];

			if (values == null) continue;
			float[] rowVals = values[row - rowMin];
			if (rowVals == null) continue;
			gl.glPointSize(3.f);
			gl.glBegin(GL.GL_POINTS);
			float x = xMin + colMin * xInc;
			for (int col = colMin; col <= colMax; col++, x += xInc)
			{
				float z = rowPointsZ[col - colMin];
				Color color = colorscale.getColor(colorscale.getIndexFromValue(rowVals[col - colMin]));
				float colorsf[] = new float[4];
				color.getRGBComponents(colorsf);
				gl.glColor4fv(colorsf, 0);
				if (z != StsParameters.nullValue)
					gl.glVertex3f(x, y, z);
			}
			gl.glEnd();
		}
	}

	public int size()
	{
		return (rowMax - rowMin + 1) * (colMax - colMin + 1);
	}

	public String toString()
	{
		String rowColString = super.toString();
		int childGridID = (childGrid == null ? -1 : childGrid.id);
		int parentGridID = (parentGrid == null ? -1 : parentGrid.id);
		return "id: " + id + " childGrid ID: " + childGridID + " parentGrid ID: " + parentGridID + " originalID: " + originalID + " nPatchPoints " + nPatchPoints;
	}

	public String toGridString()
	{
		String rowColString = super.toString();
		int childGridID = (childGrid == null ? -1 : childGrid.id);
		int parentGridID = (parentGrid == null ? -1 : parentGrid.id);
		if (patchPoints != null)
			return "id: " + id + " childGrid ID: " + childGridID + " parentGrid ID: " + parentGridID + " originalID: " + originalID + " nPatchPoints " + nPatchPoints + " " + rowColString + " zMin: " + zMin + " zMax: " + zMax;
		else
			return "id: " + id + " childGrid ID: " + childGridID + " parentGrid ID: " + parentGridID + " originalID: " + originalID + " nPatchPoints " + nPatchPoints + " " + rowColString;
	}

	public String toFamilyString()
	{
		int childGridID = (childGrid == null ? -1 : childGrid.id);
		int parentGridID = (parentGrid == null ? -1 : parentGrid.id);
		return getFamilyTypeString() + " id: " + id + " child ID: " + childGridID + " parent ID: " + parentGridID + " original ID: " + originalID;
	}

	public String parentString()
	{
		if (parentGrid == null)
			return " null ";
		else
			return new String(" " + parentGrid.id + " ");
	}

	public String childString()
	{
		if (childGrid == null)
			return " null ";
		else
			return new String(" " + childGrid.id + " ");
	}

	public StsPatchGrid getParentGrid()
	{
		if (parentGrid != null) return parentGrid;
		else return this;
	}

	public int fillHistogram(float[] data, int nValues)
	{
		for (int row = 0; row < nRows; row++)
		{
			for (int col = 0; col < nCols; col++)
			{
				float value = values[row][col];
				if (value != nullValue)
					data[nValues++] = value;
				if (nValues == data.length)
					return nValues;
			}
		}
		return nValues;
	}

	public float getXMin()
	{
		return patchVolume.xMin;
	}

	public float getXMax()
	{
		return patchVolume.xMax;
	}

	public float getYMin()
	{
		return patchVolume.yMin;
	}

	public float getYMax()
	{
		return patchVolume.yMax;
	}

	public float getXInc()
	{
		return patchVolume.xInc;
	}

	public float getYInc()
	{
		return patchVolume.yInc;
	}

	public float getRowCoor(float[] xyz)
	{
		return patchVolume.getRowCoor(xyz);
	}

	public float getColCoor(float[] xyz)
	{
		return patchVolume.getRowCoor(xyz);
	}

	public double getXOrigin()
	{
		return patchVolume.xOrigin;
	}

	public double getYOrigin()
	{
		return patchVolume.yOrigin;
	}

	public float getXSize()
	{
		return patchVolume.getXSize();
	}

	public float getYSize()
	{
		return patchVolume.getYSize();
	}

	public float getAngle()
	{
		return patchVolume.getAngle();
	}

	public float getXCoor(float rowF, float colF)
	{
		return patchVolume.getXCoor(rowF, colF);
	}

	public float getYCoor(float rowF, float colF)
	{
		return patchVolume.getYCoor(rowF, colF);
	}

	public StsPoint getPoint(int volumeRow, int volumeCol)
	{
		float[] xyz = getXYZorT(volumeRow, volumeCol);
		return new StsPoint(xyz);
	}

	public float[] getXYZorT(int volumeRow, int volumeCol)
	{
		float[] xy = patchVolume.getXYCoors(volumeRow, volumeCol);

		float z = pointsZ[volumeRow - rowMin][volumeCol - colMin];
		return new float[]{xy[0], xy[1], z};
	}

	public StsPoint getPoint(float rowF, float colF)
	{
		return null;
	} // not used

	public float[] getXYZorT(float rowF, float colF)
	{
		return null;
	} // not used

	public float[] getNormal(int row, int col)
	{
		return null;
	} // not used

	public float[] getNormal(float rowF, float colF)
	{
		return null;
	} // not used

	public void checkConstructGridNormals()
	{
		return;
	} // not used

	public float getZInc()
	{
		return 0.0f;
	} // not used

	public float getZMin()
	{
		return zMin;
	}

	public float getZMax()
	{
		return zMax;
	}

	public String getLabel()
	{
		return toString();
	}

	public float interpolateBilinearZ(StsPoint point, boolean computeIfNull, boolean setPoint)
	{
		return 0.0f;
	} // not used:  yet

	public float interpolateBilinearZ(StsGridPoint gridPoint, boolean computeIfNull, boolean setPoint)
	{
		return 0.0f;
	} // not used: yet

	public float getComputePointZ(int row, int col)
	{
		return pointsZ[row][col];
	}

	public boolean toggleSurfacePickingOn()
	{
		return false;
	}

	public void toggleSurfacePickingOff()
	{
	}

	public String getName()
	{
		return "patchGrid-" + originalID;
	}

	public StsGridPoint getSurfacePosition(StsMouse mouse, boolean display, StsGLPanel3d glPanel3d)
	{
		return null;
	}

	public void setIsVisible(boolean isVisible)
	{
		this.isVisible = isVisible;
	}

	public boolean getIsVisible()
	{
		return isVisible;
	}

	/**
	 * For this row link along the bottom of the cell whose lower-left corner is at row & col:
	 * determine the type of link to draw: NONE, ABOVE, BELOW, BOTH, or LINE
	 * if no link, return NONE;
	 * If a column link exists above but not below: return LINK_ABOVE
	 * If a column link exists below but not above: return LINK_BELOW
	 * If at least one above and one below: return LINK_BOTH
	 * If link exists, but NO col links above or below: return LINK_LINE
	 * @param row patch row of lower left of cell
	 * @param col patch col of lower rite of cell
	 * @return link type
	 */

	public byte computeRowLink(int row, int col)
	{
		if (!hasRowLink(row, col)) return LINK_NONE;
		boolean hasLinkAbove = hasColLink(row, col) || hasColLink(row, col + 1);
		boolean hasLinkBelow = hasColLink(row - 1, col) || hasColLink(row - 1, col + 1);
		if (hasLinkAbove)
		{
			if (hasLinkBelow)
				return LINK_BOTH;
			else
				return LINK_ABOVE;
		}
		else // !hasLinkAbove
		{
			if (hasLinkBelow)
				return LINK_BELOW;
			else // row link but no col links above or below
				return LINK_LINE;
		}
	}

	/**
	 * For this row link along the bottom of the cell whose lower-left corner is at row & col:
	 * determine the type of link to draw: NONE, ABOVE, BELOW, BOTH, or LINE
	 * if no link, return NONE;
	 * If a column link exists above but not below: return LINK_ABOVE
	 * If a column link exists below but not above: return LINK_BELOW
	 * If at least one above and one below: return LINK_BOTH
	 * If link exists, but NO col links above or below: return LINK_LINE
	 * @param row patch row of lower left of cell
	 * @param col patch col of lower rite of cell
	 * @return link type
	 */

	public byte computeColLink(int row, int col)
	{
		if (!hasColLink(row, col)) return LINK_NONE;
		boolean hasLinkLeft = hasRowLink(row, col - 1) || hasRowLink(row + 1, col - 1);
		boolean hasLinkRite = hasRowLink(row, col) || hasRowLink(row + 1, col);
		if (hasLinkLeft)
		{
			if (hasLinkRite)
				return LINK_BOTH;
			else
				return LINK_LEFT;
		}
		else // !hasLinkAbove
		{
			if (hasLinkRite)
				return LINK_RIGHT;
			else // row link but no col links above or below
				return LINK_LINE;
		}
	}

	public PatchPoint getVolPatchPoint(PatchPoint patchPoint)
	{
		return getVolPatchPoint(patchPoint.getRow(), patchPoint.getCol());
	}

	/**
	 * Given global row,col get patchPoint on grid at this point
	 * @param volRow global row
	 * @param volCol global col
	 * @return
	 */
	public PatchPoint getVolPatchPoint(int volRow, int volCol)
	{
		if (patchPoints != null && isInsideRowCol(volRow, volCol))
			return patchPoints[volRow - rowMin][volCol - colMin];
		return null;
	}

	/**
	 * Given local row,col get patchPoint on grid at this point
	 * @param row local row
	 * @param col local col
	 * @return
	 */
	public PatchPoint getGridPatchPoint(int row, int col)
	{
		if (patchPoints == null) return null;
		if (!isInsideGridRowCol(row, col)) return null;
		return patchPoints[row][col];
	}

	public float getVolumePointZ(int volumeRow, int volumeCol)
	{
		return getPatchPointZ(volumeRow - rowMin, volumeCol - colMin);
	}

	public float getPatchPointZ(int patchRow, int patchCol)
	{
		if (pointsZ == null)
			return nullValue;
		if (!isInsidePatchRowCol(patchRow, patchCol))
			return nullValue;
		return pointsZ[patchRow][patchCol];
	}

	public boolean isInsidePatchRowCol(int row, int col)
	{
		return row >= 0 && row < nRows && col >= 0 && col < nCols;
	}

	public boolean hasRowLink(int patchRow, int patchCol)
	{
		if (rowCorrels == null) return false;
		if (!isInsidePatchRowCol(patchRow, patchCol)) return false;
		return rowCorrels[patchRow][patchCol] > patchVolume.minLinkCorrel;
	}

	public boolean hasColLink(int patchRow, int patchCol)
	{
		if (colCorrels == null) return false;
		if (!isInsidePatchRowCol(patchRow, patchCol)) return false;
		return colCorrels[patchRow][patchCol] > patchVolume.minLinkCorrel;
	}

	public float getCurvature(int volumeRow, int volumeCol, float dataMin, float dataMax)
	{
		if (!isInsideRowCol(volumeRow, volumeCol))
		{
			StsException.systemError(this, "getValue", "volumeRow or volumeCol not inside patch");
			return 0.0f;
		}
		if (values == null) return 0.0f;
		float value = values[volumeRow - rowMin][volumeCol - colMin];
		if (value == nullValue) return 0.0f;
		if (value < dataMin) return dataMin;
		else if (value > dataMax) return dataMax;
		else return value;
	}
}
