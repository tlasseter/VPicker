package com.Sts.SeismicAttributes;

import java.util.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: Grid status flags have been encapsulated into this class to make management easier.
 * Flags are used in multiple places with the same logic, so this reduces code. It is possible to link multiple
 * instances of this class together. For example...</p>
 * <p>StsGridStatus3d parent = StsGridStatus3d.construct3dGrid(100, 100, 50);<br>
 * StsGridStatus3d child = parent.createSuccessor();</p>
 * <p>Setting a location as valid in the parent obeject will only set the flag there.
 * Setting a location as invalid in the parent will set the flag as false there and in the child object locations.</p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

abstract public class StsGridStatus
{
	static final public boolean asserts = false;
	protected StsVolumeStatusFace volumeStatusFace;
	private StsGridStatus predecessorGridStatus = null;
//	private List<StsGridStatus> successorGridStatus = new LinkedList<StsGridStatus>();
	private int nRows;
	private int nLocations;

	private int locationOKCount = 0;
	private int[] rowOKCount;
	private boolean[][] locationOK;

	abstract public StsGridStatus createSuccessor(StsVolumeStatusFace volumeStatusFace);
	abstract public boolean isColumnOK(int col);
	abstract public void setColumnOK(int col);
	abstract public void setColumnInvalid(int col);
	abstract public boolean isSliceOK(int slice);
	abstract public void setSliceOK(int slice);
	abstract public void setSliceInvalid(int slice);

	public static StsGridStatus construct3dGrid(int nRows, int nCols, int nSlices)
	{
		return construct3dGrid(null, nRows, nCols, nSlices);
	}

	public static StsGridStatus construct2dLineGrid(int nRows)
	{
		return construct2dLineGrid(null, nRows);
	}

	public static StsGridStatus construct3dGrid(StsVolumeStatusFace volumeStatusFace, int nRows, int nCols, int nSlices)
	{
		if (volumeStatusFace == null)
		{
			volumeStatusFace = new StsVolumeStatusFace()
			{
				public void setDirty()
				{
				}
			};
		}
		return new StsGridStatus3d(volumeStatusFace, nRows, nCols, nSlices);
	}

	public static StsGridStatus construct2dLineGrid(StsVolumeStatusFace volumeStatusFace, int nRows)
	{
		if (volumeStatusFace == null)
		{
			volumeStatusFace = new StsVolumeStatusFace()
			{
				public void setDirty()
				{
				}
			};
		}
		return new StsGridStatus2d(volumeStatusFace, nRows);
	}

	protected StsGridStatus(StsVolumeStatusFace volumeStatusFace, int nRows, int nCols)
	{
		this.volumeStatusFace = volumeStatusFace;
		this.nRows = nRows;
		nLocations = nRows * nCols;
		rowOKCount = new int[nRows];
		locationOK = new boolean[nRows][nCols];
	}

	protected StsGridStatus(StsVolumeStatusFace volumeStatusFace, int nRows)
	{
		this.volumeStatusFace = volumeStatusFace;
		this.nRows = nRows;
		rowOKCount = new int[nRows];
		locationOK = new boolean[nRows][];
	}

	// use this constructor to create a grid based on another grid and have them linked up.
	protected StsGridStatus(StsVolumeStatusFace volumeStatusFace, StsGridStatus predecessorGridStatus)
	{
		this.volumeStatusFace = volumeStatusFace;
		nRows = predecessorGridStatus.nRows;
		nLocations = predecessorGridStatus.nLocations;
		rowOKCount = new int[nRows];
		locationOK = new boolean[nRows][];
		for (int row = 0; row < nRows; row++)
		{
			int nCols = predecessorGridStatus.locationOK[row].length;
			locationOK[row] = new boolean[nCols];
		}

		// link up predecessor and successor
    /*
        this.predecessorGridStatus = predecessorGridStatus;
		synchronized (predecessorGridStatus.successorGridStatus)
		{
			predecessorGridStatus.successorGridStatus.add(this);
		}
    */
	}

	public boolean isGridOK()
	{
		return locationOKCount == nLocations;
	}

	public boolean isRowOK(int row)
	{
		return rowOKCount[row] == locationOK[row].length;
	}

	public void setRowOK(int row)
	{
		int maxCol = locationOK[row].length;
		for (int col = 0; col < maxCol; col++)
		{
			setLocationOK(row, col);
		}
	}

	public void setRowInvalid(int row)
	{
		volumeStatusFace.setDirty();
		int maxCol = locationOK[row].length;
		for (int col = 0; col < maxCol; col++)
		{
			setLocationInvalid(row, col);
		}
	}

	public boolean isLocationOK(int row, int col)
	{
		return locationOK[row][col];
	}

	public void setLocation(int row, int col, boolean ok)
	{
		if (ok)
		{
			setLocationOK(row, col);
		}
		else
		{
			setLocationInvalid(row, col);
		}
	}

	public void setGrid(boolean ok)
	{
		if (ok)
		{
			setGridOK();
		}
		else
		{
			setGridInvalid();
		}
	}

	public void setGridOK()
	{
		locationOKCount = nLocations;
		for (int row = 0; row < nRows; row++)
		{
			Arrays.fill(locationOK[row], true);
			rowOKCount[row] = locationOK[row].length;
		}
	}

	public void setGridInvalid()
	{
		volumeStatusFace.setDirty();
		Arrays.fill(rowOKCount, 0);
		locationOKCount = 0;
		for (int row = 0; row < nRows; row++)
		{
			Arrays.fill(locationOK[row], false);
		}
	}

	public boolean setLocationOK(int row, int col)
	{
		if (asserts)
		{
			assert row < nRows && col < locationOK[row].length;
		}
		if (locationOK[row][col])
		{
			return false;
		}

		if (asserts && predecessorGridStatus != null)
		{
			assert predecessorGridStatus.isLocationOK(row, col);
		}

		locationOK[row][col] = true;
		rowOKCount[row]++;
		locationOKCount++;

		return true;
	}

	public boolean setLocationInvalid(int row, int col)
	{
		if (! locationOK[row][col])
		{
			return false;
		}
		volumeStatusFace.setDirty();
		locationOK[row][col] = false;
		rowOKCount[row]--;
		locationOKCount--;
/*
		for (StsGridStatus gridStatus : successorGridStatus)
		{
			gridStatus.setLocationInvalid(row, col);
		}
*/
		return true;
	}

	public synchronized void initialiseRow(int row, int nCols)
	{
		if (asserts)
		{
			assert row < nRows;
		}
		if (locationOK[row] != null)
		{
			nLocations -= locationOK[row].length;
		}
		locationOK[row] = new boolean[nCols];
		nLocations += nCols;
	}
}

class StsGridStatus2d extends StsGridStatus
{
	private StsGridStatus2d(StsVolumeStatusFace volumeStatusFace, StsGridStatus2d predecessorGridStatus)
	{
		super(volumeStatusFace, predecessorGridStatus);
	}

	protected StsGridStatus2d(StsVolumeStatusFace volumeStatusFace, int nRows)
	{
		super(volumeStatusFace, nRows);
	}

	// use this to create a grid based on another grid and have them linked up.
	public StsGridStatus createSuccessor(StsVolumeStatusFace volumeStatusFace)
	{
		return new StsGridStatus2d(volumeStatusFace, this);
	}

	public boolean isColumnOK(int col)
	{
		throw new RuntimeException("isColumnOK not implemented for 2d status grid");
	}

	public void setColumnOK(int col)
	{
		throw new RuntimeException("setColumnOK not implemented for 2d status grid");
	}

	public void setColumnInvalid(int col)
	{
		throw new RuntimeException("setColumnInvalid not implemented for 2d status grid");
	}

	public boolean isSliceOK(int slice)
	{
		throw new RuntimeException("isSliceOK not implemented for 2d status grid");
	}

	public void setSliceOK(int slice)
	{
		throw new RuntimeException("setSliceOK not implemented for 2d status grid");
	}

	public void setSliceInvalid(int slice)
	{
		throw new RuntimeException("setSliceInvalid not implemented for 2d status grid");
	}
}

class StsGridStatus3d extends StsGridStatus
{
	private int nRows;
	private int nCols;
	private int nSlices;
	private int[] colOKCount;
	private boolean[] sliceOK;
	private int slicesOKCount = 0;

	private StsGridStatus3d(StsVolumeStatusFace volumeStatusFace, StsGridStatus3d predecessorGridStatus)
	{
		super(volumeStatusFace, predecessorGridStatus);
		this.nRows = predecessorGridStatus.nRows;
		this.nCols = predecessorGridStatus.nCols;
		this.nSlices = predecessorGridStatus.nSlices;
		colOKCount = new int[nCols];
		sliceOK = new boolean[nSlices];
	}

	protected StsGridStatus3d(StsVolumeStatusFace volumeStatusFace, int nRows, int nCols, int nSlices)
	{
		super(volumeStatusFace, nRows, nCols);
		this.nRows = nRows;
		this.nCols = nCols;
		this.nSlices = nSlices;
		colOKCount = new int[nCols];
		sliceOK = new boolean[nSlices];
	}

	// use this to create a grid based on another grid and have them linked up.
	public StsGridStatus createSuccessor(StsVolumeStatusFace volumeStatusFace)
	{
		return new StsGridStatus3d(volumeStatusFace, this);
	}

	public boolean setLocationOK(int row, int col)
	{
		if (super.setLocationOK(row, col))
		{
			colOKCount[col]++;
			return true;
		}
		return false;
	}

	public boolean setLocationInvalid(int row, int col)
	{
		if (super.setLocationInvalid(row, col))
		{
			colOKCount[col]--;
			return true;
		}
		slicesOKCount = 0;
		Arrays.fill(sliceOK, false);
		return false;
	}

	public boolean isColumnOK(int col)
	{
		return colOKCount[col] == nRows;
	}

	public boolean isGridOK()
	{
		return super.isGridOK() && this.nSlices == slicesOKCount;
	}

	public void setGridOK()
	{
		super.setGridOK();
		slicesOKCount = nSlices;
		for (int col = 0; col < nCols; col++)
		{
			colOKCount[col] = nRows;
		}
		Arrays.fill(sliceOK, true);
	}

	public void setGridInvalid()
	{
		super.setGridInvalid();
		slicesOKCount = 0;
		Arrays.fill(colOKCount, 0);
		Arrays.fill(sliceOK, false);
	}

	public void setColumnOK(int col)
	{
		for (int row = 0; row < nRows; row++)
		{
			setLocationOK(row, col);
		}
	}

	public void setColumnInvalid(int col)
	{
		volumeStatusFace.setDirty();
		for (int row = 0; row < nRows; row++)
		{
			setLocationInvalid(row, col);
		}
		slicesOKCount = 0;
		Arrays.fill(sliceOK, false);
	}

	public void initialiseRow(int row, int nCols)
	{
		throw new RuntimeException("initialiseRow not implemented for 3d status grid");
	}

	public boolean isSliceOK(int slice)
	{
		return sliceOK[slice];
	}

	public void setSliceOK(int slice)
	{
		slicesOKCount++;
		sliceOK[slice] = true;
	}

	public void setSliceInvalid(int slice)
	{
		slicesOKCount--;
		sliceOK[slice] = false;
		volumeStatusFace.setDirty();
	}
}

