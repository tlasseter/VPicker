package com.Sts.Interfaces;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002-2003</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.1
 */

import com.Sts.DBTypes.*;
import com.Sts.Types.*;

import javax.media.opengl.*;
import java.awt.*;

/** Interface for color display of surface thru a 3D grid of properties.
 *  The surface and the properties data are on the same XY (IJ) grid.
 *  Grid can be stratigraphic, i.e., the grid shape follows the surfaces or
 *  grid can be orthogonal in which case we must interpolate in z direction.
 */

public interface StsSurfaceDisplayable
{
	/** grid is stratigraphic and values are cell-centered */
	static final byte STRAT_CELL_CENTERED = 1;
	/** Grid is the same, but we want an approx draw to speed it up. */
	static final byte STRAT_APPROX_CELL_CENTERED = 2;
	/** Grid is XYZ orthogonal grid and evenly-spaced in z. */
	static final byte ORTHO_GRID_CENTERED = 3;

	/** Gets the display type (Cell or Grid); */
	public byte getDisplayType();
	/** This property can be isVisible on this surace. */
	public boolean isDisplayable();
	/** Call for getting a cell-centered color */
	public Color getCellColor(int row, int col, int layer);
	/** Call for getting a grid-centered color */
	public Color getGridColor(int row, int col, float z);
    /** returns colorscale for this surface displayable property */
    public StsColorscale getColorscale();
    /** return a 2D {row][col] array of Java colors */
    public Color[][] get2dColorArray(StsRotatedGridBoundingBox surfaceBoundingBox, float[][] z, float zOffset);
    /** return a 2D {row][col] array of unsigned bytes */
    public byte[] getByteArray(StsRotatedGridBoundingBox surfaceBoundingBox, float[][] z, float zOffset);
    /** returns a boolean indicating texture should use pixel (NEAREST) mode */
    public boolean getIsPixelMode();
    /** get colors GL displayList */
    public int getColorListNum(GL gl);
}
