
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Utilities;

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.Types.*;

import javax.media.opengl.*;

/** Not current used. */

public class TriangleStrip
{
    public int rowNumber = -1;
    public int firstIndex = -1;
    public int lastIndex = -1;
    public int firstSide = StsParameters.STRIP_INVALID;
    public int lastSide = StsParameters.STRIP_INVALID;
    public boolean isGap;

	// Convenience copies of useful flags
	static public final byte SURF_PNT = StsParameters.SURF_PNT;
	static public final byte SURF_BOUNDARY = StsParameters.SURF_BOUNDARY;
	static public final byte SURF_GAP = StsParameters.SURF_GAP;
	static public final byte SURF_GAP_SET = StsParameters.SURF_GAP_SET;
	static public final byte SURF_GAP_FILLED = StsParameters.SURF_GAP_FILLED;
//	static public final byte NULL_GAP_OR_BOUNDARY = StsParameters.NULL_GAP_OR_BOUNDARY;
	static public final byte SURF_GAP_NOT_FILLED = StsParameters.SURF_GAP_NOT_FILLED;

	/** Convenience copies of TStrip flags */
	static public final byte STRIP_INVALID = StsParameters.STRIP_INVALID;
	static public final byte STRIP_BOTH = StsParameters.STRIP_BOTH;
	static public final byte STRIP_BOT = StsParameters.STRIP_BOT;
	static public final byte STRIP_TOP = StsParameters.STRIP_TOP;

	static public final int nullInteger = StsParameters.nullInteger;

    public TriangleStrip() { }

    public boolean fieldsAreValid()
    {
        if (rowNumber==-1) return false;
        if (firstIndex==-1) return false;
        if (lastIndex==-1) return false;
        if (firstSide==StsParameters.STRIP_INVALID) return false;
        if (lastSide==StsParameters.STRIP_INVALID) return false;
        return true;
    }

	static public StsList computeTStrips(StsXYSurfaceGridable grid, byte[][] pointsNull)
	{
		int row, col;
		int colMin, colMax, nextColMin, nextColMax;
		int nBot, nTop;
		byte pointNull;
		boolean botNotNull, topNotNull;
		TriangleStrip t = null;
		boolean startedTStrip = false;
		int nPoints = 0;

		StsGridIterator gridIterator = new StsGridIterator(grid);
		int nGridRows = gridIterator.getNRows();
        StsList tStrips = new StsList(nGridRows, nGridRows);
		while((row = gridIterator.getNextRow()) != nullInteger)
		{
			colMin = gridIterator.colMin;
			colMax = gridIterator.colMax;

			nextColMin = gridIterator.nextColMin;
			nextColMax = gridIterator.nextColMax;

			int colStart = Math.min(colMin, nextColMin);
			int colEnd = Math.max(colMax, nextColMax);

			for (col=colStart; col <= colEnd; col++)
			{
				pointNull = pointsNull[row][col];
				botNotNull = (pointNull == SURF_PNT || pointNull == SURF_GAP_FILLED || pointNull == SURF_GAP || pointNull == SURF_GAP_NOT_FILLED);
//				botNotNull = (pointNull == NOT_NULL || pointNull == NULL_GAP_FILLED);

				pointNull = pointsNull[row+1][col];
				topNotNull = (pointNull == SURF_PNT || pointNull == SURF_GAP_FILLED || pointNull == SURF_GAP || pointNull == SURF_GAP_NOT_FILLED);
//				topNotNull = (pointNull == NOT_NULL || pointNull == NULL_GAP_FILLED);

				if (!startedTStrip)
				{
					nPoints = 0;
					if (botNotNull && topNotNull)
					{
						t = new TriangleStrip();  // create a new tstrip
						t.firstSide = STRIP_BOTH;
						t.firstIndex = col;
						t.lastIndex = col;
						nPoints += 2;
						startedTStrip = true;
					}
					else if (botNotNull) // top point null
					{
						t = new TriangleStrip();  // create a new tstrip
						t.firstSide = STRIP_BOT;
						t.firstIndex = col+1;
						nPoints++;
						startedTStrip = true;
					}
					else if (topNotNull) // bot point null
					{
						t = new TriangleStrip();  // create a new tstrip
						t.firstSide = STRIP_TOP;
						t.firstIndex = col+1;
						nPoints++;
						startedTStrip = true;
					}
				}  // !startedTStrip

				else  // tstrip started
				{
					if (!botNotNull && !topNotNull)  // both null - end tstrip
					{
						if (nPoints > 2)  // need 3 or more points
						{
							// add the triangle strip
							t.rowNumber = row;
							tStrips.add(t);
						}
						startedTStrip = false;  // reset flag
					}
					else if (botNotNull && topNotNull)  // continue tstrip
					{
						t.lastIndex = col;
						t.lastSide = STRIP_BOTH;
						nPoints += 2;
					}
					else if (botNotNull) // top point null  - end this tstrip
					{
						nPoints++;
						if (nPoints > 2)
						{
							t.lastSide = STRIP_BOT;
							t.rowNumber = row;
							tStrips.add(t);
						}
						t = new TriangleStrip();  // start a new tstrip
						t.firstSide = STRIP_BOT;
						t.firstIndex = col+1;
						nPoints = 1;
					}
					else // (topNotNull) top point not null - end this tstrip
					{
						nPoints++;
						if (nPoints > 2)
						{
							t.lastSide = STRIP_TOP;
							t.rowNumber = row;
							tStrips.add(t);
						}
						t = new TriangleStrip();  // start a new tstrip
						t.firstSide = STRIP_TOP;
						t.firstIndex = col+1;
						nPoints = 1;
					}
				}  // tstrip started
			}  // for "col" loop

			// reached end of a row
			if (startedTStrip && nPoints>2)
			{
				t.rowNumber = row;
				tStrips.add(t);
			}
			startedTStrip = false;  // reset flag
		}  // for "i" loop

		tStrips.trimToSize();  // shrink to actual allocation
        return tStrips;
	}

    static public StsList computeTStrips(StsXYSurfaceGridable grid)
	{
		int row, col;
		int colMin, colMax, nextColMin, nextColMax;
		int nBot, nTop;
		boolean botNotNull, topNotNull;
		TriangleStrip t = null;
		boolean startedTStrip = false;
		int nPoints = 0;

		StsGridIterator gridIterator = new StsGridIterator(grid);
		int nGridRows = gridIterator.getNRows();
		StsList tStrips = new StsList(nGridRows, nGridRows);
        float[][] pointsZ = grid.getPointsZ();
		while((row = gridIterator.getNextRow()) != nullInteger)
		{
			colMin = gridIterator.colMin;
			colMax = gridIterator.colMax;

			nextColMin = gridIterator.nextColMin;
			nextColMax = gridIterator.nextColMax;

			int colStart = Math.min(colMin, nextColMin);
			int colEnd = Math.max(colMax, nextColMax);

			for (col=colStart; col <= colEnd; col++)
			{
				botNotNull = pointsZ[row][col] != StsParameters.nullValue;
				topNotNull = pointsZ[row+1][col] != StsParameters.nullValue;

				if (!startedTStrip)
				{
					nPoints = 0;
					if (botNotNull && topNotNull)
					{
						t = new TriangleStrip();  // create a new tstrip
						t.firstSide = STRIP_BOTH;
						t.firstIndex = col;
						t.lastIndex = col;
						nPoints += 2;
						startedTStrip = true;
					}
					else if (botNotNull) // top point null
					{
						t = new TriangleStrip();  // create a new tstrip
						t.firstSide = STRIP_BOT;
						t.firstIndex = col+1;
						nPoints++;
						startedTStrip = true;
					}
					else if (topNotNull) // bot point null
					{
						t = new TriangleStrip();  // create a new tstrip
						t.firstSide = STRIP_TOP;
						t.firstIndex = col+1;
						nPoints++;
						startedTStrip = true;
					}
				}  // !startedTStrip

				else  // tstrip started
				{
					if (!botNotNull && !topNotNull)  // both null - end tstrip
					{
						if (nPoints > 2)  // need 3 or more points
						{
							// add the triangle strip
							t.rowNumber = row;
							tStrips.add(t);
						}
						startedTStrip = false;  // reset flag
					}
					else if (botNotNull && topNotNull)  // continue tstrip
					{
						t.lastIndex = col;
						t.lastSide = STRIP_BOTH;
						nPoints += 2;
					}
					else if (botNotNull) // top point null  - end this tstrip
					{
						nPoints++;
						if (nPoints > 2)
						{
							t.lastSide = STRIP_BOT;
							t.rowNumber = row;
							tStrips.add(t);
						}
						t = new TriangleStrip();  // start a new tstrip
						t.firstSide = STRIP_BOT;
						t.firstIndex = col+1;
						nPoints = 1;
					}
					else // (topNotNull) top point not null - end this tstrip
					{
						nPoints++;
						if (nPoints > 2)
						{
							t.lastSide = STRIP_TOP;
							t.rowNumber = row;
							tStrips.add(t);
						}
						t = new TriangleStrip();  // start a new tstrip
						t.firstSide = STRIP_TOP;
						t.firstIndex = col+1;
						nPoints = 1;
					}
				}  // tstrip started
			}  // for "col" loop

			// reached end of a row
			if (startedTStrip && nPoints>2)
			{
				t.rowNumber = row;
				tStrips.add(t);
			}
			startedTStrip = false;  // reset flag
		}  // for "i" loop

		tStrips.trimToSize();  // shrink to actual allocation
        return tStrips;
	}

    static public StsList computeTStrips(StsXYSurfaceLinkGridable grid)
	{
		int row, col;
		int colMin, colMax, nextColMin, nextColMax;
		boolean botNotNull, topNotNull;
		TriangleStrip t = null;
		boolean startedTStrip = false;
		int nPoints = 0;

        try
        {
            StsGridIterator gridIterator = new StsGridIterator(grid);
            int nGridRows = gridIterator.getNRows();
            StsList tStrips = new StsList(nGridRows, nGridRows);
            float[][] pointsZ = grid.getPointsZ();
            while((row = gridIterator.getNextRow()) != nullInteger)
            {
                colMin = gridIterator.colMin;
                colMax = gridIterator.colMax;

                nextColMin = gridIterator.nextColMin;
                nextColMax = gridIterator.nextColMax;

                int colStart = Math.min(colMin, nextColMin);
                int colEnd = Math.max(colMax, nextColMax);

                for (col=colStart; col <= colEnd; col++)
                {
                    boolean hasRowLink = grid.hasRowLink(row, col);
                    boolean hasColLink = grid.hasColLink(row, col);
                    botNotNull = pointZNotNull(pointsZ, row, col) && (hasRowLink || hasColLink);
                    topNotNull = pointZNotNull(pointsZ, row+1, col) && hasColLink;

                    if (!startedTStrip)
                    {
                        nPoints = 0;
                        if (botNotNull && topNotNull)
                        {
                            t = new TriangleStrip();  // create a new tstrip
                            t.firstSide = STRIP_BOTH;
                            t.firstIndex = col;
                            t.lastIndex = col;
                            nPoints += 2;
                            startedTStrip = true;
                        }
                        else if (botNotNull) // top point null
                        {
                            t = new TriangleStrip();  // create a new tstrip
                            t.firstSide = STRIP_BOT;
                            t.firstIndex = col+1;
                            nPoints++;
                            startedTStrip = true;
                        }
                        else if (topNotNull) // bot point null
                        {
                            t = new TriangleStrip();  // create a new tstrip
                            t.firstSide = STRIP_TOP;
                            t.firstIndex = col+1;
                            nPoints++;
                            startedTStrip = true;
                        }
                    }  // !startedTStrip

                    else  // tstrip started
                    {
                        if (!botNotNull && !topNotNull)  // both null - end tstrip
                        {
                            if (nPoints > 2)  // need 3 or more points
                            {
                                // add the triangle strip
                                t.rowNumber = row;
                                tStrips.add(t);
                            }
                            startedTStrip = false;  // reset flag
                        }
                        else if (botNotNull && topNotNull)  // continue tstrip
                        {
                            t.lastIndex = col;
                            t.lastSide = STRIP_BOTH;
                            nPoints += 2;
                        }
                        else if (botNotNull) // top point null  - end this tstrip
                        {
                            nPoints++;
                            if (nPoints > 2)
                            {
                                t.lastSide = STRIP_BOT;
                                t.rowNumber = row;
                                tStrips.add(t);
                            }
                            t = new TriangleStrip();  // start a new tstrip
                            t.firstSide = STRIP_BOT;
                            t.firstIndex = col+1;
                            nPoints = 1;
                        }
                        else // (topNotNull) top point not null - end this tstrip
                        {
                            nPoints++;
                            if (nPoints > 2)
                            {
                                t.lastSide = STRIP_TOP;
                                t.rowNumber = row;
                                tStrips.add(t);
                            }
                            t = new TriangleStrip();  // start a new tstrip
                            t.firstSide = STRIP_TOP;
                            t.firstIndex = col+1;
                            nPoints = 1;
                        }
                    }  // tstrip started
                }  // for "col" loop

                // reached end of a row
                if (startedTStrip && nPoints>2)
                {
                    t.rowNumber = row;
                    tStrips.add(t);
                }
                startedTStrip = false;  // reset flag
            }  // for "i" loop

            tStrips.trimToSize();  // shrink to actual allocation
            return tStrips;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(TriangleStrip.class, "computeTStrips", e);
            return null;
        }
    }

    static private boolean pointZNotNull(float[][] pointsZ, int row, int col)
    {
        if(pointsZ == null) return false;
        if(pointsZ[row] == null) return false;
        return pointsZ[row][col] != StsParameters.nullValue;
    }

   /** construct normals from tstrips */
	static public float[][][] computeTStripNormals(StsXYSurfaceGridable grid, StsList tStrips, float[][] pointsZ)
	{
		int i, j, n;

        int nRows = grid.getNRows();
        int nCols = grid.getNCols();
        float xInc = grid.getXInc();
        float yInc = grid.getYInc();

        float[][][] normals = new float[nRows][nCols][3];

		float[] idif = new float[3];
		float[] jdif = new float[3];
		float[] normal;

		idif[0] = idif[1] = yInc;
		jdif[0] = xInc;

		int nTStrips = tStrips.getSize();
		for (n = 0; n < nTStrips; n++)
		{
			TriangleStrip t = (TriangleStrip)tStrips.getElement(n);
			i = t.rowNumber;
			j = t.firstIndex;

			// first side
			if (t.firstSide == STRIP_TOP)
			{
				idif[2] = pointsZ[i+1][j] - pointsZ[i][j];
				jdif[2] = pointsZ[i+1][j] - pointsZ[i+1][j-1];
				normals[i+1][j-1] = StsMath.leftCrossProduct(idif, jdif);
		   }
			else if (t.firstSide == STRIP_BOT)
			{
				idif[2] = pointsZ[i+1][j] - pointsZ[i][j];
				jdif[2] = pointsZ[i][j] - pointsZ[i][j-1];
				normals[i][j-1] = StsMath.leftCrossProduct(idif, jdif);
			}

			// middle sides excluding last middle side
			for (; j<t.lastIndex; j++)
			{
				idif[2] = pointsZ[i+1][j] - pointsZ[i][j];

				jdif[2] = pointsZ[i+1][j+1] - pointsZ[i+1][j];
				normals[i+1][j] = StsMath.leftCrossProduct(idif, jdif);

				jdif[2] = pointsZ[i][j+1] - pointsZ[i][j];
				normals[i][j] = StsMath.leftCrossProduct(idif, jdif);
			}
			// j has been bumped to next j by for loop above
			// last side: last idif is from bottom row if available; otherwise top row

			idif[2] = pointsZ[i+1][j] - pointsZ[i][j];
			normal = StsMath.leftCrossProduct(idif, jdif);
			normals[i][j] = normal;
			normals[i+1][j] = normal;

			if (t.lastSide == STRIP_BOT)
			{
				jdif[2] = pointsZ[i][j+1] -  pointsZ[i][j];
				normals[i][j+1] = StsMath.leftCrossProduct(idif, jdif);
			}
			else if (t.lastSide == STRIP_TOP)
			{
				jdif[2] = pointsZ[i+1][j+1] - pointsZ[i+1][j];
				normals[i+1][j+1] = StsMath.leftCrossProduct(idif, jdif);
			}
		}
        return normals;
	}

	/* display surface containing null Z values by drawing triangle strips */
	static public void drawSurfaceFillWithNulls(GL gl, StsXYSurfaceGridable grid, StsList tStrips, float[][] pointsZ, float[][][] normals, boolean isLighted)
	{
		TriangleStrip t = null;
		float[] point = new float[3];
		int i, j = 0;
		int n;


		float startX = grid.getXMin();
		float startY = grid.getYMin();
        float xInc = grid.getXInc();
        float yInc = grid.getYInc();
        int rowMin = grid.getRowMin();
        int colMin = grid.getColMin();
        try
		{
//			StsColor.setGLColor(gl, this.stsColor);

			int nTStrips = tStrips.getSize();
			for (int nt = 0; nt < nTStrips; nt++)
			{
				t = (TriangleStrip)tStrips.getElement(nt);
				i = t.rowNumber;
				j = t.firstIndex;

				gl.glBegin(GL.GL_TRIANGLE_STRIP);

				// first side
				if (t.firstSide == STRIP_BOT)
				{
					point[0] = startX + (j-1+colMin)*xInc;
					point[1] = startY + (i+rowMin)*yInc;
					point[2] = pointsZ[i][j-1];

					if(isLighted) gl.glNormal3fv(normals[i][j-1], 0);
					gl.glVertex3fv(point, 0);
					if(isLighted) gl.glNormal3fv(normals[i][j-1], 0);
					gl.glVertex3fv(point, 0);

					point[0] += xInc;
				}
				else if (t.firstSide == STRIP_TOP)
				{
					point[0] = startX + (j-1+colMin)*xInc;
					point[1] = startY + (i+1+rowMin)*yInc;
					point[2] = pointsZ[i+1][j-1];

					if(isLighted) gl.glNormal3fv(normals[i+1][j-1], 0);
					gl.glVertex3fv(point, 0);
					if(isLighted) gl.glNormal3fv(normals[i+1][j-1], 0);
					gl.glVertex3fv(point, 0);

					point[0] += xInc;
					point[1] -= yInc;
				}
				else
				{
					point[0] = startX + (j+colMin)*xInc;
					point[1] = startY + (i+rowMin)*yInc;
				}
				// middle sides

				for (; j<=t.lastIndex; j++)
				{
					point[2] = pointsZ[i][j];

					if(isLighted) gl.glNormal3fv(normals[i][j], 0);
					gl.glVertex3fv(point, 0);

					point[1] += yInc;
					point[2] = pointsZ[i+1][j];

					if(isLighted) gl.glNormal3fv(normals[i+1][j], 0);
					gl.glVertex3fv(point, 0);

					point[0] += xInc;
					point[1] -= yInc;
				}

				// last side: j has been bumped to next j in for-loop above
				if (t.lastSide == STRIP_BOT)
				{
					point[2] = pointsZ[i][j];

					if(isLighted) gl.glNormal3fv(normals[i][j], 0);
					gl.glVertex3fv(point, 0);
				}
				else if (t.lastSide == STRIP_TOP)
				{
					point[1] += yInc;
					point[2] = pointsZ[i+1][j];

					if(isLighted) gl.glNormal3fv(normals[i+1][j], 0);
					gl.glVertex3fv(point, 0);
				}

				gl.glEnd();
			}
		}
		catch(Exception e)
		{
			StsException.outputException("Exception in display. " +
				"row: " + t.rowNumber + " col: " + j, e, StsException.WARNING);
		}
	}
	static public void drawSurfaceFillWithNullsAndColor(GL gl, StsXYSurfaceGridable grid, StsList tStrips,
                            float[][] pointsZ, float[][][] normals, float[][] values, StsColorscale colorscale, boolean isLighted)
	{
		TriangleStrip t = null;
		float[] point = new float[3];
		int i, j = 0;
		int n;


		float startX = grid.getXMin();
		float startY = grid.getYMin();
        float xInc = grid.getXInc();
        float yInc = grid.getYInc();
        int rowMin = grid.getRowMin();
        int colMin = grid.getColMin();
        try
		{
//			StsColor.setGLColor(gl, this.stsColor);

			int nTStrips = tStrips.getSize();
			for (int nt = 0; nt < nTStrips; nt++)
			{
				t = (TriangleStrip)tStrips.getElement(nt);
				i = t.rowNumber;
				j = t.firstIndex;

				gl.glBegin(GL.GL_TRIANGLE_STRIP);

				// first side
				if (t.firstSide == STRIP_BOT)
				{
					point[0] = startX + (j-1+colMin)*xInc;
					point[1] = startY + (i+rowMin)*yInc;
					point[2] = pointsZ[i][j-1];
                    setColor(gl, values, i, j-1, colorscale);
                    if(isLighted) gl.glNormal3fv(normals[i][j-1], 0);
					gl.glVertex3fv(point, 0);
                    // setBeachballColors(gl, values, i, j-1, colorscale);
                    // if(isLighted) gl.glNormal3fv(normals[i][j-1], 0);
					// gl.glVertex3fv(point, 0);

					point[0] += xInc;
				}
				else if (t.firstSide == STRIP_TOP)
				{
					point[0] = startX + (j-1+colMin)*xInc;
					point[1] = startY + (i+1+rowMin)*yInc;
					point[2] = pointsZ[i+1][j-1];
                    setColor(gl, values, i+1, j-1, colorscale);  
					if(isLighted) gl.glNormal3fv(normals[i+1][j-1], 0);
					gl.glVertex3fv(point, 0);
                    // setBeachballColors(gl, values, i+1, j-1, colorscale);
                    // if(isLighted) gl.glNormal3fv(normals[i+1][j-1], 0);
					// gl.glVertex3fv(point, 0);

					point[0] += xInc;
					point[1] -= yInc;
				}
				else
				{
					point[0] = startX + (j+colMin)*xInc;
					point[1] = startY + (i+rowMin)*yInc;
				}
				// middle sides

				for (; j<=t.lastIndex; j++)
				{
					point[2] = pointsZ[i][j];
                    setColor(gl, values, i, j, colorscale);
					if(isLighted) gl.glNormal3fv(normals[i][j], 0);
					gl.glVertex3fv(point, 0);

					point[1] += yInc;
					point[2] = pointsZ[i+1][j];
                    setColor(gl, values, i+1, j, colorscale);
					if(isLighted) gl.glNormal3fv(normals[i+1][j], 0);
					gl.glVertex3fv(point, 0);

					point[0] += xInc;
					point[1] -= yInc;
				}

				// last side: j has been bumped to next j in for-loop above
				if (t.lastSide == STRIP_BOT)
				{
					point[2] = pointsZ[i][j];
                    setColor(gl, values, i, j, colorscale);
					if(isLighted) gl.glNormal3fv(normals[i][j], 0);
					gl.glVertex3fv(point, 0);
				}
				else if (t.lastSide == STRIP_TOP)
				{
					point[1] += yInc;
					point[2] = pointsZ[i+1][j];
                    setColor(gl, values, i+1, j, colorscale);
					if(isLighted) gl.glNormal3fv(normals[i+1][j], 0);
					gl.glVertex3fv(point, 0);
				}
				gl.glEnd();
			}
		}
		catch(Exception e)
		{
			StsException.outputException("Exception in display. " +
				"row: " + t.rowNumber + " col: " + j, e, StsException.WARNING);
		}
	}

    static private void setColor(GL gl, float[][] values, int row, int col, StsColorscale colorscale)
    {
        StsColor color = colorscale.getStsColor(colorscale.getIndexFromValue(values[row][col]));
        color.setGLColor(gl);
    }


    public void print() { System.out.println(toString()); }

    public String toString()
    {
        return new String("triStrip row: " + rowNumber + " firstIndex: " + firstIndex +
                           " firstSide: " + firstSide + " lastIndex: " + lastIndex +
                           " lastSide: " + lastSide);
    }
}

