package com.Sts.Utilities;

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;

import javax.media.opengl.*;
import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Oct 15, 2009
 * Time: 1:26:36 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsDiamondStrips
{
    StsXYSurfaceLinkGridable grid;
    int nRows;
    int nCols;
    float xInc;
    float yInc;
    float xMin;
    float yMin;
    int rowMin;
    int colMin;
    float[][] pointsZ;
    float[][][] normals;
	float[][] centerPointsZ;
    float[][][] centerNormals;
	StsColor[][] centerColors;
	/** rowLinks are types: NONE, LINE, ABOVE, BELOW, BOTH. Array is one row larger than grid and is filled with LINK_NONE */
    byte[][] rowLinks;
	/** colLinks are types: NONE, LINE, LEFT, RIGHT, BOTH. Array is one col larger than grid and is filled with LINK_NONE */
    byte[][] colLinks;
	/** display property on grid */
	boolean displayProperty = false;

    private float[][] values;
    float[][] rowCenterValues;
    float[][] colCenterValues;

    static final float nullValue = StsParameters.nullValue;
    static final float[] verticalNormal = new float[]{0.0f, 0.0f, -1.0f};

    public StsDiamondStrips(StsXYSurfaceLinkGridable grid)
    {
        this.grid = grid;
        nRows = grid.getNRows();
        nCols = grid.getNCols();
        xInc = grid.getXInc();
        yInc = grid.getYInc();
        xMin = grid.getXMin();
        yMin = grid.getYMin();
        rowMin = grid.getRowMin();
        colMin = grid.getColMin();
        pointsZ = grid.getPointsZ();
        normals = StsToolkit.computeSmoothNormals(pointsZ, nRows, nCols, xInc, yInc);
        computeHasRowLinks();
        computeHasColLinks();
        computeGridCenterValues();
    }

    private void computeGridCenterValues()
    {
        int row = -1, col = -1;

        centerPointsZ = new float[nRows][nCols];
        centerNormals = new float[nRows][nCols][];

        try
        {
            for (row = 0; row < nRows; row++)
                Arrays.fill(centerPointsZ[row], nullValue);
            // iterate over grid cells and assign values for each
            for (row = 0; row < nRows; row++)
            {
                for (col = 0; col < nCols; col++)
                {
                    boolean hasRowLink = hasRowLink(row, col);
                    boolean hasColLink = hasColLink(row, col);
                    boolean hasNextRowLink = hasRowLink(row + 1, col);
                    boolean hasNextColLink = hasColLink(row, col + 1);
                    int nLinks = 0;
                    if (hasRowLink) nLinks++;
                    if (hasColLink) nLinks++;
                    if (hasNextRowLink) nLinks++;
                    if (hasNextColLink) nLinks++;
                    boolean[][] usePoints = getUsePointsFromLinks(hasRowLink, hasColLink, hasNextRowLink, hasNextColLink);
					averageGridCellPointsAndNormals(row, col, usePoints, 1);
				/*
                    if (nLinks > 2) // we have 3 or 4 links around grid, so use common center and normal
                    {
                        // averages 3 or 4 points and normals, ignoring the null one (all 4 should be present, though)
                        averageGridCellPointsAndNormals(row, col, usePoints, 4);
                    }
                    else if (nLinks == 2) // could be two links on opposite sides or 2 forming an L; if former than compute two centers, otherwise compute shared center
                    {
                        if (computeRowLink)
                        {
                            if (hasNextRowLink) // compute two column centers and normals
                            {
                                if(col < nCols+1) botCenterNormals[row][col] = addVectorsNormalize(normals[row][col], normals[row][col + 1]);
                                if(row < nRows+1) topCenterNormals[row][col] = addVectorsNormalize(normals[row + 1][col], normals[row + 1][col + 1]);
                                if(col < nCols+1) botCenterPointsZ[row][col] = averageValues(pointsZ[row][col], pointsZ[row][col + 1]);
                                if(row < nRows+1) topCenterPointsZ[row][col] = averageValues(pointsZ[row + 1][col], pointsZ[row + 1][col + 1]);
                            }
                            // with two links we will have either 3 or 4 active points, so min argument below is 3
                            else if (hasColLink) // 1
                            {
                                averageGridCellPointsAndNormals(row, col, usePoints, 3);
                            }
                            else // hasNextColLink 2
                            {
                                averageGridCellPointsAndNormals(row, col, usePoints, 3);
                            }
                        }
                        else if (hasColLink)
                        {
                            if (hasNextColLink) // compute two row centers and normals   5
                            {
                                if(row < nRows+1) centerNormals[row][col] = addVectorsNormalize(normals[row][col], normals[row + 1][col]);
                                if(col < nCols+1) riteCenterNormals[row][col] = addVectorsNormalize(normals[row][col + 1], normals[row + 1][col + 1]);
                                if(row < nRows+1) centerPointsZ[row][col] = averageValues(pointsZ[row][col], pointsZ[row + 1][col]);
                                if(col < nCols+1) riteCenterPointsZ[row][col] = averageValues(pointsZ[row][col + 1], pointsZ[row + 1][col + 1]);
                            }
                            else //  no rowLink or nextColLink so must have nextRowlink 4
                            {
                                averageGridCellPointsAndNormals(row, col, usePoints, 3);
                            }
                        }
                        else if (hasNextRowLink)// doesn't have rowLink or colLink so can only have nextColLink 3
                        {
                            averageGridCellPointsAndNormals(row, col, usePoints, 3);
                        }
                    }
                    else if (nLinks == 1)
                    {
                        if (computeRowLink)
                        {
                            botCenterNormals[row][col] = addVectorsNormalize(normals[row][col], normals[row][col + 1]);
                            botCenterPointsZ[row][col] = averageValues(pointsZ[row][col], pointsZ[row][col + 1]);
                        }
                        else if (hasNextRowLink)
                        {
                            topCenterNormals[row][col] = addVectorsNormalize(normals[row + 1][col], normals[row + 1][col + 1]);
                            topCenterPointsZ[row][col] = averageValues(pointsZ[row + 1][col], pointsZ[row + 1][col + 1]);
                        }
                        else if (hasColLink)
                        {
                            centerNormals[row][col] = addVectorsNormalize(normals[row][col], normals[row + 1][col]);
                            centerPointsZ[row][col] = averageValues(pointsZ[row][col], pointsZ[row + 1][col]);
                        }
                        else // hasNextColLink
                        {
                            riteCenterNormals[row][col] = addVectorsNormalize(normals[row][col + 1], normals[row + 1][col + 1]);
                            riteCenterPointsZ[row][col] = averageValues(pointsZ[row][col + 1], pointsZ[row + 1][col + 1]);
                        }
                    }
                    else // nLinks == 0
                    {

                    }
                */
                }
            }
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "computeGridCenterValues", "Failed at row: " + row + " col: " + col, e);
        }
    }

    public static final float[] addVectorsNormalize(float[] a, float[] b)
    {
        float[] vector = StsMath.addVectorsNormalize(a, b);
        if(vector != null) return vector;
        else               return verticalNormal;
    }

    private boolean[][] getUsePointsFromLinks(boolean hasRowLink, boolean hasColLink, boolean hasNextRowLink, boolean hasNextColLink)
    {
        boolean[][] usePoints = new boolean[2][2];
        if (hasRowLink || hasColLink)
            usePoints[0][0] = true;
        if (hasRowLink || hasNextColLink)
            usePoints[0][1] = true;
        if (hasNextColLink || hasNextRowLink)
            usePoints[1][1] = true;
        if (hasNextRowLink || hasColLink)
            usePoints[1][0] = true;
        return usePoints;
    }

    private void averageGridCellPointsAndNormals(int row, int col, boolean[][] usePoints, int min)
    {
		centerNormals[row][col] = averageGridNormals(row, col, usePoints, min);
		centerPointsZ[row][col] = averageGridPoints(row, col, usePoints, min);
    }

    private float[] averageGridNormals(int row, int col, boolean[][] usePoints, int min)
    {
        try
        {
            float[][] gridNormals = new float[4][];
            if (usePoints[0][0])
                gridNormals[0] = normals[row][col];
            if (usePoints[0][1] && col < nCols)
                gridNormals[1] = normals[row][col + 1];
            if (usePoints[1][1] && row < nRows && col < nCols)
                gridNormals[2] = normals[row + 1][col + 1];
            if (usePoints[1][0] && row < nRows)
                gridNormals[3] = normals[row + 1][col];

            float[] normal = StsMath.addVectorsNormalize(gridNormals, 3, min);
            return normal;
        }
        catch(Exception e)
        {
            return verticalNormal;
        }
    }

    private float averageGridPoints(int row, int col, boolean[][] usePoints, int min)
    {
        try
        {
            float[] gridPoints = new float[4];
            Arrays.fill(gridPoints, nullValue);
            if (usePoints[0][0])
                gridPoints[0] = pointsZ[row][col];
            if (usePoints[0][1])
                gridPoints[1] = pointsZ[row][col + 1];
            if (usePoints[1][1])
                gridPoints[2] = pointsZ[row + 1][col + 1];
            if (usePoints[1][0])
                gridPoints[3] = pointsZ[row + 1][col];

            return StsMath.average(gridPoints, nullValue, min);
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "averageGridPoints", e);
            return nullValue;
        }
    }

    private void computeGridRowCenterValues()
    {
        centerPointsZ = new float[nRows][nCols - 1];

        for (int row = 0; row < nRows; row++)
        {
            for (int col = 0; col < nCols - 1; col++)
                centerPointsZ[row][col] = averageValues(pointsZ[row][col], pointsZ[row][col + 1]);
        }
    }

    private float averageValues(float v1, float v2)
    {
        if (v1 == nullValue || v2 == nullValue) return nullValue;
        return (v1 + v2) / 2;
    }

    private void computeHasRowLinks()
    {
        rowLinks = new byte[nRows+1][nCols];

        for (int row = 0; row < nRows; row++)
        {
            for (int col = 0; col < nCols - 1; col++)
                rowLinks[row][col] = grid.computeRowLink(row, col);
            // rowLinks[row][nCols - 1] = false;
        }
        // Arrays.fill(rowLinks[nRows - 1], false);
    }

    private void computeHasColLinks()
    {
        colLinks = new byte[nRows][nCols+1];

        for (int col = 0; col < nCols; col++)
        {
            for (int row = 0; row < nRows - 1; row++)
                colLinks[row][col] = grid.computeColLink(row, col);
            // colLinks[nRows-1][col] = false;
        }
        //for(int row = 0; row < nRows - 1; row++)
        //    colLinks[row][nCols-1] = grid.hasColLink(row, nCols);
    }

    final boolean hasRowLink(int row, int col)
    {
		return rowLinks[row][col] != StsPatchGrid.LINK_NONE;
    }

    final boolean hasColLink(int row, int col)
    {
		return colLinks[row][col] != StsPatchGrid.LINK_NONE;
    }

    final float getPointZ(int row, int col)
    {
        float pointZ = pointsZ[row][col];
        if(pointZ != nullValue) return pointZ;
        StsException.systemError(this, "getPointZ", "NULL Z FOR grid: " + ((StsPatchGrid)grid).id + " at volRow: " + (row + grid.getRowMin()) +
                " volCol: " + (col + grid.getColMin()));
        return 0.0f;
    }

    final float[] getNormal(int row, int col)
    {
        float[] normal = normals[row][col];
        if (normal != null) return normal;
        StsException.systemDebug(this, "getNormal", "Failed for row: " + row + " col: " + col);
        return verticalNormal;
    }

    final float[] getLeftCenterNormal(int row, int col)
    {
        if(centerNormals[row][col-1] != null) return centerNormals[row][col-1];
        StsException.systemDebug(this, "getLeftCenterNormal", "Failed for row: " + row + " col: " + col);
        return verticalNormal;
     }

    final float[] getRiteCenterNormal(int row, int col)
    {
        if(centerNormals[row][col] != null) return centerNormals[row][col];
        StsException.systemDebug(this, "getRiteCenterNormal", "Failed for row: " + row + " col: " + col);
        return verticalNormal;    }

    final float[] getBotCenterNormal(int row, int col)
    {
        if(centerNormals[row-1][col] != null) return centerNormals[row-1][col];
        StsException.systemDebug(this, "getBotCenterNormal", "Failed for row: " + row + " col: " + col);
        return verticalNormal;
    }

    final float[] getTopCenterNormal(int row, int col)
    {
        if(centerNormals[row][col] != null) return centerNormals[row][col];
        StsException.systemDebug(this, "getTopCenterNormal", "Failed for row: " + row + " col: " + col);
        return verticalNormal;
    }

    public void drawSurfaceFillWithNulls(GL gl)
    {
        try
        {
            drawRowDiamondsWithNulls(gl);
            drawColDiamondsWithNulls(gl);
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "drawSurfaceFillWithNulls", e);
        }
    }

    private void drawRowDiamondsWithNulls(GL gl)
    {
        float rowY = yMin + rowMin * yInc;
        for (int row = 0; row < nRows; row++, rowY += yInc)
        {
            float colX = xMin + colMin * xInc;
            for (int col = 0; col < nCols - 1; col++, colX += xInc)
                drawRowDiamond(rowLinks[row][col], row, col, rowY, colX, gl);
        }
    }

    private void drawColDiamondsWithNulls(GL gl)
    {
        float colX = xMin + colMin * xInc;
        for (int col = 0; col < nCols; col++, colX += xInc)
        {
            float rowY = yMin + rowMin * yInc;
            for (int row = 0; row < nRows - 1; row++, rowY += yInc)
                drawColDiamond(colLinks[row][col], row, col, rowY, colX, gl);
        }
    }

    /** row diamond is from row,col (left) to row,col+1 (right) with botCenter at row,col (above) and topCenter at row-1,col (below) */
    private void drawRowDiamond(byte linkType, int row, int col, float rowY, float colX, GL gl)
	{
		try
		{
			if (linkType == StsPatchGrid.LINK_NONE)
				return;

            if(StsPatchVolume.debug)
            {
                float rowDZ = rowDZ(row, col);
                if(rowDZ > maxDZ)
                    StsException.systemDebug(this, "DrawRowDiamond", "rowDZ " + rowDZ + " exceeds maxDZ");
            }

		//	if(linkType == StsPatchGrid.LINK_NONE)
		//		drawRowLine(gl, row, col, rowY, colX);
		//	else // draw above and/or below triangles
			{
				if (linkType == StsPatchGrid.LINK_ABOVE || linkType == StsPatchGrid.LINK_BOTH)
					drawTopTriangle(gl, row, col, rowY, colX);
				if (linkType == StsPatchGrid.LINK_BELOW || linkType == StsPatchGrid.LINK_BOTH)
					drawBotTriangle(gl, row, col, rowY, colX);
			}
		}
		catch(Exception e)
		{
			StsException.outputWarningException(this, "drawRowDiamond", e);
		}
		finally
		{
			gl.glEnd();
		}
	}

    private void drawBotTriangle(GL gl, int row, int col, float rowY, float colX)
    {
        try
        {
            gl.glBegin(GL.GL_TRIANGLES);
            drawLeftPoint(gl, row, col, rowY, colX);
            drawBotCenterPoint(gl, row, col, rowY, colX);
            drawRitePoint(gl, row, col, rowY, colX);
            gl.glEnd();
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "drawRowDiamond", e);
        }
        finally
        {
            gl.glEnd();
        }
    }

    private void drawTopTriangle(GL gl, int row, int col, float rowY, float colX)
    {
        try
        {
            gl.glBegin(GL.GL_TRIANGLES);
            drawLeftPoint(gl, row, col, rowY, colX);
            drawRitePoint(gl, row, col, rowY, colX);
            drawTopCenterPoint(gl, row, col, rowY, colX);
            gl.glEnd();
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "drawRowDiamond", e);
        }
        finally
        {
            gl.glEnd();
        }
    }

	private void drawRowLine(GL gl, int row, int col, float rowY, float colX)
	{
		try
		{
			gl.glBegin(GL.GL_LINES);
			drawLeftPoint(gl, row, col, rowY, colX);
			drawRitePoint(gl, row, col, rowY, colX);
			gl.glEnd();
		}
		catch(Exception e)
		{
			StsException.outputWarningException(this, "drawRowDiamond", e);
		}
		finally
		{
			gl.glEnd();
		}
	}
    /** draw the center below this row in this col which is the topCenter of the row below.
     *  coordinates are at the given row,col, so add xInc/2 and subtract yInc/2 */
    private void drawBotCenterPoint(GL gl, int row, int col, float rowY, float colX)
    {
        gl.glNormal3fv(getBotCenterNormal(row, col), 0);
        gl.glVertex3f(colX + xInc / 2, rowY - yInc / 2, centerPointsZ[row-1][col]);
    }

    /** draw the center above this row in this col which is the botCenter of this row.
     *  coordinates are at the given row,col, so add xInc/2 and add yInc/2 */
    private void drawTopCenterPoint(GL gl, int row, int col, float rowY, float colX)
    {
        gl.glNormal3fv(getTopCenterNormal(row, col), 0);
        gl.glVertex3f(colX + xInc / 2, rowY + yInc / 2, centerPointsZ[row][col]);
    }

    /** col diamond is from row,col (bot) to row+1,col (top) with centers at row-1,col (left) and row,col (rite) */
    private void drawColDiamond(byte linkType, int row, int col, float rowY, float colX, GL gl)
    {
		try
		{
			if (linkType == StsPatchGrid.LINK_NONE) return;

            if(StsPatchVolume.debug)
            {
                float colDZ = colDZ(row, col);
                if(colDZ > maxDZ)
                    StsException.systemDebug(this, "drawColDiamond", "colDZ " + colDZ + " exceeds maxDZ");
            }

        //    if(linkType == StsPatchGrid.LINK_LINE)
		//		drawColLine(gl, row, col, rowY, colX);
		//	else // draw left and or right triangles
			{
				if (linkType == StsPatchGrid.LINK_LEFT || linkType == StsPatchGrid.LINK_BOTH)
					drawLeftTriangle(gl, row, col, rowY, colX);
				if (linkType == StsPatchGrid.LINK_RIGHT || linkType == StsPatchGrid.LINK_BOTH)
					drawRiteTriangle(gl, row, col, rowY, colX);
			}
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "drawColDiamond", e);
		}
		finally
		{
			gl.glEnd();
		}
    }

    static final float maxDZ = 20.0f;

    private float colDZ(int row, int col)
    {
        return Math.abs(getPointZ(row, col) - getPointZ(row+1, col));
    }

    private float rowDZ(int row, int col)
    {
        return Math.abs(getPointZ(row, col) - getPointZ(row, col+1));
    }

    private void drawLeftTriangle(GL gl, int row, int col, float rowY, float colX)
    {
        try
        {
            gl.glBegin(GL.GL_TRIANGLES);
            drawBotPoint(gl, row, col, rowY, colX);
            drawTopPoint(gl, row, col, rowY, colX);
            drawLeftCenterPoint(gl, row, col, rowY, colX);
            gl.glEnd();
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "drawRowDiamond", e);
        }
        finally
        {
            gl.glEnd();
        }
    }

    private void drawRiteTriangle(GL gl, int row, int col, float rowY, float colX)
    {
        try
        {
            gl.glBegin(GL.GL_TRIANGLES);
            drawBotPoint(gl, row, col, rowY, colX);
            drawRiteCenterPoint(gl, row, col, rowY, colX);
            drawTopPoint(gl, row, col, rowY, colX);
            gl.glEnd();
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "drawRowDiamond", e);
        }
        finally
        {
            gl.glEnd();
        }
    }

	private void drawColLine(GL gl, int row, int col, float rowY, float colX)
	{
		try
		{
			gl.glBegin(GL.GL_LINES);
			drawBotPoint(gl, row, col, rowY, colX);
			drawTopPoint(gl, row, col, rowY, colX);
			gl.glEnd();
		}
		catch(Exception e)
		{
			StsException.outputWarningException(this, "drawRowDiamond", e);
		}
		finally
		{
			gl.glEnd();
		}
	}

    /** draw the center left of this col in this row which is the riteCenter of the cell to the left.
     *  coordinates are at the given row,col, so subtract xInc/2 and add yInc/2 */
    private void drawLeftCenterPoint(GL gl, int row, int col, float rowY, float colX)
    {
        gl.glNormal3fv(getLeftCenterNormal(row, col), 0);
        gl.glVertex3f(colX - xInc / 2, rowY + yInc / 2, centerPointsZ[row][col-1]);
    }

    /** draw the center rite of this col in this row which is the leftCenter of this cell.
     *  coordinates are at the given row,col, so add xInc/2 and add yInc/2 */
    private void drawRiteCenterPoint(GL gl, int row, int col, float rowY, float colX)
    {

        gl.glNormal3fv(getRiteCenterNormal(row, col), 0);
        gl.glVertex3f(colX + xInc / 2, rowY + yInc / 2, centerPointsZ[row][col]);
    }

    private void drawLeftPoint(GL gl, int row, int col, float rowY, float colX)
    {
        drawPoint(gl, row, col, rowY, colX);
    }

    private void drawRitePoint(GL gl, int row, int col, float rowY, float colX)
    {
        drawPoint(gl, row, col+1, rowY, colX + xInc);
    }

    private void drawBotPoint(GL gl, int row, int col, float rowY, float colX)
    {
        drawPoint(gl, row, col, rowY, colX);
    }

    private void drawTopPoint(GL gl, int row, int col, float rowY, float colX)
    {
        drawPoint(gl, row + 1, col, rowY + yInc, colX);
    }

    private void drawPoint(GL gl, int row, int col, float rowY, float colX)
    {
        gl.glNormal3fv(getNormal(row, col), 0);
        gl.glVertex3f(colX, rowY, getPointZ(row, col));
    }
    public void drawGridLines(GL gl)
    {
        try
        {
            drawRowGridLines(gl);
            drawColGridLines(gl);
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "drawSurfaceFillWithNulls", e);
        }
    }

    private void drawRowGridLines(GL gl)
    {
        float rowY = yMin + rowMin * yInc;
        for (int row = 0; row < nRows; row++, rowY += yInc)
        {
            float colX = xMin + colMin * xInc;
            for (int col = 0; col < nCols - 1; col++, colX += xInc)
                if (rowLinks[row][col] != StsPatchGrid.LINK_NONE) drawRowGridLine(gl, row, col, rowY, colX);
        }
    }

    private void drawRowGridLine(GL gl, int row, int col, float rowY, float colX)
    {
        try
        {
            gl.glBegin(GL.GL_LINES);
            drawLeftPoint(gl, row, col, rowY, colX); // left
            drawRitePoint(gl, row, col, rowY, colX + xInc / 2); // rite
            gl.glEnd();
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "drawRowDiamond", e);
        }
        finally
        {
            gl.glEnd();
        }
    }

    private void drawColGridLines(GL gl)
    {
        float colX = xMin + colMin * xInc;
        for (int col = 0; col < nCols; col++, colX += xInc)
        {
            float rowY = yMin + rowMin * yInc;
            for (int row = 0; row < nRows - 1; row++, rowY += yInc)
                if (rowLinks[row][col] != StsPatchGrid.LINK_NONE) drawColGridLine(gl, row, col, rowY, colX);
        }
    }

    private void drawColGridLine(GL gl, int row, int col, float rowY, float colX)
    {
        try
        {
            gl.glBegin(GL.GL_LINES);
            drawBotPoint(gl, row, col, rowY, colX);
            drawTopPoint(gl, row, col, rowY, colX);
            gl.glEnd();
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "drawRowDiamond", e);
        }
        finally
        {
            gl.glEnd();
        }
    }

	public float[][] getValues()
	{
		return values;
	}

	public void setValues(float[][] values)
	{
		this.values = values;
	}
}
