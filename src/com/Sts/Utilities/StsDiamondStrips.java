package com.Sts.Utilities;

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
    float[][] rowLeftCenterPointsZ;
    float[][] colBotCenterPointsZ;
    float[][] rowRiteCenterPointsZ;
    float[][] colTopCenterPointsZ;
    float[][] pointsZ;
    float[][][] normals;
    float[][][] rowLeftCenterNormals;
    float[][][] colBotCenterNormals;
    float[][][] rowRiteCenterNormals;
    float[][][] colTopCenterNormals;
    boolean[][] hasRowLinks;
    boolean[][] hasColLinks;

    float[][] values;
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

        rowLeftCenterPointsZ = new float[nRows][nCols];
        rowRiteCenterPointsZ = new float[nRows][nCols];
        colBotCenterPointsZ = new float[nRows][nCols];
        colTopCenterPointsZ = new float[nRows][nCols];
        rowLeftCenterNormals = new float[nRows][nCols][];
        rowRiteCenterNormals = new float[nRows][nCols][];
        colBotCenterNormals = new float[nRows][nCols][];
        colTopCenterNormals = new float[nRows][nCols][];

        try
        {
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

                    if (nLinks > 2) // we have 3 or 4 links around grid, so use common center and normal
                    {
                        // averages 3 or 4 points and normals, ignoring the null one (all 4 should be present, though)
                        averageGridCellPointsAndNormals(row, col, usePoints, 4);
                    }
                    else if (nLinks == 2) // could be two links on opposite sides or 2 forming an L; if former than compute two centers, otherwise compute shared center
                    {
                        if (hasRowLink)
                        {
                            if (hasNextRowLink) // compute two column centers and normals 6
                            {
                                colBotCenterNormals[row][col] = StsMath.addVectorsNormalize(normals[row][col], normals[row][col + 1]);
                                colTopCenterNormals[row][col] = StsMath.addVectorsNormalize(normals[row + 1][col], normals[row + 1][col + 1]);
                                colBotCenterPointsZ[row][col] = averageValues(pointsZ[row][col], pointsZ[row][col + 1]);
                                colTopCenterPointsZ[row][col] = averageValues(pointsZ[row + 1][col], pointsZ[row + 1][col + 1]);
                            }
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
                                rowLeftCenterNormals[row][col] = StsMath.addVectorsNormalize(normals[row][col], normals[row + 1][col]);
                                rowRiteCenterNormals[row][col] = StsMath.addVectorsNormalize(normals[row][col + 1], normals[row + 1][col + 1]);
                                rowLeftCenterPointsZ[row][col] = averageValues(pointsZ[row][col], pointsZ[row + 1][col]);
                                rowRiteCenterPointsZ[row][col] = averageValues(pointsZ[row][col + 1], pointsZ[row + 1][col + 1]);
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
                        if (hasRowLink)
                        {
                            colBotCenterNormals[row][col] = StsMath.addVectorsNormalize(normals[row][col], normals[row][col + 1]);
                            colBotCenterPointsZ[row][col] = averageValues(pointsZ[row][col], pointsZ[row][col + 1]);
                        }
                        else if (hasNextRowLink)
                        {
                            colTopCenterNormals[row][col] = StsMath.addVectorsNormalize(normals[row + 1][col], normals[row + 1][col + 1]);
                            colTopCenterPointsZ[row][col] = averageValues(pointsZ[row + 1][col], pointsZ[row + 1][col + 1]);
                        }
                        else if (hasColLink)
                        {
                            rowLeftCenterNormals[row][col] = StsMath.addVectorsNormalize(normals[row][col], normals[row + 1][col]);
                            rowLeftCenterPointsZ[row][col] = averageValues(pointsZ[row][col], pointsZ[row + 1][col]);
                        }
                        else // hasNextColLink
                        {
                            rowRiteCenterNormals[row][col] = StsMath.addVectorsNormalize(normals[row][col + 1], normals[row + 1][col + 1]);
                            rowRiteCenterPointsZ[row][col] = averageValues(pointsZ[row][col + 1], pointsZ[row + 1][col + 1]);
                        }
                    }
                    else // nLinks == 0
                    {

                    }
                }
            }
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "computeGridCenterValues", "Failed at row: " + row + " col: " + col, e);
        }
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
        float[] normal = averageGridNormals(row, col, usePoints, min);
        rowLeftCenterNormals[row][col] = normal;
        colBotCenterNormals[row][col] = normal;
        rowRiteCenterNormals[row][col] = normal;
        colTopCenterNormals[row][col] = normal;
        float centerPointZ = averageGridPoints(row, col, usePoints, min);
        rowLeftCenterPointsZ[row][col] = centerPointZ;
        colBotCenterPointsZ[row][col] = centerPointZ;
        rowRiteCenterPointsZ[row][col] = centerPointZ;
        colTopCenterPointsZ[row][col] = centerPointZ;
    }

    private float[] averageGridNormals(int row, int col, boolean[][] usePoints, int min)
    {
        float[][] gridNormals = new float[4][];
        if (usePoints[0][0])
            gridNormals[0] = normals[row][col];
        if (usePoints[0][1])
            gridNormals[1] = normals[row][col + 1];
        if (usePoints[1][1])
            gridNormals[2] = normals[row + 1][col + 1];
        if (usePoints[1][0])
            gridNormals[3] = normals[row + 1][col];

        float[] normal = StsMath.addVectorsNormalize(gridNormals, 3, min);
        return normal;
    }

    private float averageGridPoints(int row, int col, boolean[][] usePoints, int min)
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

    private void computeGridRowCenterValues()
    {
        rowLeftCenterPointsZ = new float[nRows][nCols - 1];

        for (int row = 0; row < nRows; row++)
        {
            for (int col = 0; col < nCols - 1; col++)
                rowLeftCenterPointsZ[row][col] = averageValues(pointsZ[row][col], pointsZ[row][col + 1]);
        }
    }

    private float averageValues(float v1, float v2)
    {
        if (v1 == nullValue || v2 == nullValue) return nullValue;
        return (v1 + v2) / 2;
    }

    private void computeGridColCenterValues()
    {
        colBotCenterPointsZ = new float[nRows - 1][nCols];

        for (int col = 0; col < nCols; col++)
        {
            for (int row = 0; row < nRows - 1; row++)
                colBotCenterPointsZ[row][col] = averageValues(pointsZ[row][col], pointsZ[row + 1][col]);
        }
    }

    private void computeRowCenterNormals()
    {
        rowLeftCenterNormals = new float[nRows][nCols - 1][3];

        for (int row = 0; row < nRows; row++)
            for (int col = 0; col < nCols - 1; col++)
                rowLeftCenterNormals[row][col] = StsMath.addVectorsNormalize(normals[row][col], normals[row][col + 1]);
    }

    private void computeColCenterNormals()
    {
        colBotCenterNormals = new float[nRows - 1][nCols][3];

        for (int col = 0; col < nCols; col++)
            for (int row = 0; row < nRows - 1; row++)
                colBotCenterNormals[row][col] = StsMath.addVectorsNormalize(normals[row][col], normals[row + 1][col]);
    }

    private void computeHasRowLinks()
    {
        hasRowLinks = new boolean[nRows][nCols];

        for (int row = 0; row < nRows; row++)
        {
            for (int col = 0; col < nCols - 1; col++)
                hasRowLinks[row][col] = grid.hasRowLink(row, col);
            // hasRowLinks[row][nCols - 1] = false;
        }
        // Arrays.fill(hasRowLinks[nRows - 1], false);
    }

    private void computeHasColLinks()
    {
        hasColLinks = new boolean[nRows][nCols];

        for (int col = 0; col < nCols; col++)
        {
            for (int row = 0; row < nRows - 1; row++)
                hasColLinks[row][col] = grid.hasColLink(row, col);
            // hasColLinks[nRows-1][col] = false;
        }
        //for(int row = 0; row < nRows - 1; row++)
        //    hasColLinks[row][nCols-1] = grid.hasColLink(row, nCols);
    }

    final boolean hasRowLink(int row, int col)
    {
        if(row < nRows && col < nCols)
            return grid.hasRowLink(row, col);
        else
            return false;
    }

    final boolean hasColLink(int row, int col)
    {
        if(row < nRows && col < nCols)
            return grid.hasColLink(row, col);
        else
            return false;
    }

    final float getPointZ(int row, int col)
    {
        return pointsZ[row][col];
    }

    final float getRowCenterPointZ(int row, int col)
    {
        return rowLeftCenterPointsZ[row][col];
    }

    final float getColCenterPointZ(int row, int col)
    {
        return colBotCenterPointsZ[row][col];
    }

    final float[] getNormal(int row, int col)
    {
        float[] normal = normals[row][col];
        if (normal != null) return normal;
        else return verticalNormal;
    }

    final float[] getRowLeftCenterNormal(int row, int col)
    {
        float[] normal = rowLeftCenterNormals[row][col];
        if (normal != null) return normal;
        else return verticalNormal;
    }

    final float[] getColCenterNormal(int row, int col)
    {
        float[] normal = colBotCenterNormals[row][col];
        if (normal != null) return normal;
        else return verticalNormal;
    }

    final float getValue(int row, int col)
    {
        return values[row][col];
    }

    final float getRowCenterValue(int row, int col)
    {
        return rowCenterValues[row][col];
    }

    final float getColCenterValue(int row, int col)
    {
        return colCenterValues[row][col];
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
                if (hasRowLinks[row][col]) drawRowDiamond(gl, row, col, rowY, colX);
        }
    }

    private void drawColDiamondsWithNulls(GL gl)
    {
        float colX = xMin + colMin * xInc;
        for (int col = 0; col < nCols; col++, colX += xInc)
        {
            float rowY = yMin + rowMin * yInc;
            for (int row = 0; row < nRows - 1; row++, rowY += yInc)
                if (hasColLinks[row][col]) drawColDiamond(gl, row, col, rowY, colX);
        }
    }

    private void drawRowDiamond(GL gl, int row, int col, float rowY, float colX)
    {
        if (row == 0)
            drawTopRowTriangle(gl, row, col, rowY, colX);
        else if (row == nRows - 1)
            drawBotRowTriangle(gl, row, col, rowY, colX);
        else
        {
            try
            {
                gl.glBegin(GL.GL_QUADS);
                drawBotRowCenterPoint(gl, row, col, rowY, colX);
                drawRiteRowPoint(gl, row, col, rowY, colX);
                drawTopRowCenterPoint(gl, row, col, rowY, colX);
                drawLeftRowPoint(gl, row, col, rowY, colX);
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
    }

    private void drawBotRowTriangle(GL gl, int row, int col, float rowY, float colX)
    {
        try
        {
            gl.glBegin(GL.GL_TRIANGLES);
            drawBotRowCenterPoint(gl, row, col, rowY, colX);
            drawRiteRowPoint(gl, row, col, rowY, colX);
            drawLeftRowPoint(gl, row, col, rowY, colX);
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

    private void drawTopRowTriangle(GL gl, int row, int col, float rowY, float colX)
    {
        try
        {
            gl.glBegin(GL.GL_TRIANGLES);
            drawTopRowCenterPoint(gl, row, col, rowY, colX);
            drawLeftRowPoint(gl, row, col, rowY, colX);
            drawRiteRowPoint(gl, row, col, rowY, colX);
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

    private void drawBotRowCenterPoint(GL gl, int row, int col, float rowY, float colX)
    {
        gl.glNormal3fv(colTopCenterNormals[row - 1][col], 0);
        gl.glVertex3f(colX + xInc / 2, rowY - yInc / 2, colTopCenterPointsZ[row - 1][col]);
    }

    private void drawTopRowCenterPoint(GL gl, int row, int col, float rowY, float colX)
    {
        gl.glNormal3fv(colBotCenterNormals[row][col], 0);
        gl.glVertex3f(colX + xInc / 2, rowY + yInc / 2, colBotCenterPointsZ[row][col]);
    }

    private void drawColDiamond(GL gl, int row, int col, float rowY, float colX)
    {
        if (col == 0)
            drawRiteColTriangle(gl, row, col, rowY, colX);
        else if (col == nCols - 1)
            drawLeftColTriangle(gl, row, col, rowY, colX);
        else
        {
            try
            {
                gl.glBegin(GL.GL_QUADS);
                drawLeftColCenterPoint(gl, row, col, rowY, colX);
                drawBotColPoint(gl, row, col, rowY, colX);
                drawRiteColCenterPoint(gl, row, col, rowY, colX);
                drawTopColPoint(gl, row, col, rowY, colX);
                gl.glEnd();
            } catch (Exception e)
            {
                StsException.outputWarningException(this, "drawRowDiamond", e);
            } finally
            {
                gl.glEnd();
            }
        }
    }

    private void drawLeftColTriangle(GL gl, int row, int col, float rowY, float colX)
    {
        try
        {
            gl.glBegin(GL.GL_TRIANGLES);
            drawLeftColCenterPoint(gl, row, col, rowY, colX);
            drawBotColPoint(gl, row, col, rowY, colX);
            drawTopColPoint(gl, row, col, rowY, colX);
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

    private void drawRiteColTriangle(GL gl, int row, int col, float rowY, float colX)
    {
        try
        {
            gl.glBegin(GL.GL_TRIANGLES);
            drawRiteColCenterPoint(gl, row, col, rowY, colX);
            drawTopColPoint(gl, row, col, rowY, colX);
            drawBotColPoint(gl, row, col, rowY, colX);
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

    private void drawLeftColCenterPoint(GL gl, int row, int col, float rowY, float colX)
    {
        gl.glNormal3fv(rowRiteCenterNormals[row][col - 1], 0);
        gl.glVertex3f(colX - xInc / 2, rowY + yInc / 2, rowRiteCenterPointsZ[row][col - 1]);
    }

    private void drawRiteColCenterPoint(GL gl, int row, int col, float rowY, float colX)
    {
        gl.glNormal3fv(rowLeftCenterNormals[row][col], 0);
        gl.glVertex3f(colX + xInc / 2, rowY + yInc / 2, rowLeftCenterPointsZ[row][col]);
    }

    private void drawLeftRowPoint(GL gl, int row, int col, float rowY, float colX)
    {
        drawPoint(gl, row, col, rowY, colX);
    }

    private void drawRiteRowPoint(GL gl, int row, int col, float rowY, float colX)
    {
        drawPoint(gl, row, col+1, rowY, colX + xInc);
    }

    private void drawBotColPoint(GL gl, int row, int col, float rowY, float colX)
    {
        drawPoint(gl, row, col, rowY, colX);
    }

    private void drawTopColPoint(GL gl, int row, int col, float rowY, float colX)
    {
        drawPoint(gl, row+1, col, rowY + yInc, colX);
    }

    private void drawPoint(GL gl, int row, int col, float rowY, float colX)
    {
        gl.glNormal3fv(getNormal(row, col), 0);
        if(pointsZ[row][col] == nullValue)
            StsException.systemError(this, "drawPoint", "Null Z at row: " + row + " col: " + col);
        gl.glVertex3f(colX, rowY, pointsZ[row][col]);
    }
}
