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
    float[][] leftCenterPointsZ;
    float[][] botCenterPointsZ;
    float[][] riteCenterPointsZ;
    float[][] topCenterPointsZ;
    float[][] pointsZ;
    float[][][] normals;
    float[][][] leftCenterNormals;
    float[][][] botCenterNormals;
    float[][][] riteCenterNormals;
    float[][][] topCenterNormals;
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

        leftCenterPointsZ = new float[nRows][nCols];
        riteCenterPointsZ = new float[nRows][nCols];
        botCenterPointsZ = new float[nRows][nCols];
        topCenterPointsZ = new float[nRows][nCols];
        leftCenterNormals = new float[nRows][nCols][];
        riteCenterNormals = new float[nRows][nCols][];
        botCenterNormals = new float[nRows][nCols][];
        topCenterNormals = new float[nRows][nCols][];

        try
        {
            for (row = 0; row < nRows; row++)
            {
                Arrays.fill(leftCenterPointsZ[row], nullValue);
                Arrays.fill(riteCenterPointsZ[row], nullValue);
                Arrays.fill(botCenterPointsZ[row], nullValue);
                Arrays.fill(topCenterPointsZ[row], nullValue);
            }
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
                                if(row < nRows+1) leftCenterNormals[row][col] = addVectorsNormalize(normals[row][col], normals[row + 1][col]);
                                if(col < nCols+1) riteCenterNormals[row][col] = addVectorsNormalize(normals[row][col + 1], normals[row + 1][col + 1]);
                                if(row < nRows+1) leftCenterPointsZ[row][col] = averageValues(pointsZ[row][col], pointsZ[row + 1][col]);
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
                        if (hasRowLink)
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
                            leftCenterNormals[row][col] = addVectorsNormalize(normals[row][col], normals[row + 1][col]);
                            leftCenterPointsZ[row][col] = averageValues(pointsZ[row][col], pointsZ[row + 1][col]);
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
        float[] normal = averageGridNormals(row, col, usePoints, min);
        leftCenterNormals[row][col] = normal;
        botCenterNormals[row][col] = normal;
        riteCenterNormals[row][col] = normal;
        topCenterNormals[row][col] = normal;
        float centerPointZ = averageGridPoints(row, col, usePoints, min);
        leftCenterPointsZ[row][col] = centerPointZ;
        botCenterPointsZ[row][col] = centerPointZ;
        riteCenterPointsZ[row][col] = centerPointZ;
        topCenterPointsZ[row][col] = centerPointZ;
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
        leftCenterPointsZ = new float[nRows][nCols - 1];

        for (int row = 0; row < nRows; row++)
        {
            for (int col = 0; col < nCols - 1; col++)
                leftCenterPointsZ[row][col] = averageValues(pointsZ[row][col], pointsZ[row][col + 1]);
        }
    }

    private float averageValues(float v1, float v2)
    {
        if (v1 == nullValue || v2 == nullValue) return nullValue;
        return (v1 + v2) / 2;
    }

    private void computeGridColCenterValues()
    {
        botCenterPointsZ = new float[nRows - 1][nCols];

        for (int col = 0; col < nCols; col++)
        {
            for (int row = 0; row < nRows - 1; row++)
                botCenterPointsZ[row][col] = averageValues(pointsZ[row][col], pointsZ[row + 1][col]);
        }
    }

    private void computeRowCenterNormals()
    {
        leftCenterNormals = new float[nRows][nCols - 1][3];

        for (int row = 0; row < nRows; row++)
            for (int col = 0; col < nCols - 1; col++)
                leftCenterNormals[row][col] = StsMath.addVectorsNormalize(normals[row][col], normals[row][col + 1]);
    }

    private void computeColCenterNormals()
    {
        botCenterNormals = new float[nRows - 1][nCols][3];

        for (int col = 0; col < nCols; col++)
            for (int row = 0; row < nRows - 1; row++)
                botCenterNormals[row][col] = StsMath.addVectorsNormalize(normals[row][col], normals[row + 1][col]);
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
        float pointZ = pointsZ[row][col];
        if(pointZ != nullValue) return pointZ;
        StsException.systemError(this, "getPointZ", "Null Z at row: " + row + " col: " + col);
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
        float[] normal = leftCenterNormals[row][col];
        if(normal != null) return normal;
        StsException.systemDebug(this, "getLeftCenterNormal", "Failed for row: " + row + " col: " + col);
        return verticalNormal;
     }

    final float[] getRiteCenterNormal(int row, int col)
    {
        float[] normal = riteCenterNormals[row][col];
        if(normal != null) return normal;
        StsException.systemDebug(this, "getRiteCenterNormal", "Failed for row: " + row + " col: " + col);
        return verticalNormal;    }

    final float[] getBotCenterNormal(int row, int col)
    {
        float[] normal = botCenterNormals[row][col];
        if(normal != null) return normal;
        StsException.systemDebug(this, "getBotCenterNormal", "Failed for row: " + row + " col: " + col);
        return verticalNormal;
    }

    final float[] getTopCenterNormal(int row, int col)
    {
        float[] normal = topCenterNormals[row][col];
        if(normal != null) return normal;
        StsException.systemDebug(this, "getTopCenterNormal", "Failed for row: " + row + " col: " + col);
        return verticalNormal;
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

    /** row diamond is from row,col (left) to row,col+1 (right) with botCenter at row,col (above) and topCenter at row-1,col (below) */
    private void drawRowDiamond(GL gl, int row, int col, float rowY, float colX)
    {
        if (row == 0) // only draw triangle above row axis */
            drawTopTriangle(gl, row, col, rowY, colX);
        else if (row == nRows - 1)
            drawBotTriangle(gl, row, col, rowY, colX);
        else
        {
            try
            {
                gl.glBegin(GL.GL_QUADS);
                drawLeftPoint(gl, row, col, rowY, colX); // left
                drawBotCenterPoint(gl, row, col, rowY, colX); // below
                drawRitePoint(gl, row, col, rowY, colX + xInc / 2); // rite
                drawTopCenterPoint(gl, row, col, rowY, colX); // above

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

    /** draw the center below this row in this col which is the topCenter of the row below.
     *  coordinates are at the given row,col, so add xInc/2 and subtract yInc/2 */
    private void drawBotCenterPoint(GL gl, int row, int col, float rowY, float colX)
    {
        gl.glNormal3fv(getTopCenterNormal(row-1, col), 0);
        gl.glVertex3f(colX + xInc / 2, rowY - yInc / 2, topCenterPointsZ[row-1][col]);
    }

    /** draw the center above this row in this col which is the botCenter of this row.
     *  coordinates are at the given row,col, so add xInc/2 and add yInc/2 */
    private void drawTopCenterPoint(GL gl, int row, int col, float rowY, float colX)
    {
        gl.glNormal3fv(getBotCenterNormal(row, col), 0);
        gl.glVertex3f(colX + xInc / 2, rowY + yInc / 2, botCenterPointsZ[row][col]);
    }

    /** col diamond is from row,col (bot) to row+1,col (top) with centers at row-1,col (left) and row,col (rite) */
    private void drawColDiamond(GL gl, int row, int col, float rowY, float colX)
    {
        if (col == 0)
            drawRiteTriangle(gl, row, col, rowY, colX);
        else if (col == nCols - 1)
            drawLeftTriangle(gl, row, col, rowY, colX);
        else
        {
            try
            {
                gl.glBegin(GL.GL_QUADS);
                drawLeftCenterPoint(gl, row, col, rowY, colX);
                drawBotPoint(gl, row, col, rowY, colX);
                drawRiteCenterPoint(gl, row, col, rowY, colX);
                drawTopPoint(gl, row, col, rowY, colX);
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

    /** draw the center left of this col in this row which is the riteCenter of the cell to the left.
     *  coordinates are at the given row,col, so subtract xInc/2 and add yInc/2 */
    private void drawLeftCenterPoint(GL gl, int row, int col, float rowY, float colX)
    {
        gl.glNormal3fv(getRiteCenterNormal(row, col - 1), 0);
        gl.glVertex3f(colX - xInc / 2, rowY + yInc / 2, leftCenterPointsZ[row][col-1]);
    }

    /** draw the center rite of this col in this row which is the leftCenter of this cell.
     *  coordinates are at the given row,col, so add xInc/2 and add yInc/2 */
    private void drawRiteCenterPoint(GL gl, int row, int col, float rowY, float colX)
    {

        gl.glNormal3fv(getLeftCenterNormal(row, col), 0);
        gl.glVertex3f(colX + xInc / 2, rowY + yInc / 2, leftCenterPointsZ[row][col]);
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
        drawPoint(gl, row+1, col, rowY + yInc, colX);
    }

    private void drawPoint(GL gl, int row, int col, float rowY, float colX)
    {
        gl.glNormal3fv(getNormal(row, col), 0);
        gl.glVertex3f(colX, rowY, getPointZ(row, col));
    }
}
