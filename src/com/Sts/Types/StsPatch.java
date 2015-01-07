
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;


import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Reflect.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;

public class StsPatch
{
   // arrays of StsPoint for points & normals
    transient private StsPoint[][] points = null;
    transient private StsPoint[][] normals = null;
    transient private float[] planarNormal = null;
    transient boolean isPlanar;
    transient int nRows, nCols;
    transient float zTop, zBot, dZ, rdZ;
    transient StsColor stsColor;
    transient private int listNum = 0; // display list number (>0)
    transient private StsSection section = null;

    public static final int RIGHT = StsParameters.RIGHT;
    public static final int LEFT = StsParameters.LEFT;

    private static final int STRIP_LEFT = StsParameters.STRIP_BOT;
    private static final int STRIP_RIGHT = StsParameters.STRIP_TOP;

    private static final float[] verticalNormal = new float[] { 0.0f, 0.0f, 1.0f };

    static final boolean debug = false;

    public StsPatch(int nRows, int nCols, StsSection section)
    {
        this.nRows = nRows;
        this.nCols = nCols;
        this.section = section;
        isPlanar = section.isPlanar();

        points = new StsPoint[nRows][nCols];
        if(!isPlanar) normals = new StsPoint[nRows][nCols];
    }

    public StsPatch(int nRows, int nCols, StsColor color, StsSection section)
    {
        this(nRows, nCols, section);
        stsColor = color;
    }

    public void initialize(float zTop, float zBot, float dZ)
    {
        this.zTop = zTop;
        this.zBot = zBot;
        this.dZ = dZ;
        this.rdZ = 1.0f/dZ;
    }
    // Accessors
    public void setStsColor(StsColor color) { stsColor = color; }
    public StsColor getStsColor() { return stsColor; }

    private float[] getValues(int dim)
    {
        float[] values = new float[nCols*nRows];
        int n=0;
        for( int i=0; i<nCols; i++ )
            for( int j=0; j<nRows; j++, n++ )
                values[n] = points[j][i].getPointValues()[dim];
        return values;
    }

    public float[] getXValues() { return getValues(0); }
    public float[] getYValues() { return getValues(1); }
    public float[] getZValues() { return getValues(2); }
    public StsPoint[][] getPoints() { return points; }
    public StsPoint[][] getNormals() { return normals; }
    public int getNRows() { return nRows; }
    public int getNCols() { return nCols; }
    public float getZTop() { return points[0][0].getZ(); }
    public float getZBot() { return points[nRows-1][0].getZ(); }
    public StsSection getSection() { return section; }

    public StsPoint getStsPoint(int row, int col)
    {
        if(!inside(row, col)) return null;
        return points[row][col];
    }

    public float[] getXYZorT(int row, int col)
    {
        if(!inside(row, col)) return null;
        return points[row][col].getXYZorT();
    }
/*
    public boolean isFinalPointNull(int row, int col)
    {
        return false;
    }
*/

    /*
    public StsPatch(int nRows, int nCols, StsSurfaceEdge botEdge, StsSurfaceEdge topEdge, StsSection section)
        throws StsException
    {
        int nc, nr;
        float z;
        StsWell well;
        float fRow;
        StsPoint[] botTranslatedPoints, topTranslatedPoints;
        StsLineCoordinates botCoor, topCoor;
        this.nRows = nRows;
        this.nCols = nCols;

        int noPoints = nRows*nCols;

        points = new StsPoint[nRows][nCols];
        normals = new StsPoint[nRows][nCols];
        for (int i=0; i<nRows; i++)
        {
            for (int j=0; j<nCols; j++)
            {
                points[i][j] = null;
                normals[i][j] = null;
            }
        }

        // set points along top and bottom patch edges

        StsPoint[] botPoints = botEdge.getPoints();
        StsPoint[] topPoints = topEdge.getPoints();

        zBot = botPoints[0].getZ();
        zTop = topPoints[0].getZ();
        dZ = (zTop - zBot)/(nRows-1);
        rdZ = 1.0f/dZ;

        for(nc = 0; nc < nCols; nc++)
        {
            points[0][nc] = new StsPoint(botPoints[nc]);
            points[nRows-1][nc] = new StsPoint(topPoints[nc]);
        }

        // set points up left side

        well = section.getFirstWell();
        z = zBot;
        for(nr = 1; nr < nRows-1; nr++)
        {
            z += dZ;
            points[nr][0] = well.getPointAtZ(z, true);
        }

        // set points up rite side

        well = section.getLastWell();
        z = zBot;
        for(nr = 1; nr < nRows-1; nr++)
        {
            z += dZ;
            points[nr][nCols-1] = well.getPointAtZ(z, true);
        }

        // now interpolate interior points

        float dfCol = 1.0f/(nCols - 1);
        float dfRow = 1.0f/(nRows - 1);

	    botCoor = botEdge.computeCoordinates();
        topCoor = topEdge.computeCoordinates();

        fRow = dfRow;

        for(nr = 1; nr < nRows-1; nr++)
        {
            StsPoint leftPoint = points[nr][0];
            StsPoint ritePoint = points[nr][nCols-1];

            botTranslatedPoints = StsLineCoordinates.computeTranslatedPoints(botPoints, botCoor, leftPoint, ritePoint);
            topTranslatedPoints = StsLineCoordinates.computeTranslatedPoints(topPoints, topCoor, leftPoint, ritePoint);

            for(nc = 1; nc < nCols-1; nc++)
            {
                points[nr][nc] = new StsPoint();
                points[nr][nc].interpolatePoints(botTranslatedPoints[nc], topTranslatedPoints[nc], fRow);
            }

            fRow += dfRow;
        }

        makeNormals();
    }
    */

    public boolean inside(int row, int col)
    {
        return row >= 0 && row < nRows && col >= 0 && col < nCols;
    }

 	public void checkMakeNormals()
    {
        if(normals != null) return;
        makeNormals();
    }

    public void makeNormals()
    {
        if(isPlanar) makePlanarSurfaceNormal();
        else makeCurvedSurfaceNormals();
	}

    public void makePlanarSurfaceNormal()
    {
        StsPoint idif = StsPoint.subPointsStatic(points[0][nCols-1], points[0][0]);
        planarNormal = new float[] { idif.v[1], -idif.v[0], 0.0f };
	}

    public void makeCurvedSurfaceNormals()
    {
		int i, j;
		int im1, ip1, jm1, jp1;

        try
        {
            StsPoint idif = new StsPoint(3);
            StsPoint jdif = new StsPoint(3);

            for (i=0; i<nRows; i++)
            {
                im1 = i==0 ? i : i-1;
                ip1 = i==nRows-1 ? i : i+1;

                for (j=0; j<nCols; j++)
                {
                    jm1 = j==0 ? j : j-1;
                    jp1 = j==nCols-1 ? j : j+1;

                    idif.subPoints(points[ip1][j], points[im1][j]);
                    jdif.subPoints(points[i][jp1], points[i][jm1]);
                    if (normals[i][j]==null) normals[i][j] = new StsPoint(3);
                    normals[i][j].leftCrossProduct(idif, jdif);
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsPatch.makeNormals() failed.",
                e, StsException.WARNING);
        }
	}
/*
    public float[] getNormalZF(float z, float f)
    {
        float rowF = (zBot - z)*rdZ;

        f = StsMath.minMax(f, 0.0f, 1.0f);
        float colF = f*(nCols-1);

        return getNormal(rowF, colF);
    }
*/
    public float getZAtRow(int row)
    {
        return zTop + row*dZ;
    }
/*
    public float getFAtCol(int col)
    {
        return ((float)col)/(nCols-1);
    }
*/

    public float[] getXYZorT(float rowF, float colF)
    {
        StsPoint point = getStsPoint(rowF, colF);
        return point.getXYZorT();
    }
    public StsPoint getStsPoint(float rowF, float colF)
    {
        StsPoint point0, point1, point2;

        try
        {
            int row = (int)rowF;
			row = StsMath.minMax(row, 0, nRows-2);
            float dR = rowF - row;

            int col = (int)colF;
			col = StsMath.minMax(col, 0, nCols-2);
            float dC = colF - col;

            if(dR == 0.0f)
            {
                if(dC == 0.0f)
                    return points[row][col];
                else
                {
                    point0 = StsPoint.staticInterpolatePoints(points[row][col], points[row][col+1], dC);
                    return point0;
                }
            }
            else
            {
                point0 = StsPoint.staticInterpolatePoints(points[row][col], points[row+1][col], dR);
                if(dC == 0.0f)
                    return point0;
                else
                {
                    point1 = StsPoint.staticInterpolatePoints(points[row][col+1], points[row+1][col+1], dR);
                    point2 = StsPoint.staticInterpolatePoints(point0, point1, dC);
                    return point2;
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputException(e, StsException.WARNING);
            return null;
        }
    }

    public float[] getPoint(float rowF, float colF)
    {
		return getStsPoint(rowF, colF).v;
    }

    public float[] getNormal(float rowF, float colF)
    {
        StsPoint norm0, norm1, norm2;

        try
        {
            if(isPlanar) return planarNormal;

            rowF = StsMath.minMax(rowF, 0.0f, (float)(nRows-1));
            int row = (int)rowF;
            float dR = rowF - row;

            colF = StsMath.minMax(colF, 0.0f, (float)(nCols-1));
            int col = (int)colF;
            float dC = colF - col;

            if(dR == 0.0f)
            {
                if(dC == 0.0f)
                    return normals[row][col].v;
                else
                {
                    norm0 = StsPoint.staticInterpolatePoints(normals[row][col], normals[row][col+1], dC);
                    return norm0.v;
                }
            }
            else
            {
                norm0 = StsPoint.staticInterpolatePoints(normals[row][col], normals[row+1][col], dR);
                if(dC == 0.0f)
                    return norm0.v;
                else
                {
                    norm1 = StsPoint.staticInterpolatePoints(normals[row][col+1], normals[row+1][col+1], dR);
                    norm2 = StsPoint.staticInterpolatePoints(norm0, norm1, dC);
                    return norm2.v;
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsPatch.getNormal() failed for rowF " + rowF + " colF " + colF +
                " nRows " + nRows + " nCols " + nCols,
                e, StsException.WARNING);
            return verticalNormal;
        }
    }

    public float getRowF(float z)
    {
        return (float)(((double)z - (double)zTop)/(double)dZ);
    }

	public StsPoint[] getPointsAtZ(float z)
	{
		int row = -1, col = -1;
        try
        {
		    float rowF = getRowF(z);

            row = (int)rowF;
			row = StsMath.minMax(row, 0, nRows-2);
            float dR = rowF - row;

            if(dR == 0.0f)
                return points[row];
			else
			{
				StsPoint[] rowPoints = new StsPoint[nCols];
				for(col = 0; col < nCols; col++)
				{
                    StsPoint point = StsPoint.staticInterpolatePoints(points[row][col], points[row+1][col], dR);
                    rowPoints[col] = new StsPoint(point.getX(), point.getY(), z);
                }
				return rowPoints;
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsPatch.getPointsAtZ() failed at row " + row + " col " + col + ".",
				    e, StsException.WARNING);
            return null;
        }
    }

    public float[] getPlanarNormal() { return planarNormal; }

    public float[] getNormal(int row, int col)
    {
        if(!inside(row, col)) return null;
        return normals[row][col].v;
    }

    public void checkAddTimeOrDepth(StsSeismicVelocityModel velocityModel, boolean isDepth)
    {
        int row = -1, col = -1;
        try
        {
            if (isDepth)
            {
                for (row = 0; row < nRows; row++)
                {
                    for (col = 0; col < nCols; col++)
                    {
                        StsPoint point = points[row][col];
                        float z = point.getZ();
                        float t = (float)velocityModel.getT(point.v);
                        point.setT(t);
                    }
                }
            }
            else
            {
                for (row = 0; row < nRows; row++)
                {
                    for (col = 0; col < nCols; col++)
                    {
                        StsPoint point = points[row][col];
                        float t = point.getT();
                        float z = (float)velocityModel.getZ(point.getX(), point.getY(), point.getT());
                        point.setZ(z);
                    }
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsPatch.checkAddTimeOrDepth() failed for row " + row + " col " + col,
                                         e, StsException.WARNING);
        }
    }

	public void displaySurface(StsModel model, int side, boolean isDepth)
	{
        StsSeismicVolume sv;

        if (points == null) return;

        GL gl = model.getWin3dGL();
        if(gl == null) return;

        StsSeismicVolumeClass volumeClass = (StsSeismicVolumeClass)model.getStsClass(StsSeismicVolume.class);
        if((sv = volumeClass.getSeismicVolumeDisplayableOnSection()) != null)
            displaySeismicOnSection(model, gl, sv, isDepth);
        else
            displaySurface1(model, gl, side, isDepth);
    }

    private void displaySeismicOnSection(StsModel model, GL gl, StsSeismicVolume sv, boolean isDepth)
    {
        boolean useDisplayLists = model.useDisplayLists;

        if (useDisplayLists)
		{
            if (listNum == 0)  // build a new display lst
            {
                listNum = gl.glGenLists(1);
                if(listNum == 0)
                {
                    StsMessageFiles.logMessage("System Error in StsPatch.displaySeismicOnSection(): " +
                                           "Failed to allocate a display list");
                    return;
                }

                gl.glNewList(listNum, GL.GL_COMPILE_AND_EXECUTE);
                drawSeismicOnSection(gl, sv, isDepth);
                gl.glEndList();
            }
            else
			    gl.glCallList(listNum);

			// timer.stop("display list setup  1D: ");
		}
        else
        {
            if (listNum != 0)  // delete existing display list
            {
                gl.glDeleteLists(listNum, 1);
                listNum = 0;
            }
            drawSeismicOnSection(gl, sv, isDepth);
        }
    }

    private void drawSeismicOnSection(GL gl, StsSeismicVolume sv, boolean isDepth)
    {
        try
        {
            float zMin = sv.getZMin();
            float zMax = sv.getZMax();
            float dZ = zMax - zMin;

            float rowMinF = section.getRowF(zMin);
            float rowMaxF = section.getRowF(zMax);
            float dRow = dZ/section.getDZRow();

            int minPlane = 0;
            int maxPlane = sv.getNSlices()-1;

            if(rowMinF < 0.0)
            {
                float nRowsF = -rowMinF/dRow;
                int nRows = StsMath.ceiling(nRowsF);
                rowMinF += dRow*nRows;
                minPlane = nRows;
            }
            if(rowMaxF > nRows-1)
            {
                float nRowsF = (rowMaxF - (nRows-1))/dRow;
                int nRows = StsMath.ceiling(nRowsF);
                rowMaxF -= dRow*nRows;
                maxPlane -= nRows;
            }

            SeismicColPoint[] colPoints0, colPoints1;

            colPoints1 = getSeismicColPoints(0.0f, rowMinF, dRow, sv, minPlane, maxPlane);
            for(int col = 1; col < nCols; col++)
            {
                colPoints0 = colPoints1;
                colPoints1 = getSeismicColPoints((float)col, rowMinF, dRow, sv, minPlane, maxPlane);

                int nIncs = colPoints0[0].getNIncs(colPoints1[0]);

                if(nIncs <= 1)
                    drawSeismicOnSectionCol(colPoints0, colPoints1, gl, isDepth);
                else
                {
                    float df = 1.0f/nIncs;
                    float f = 0.0f;
                    SeismicColPoint[] colPointsA, colPointsB;
                    colPointsB = colPoints0;
                    for(int n = 1; n < nIncs; n++)
                    {
                        colPointsA = colPointsB;
                        f += df;
                        colPointsB = interpolateColPoints(colPoints0, colPoints1, f, sv);
                        drawSeismicOnSectionCol(colPointsA, colPointsB, gl, isDepth);
                    }
                    drawSeismicOnSectionCol(colPointsB, colPoints1, gl, isDepth);
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsPatch.drawSeismicOnSection() failed.",
                e, StsException.WARNING);
        }
    }

    private void drawSeismicOnSectionCol(SeismicColPoint[] colPoints0, SeismicColPoint[] colPoints1, GL gl, boolean isDepth)
    {
        SeismicColPoint colPoint;

        try
        {
            int nPoints = colPoints0.length;
            gl.glBegin(GL.GL_QUAD_STRIP);
            if(isDepth)
            {
                for (int n = nPoints - 1; n >= 0; n--)
                {
                    colPoint = colPoints0[n];
                    colPoint.stsColor.setGLColor(gl);
                    gl.glNormal3fv(colPoint.normal, 0);
                    gl.glVertex3fv(colPoint.point, 0);
                    colPoint = colPoints1[n];
                    colPoint.stsColor.setGLColor(gl);
                    gl.glNormal3fv(colPoint.normal, 0);
                    gl.glVertex3fv(colPoint.point, 0);
                }
            }
            else
            {
                for (int n = nPoints - 1; n >= 0; n--)
                {
                    colPoint = colPoints0[n];
                    colPoint.stsColor.setGLColor(gl);
                    gl.glNormal3fv(colPoint.normal, 0);
                    gl.glVertex3f(colPoint.point[0], colPoint.point[1], colPoint.point[4]);
                    colPoint = colPoints1[n];
                    colPoint.stsColor.setGLColor(gl);
                    gl.glNormal3fv(colPoint.normal, 0);
                    gl.glVertex3f(colPoint.point[0], colPoint.point[1], colPoint.point[4]);
                }
            }
            gl.glEnd();
        }
        catch(Exception e)
        {
            StsException.outputException("StsPatch.drawSeismicOnSectionCol() failed.",
                e, StsException.WARNING);
        }
    }

    private SeismicColPoint[] getSeismicColPoints(float colF, float rowMinF, float dRow,
                                                   StsSeismicVolume sv, int minPlane, int maxPlane)
    {
        float[] normal, point;

        int nPlanes = maxPlane - minPlane + 1;
        SeismicColPoint[] colPoints = new SeismicColPoint[nPlanes];

        float rowF = rowMinF;
        int n = 0;
        for(int plane = minPlane; plane <= maxPlane; plane++, n++)
        {
            SeismicColPoint colPoint = new SeismicColPoint();
            colPoint.normal = getNormal(rowF, colF);
            colPoint.point = getPoint(rowF, colF);

            colPoint.computeSeismicRowAndCol(sv);
            colPoint.plane = plane;
            colPoint.computeStsColor(sv);
            colPoints[n] = colPoint;
            rowF += dRow;
        }
        return colPoints;
    }

    private SeismicColPoint[] interpolateColPoints(SeismicColPoint[] colPoints0, SeismicColPoint[] colPoints1,
                                float f, StsSeismicVolume sv)
    {
        if(colPoints0 == null || colPoints1 == null) return null;
        int nPoints = colPoints0.length;
        SeismicColPoint[] colPoints = new SeismicColPoint[nPoints];
        for(int n = 0; n < nPoints; n++)
            colPoints[n] = interpolateColPoint(colPoints0[n], colPoints1[n], f, sv);
        return colPoints;
    }

    private SeismicColPoint interpolateColPoint(SeismicColPoint colPoint0, SeismicColPoint colPoint1,
                                float f, StsSeismicVolume sv)
    {
        SeismicColPoint colPoint = new SeismicColPoint();
        colPoint.point = StsMath.interpolate(colPoint0.point, colPoint1.point, f);
        colPoint.normal = StsMath.interpolate(colPoint0.normal, colPoint1.normal, f);
        colPoint.seismicRowF = colPoint0.seismicRowF + f*(colPoint1.seismicRowF - colPoint0.seismicRowF);
        colPoint.seismicColF = colPoint0.seismicColF + f*(colPoint1.seismicColF - colPoint0.seismicColF);
        colPoint.plane = colPoint0.plane;
        colPoint.computeStsColor(sv);
        return colPoint;
    }

    class SeismicColPoint
    {
        float[] point, normal;
        float seismicRowF, seismicColF;
        int plane;
        StsColor stsColor;

        SeismicColPoint()
        {
        }

        int getNIncs(SeismicColPoint otherPoint)
        {
            float dRowF = otherPoint.seismicRowF - seismicRowF;
            float dColF = otherPoint.seismicColF - seismicColF;
            float maxD = Math.max(Math.abs(dRowF), Math.abs(dColF));
            return StsMath.ceiling(maxD);
        }

        void computeSeismicRowAndCol(StsSeismicVolume sv)
        {
            seismicRowF = sv.getRowCoor(point[1]);
            seismicColF = sv.getColCoor(point[0]);
        }

        void computeStsColor(StsSeismicVolume sv)
        {
            stsColor = sv.getStsColor(seismicRowF, seismicColF, plane);
        }
    }
/*
    private void drawSeismicOnSection(GL gl, StsSeismicVolume sv)
    {
        float zMin = sv.getZMin();
        float zMax = sv.getZMax();
        float dZ = sv.getDZ();

        float rowMinF = section.getRowF(zMin);
        float rowMaxF = section.getRowF(zMax);
        float dRow = dZ/section.getDZRow();

        int minPlane = 0;
        int maxPlane = sv.getNSamples()-1;

        if(rowMinF < 0.0)
        {
            float nRowsF = -rowMinF/dRow;
            int nRows = StsMath.ceiling(nRowsF);
            rowMinF += dRow*nRows;
            minPlane = nRows;
        }
        if(rowMaxF > nRows-1)
        {
            float nRowsF = (rowMaxF - (nRows-1))/dRow;
            int nRows = StsMath.ceiling(nRowsF);
            rowMaxF -= dRow*nRows;
            maxPlane -= nRows;
        }

        int nPlanes = maxPlane - minPlane + 1;

        float colF0, colF1;
        float[] normal, point;
        int value;

        colF1 = 0.0f;
        for(int col = 1; col < nCols; col++)
        {
            colF0 = colF1;
            colF1 = (float)col;
            float rowF = rowMaxF;
//            float rowF = rowMinF;

            gl.glBegin(GL.GL_QUAD_STRIP);

            for(int n = 0; n < nPlanes; n++)
            {
                normal = getNormal(rowF, colF0);
                point = getPoint(rowF, colF0);
		        StsColor.setGLColor(gl, sv.getStsColor(point));
                gl.glNormal3fv(normal);
                gl.glVertex3fv(point);

                normal = getNormal(rowF, colF1);
                point = getPoint(rowF, colF1);
		        StsColor.setGLColor(gl, sv.getStsColor(point));
                gl.glNormal3fv(normal);
                gl.glVertex3fv(point);

                rowF -= dRow;
            }

            gl.glEnd();
        }
    }
*/
    private void displaySurface1(StsModel model, GL gl, int side, boolean isDepth)
    {
        if(stsColor != null)
		    stsColor.setGLColor(gl);
        else
		    StsColor.PURPLE.setGLColor(gl);

        boolean useDisplayLists = model.useDisplayLists;

        if (useDisplayLists)
		{
            if (listNum == 0)  // build a new display lst
            {
                listNum = gl.glGenLists(1);
                if(listNum == 0)
                {
                    StsMessageFiles.logMessage("System Error in StsPatch.display: " +
                                           "Failed to allocate a display list");
                    return;
                }

                gl.glNewList(listNum, GL.GL_COMPILE);
                drawSurface(gl, isDepth);
                gl.glEndList();
            }
            drawDisplayList(gl, listNum, side, isDepth);

			// timer.stop("display list setup  1D: ");
		}
        else
        {
            if (listNum != 0)  // delete existing display list
            {
                gl.glDeleteLists(listNum, 1);
                listNum = 0;
            }
            drawSurface(gl, side, isDepth);
        }
    }

    private void drawSurface(GL gl, int side, boolean isDepth)
    {
        if(side == 0)
        {
            drawSurface(gl, isDepth);
        }
        else
        {
            gl.glEnable(GL.GL_CULL_FACE);

            if(side == RIGHT) gl.glCullFace(GL.GL_BACK);
            else          gl.glCullFace(GL.GL_FRONT);

            drawSurface(gl, isDepth);

            gl.glEnable(GL.GL_POLYGON_STIPPLE);
            gl.glPolygonStipple(StsGraphicParameters.halftone, 0);

            if(side == RIGHT) gl.glCullFace(GL.GL_FRONT);
            else          gl.glCullFace(GL.GL_BACK);

            drawSurface(gl, isDepth);

            gl.glDisable(GL.GL_CULL_FACE);
            gl.glDisable(GL.GL_POLYGON_STIPPLE);
        }
    }

    private void drawDisplayList(GL gl, int listNum, int side, boolean isDepth)
    {
//    	System.out.println("drawDisplayList for section side " + side +
//        					"\nright = " + RIGHT +
//        					"\tleft = " + LEFT);
        if(side == 0)
        {
            drawSurface(gl, isDepth);
        }
        else
        {
            gl.glEnable(GL.GL_CULL_FACE);

            if(side == RIGHT) gl.glCullFace(GL.GL_BACK);
            else          gl.glCullFace(GL.GL_FRONT);

			gl.glCallList( listNum );

            gl.glEnable(GL.GL_POLYGON_STIPPLE);
            gl.glPolygonStipple(StsGraphicParameters.halftone, 0);

            if(side == RIGHT) gl.glCullFace(GL.GL_FRONT);
            else          gl.glCullFace(GL.GL_BACK);

			gl.glCallList( listNum );

            gl.glDisable(GL.GL_CULL_FACE);
            gl.glDisable(GL.GL_POLYGON_STIPPLE);
        }
    }

    private void drawSurface(GL gl, boolean isDepth)
    {
        if (section == null) return;
        for (int i = 0; i < nRows-1; i++)
            drawRowQuadStrip(gl, i, 0, nCols, isDepth);
    }

    public void drawRowQuadStrip(GL gl, int irow, int beginCol, int endCol, boolean isDepth)
    {
        gl.glBegin(GL.GL_QUAD_STRIP);

        if(!isPlanar)
        {
            for (int j = beginCol; j < endCol; j++)
            {
                gl.glNormal3fv(normals[irow][j].v, 0);
                gl.glVertex3fv(points[irow][j].v, 0);
                gl.glVertex3fv(points[irow + 1][j].v, 0);
            }
        }
        else
        {
            for (int j = beginCol; j < endCol; j++)
            {
                gl.glNormal3fv(planarNormal, 0);
                gl.glVertex3fv(points[irow][j].v, 0);
                gl.glVertex3fv(points[irow + 1][j].v, 0);
            }
        }
        gl.glEnd();
    }

    public void drawRowQuadCell(GL gl, int irow, int jcol)
    {
        gl.glBegin(GL.GL_QUAD_STRIP);

        if(!isPlanar)
        {
            gl.glNormal3fv(normals[irow][jcol].v, 0);
            gl.glVertex3fv(points[irow][jcol].v, 0);
            gl.glNormal3fv(normals[irow+1][jcol].v, 0);
            gl.glVertex3fv(points[irow+1][jcol].v, 0);

            gl.glNormal3fv(normals[irow][jcol+1].v, 0);
            gl.glVertex3fv(points[irow][jcol+1].v, 0);
            gl.glNormal3fv(normals[irow+1][jcol+1].v, 0);
            gl.glVertex3fv(points[irow+1][jcol+1].v, 0);
        }
        else
        {
            gl.glNormal3fv(planarNormal, 0);

            gl.glVertex3fv(points[irow][jcol].v, 0);
            gl.glVertex3fv(points[irow+1][jcol].v, 0);

            gl.glVertex3fv(points[irow][jcol+1].v, 0);
            gl.glVertex3fv(points[irow+1][jcol+1].v, 0);
        }

        gl.glEnd();
    }

	public void displayGrid(StsGLPanel3d glPanel3d)
    {
		GL gl = glPanel3d.getGL();
        if(gl == null) return;

        gl.glDisable(GL.GL_LIGHTING);
        glPanel3d.setViewShift(gl, StsGraphicParameters.gridShift);

        gl.glLineWidth(StsGraphicParameters.gridLineWidth);

        StsColor gridColor = StsColor.BLACK;
        if(!section.getDrawSurface()) gridColor = stsColor;
        StsGLDraw.drawGridLines(gl, gridColor, StsGraphicParameters.gridLineWidth,
                      points, nRows, nCols); // , 4, 4);

        gl.glEnable(GL.GL_LIGHTING);
        glPanel3d.resetViewShift(gl);
    }

    public void deleteDisplayLists(GL gl)
    {
        if(listNum != 0)
        {
            gl.glDeleteLists(listNum, 1);
            listNum = 0;
        }
    }

    public void insertPatchInPatch(StsPatch childPatch,
                    int parentRowNo, int parentColNo) throws StsException
    {
        int cr, cc, pr, pc;

        if(childPatch == null)
            throw new StsException(StsException.WARNING, "StsPatch.insertPatchInPatch(...)",
                                   "ChildPatch is NULL");

        StsPoint[][] childPoints = childPatch.getPoints();
        int childNRows = childPatch.getNRows();
        int childNCols = childPatch.getNCols();

        for(cr = 0, pr = parentRowNo; cr < childNRows; cr++, pr++)
            for(cc = 0, pc = parentColNo; cc < childNCols; cc++, pc++)
                this.points[pr][pc] = childPoints[cr][cc];
    }

// Picking routines

    /** This routine draws all grid quads in picking */

    public double[] getPointOnPatch(StsGLPanel3d glPanel3d, StsMousePoint mousePoint)
    {
        try
        {
    	    StsMethod pickMethod = new StsMethod(this, "pickOnPatch", glPanel3d);
            if(StsJOGLPick.pick3d(glPanel3d, pickMethod, StsMethodPick.PICKSIZE_SMALL, StsMethodPick.PICK_CLOSEST))
            {
                StsPickItem pickItem = StsJOGLPick.pickItems[0];
                int rowNo = pickItem.names[0];
                int colNo = pickItem.names[1];
                int triNo = pickItem.names[2];
                return getMouseGridIntersect(glPanel3d, mousePoint, rowNo, colNo, triNo);
            }
            else
                return null;
        }
        catch(Exception e)
        {
            StsException.outputException("StsPatch.getPointOnPatch() failed.", e, StsException.WARNING);
            return null;
        }

	}

    public void pickOnPatch(StsGLPanel3d glPanel3d)
    {
    	StsPoint p00, p01, p10, p11;

        GL gl = glPanel3d.getGL();
        if(gl == null) return;

        /** draw two triangles for each quad: nameList has row & col of lower left
          *	of quad. Lower left triangle is 0 and upper right is 1.
          */

        for(int i = 0; i < nRows-1; i++)
        {
            p01 = points[i][0];
            p11 = points[i+1][0];

            for(int j = 0; j < nCols-1; j++)
            {
                p00 = p01;
                p10 = p11;

                p01 = points[i][j+1];
                p11 = points[i+1][j+1];

                gl.glInitNames();
                gl.glPushName(i);
                gl.glPushName(j);
                gl.glPushName(0);

                gl.glBegin(GL.GL_TRIANGLE_STRIP);
                gl.glVertex3fv(p00.v, 0);
                gl.glVertex3fv(p10.v, 0);
                gl.glVertex3fv(p01.v, 0);
                gl.glEnd();

                gl.glPopName();
                gl.glPushName(1);

                gl.glBegin(GL.GL_TRIANGLE_STRIP);
                gl.glVertex3fv(p10.v, 0);
                gl.glVertex3fv(p01.v, 0);
                gl.glVertex3fv(p11.v, 0);
                gl.glEnd();

                gl.glPopName();
                gl.glPopName();
            }
        }
    }

	private double[] getMouseGridIntersect(StsGLPanel3d glPanel3d, StsMousePoint mousePoint, int rowNo, int colNo, int triNo)
    {
        if(glPanel3d == null) return null;

    	StsPoint[] triPoints = new StsPoint[3];

        if(triNo == 0)
        {
            triPoints[0] = points[rowNo][colNo];
            triPoints[1] = points[rowNo+1][colNo];
            triPoints[2] = points[rowNo][colNo+1];
        }
        else
        {
            triPoints[0] = points[rowNo+1][colNo];
            triPoints[1] = points[rowNo][colNo+1];
            triPoints[2] = points[rowNo+1][colNo+1];
        }

        return glPanel3d.getMouseTriangleIntersect(mousePoint, triPoints);
   }

    /** This method divides area to be searched into 8  by 8 sub-quads and continues to
      * subdivide until single grid cell is found */

    /** Pick limits on patch */

    int rowMin, rowMax, colMin, colMax, rowInc, colInc;

    public StsPoint getPointOnPatchSearch(StsGLPanel3d glPanel3d, StsMousePoint mousePoint)
    {
        StsPickItem[] pickItems;

        try
        {
            int rowNo = 0, colNo = 0, triNo = 0;
            int i;

            StsMethod pickMethod = new StsMethod(this, "pickDecimatedPatch", glPanel3d);
            StsJOGLPick picker = new StsJOGLPick(glPanel3d, pickMethod, StsJOGLPick.PICKSIZE_SMALL, StsJOGLPick.PICK_CLOSEST);

            int nrm1 = nRows - 1;
            int ncm1 = nCols - 1;

            rowMin = 0;
            rowMax = nrm1;
            colMin = 0;
            colMax = ncm1;

            while(true)
            {
                if(debug)
                {
                    System.out.println("rows: " + rowMin + " " + rowMax + " " + rowInc);
                    System.out.println("cols: " + colMin + " " + colMax + " " + colInc);
                }
                if(!picker.methodPick3d()) return null;

                if(rowInc > 1 || colInc > 1)
                {
                    rowMin = nrm1 - 1;
                    colMin = ncm1 - 1;
                    rowMax = 0;
                    colMax = 0;

                    pickItems = picker.pickItems;

                    for(i = 0; i < pickItems.length; i++)
                    {
                        int row = pickItems[i].names[0];
                        int col = pickItems[i].names[1];

                        if(debug) System.out.println("hit: " + i + " row: " + row + " col: " + col);

                        rowMin = Math.min(row, rowMin);
                        colMin = Math.min(col, colMin);
                        rowMax = Math.max(row + rowInc, rowMax);
                        colMax = Math.max(col + colInc, colMax);
                    }
                    rowMax = Math.min(rowMax, nRows-1);
                    colMax = Math.min(colMax, nCols-1);
                }
                else
                {
                    pickItems = picker.pickItems;

                    float zMin = StsParameters.largeFloat;

                    for(i = 0; i < pickItems.length; i++)
                    {
                        if(pickItems[i].zMin < zMin)
                        {
                            rowNo = pickItems[i].names[0];
                            colNo = pickItems[i].names[1];
                            triNo = pickItems[i].names[2];
                        }
                    }

                    if(debug) System.out.println("Picked row: " + rowNo + " col: " + colNo + " tri: " + triNo);

                    double[] xyz = getMouseGridIntersect(glPanel3d, mousePoint, rowNo, colNo, triNo);
                    StsPoint point = new StsPoint(5);
                    point.setX((float)xyz[0]);
                    point.setY((float)xyz[1]);
                    point.setZorT((float)xyz[2]);
                    return point;
                }
            }
        }
        catch (Exception e) { return null; }
    }

    public void pickDecimatedPatch(StsGLPanel3d glPanel3d)
    {
        int dRows, dCols;
        int i0, i1, j0, j1;
        StsPoint p00, p01, p10, p11;
        int i, j;
        int nDecRows = 8, nDecCols = 8;

        dRows = rowMax - rowMin;
        rowInc = Math.max(1, dRows/nDecRows);
        nDecRows = dRows/rowInc;
        if(nDecRows*rowInc < dRows) nDecRows++;

        dCols = colMax - colMin;
        colInc = Math.max(1, dCols/nDecCols);
        nDecCols = dCols/colInc;
        if(nDecCols*colInc < dCols) nDecCols++;

		GL gl = glPanel3d.getGL();
        if(gl == null) return;

        i1 = rowMin;

        for(i = 0; i < nDecRows; i++)
        {
            i0 = i1;
            i1 = Math.min(rowMax, i0 + rowInc);

            gl.glInitNames();
            gl.glPushName(i0);

            p01 = points[i0][colMin];
            p11 = points[i1][colMin];

            j1 = colMin;

            for(j = 0; j < nDecCols; j++)
            {
                p00 = p01;
                p10 = p11;

                j0 = j1;
                j1 = Math.min(colMax, j0 + colInc);

                p01 = points[i0][j1];
                p11 = points[i1][j1];

                gl.glPushName(j0);
                gl.glPushName(0);

                gl.glBegin(GL.GL_TRIANGLE_STRIP);
                gl.glVertex3fv(p00.v, 0);
                gl.glVertex3fv(p10.v, 0);
                gl.glVertex3fv(p01.v, 0);
                gl.glEnd();


                gl.glPopName();
                gl.glPushName(1);

                gl.glBegin(GL.GL_TRIANGLE_STRIP);
                gl.glVertex3fv(p10.v, 0);
                gl.glVertex3fv(p01.v, 0);
                gl.glVertex3fv(p11.v, 0);
                gl.glEnd();

                gl.glPopName();
                gl.glPopName();
            }
        }
    }
}
