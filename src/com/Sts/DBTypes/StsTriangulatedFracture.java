package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.StsCursor3d;
import com.Sts.Types.StsBoundingBox;
import com.Sts.Types.StsGolderFracture;
import com.Sts.Types.StsPoint;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.Utilities.*;
import com.Sts.Utilities.Triangulation.*;

import javax.media.opengl.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jun 4, 2008
 * Time: 8:04:29 PM
 *
 * This class takes a set of 3D points, constructs a least-squares plane thru them, projects the points to the plane,
 * and Delaunay triangulates these points.  The convex hull of this triangulation is a polygon which can then be
 * drawn to display the plane.
 */
public class StsTriangulatedFracture extends StsMainObject implements StsTreeObjectI
{
    transient private final StsColor stsColor = StsColor.RED;
    transient private double zMin;
    transient private double zMax;
    transient private double dipAngle;
    transient private double azimuth;

	transient double[][] points;
    transient int[][] triangleIndexes;
    transient double[] normal;
	transient double[][] polygonPoints;
	transient StsLstSqPointsPlane plane;
	transient float sliceZ;
	transient double[][] slicePoints = null;
	transient StsBoundingBox polygonBoundingBox;

    static final double[] vertical = new double[]{0.0, 0.0, -1.0};

    static public StsFieldBean[] displayFields = null;
    static StsDoubleFieldBean zMinBean;
    static StsDoubleFieldBean zMaxBean;
    static StsDoubleFieldBean dipAngleBean;
    static StsDoubleFieldBean azimuthBean;

    static protected StsObjectPanel objectPanel = null;
	static StsTriangulatedFractureClass triangulatedFractureClass = null;

	static final boolean debug = false;

    public StsTriangulatedFracture()
    {
    }

    private StsTriangulatedFracture(StsLstSqPointsPlane plane, String name, StsColor color)
    {
        super(false);
        setName(name);
        // stsColor = color;
		initializePlane(plane);
        addToModel();
    }

	private StsTriangulatedFracture(double[][] xyzPoints, String name, StsColor color)
	{
		super(false);
		setName(name);
		// stsColor = color;
		initializePoints(xyzPoints);
		addToModel();
	}

    static public StsTriangulatedFracture constructor(StsLstSqPointsPlane plane, String name, StsColor color)
    {
        try
        {
            return new StsTriangulatedFracture(plane, name, color);
        }
        catch (Exception e)
        {
            StsException.systemError(StsTriangulatedFracture.class, "constructor");
            return null;
        }
    }

	static public StsTriangulatedFracture constructor(double[][] xyzPoints, String name, StsColor color)
	{
		try
		{
			return new StsTriangulatedFracture(xyzPoints, name, color);
		}
		catch (Exception e)
		{
			StsException.systemError(StsTriangulatedFracture.class, "constructor");
			return null;
		}
	}
	public void initializePoints(double[][] xyzPoints)
	{
		initializePlane(new StsLstSqPointsPlane(xyzPoints));
	}
	public void initializePlane(StsLstSqPointsPlane plane)
	{
		this.plane = plane;
        triangulate(plane);
        computeRotatedPoints();
	}

    public void checkCreateFractureSetClass()
    {
        if(triangulatedFractureClass == null)
            triangulatedFractureClass = (StsTriangulatedFractureClass)currentModel.getCreateStsClass(StsTriangulatedFracture.class);
    }

	/** Given a least-squares plane thru the points, compute the normal axis direction and the angle to rotate this
	 *  axis to the vertical.  Apply this rotation to all the points which rotates the plane to the horizontal where
	 *  the x and y coordinates of the rotated points can be used as the coordinates for the Delaunay triangulation.
	 *
	 * @param plane a least-squares plane thru the points
	 */
    private void triangulate(StsLstSqPointsPlane plane)
    {
        normal = plane.getNormal();
        double[] axis = StsMath.cross(normal, vertical);
        setAzimuth(StsMath.atan2d(axis[1], axis[0]));
        if(debug) StsToolkit.print(axis, "Axis: ");
        double axisLen = StsMath.length(axis);
        float angleRad = (float) Math.asin(axisLen);
        float angle = (float) StsMath.DEGperRAD * angleRad;
        if(debug) System.out.println("Rotation angle: " + angle);
        dipAngle = Math.abs(90 - angle);
        StsQnion q = new StsQnion(axis, angleRad);
        points = plane.getPoints();
        int nPoints = points.length;
		// compute the points on the plane rotated to the horizontal
        double[][] rotatedPoints = new double[nPoints][];
        for (int n = 0; n < nPoints; n++)
            rotatedPoints[n] = q.leftRotateVec(points[n]);
		// tiangulate these points; this triangulation can then be applied to the points on the plane (plane.points)
		// or the triangulation could be applied to the original points (which are in the calling class)
        StsDelaunayTriangulation triangulation = new StsDelaunayTriangulation(rotatedPoints);
        triangleIndexes = triangulation.triangleIndexes;
		//polygonPoints = getHullPoints(triangulation.hullPointIndexes);      // Fails to produce a plane in most cases (Demo dataset fails on 3 of 4 stages)
        polygonPoints = getHullPoints(triangulation.nextHullPointIndexes);    // Produces a reasonable result so using it until statement above can be fixed.
		normal = triangulation.normal;
        computeZRange();
		// for debug purposes rotate all the points to the horizontal plare and reOrigin with respect to first point
		if(debug)
		{
			//StsLstSqPointsPlane rotatedPlane = StsLstSqPointsPlane.rotatePointsToHorizPlane(points);
			//double[] centroid = rotatedPlane.centroid;
			//double[][] points = rotatedPlane.points;
			for(double[] point : points)
			{
			//	point[0] -= centroid[0];
			//	point[1] -= centroid[1];
			//	point[2] -= centroid[2];
				System.out.println(point[0] + " " + point[1] + " " + point[2]);
			}
		}
    }

	private double[][] getHullPoints(int[] hullPointIndexes)
	{
		int nHullPoints = hullPointIndexes.length;
		double[][] hullPoints = new double[nHullPoints][];
		for(int n = 0; n < nHullPoints; n++)
			hullPoints[n] = points[hullPointIndexes[n]];
		return hullPoints;
	}

    private void computeRotatedPoints()
    {
        StsProject project = currentModel.getProject();
        for (int n = 0; n < points.length; n++)
        {
            float[] xy = project.getRotatedRelativeXYFromUnrotatedRelativeXY((float)points[n][0], (float)points[n][1]);
            points[n][0] = xy[0];
            points[n][1] = xy[1];
        }
    }

    private void computeZRange()
    {
        if (points == null) return;
        zMin = points[0][2];
        zMax = zMin;

        for (int n = 1; n < points.length; n++)
        {
            double z = points[n][2];
            if (z < zMin) zMin = z;
            else if (z > zMax) zMax = z;
        }
    }

    public void displayPolygons(StsGLPanel glPanel)
    {
		GL gl = glPanel.getGL();
		drawBestFitPolygon(gl);
		// drawStimulatedArea(gl);
	}

	private void drawBestFitPolygon(GL gl)
	{
		if(polygonPoints == null) return;
		try
		{
			gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL);
			StsColor.RED.setGLColor(gl);
			gl.glBegin(GL.GL_POLYGON);
			gl.glNormal3dv(normal, 0);
			for (int n = 0; n < polygonPoints.length; n++)
				gl.glVertex3dv(polygonPoints[n], 0);

				//gl.glLightModeli(GL.GL_LIGHT_MODEL_TWO_SIDE, GL.GL_FALSE);
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "display.", e);
        }
        finally
        {
            gl.glEnd();
        }
    }
/*
	private void drawStimulatedArea(GL gl)
	{
	    float zCoor = currentModel.win3d.getCursor3d().getCurrentDirCoordinate(StsCursor3d.ZDIR);
		if(zCoor != sliceZ || slicePoints == null)
		{
			computeSliceIntersects(zCoor);
	    	if (slicePoints == null) return;
			sliceZ = zCoor;
		}
		try
		{
			gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL);
			StsColor.BLACK.setGLColor(gl);
			gl.glBegin(GL.GL_LINES);
			int nPoints = slicePoints.length;
			for (int n = 0; n < nPoints; n++)
				gl.glVertex3dv(slicePoints[n], 0);
			gl.glEnd();
			int scale = 100;
			for(int i=0; i<nPoints; i++)
			{
				if(points[i] == null) continue;
				// scale = getAreaScale();
				float[] xyz = StsMath.convertDoubleToFloatArray(slicePoints[i]);
				StsGLDraw.drawFilledEllipse(xyz, StsColor.BLUE, currentModel.win3d.getGlPanel3d(), scale, scale, 0.0f, 2.0f);
			}
			int nLines = nPoints/2;
			int i = 0;
			for(int n = 0; n < nLines; n++)
			{
				double[] p0 = slicePoints[i++];
				double[] p1 = slicePoints[i++];
				double[] dp = StsMath.subtract(p1, p0);
				StsMath.normalize(dp);
				double[] xyOffset = new double[] { scale*dp[1], -scale*dp[0] };

			}
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "display.", e);
        }
        finally
        {
            gl.glEnd();
        }
	}
*/
	public double[][] computeLineIntersectAtZ(float z)
	{
		int nPolygonPoints = polygonPoints.length;
		ArrayList<double[]> slicePointsList = new ArrayList<double[]>();
		double[] nextPoint = polygonPoints[nPolygonPoints-1];
		double nextZ =nextPoint[2];
		for(int n = 0; n < nPolygonPoints; n++)
		{
			double[] prevPoint = nextPoint;
			double prevZ = nextZ;
			nextPoint = polygonPoints[n];
			nextZ = nextPoint[2];
			double f;
			if(nextZ == prevZ)
				f = 0.5;
			else
			{
				f =  (z - prevZ)/(nextZ - prevZ);
				if(f < 0.0 || f > 1.0) continue;
				if(f == 0.0) f = 0.001;
				else if(f == 1.0) f = 0.999;
			}
			double[] slicePoint = StsMath.interpolate(prevPoint, nextPoint, f);
			slicePointsList.add(slicePoint);
		}
		if(slicePointsList.size() > 2) orderSlicePoints(slicePointsList);
		return slicePointsList.toArray(new double[0][]);
	}

	public StsPoint[] getLineIntersectionAtZ(float z)
	{
		double[][] slicePoints = computeLineIntersectAtZ(z);
		if(slicePoints.length < 2) return null;
		StsPoint firstPoint = new StsPoint(slicePoints[0]);
		StsPoint lastPoint = new StsPoint(slicePoints[slicePoints.length-1]);
		return new StsPoint[] { firstPoint, lastPoint };
	}
	private void orderSlicePoints(ArrayList<double[]> pointsList)
	{
		StsBoundingBox box = new StsBoundingBox(pointsList.toArray(new double[0][]));
		if(box.getXSize() > box.getYSize())
			orderPoints(pointsList, 0);
		else
			orderPoints(pointsList, 1);
		int nPoints = pointsList.size();
		if(!StsMath.isEven(nPoints))
			pointsList.remove(nPoints/2);
	}

	private void orderPoints(ArrayList<double[]> pointsList, int nCoor)
	{
		Comparator<double[]> comparator = new ArrayComparatorDouble(nCoor);
		Collections.sort(pointsList, comparator);
	}

	class ArrayComparatorDouble implements Comparator<double[]>
	{
		int coor;

		ArrayComparatorDouble(int coor) { this.coor = coor; }

		public int compare(double[] point, double[] other)
		{
			return Double.compare(point[coor], other[coor]);
		}
	}

	private void drawStimulatedPolygon(GL gl)
	{
	    float zCoor = currentModel.win3d.getCursor3d().getCurrentDirCoordinate(StsCursor3d.ZDIR);
		if(zCoor != sliceZ || slicePoints == null)
		{
			computeSlicePolygonIntersects(zCoor);
	    	if (slicePoints == null) return;
			sliceZ = zCoor;
		}
		try
		{
			gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL);
			StsColor.RED.setGLColor(gl);
			gl.glBegin(GL.GL_POLYGON);
			gl.glNormal3dv(vertical, 0);
			for (int n = 0; n < slicePoints.length; n++)
				gl.glVertex3dv(slicePoints[n], 0);

				//gl.glLightModeli(GL.GL_LIGHT_MODEL_TWO_SIDE, GL.GL_FALSE);
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "display.", e);
        }
        finally
        {
            gl.glEnd();
        }
	}

	private void computeSlicePolygonIntersects(float z)
	{
		slicePoints = new double[3*triangleIndexes.length][];
		int[][] intersectEdges = new int[3*triangleIndexes.length][2];
		int nSlicePoints = 0;
		for(int[] triangleIndex : triangleIndexes)
		{
			int nextIndex = triangleIndex[2];
			double[] nextPoint = points[nextIndex];
			double nextZ = nextPoint[2];
			for(int i = 0; i < 3; i++)
			{
				int prevIndex = nextIndex;
				double[] prevPoint = nextPoint;
				double prevZ = nextZ;
				nextIndex = triangleIndex[i];
				nextPoint = points[nextIndex];
				nextZ = nextPoint[2];
				if(nextIndex < prevIndex) continue;
				double f;
				if(nextZ == prevZ)
					f = 0.5;
				else
				{
					f =  (z - prevZ)/(nextZ - prevZ);
					if(f < 0.0 || f > 1.0) continue;
				}
				int[] intersectEdge = new int[] { prevIndex, nextIndex };
				if(hasIntersectEdge(intersectEdge, intersectEdges, nSlicePoints)) continue;
				double[] slicePoint = StsMath.interpolate(prevPoint, nextPoint, f);
				slicePoints[nSlicePoints] = slicePoint;
				intersectEdges[nSlicePoints++] = intersectEdge;
			}
		}
		if(nSlicePoints > 2)
		{
			slicePoints = (double[][])StsMath.trimArray(slicePoints, nSlicePoints);
			StsDelaunayTriangulation sliceTriangulation = new StsDelaunayTriangulation(slicePoints);
			slicePoints = getHullPoints(sliceTriangulation.hullPointIndexes);
		}
	}

	public boolean hasIntersectEdge(int[] intersectEdge, int[][] intersectEdges, int nIntersectEdges)
	{
		int i0 = intersectEdge[0];
		int i1 = intersectEdge[1];
		for(int[] edge : intersectEdges)
		{
			if(i0 == edge[0] && i1 == edge[1]) return true;
			if(i1 == edge[0] && i0 == edge[1]) return true;
		}
		return false;
	}

    public StsPoint[] getSectionIntersectionAtZ(float zCoor)
    {
    	return null;
    }

    public boolean initialize(StsModel model)
    {
        return true;
    }

    public StsFieldBean[] getDisplayFields()
    {
        try
        {
            if (displayFields == null)
            {
                displayFields = new StsFieldBean[]
                    {
                        zMinBean = new StsDoubleFieldBean(StsTriangulatedFracture.class, "zMin", true, "Depth Min"),
                        zMaxBean = new StsDoubleFieldBean(StsTriangulatedFracture.class, "zMax", true, "Depth Max"),
                        dipAngleBean = new StsDoubleFieldBean(StsTriangulatedFracture.class, "dipAngle", true, "Dip Angle"),
                        azimuthBean = new StsDoubleFieldBean(StsTriangulatedFracture.class, "azimuth", true, "Azimuth")
                    };
            }
            return displayFields;
        }
        catch (Exception e)
        {
            StsException.outputException("StsSensor.getDisplayFields() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public StsFieldBean[] getPropertyFields()
    {
        return null;
    }

    public Object[] getChildren()
    {
        return new Object[0];
    }

    public boolean anyDependencies()
    {
        return false;
    }

    public StsObjectPanel getObjectPanel()
    {
        if (objectPanel == null)
        {
            objectPanel = StsObjectPanel.constructor(this, true);
        }
        return objectPanel;
    }

    public void treeObjectSelected()
    {
    }

    public double getZMin()
    {
        return zMin;
    }

    public void setZMin(double zMin)
    {
        this.zMin = zMin;
    }

    public double getZMax()
    {
        return zMax;
    }

    public void setZMax(double zMax)
    {
        this.zMax = zMax;
    }

    public double getDipAngle()
    {
        return dipAngle;
    }

    public void setDipAngle(double dipAngle)
    {
        this.dipAngle = dipAngle;
    }

    public double getAzimuth()
    {
        return azimuth;
    }

    public void setAzimuth(double azimuth)
    {
        this.azimuth = azimuth;
    }
	public boolean intersects(StsGolderFracture fracture)
	{
		StsPoint[][] fractureLines = fracture.getTopAndBottomLines();
		StsPoint[] topLine = fractureLines[0];
		StsPoint[] botLine = fractureLines[1];
		if(intersectsLine(topLine)) return true;
		return intersectsLine(botLine);
	}

	public boolean intersectsLine(StsPoint[] line)
	{
		StsPoint[] sectionEdgePoints = getLineIntersectionAtZ(line[0].getZ());
		if(sectionEdgePoints == null) return false;
		StsPoint point1 = sectionEdgePoints[0];
		for(int n = 1; n < sectionEdgePoints.length; n++)
		{
			StsPoint point0 = point1;
			point1 = sectionEdgePoints[n];
			if(StsMath.lineIntersectXY(line[0], line[1], point0, point1)) return true;
		}
		return false;
	}
}
