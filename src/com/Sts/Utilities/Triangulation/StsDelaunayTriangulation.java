package com.Sts.Utilities.Triangulation;

import com.Sts.IO.StsAsciiFile;
import com.Sts.IO.StsFile;
import com.Sts.Types.StsBoundingBox;
import com.Sts.Utilities.*;

import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
  * User: Tom Lasseter
  * Date: Jun 5, 2008
  * Time: 6:07:49 AM
  * <p/>
  * ported from p bourke's triangulate.c
  * http://astronomy.swin.edu.au/~pbourke/terrain/triangulate/triangulate.c
  * <p/>
  * fjenett, 20th february 2005, offenbach-germany.
  * contact: http://www.florianjenett.de/
  * <p/>
  */


 public class StsDelaunayTriangulation
 {
     double[][] points;
     int nPoints;
	 int nTotalPoints;
     public int[][] triangleIndexes;
     public int[][] edgeIndexes;
	 public int nHullPoints = 0;
     private int nTriangles;
	 public int[] nextHullPointIndexes;
	 public int[] hullPointIndexes;
	 // public double[][] hullPoints;
	 public double[] normal;

     static final double EPSILON = 0.000001;
	 static final boolean debug = false;

     /*
         Return TRUE if a point (xp,yp) is inside the circumcircle made up
         of the points (x1,y1), (x2,y2), (x3,y3)
         The circumcircle centre is returned in (xc,yc) and the radius r
         NOTE: A point on the edge is inside the circumcircle
     */

     public StsDelaunayTriangulation(double[][] initialPoints)
     {
         nPoints = initialPoints.length;
         nTotalPoints = nPoints + 3;
         points = new double[nTotalPoints][3];
         System.arraycopy(initialPoints, 0, points, 0, nPoints);
         triangulate();
         triangleIndexes = (int[][]) StsMath.trimArray(triangleIndexes, nTriangles);
     }

     boolean circumCircle(double xp, double yp, double x1, double y1, double x2, double y2, double x3, double y3, double[] circle)
     {
         double m1, m2, mx1, mx2, my1, my2;
         double dx, dy, rsqr, drsqr;
         double xc, yc, r;

         /* Check for coincident points */

         if (Math.abs(y1 - y2) < EPSILON && Math.abs(y2 - y3) < EPSILON)
         {
             System.out.println("CircumCircle: Points are coincident.");
             return false;
         }

         if (Math.abs(y2 - y1) < EPSILON)
         {
             m2 = -(x3 - x2) / (y3 - y2);
             mx2 = (x2 + x3) / 2.0;
             my2 = (y2 + y3) / 2.0;
             xc = (x2 + x1) / 2.0;
             yc = m2 * (xc - mx2) + my2;
         }
         else if (Math.abs(y3 - y2) < EPSILON)
         {
             m1 = -(x2 - x1) / (y2 - y1);
             mx1 = (x1 + x2) / 2.0;
             my1 = (y1 + y2) / 2.0;
             xc = (x3 + x2) / 2.0;
             yc = m1 * (xc - mx1) + my1;
         }
         else
         {
             m1 = -(x2 - x1) / (y2 - y1);
             m2 = -(x3 - x2) / (y3 - y2);
             mx1 = (x1 + x2) / 2.0;
             mx2 = (x2 + x3) / 2.0;
             my1 = (y1 + y2) / 2.0;
             my2 = (y2 + y3) / 2.0;
             xc = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2);
             yc = m1 * (xc - mx1) + my1;
         }

         dx = x2 - xc;
         dy = y2 - yc;
         rsqr = dx * dx + dy * dy;
         r = Math.sqrt(rsqr);

         dx = xp - xc;
         dy = yp - yc;
         drsqr = dx * dx + dy * dy;

         circle[0] = xc;
         circle[1] = yc;
         circle[2] = r;

         return (drsqr <= rsqr ? true : false);
     }

     /*
         Triangulation subroutine
         Takes as input NV vertices in array pxyz
         Returned is a list of ntri triangular faces in the array v
         These triangles are arranged in a consistent clockwise order.
         The triangle array 'v' should be malloced to 3 * nv
         The vertex array pxyz must be big enough to hold 3 more points
         The vertex array must be sorted in increasing x values say

         qsort(p,nv,sizeof(XYZ),XYZCompare);

         int XYZCompare(void *v1,void *v2)
         {
             XYZ *p1,*p2;
             p1 = v1;
             p2 = v2;
             if (p1->x < p2->x)
                 return(-1);
             else if (p1->x > p2->x)
                 return(1);
             else
                 return(0);
         }
     */

     int triangulate()
     {
         triangleIndexes = new int[nPoints * 3][3];
         
         boolean complete[];
         int nedge = 0;
         int trimax, emax = 200;
         int status = 0;

         boolean inside;
         //int 	i, j, k;
         double xp, yp, x1, y1, x2, y2, x3, y3, xc, yc, r;
         double xmin, xmax, ymin, ymax, xmid, ymid;
         double dx, dy, dmax;

         nTriangles = 0;

         /* Allocate memory for the completeness list, flag for each triangle */
         trimax = 4 * nPoints;
         complete = new boolean[trimax];
         // for (int ic = 0; ic < trimax; ic++) complete[ic] = false;

         /* Allocate memory for the edge list */
         edgeIndexes = new int[emax][2];

         /*
 		 Find the maximum and minimum vertex bounds.
         This is to allow calculation of the bounding triangle
         */
         xmin = points[0][0];
         ymin = points[0][1];
         xmax = xmin;
         ymax = ymin;
         for (int i = 1; i < nPoints; i++)
         {
             if (points[i][0] < xmin) xmin = points[i][0];
             if (points[i][0] > xmax) xmax = points[i][0];
             if (points[i][1] < ymin) ymin = points[i][1];
             if (points[i][1] > ymax) ymax = points[i][1];
         }
         dx = xmax - xmin;
         dy = ymax - ymin;
         dmax = (dx > dy) ? dx : dy;
         xmid = (xmax + xmin) / 2.0;
         ymid = (ymax + ymin) / 2.0;

         /*
		 Set up the supertriangle
		 This is a triangle which encompasses all the sample points.
		 The supertriangle coordinates are added to the end of the
		 vertex list. The supertriangle is the first triangle in
		 the triangle list.
         */
         points[nPoints + 0][0] = xmid - 2.0 * dmax;
         points[nPoints + 0][1] = ymid - dmax;
         points[nPoints + 0][2] = 0.0;
         points[nPoints + 1][0] = xmid;
         points[nPoints + 1][1] = ymid + 2.0 * dmax;
         points[nPoints + 1][2] = 0.0;
         points[nPoints + 2][0] = xmid + 2.0 * dmax;
         points[nPoints + 2][1] = ymid - dmax;
         points[nPoints + 2][2] = 0.0;
         triangleIndexes[0][0] = nPoints;
         triangleIndexes[0][1] = nPoints + 1;
         triangleIndexes[0][2] = nPoints + 2;
         complete[0] = false;
         nTriangles = 1;
         // Include each point one at a time into the existing mesh
         for (int i = 0; i < nPoints; i++)
         {

             xp = points[i][0];
             yp = points[i][1];
             nedge = 0;

             /*
                 Set up the edge buffer.
                 If the point (xp,yp) lies inside the circumcircle then the
                 three edges of that triangle are added to the edge buffer
                 and that triangle is removed.
             */
             double[] center = new double[3];
             for (int j = 0; j < nTriangles; j++)
             {
                 if (complete[j]) continue;
                 x1 = points[triangleIndexes[j][0]][0];
                 y1 = points[triangleIndexes[j][0]][1];
                 x2 = points[triangleIndexes[j][1]][0];
                 y2 = points[triangleIndexes[j][1]][1];
                 x3 = points[triangleIndexes[j][2]][0];
                 y3 = points[triangleIndexes[j][2]][1];
                 inside = circumCircle(xp, yp, x1, y1, x2, y2, x3, y3, center);
                 if (inside)
                 {
                     /* Check that we haven't exceeded the edge list size */
                     if (nedge + 3 >= emax)
                     {
                         emax += 100;
                         int[][] newEdgeIndexes = new int[emax][2];
                         System.arraycopy(edgeIndexes, 0, newEdgeIndexes, 0, edgeIndexes.length);
                         edgeIndexes = newEdgeIndexes;
                     }
                     edgeIndexes[nedge + 0][0] = triangleIndexes[j][0];
                     edgeIndexes[nedge + 0][1] = triangleIndexes[j][1];
                     edgeIndexes[nedge + 1][0] = triangleIndexes[j][1];
                     edgeIndexes[nedge + 1][1] = triangleIndexes[j][2];
                     edgeIndexes[nedge + 2][0] = triangleIndexes[j][2];
                     edgeIndexes[nedge + 2][1] = triangleIndexes[j][0];
                     nedge += 3;
                     triangleIndexes[j][0] = triangleIndexes[nTriangles - 1][0];
                     triangleIndexes[j][1] = triangleIndexes[nTriangles - 1][1];
                     triangleIndexes[j][2] = triangleIndexes[nTriangles - 1][2];
                     complete[j] = complete[nTriangles - 1];
                     nTriangles--;
                     j--;
                 }
             }

             /*
                 Tag multiple edges
                 Note: if all triangles are specified anticlockwise then all
                 interior edges are opposite pointing in direction.
             */
             for (int j = 0; j < nedge - 1; j++)
             {
                 //if ( !(edges[j][0] < 0 && edges[j][1] < 0) )
                 for (int k = j + 1; k < nedge; k++)
                 {
                     if ((edgeIndexes[j][0] == edgeIndexes[k][1]) && (edgeIndexes[j][1] == edgeIndexes[k][0]))
                     {
                         edgeIndexes[j][0] = -1;
                         edgeIndexes[j][1] = -1;
                         edgeIndexes[k][0] = -1;
                         edgeIndexes[k][1] = -1;
                     }
                     /* Shouldn't need the following, see note above */
                     if ((edgeIndexes[j][0] == edgeIndexes[k][0]) && (edgeIndexes[j][1] == edgeIndexes[k][1]))
                     {
                         edgeIndexes[j][0] = -1;
                         edgeIndexes[j][1] = -1;
                         edgeIndexes[k][0] = -1;
                         edgeIndexes[k][1] = -1;
                     }
                 }
             }

             /*
                 Form new triangles for the current point
                 Skipping over any tagged edges.
                 All edges are arranged in clockwise order.
             */
             for (int j = 0; j < nedge; j++)
             {
                 if (edgeIndexes[j][0] == -1 || edgeIndexes[j][1] == -1)
                     continue;
                 if (nTriangles >= trimax) return -1;
                 triangleIndexes[nTriangles][0] = edgeIndexes[j][0];
                 triangleIndexes[nTriangles][1] = edgeIndexes[j][1];
                 triangleIndexes[nTriangles][2] = i;
                 complete[nTriangles] = false;
                 nTriangles++;
             }
         }

         /*
 			Remove triangles with supertriangle vertices
            These are triangles which have a vertex number greater than nPoints
         */

		 if(debug) print();
		 nextHullPointIndexes = new int[nPoints];
		 hullPointIndexes = new int[nPoints];
         for (int i = 0; i < nTriangles; i++)
         {
             if (triangleIndexes[i][0] >= nPoints || triangleIndexes[i][1] >= nPoints || triangleIndexes[i][2] >= nPoints)
             {
				 addNextPointIndex(i);
                 triangleIndexes[i] = triangleIndexes[nTriangles - 1];
                 nTriangles--;
                 i--;
             }
         }
		 orderEdgePoints();
		 computeNormal();
         return nTriangles;
     }

	 private void addNextPointIndex(int nTriangle)
	 {
		int[] triangle = triangleIndexes[nTriangle];
		int i1 = triangle[2];
		for(int i = 0; i < 3; i++)
		{
			int i0 = i1;
			i1 = triangle[i];
			if(i0 < nPoints && i1 < nPoints)
			{
				nextHullPointIndexes[i0] = i1;
			}
		}
	 }

	 private void computeNormal()
	 {
		 int n00, n01, n10, n11;
		 n00 = 0;
		 n01 = nPoints/2;
		 n10 = n01/2;
		 n11 = (n01 + nPoints)/2;
		 double[] v0 = StsMath.subtract(points[n01],  points[n00]);
		 double[] v1 = StsMath.subtract(points[n11],  points[n10]);
		 normal = StsMath.leftCrossProductDouble(v0, v1);
		 StsMath.normalize(normal);
	 }

	 private void orderEdgePoints()
	 {
		 int hullPointIndex;
		 int nextHullPointIndex = -1;
		 for(hullPointIndex = 0; hullPointIndex < nPoints; hullPointIndex++)
		 {
			 nextHullPointIndex = nextHullPointIndexes[hullPointIndex];
			 if(nextHullPointIndex != -1) break;
		 }
		 if(nextHullPointIndex == -1) return;

		 hullPointIndexes[nHullPoints++] = hullPointIndex;
		 int startIndex = hullPointIndex;
		 while(nHullPoints < nPoints && nextHullPointIndex != startIndex && nextHullPointIndex != -1)
		 {
			 hullPointIndexes[nHullPoints++] = nextHullPointIndex;
			 nextHullPointIndex = nextHullPointIndexes[nextHullPointIndex];
		 }
		 hullPointIndexes = (int[])StsMath.trimArray(hullPointIndexes, nHullPoints);
	 }

/*
	 private void orderEdgePoints()
	 {
		hullPointIndexes = new int[nPoints];
		int firstPoint, nextPoint = 0;
		int nOrderedPoints = 0;
		for(firstPoint = 0; firstPoint < nPoints; firstPoint++)
		{
			int[] connectedEdgePoint = connectedEdgePoints[firstPoint];
			if(connectedEdgePoint != null)
			{
				nextPoint = connectedEdgePoint[1];
				hullPointIndexes[nOrderedPoints++] = firstPoint;
				break;
			}
		}
		int prevPoint = firstPoint;
		while(nextPoint != firstPoint)
		{
	    	int[] connectedEdgePoint = connectedEdgePoints[nextPoint];
			if(connectedEdgePoint[0] == prevPoint)
			{
				hullPointIndexes[nOrderedPoints++] = nextPoint;
				prevPoint = nextPoint;
				nextPoint = connectedEdgePoint[1];
			}
			else
			if(connectedEdgePoint[1] == prevPoint)
			{
				hullPointIndexes[nOrderedPoints++] = nextPoint;
				prevPoint = nextPoint;
				nextPoint = connectedEdgePoint[0];
			}
		}
		hullPointIndexes = (int[])StsMath.trimArray(hullPointIndexes, nOrderedPoints);
		hullPoints = new double[nOrderedPoints][];
		for(int n = 0; n < nOrderedPoints; n++)
			hullPoints[n] = points[hullPointIndexes[n]];
	 }
*/
     public void print()
     {
	  	System.out.println("Points");
         for (int tt = 0; tt < points.length; tt++)
             System.out.println(tt + " " + points[tt][0] + " " + points[tt][1] +  " " + points[tt][2]);
		 System.out.println("Triangle indexes");
         for (int tt = 0; tt < nTriangles; tt++)
             System.out.println(tt + " " + triangleIndexes[tt][0] + " " + triangleIndexes[tt][1] + " " + triangleIndexes[tt][2]);
		 for (int i = 0; i < nTriangles; i++)
         {
             if (triangleIndexes[i][0] >= nPoints || triangleIndexes[i][1] >= nPoints || triangleIndexes[i][2] >= nPoints)
				 System.out.println("Will remove triangle " + i);
     	}
	 }

	 public static void main(String[] args)
     {
		 StsAsciiFile file = new StsAsciiFile(StsFile.constructor(args[0]));
		 file.openRead();
		 String line = null;
		 double[][] points = null;
		 double[] point;
		 int nPoints = 0;
		 while(true)
		 {
			 try
			 {
				line = file.readLine();
			 }
			 catch(Exception e)
			 {
				 StsException.outputWarningException(StsDelaunayTriangulation.class, "main", "Failed to read line from file " + args[0], e);
			 }
			 if(line == null) break;
			 String[] tokens = StsStringUtils.getTokens(line);
			 int nTokens = tokens.length;
			 point = new double[nTokens];
			 for(int n = 0; n < nTokens; n++)
				point[n] = Double.parseDouble(tokens[n]);
			 points = (double[][])StsMath.arrayAddElement(points, point);
			 nPoints++;
		 }
		 StsBoundingBox box = new StsBoundingBox(points);
		 double[][] normalizedPoints = box.normalizePoints(points);
         StsDelaunayTriangulation triangulation = new StsDelaunayTriangulation(normalizedPoints);
         //triangulation.print();
     }

     public static void main0(String[] args)
     {
         int nv = 20;
         if (args.length > 0 && args[0] != null) nv = new Integer(args[0]).intValue();
         if (nv <= 0 || nv > 1000) nv = 20;

         //System.out.println("Creating " + nv + " random points.");

         double[][] points = new double[nv + 3][3];

         for (int i = 0; i < points.length; i++)
             points[i] = new double[]{i * 4.0, 400.0 * Math.random(), 0.0};

         StsDelaunayTriangulation triangulation = new StsDelaunayTriangulation(points);
         triangulation.print();
     }
 }
