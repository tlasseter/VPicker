package com.Sts.Utilities;

import Jama.*;

/**
 * Created by IntelliJ IDEA.
  * User: Tom Lasseter
  * Date: Jun 5, 2008
  * Time: 4:09:51 PM
  * To change this template use File | Settings | File Templates.
  */
 public class StsLstSqPointsPlane
 {
     public double[][] points;
     int nPoints;
     public double[] centroid;
     public double[] normal;

     public StsLstSqPointsPlane()
     {
     }

     public StsLstSqPointsPlane(double[][] points)
     {
         this.points = points;
         nPoints = points.length;
         Matrix a = new Matrix(points);
         //        System.out.println("Points");
         //        a.print(10, 2);
         centroid = subtractCentroid(points);
         //        print(centroid, "centroid: ");
         //        System.out.println("A");
         //        a.print(10, 2);
         SingularValueDecomposition svd = a.svd();
         Matrix v = svd.getV();
         //        System.out.println("v");
         //        v.print(10, 2);
         double[] singularValues = svd.getSingularValues();
         //        print(singularValues, "Singular Values: ");
         normal = v.getCol(2);
         //        print(normal, "Normal: ");
         projectPoints(points, normal, centroid);
         //        System.out.println("Projected points");
         //        a.print(10, 2);
     }

     void print(double[] values, String message)
     {
         System.out.print(message + ": ");
         for (int n = 0; n < values.length; n++)
             System.out.print(" " + values[n]);
         System.out.println();
     }

     public double[] subtractCentroid(double[][] points)
     {
         int nCoors = points[0].length;
         double[] centroid = new double[nCoors];
         for (int n = 0; n < nPoints; n++)
             for (int i = 0; i < nCoors; i++)
                 centroid[i] += points[n][i];

         for (int i = 0; i < nCoors; i++)
             centroid[i] /= nPoints;

         for (int n = 0; n < nPoints; n++)
         {
             points[n][0] -= centroid[0];
             points[n][1] -= centroid[1];
             points[n][2] -= centroid[2];
         }
         return centroid;
     }

	 /** project the points to the plane */
     public void projectPoints(double[][] points, double[] normal, double[] centroid)
     {
         int nCoors = points[0].length;
         for (int n = 0; n < nPoints; n++)
         {
             double dot = StsMath.dot(points[n], normal);
             for (int i = 0; i < nCoors; i++)
                 points[n][i] += (-dot * normal[i] + centroid[i]);
         }
     }

     public double[][] getPoints() { return points; }
     public double[] getNormal() { return normal; }

	 static public StsLstSqPointsPlane rotatePointsToHorizPlane(double[][] points)
	 {
		double[] vertical = new double[] { 0.0, 0.0, 1.0 };
		double[][] rotatedPoints = StsMath.copy(points);
		StsLstSqPointsPlane plane = new StsLstSqPointsPlane(rotatedPoints);
		StsMath.normalize(plane.normal);
		double[] rotationAxis = StsMath.cross(plane.normal, vertical);
		double angle = Math.acos(StsMath.dot(plane.normal, vertical));
		StsQnion qVertical = new StsQnion(rotationAxis, angle);
		int nPoints = points.length;
		for (int n = 0; n < nPoints; n++)
			rotatedPoints[n] = qVertical.leftRotateVec(rotatedPoints[n]);
		plane.centroid = qVertical.leftRotateVec(plane.centroid);
		return plane;
	 }

     static public void main(String[] args)
     {
		 double[][] points = new double[][] { { 1.0, 1.0, 1.0 }, { 1.0, 0.0, 0.0 }, { 1.0, -1.0, 1.0 } };
		 StsLstSqPointsPlane plane = rotatePointsToHorizPlane(points);
		 int nPoints = points.length;
		 double[][] rotatedPoints = plane.points;
		 double[] centroid = plane.centroid;
		 for(double[] rotatedPoint : rotatedPoints)
		 {
				rotatedPoint = StsMath.subtract(rotatedPoint, centroid);
			 System.out.println(rotatedPoint[0] + " " + rotatedPoint[1] + " " + rotatedPoint[2]);
		 }



	 }
 }
