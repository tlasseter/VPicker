
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Utilities;

public class StsRotationMatrix
{
    double[][] matrix = { {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1} };
    double[] originPoint;
    double angleRad;

    static final byte X = 1;
    static final byte Y = 2;
    static final byte Z = 3;

    public StsRotationMatrix()
    {
    }

    static public StsRotationMatrix constructRotationMatrix(double[] originPoint, double[] oldPoint, double[] newPoint)
    {

        StsRotationMatrix rotMatrix = new StsRotationMatrix();
        if(!rotMatrix.computeRotationMatrix(originPoint, oldPoint, newPoint)) return null;
        return rotMatrix;
    }


    static public StsRotationMatrix constructRotationMatrix(float[] cornerFloats, float[] oldFloats, float[] newFloats)
    {
        StsRotationMatrix rotMatrix = new StsRotationMatrix();
        double[] cornerDoubles = new double[3];
        double[] oldDoubles = new double[3];
        double[] newDoubles = new double[3];
        for(int n = 0; n < 3; n++) cornerDoubles[n] = cornerFloats[n];
        for(int n = 0; n < 3; n++) oldDoubles[n] = oldFloats[n];
        for(int n = 0; n < 3; n++) newDoubles[n] = newFloats[n];

        if(!rotMatrix.computeRotationMatrix(cornerDoubles, oldDoubles, newDoubles)) return null;
        return rotMatrix;
    }

    public boolean computeRotationMatrix(double[] originPoint, double[] oldPoint, double[] newPoint)
    {
        double[] oldVector, newVector, axis;
        double axisLen, zAxisLen, rAxisLen;
        double sinx, cosx, siny, cosy, sinz, cosz;

        this.originPoint = StsMath.copy(originPoint);

    /*	compute rotation vectors */

        oldVector = StsMath.subtract(oldPoint, this.originPoint);
        newVector = StsMath.subtract(newPoint, this.originPoint);

    /* reverse z since its not right-handed */
    /* and scale it						    */

        this.originPoint[2] *= -1;
        oldVector[2]  *= -1;
        newVector[2]  *= -1;

    /* normalize vectors */

        StsMath.normalize(oldVector);
        StsMath.normalize(newVector);

    /*	compute axis of rotation */

        axis = StsMath.cross(oldVector, newVector);

    /*	compute angle of rotation around axis */

        axisLen = StsMath.length(axis);

        sinx = axisLen;
        cosx = StsMath.dot(oldVector, newVector);
        angleRad = Math.asin(axisLen);

    /* 	if rotation angle small, return FALSE */

        if(sinx < 0.0001) return(true); // rotation will use identity matrix

    /* 	compute rotation around z and y axes */

        zAxisLen = Math.sqrt(axis[0]*axis[0] + axis[1]*axis[1]);

    /* if zAxisLen == 0.0, rotation axis is z-axis	*/

        if(zAxisLen != 0.0)
        {
            rAxisLen = 1.0/zAxisLen;
            cosz = axis[0]*rAxisLen;
            sinz = -axis[1]*rAxisLen;
        }
        else
        {
            cosz = cosx;
            sinz = sinx;
        }

        rAxisLen = 1.0/axisLen;
        cosy = zAxisLen*rAxisLen;
        siny = axis[2]*rAxisLen;

    /* 	build composite rotation matrix */
    /* matrix is initialized to the identity matrix to start */

    /*	unscale for reverse z */

        matrix[2][2] = -1;

    /*	untranslate axis to origin */

        matrixTranslate(1.0, this.originPoint);

        matrixRotate(Z, -sinz, cosz);
        matrixRotate(Y, -siny, cosy);
        matrixRotate(X, sinx, cosx);
        matrixRotate(Y, siny, cosy);
        matrixRotate(Z, sinz, cosz);

    /*	untranslate axis to origin */

        matrixTranslate(-1.0, this.originPoint);

    /*	scale for reverse z */

        matrix[2][0] *= -1;
        matrix[2][1] *= -1;
        matrix[2][2] *= -1;

    /*  Debug for trashed matrix					*/
    /*	assert( fabs(rotMatrix[0][0]) <= 1.e10 );	*/

        return true;
    }

    static public StsRotationMatrix constructRotationMatrix(double[] axisPoint, double[] axis, double angle)
    {

        StsRotationMatrix rotMatrix = new StsRotationMatrix();
        if(!rotMatrix.computeRotationMatrix(axisPoint, axis, angle)) return null;
        return rotMatrix;
    }

    static public StsRotationMatrix constructRotationMatrix(float[] originPointF, float[] axisF, float angleF)
    {
        double[] originPoint = StsMath.copyDouble(originPointF);
        double[] axis = StsMath.copyDouble(axisF);
        return constructRotationMatrix(originPoint, axis, (double)angleF);
    }

    public boolean computeRotationMatrix(double[] originPoint, double[] axis, double angle)
    {
        double zAxisLen, rAxisLen;
        double sinx, cosx, siny, cosy, sinz, cosz;

        this.originPoint = StsMath.copy(originPoint);

    /* reverse z since its not right-handed */
    /* and scale it						    */

        this.originPoint[2] *= -1;
        axis[2] *= -1;

    /* normalize vectors */

        StsMath.normalize(axis);

    /*	compute angle of rotation around axis */

        sinx = StsMath.sind(angle);
        cosx = StsMath.cosd(angle);
        angleRad = angle/StsMath.DEGperRAD;

    /* 	if rotation angle small, return FALSE */

        if(Math.abs(sinx) < 0.0001) return(true); // rotation will use identity matrix

    /* 	compute rotation around z and y axes */

        zAxisLen = Math.sqrt(axis[0]*axis[0] + axis[1]*axis[1]);

    /* if zAxisLen == 0.0, rotation axis is z-axis	*/

        if(zAxisLen != 0.0)
        {
            rAxisLen = 1.0/zAxisLen;
            cosz = axis[0]*rAxisLen;
            sinz = -axis[1]*rAxisLen;
        }
        else
        {
            cosz = cosx;
            sinz = sinx;
        }

        cosy = zAxisLen;
        siny = axis[2];

    /* 	build composite rotation matrix */
    /* matrix is initialized to the identity matrix to start */

    /*	unscale for reverse z */

        matrix[2][2] = -1;

    /*	untranslate axis to origin */

        matrixTranslate(1.0, this.originPoint);

        matrixRotate(Z, -sinz, cosz);
        matrixRotate(Y, -siny, cosy);
        matrixRotate(X, sinx, cosx);
        matrixRotate(Y, siny, cosy);
        matrixRotate(Z, sinz, cosz);

    /*	untranslate axis to origin */

        matrixTranslate(-1.0, this.originPoint);

    /*	scale for reverse z */

        matrix[2][0] *= -1;
        matrix[2][1] *= -1;
        matrix[2][2] *= -1;

    /*  Debug for trashed matrix					*/
    /*	assert( fabs(rotMatrix[0][0]) <= 1.e10 );	*/

        return true;
    }

    static private double[][] getIdentityMatrix()
    {
       double[][] iMatrix = new double[4][4];
       iMatrix[0][0] = 1; iMatrix[1][1] = 1; iMatrix[2][2] = 1; iMatrix[3][3] = 1;
       return iMatrix;
    }

    private void matrixTranslate(double sign, double[] vect)
    {
        for(int j = 0; j < 3; j++)
            for(int i = 0; i < 3; i++)
                matrix[3][j] += sign*vect[i]*matrix[i][j];
    }

    private void matrixRotate(int axis, double sin, double cos)
    {
        double[][] rotMatrix = getIdentityMatrix();

        switch(axis)
        {
            case X:

            rotMatrix[1][1] = cos;
            rotMatrix[1][2] = sin;
            rotMatrix[2][1] = -sin;
            rotMatrix[2][2] = cos;
            rowMultMatrix(1, 2, rotMatrix, matrix);
            break;

            case Y:

            rotMatrix[0][0] = cos;
            rotMatrix[0][2] = -sin;
            rotMatrix[2][0] = sin;
            rotMatrix[2][2] = cos;
            rowMultMatrix(0, 2, rotMatrix, matrix);
            break;

            case Z:

            rotMatrix[0][0] = cos;
            rotMatrix[0][1] = sin;
            rotMatrix[1][0] = -sin;
            rotMatrix[1][1] = cos;
            rowMultMatrix(0, 1, rotMatrix, matrix);
        }
    }

    static void rowMultMatrix(int imin, int imax, double[][] a, double[][] b)
    {
        int i, j, k;
        double[][] c = getIdentityMatrix();

        for(i = imin; i <= imax; i++)
            for(j = 0; j < 4; j++)
            {
                c[i][j] = 0.0;
                for(k = 0; k < 4; k++)
                    c[i][j] += a[i][k]*b[k][j];
            }

        for(i = imin; i <= imax; i++)
            for(j = 0; j < 4; j++)
                b[i][j] = c[i][j];
    }

    public double[] pointRotate(double[] vIn)
    {
        int i, j;

        double[] vOut = new double[3];

        for(j = 0; j < 3; j++)
        {
            vOut[j] = matrix[3][j];
            for(i = 0; i < 3; i++)
                vOut[j] += vIn[i]*matrix[i][j];
        }

        return vOut;
    }

    public void pointRotate(float[] vIn)
    {
        int i, j;

        double[] vOut = new double[3];

        for(j = 0; j < 3; j++)
        {
            vOut[j] = matrix[3][j];
            for(i = 0; i < 3; i++)
                vOut[j] += vIn[i]*matrix[i][j];
        }

        for(i = 0; i < 3; i++)
            vIn[i] = (float)vOut[i];
    }

    public double getAngleRad() { return angleRad; }

    public static void main(String[] args)
    {
        double[] originPoint, oldPoint, newPoint, testPoint, rotatedPoint, axis;
        double angle;
        StsRotationMatrix rotMatrix;

        originPoint = new double[] { 1.0, 0.0, 0.0 };
        oldPoint = new double[] { 0.0, 0.0, 1.0 };
        newPoint = new double[] { 2.0, 0.0, 1.0 };

//        originPoint = new double[] { 1.0, 1.0, 1.0 };
//        oldPoint = new double[] { 2.0, 1.0, 2.0 };
//        newPoint = new double[] { 0.0, 1.0, 2.0 };

        rotMatrix = constructRotationMatrix(originPoint, oldPoint, newPoint);

        testPoint = new double[] { 1.0, 0.0, 1.0 };
//        testPoint = new double[] { 1.5, 1.0, 1.5 };
        System.out.println("point: " + StsMath.toString(testPoint));
        rotatedPoint = rotMatrix.pointRotate(testPoint);
        System.out.println("rotatedPoint: " + StsMath.toString(rotatedPoint));

        originPoint = new double[] { 10.0, 100.0, 0.0 };
        axis = new double[] { 0.0, -1.0, 0.0 };
        angle = 45;
        testPoint = new double[] { 20.0, 100.0,  };
        rotMatrix = constructRotationMatrix(originPoint, axis, angle);
        rotatedPoint = rotMatrix.pointRotate(testPoint);
        System.out.println("rotatedPoint: " + StsMath.toString(rotatedPoint));

        axis = new double[] { 0, 0, -1.0 };
        originPoint = new double[] { 10.0, 10.0, 0.0 };
        double sqrt2 = Math.sqrt(2);
        double length = 10*sqrt2;
        angle = 90;

        double angDeg = angle*StsMath.DEGperRAD;
        rotMatrix = constructRotationMatrix(originPoint, axis, angle);
        testPoint = new double[] {9, 9, 0 };
        rotatedPoint = rotMatrix.pointRotate(testPoint);
        System.out.println("rotatedPoint: " + StsMath.toString(rotatedPoint));
    }
}
