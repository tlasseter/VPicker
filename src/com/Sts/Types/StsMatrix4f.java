

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.Utilities.*;

public class StsMatrix4f
{
    float m[] = new float[16];

    public StsMatrix4f()
    {
    }

    public void setIdentity()
    {

        m[ 0] = 1.0f; m[ 1] = 0.0f; m[ 2] = 0.0f; m[ 3] = 0.0f;
        m[ 4] = 0.0f; m[ 5] = 1.0f; m[ 6] = 0.0f; m[ 7] = 0.0f;
        m[ 8] = 0.0f; m[ 9] = 0.0f; m[10] = 1.0f; m[11] = 0.0f;
        m[12] = 0.0f; m[13] = 0.0f; m[14] = 0.0f; m[15] = 1.0f;
    }

    public void translate(float x, float y, float z)
    {
        int n;

        for(n = 0; n < 4; n++)
            m[n] = m[n] + x*m[n+12];
        for(n = 4; n < 8; n++)
           m[n] = m[n] + y*m[n+8];
        for(n = 8; n < 12; n++)
           m[n] = m[n] + z*m[n+4];
    }

    public void scale(float x, float y, float z)
    {
        int n;

        for(n = 0; n < 4; n++)
            m[n] = x*m[n];
        for(n = 4; n < 8; n++)
           m[n] = y*m[n];
        for(n = 8; n < 12; n++)
           m[n] = z*m[n];
    }

    public void rotX(float a)
    {
        int n;
        float temp[] = new float[12];

        float cos = (float)StsMath.cosd(a);
        float sin = (float)StsMath.sind(a);

        for(n = 4; n < 8; n++)
        {
            temp[n] = m[n];
            m[n] = cos*m[n] - sin*m[n+4];
        }
        for(n = 8; n < 12; n++)
            m[n] = sin*temp[n-4] + cos*m[n];
    }

    public void rotY(float a)
    {
        int n;
        float temp[] = new float[4];

        float cos = (float)StsMath.cosd(a);
        float sin = (float)StsMath.sind(a);

        for(n = 0; n < 4; n++)
        {
            temp[n] = m[n];
            m[n] = cos*m[n] + sin*m[n+8];
        }
        for(n = 8; n < 12; n++)
            m[n] = -sin*temp[n-8] + cos*m[n];
    }

   public void rotZ(float a)
    {
        int n;
        float temp[] = new float[4];

        float cos = (float)StsMath.cosd(a);
        float sin = (float)StsMath.sind(a);

       for(n = 0; n < 4; n++)
       {
            temp[n] = m[n];
            m[n] = cos*m[n] - sin*m[n+4];
       }
        for(n = 4; n < 8; n++)
            m[n] = sin*temp[n-4] + cos*m[n];
    }

    public void translateInverse(float x, float y, float z)
    {
        translate(-x, -y, -z);
    }

    public void scaleInverse(float x, float y, float z)
    {
        scale(1.0f/x, 1.0f/y, 1.0f/z);
    }

    public void rotXInverse(float a)
    {
        rotX(-a);
    }

    public void rotYInverse(float a)
    {
        rotY(-a);
    }
    public void rotZInverse(float a)
    {
        rotZ(-a);
    }

    public float[] vectorMult3d(float[] vIn)
    {
        int c, n;

        float[] vOut = new float[3];

        for(c = 0, n = 0; c < 3; c++)
        {
            vOut[c]  = m[n++]*vIn[0];
            vOut[c] += m[n++]*vIn[1];
            vOut[c] += m[n++]*vIn[2];
            vOut[c] += m[n++];
        }
        return vOut;
    }

    public void vectorMult3d(StsPoint vIn, StsPoint vOut)
    {
        int c, n;

        for(c = 0, n = 0; c < 3; c++)
        {
            vOut.v[c]  = m[n++]*vIn.v[0];
            vOut.v[c] += m[n++]*vIn.v[1];
            vOut.v[c] += m[n++]*vIn.v[2];
            vOut.v[c] += m[n++];
        }
    }

    public static void main(String[] args)
    {
        StsMatrix4f matrix = new StsMatrix4f();
        matrix.setIdentity();
        matrix.translate(0.0f, 0.0f, 0.0f);
        matrix.rotZ(-30.0f);

        StsMatrix4f inverseMatrix = new StsMatrix4f();
        inverseMatrix.rotZ(30.0f);
        inverseMatrix.translate(0.0f, 0.0f, 0.0f);

        float[] vIn = new float[] { 20.0f, 15.0f, 0.0f };
        float[] vOut = matrix.vectorMult3d(vIn);
        System.out.println("Translate: vector should be (25, 0, 0): " +
                           vOut[0] + " " + vOut[1] + " " + vOut[2]);
        float[] vInverse = inverseMatrix.vectorMult3d(vOut);
        System.out.println("Translate: inverse vector should be (20, 15, 0): " +
                           vInverse[0] + " " + vInverse[1] + " " + vInverse[2]);


/*
        matrix.translate(10.0f, 20.0f, 30.0f);
        matrix.rotX(45.0f);
        matrix.rotY(55.0f);
        matrix.rotZ(65.0f);
        matrix.scale(0.1f, 0.1f, 0.1f);
        matrix.scaleInverse(0.1f, 0.1f, 0.1f);
        matrix.rotZInverse(65.0f);
        matrix.rotYInverse(55.0f);
        matrix.rotXInverse(45.0f);
        matrix.translateInverse(10.0f, 20.0f, 30.0f);

        System.out.println("Should be identity matrix: ");
        System.out.println("     " + matrix.m[ 0] + " " + matrix.m[ 1] + " " + matrix.m[ 2] + " " + matrix.m[ 3]);
        System.out.println("     " + matrix.m[ 4] + " " + matrix.m[ 5] + " " + matrix.m[ 6] + " " + matrix.m[ 7]);
        System.out.println("     " + matrix.m[ 8] + " " + matrix.m[ 9] + " " + matrix.m[10] + " " + matrix.m[11]);
        System.out.println("     " + matrix.m[12] + " " + matrix.m[13] + " " + matrix.m[14] + " " + matrix.m[15] );

        StsPoint vIn = new StsPoint(100.0f, 200.0f, 300.0f);
        StsPoint vOut = new StsPoint();

        matrix.setIdentity();
        matrix.translate(10.0f, 20.0f, 30.0f);
        matrix.vectorMult3d(vIn, vOut);
        System.out.println("Translate: vector should be (110, 220, 330): " +
        vOut.v[0] + " " + vOut.v[1] + " " + vOut.v[2]);

        matrix.setIdentity();
        matrix.rotZ(45.0f);

        vIn.v[0] = 100.f; vIn.v[1] = 0.f; vIn.v[2] = 0.f;
        matrix.vectorMult3d(vIn, vOut);
        System.out.println("RotZ: vector should be (70.7, 70.7, 0.0): " +
        vOut.v[0] + " " + vOut.v[1] + " " + vOut.v[2]);
*/
      }
}
