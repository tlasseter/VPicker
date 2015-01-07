package com.Sts.Utilities;

// Some math routines borrowed from the VRML implementation.
// The original file was VbRotation.java

import com.Sts.Types.*;

public class StsQnion
 {
     double[] quat;

     public StsQnion()
     {
         quat = new double[4];
     }

     public StsQnion(float q0, float q1, float q2, float q3)
     {
         quat = new double[4];
         quat[0] = q0;
         quat[1] = q1;
         quat[2] = q2;
         quat[3] = q3;
         //       normalize();
     }

     public StsQnion(double q0, double q1, double q2, double q3)
     {
         quat = new double[4];
         quat[0] = q0;
         quat[1] = q1;
         quat[2] = q2;
         quat[3] = q3;
         //       normalize();
     }

     public StsQnion(StsPoint axis, float radians)
     {
         quat = new double[4];
         setValue(axis, radians);
     }

     public StsQnion(double[] axis, float radians)
     {
         quat = new double[4];
         setValue(axis, radians);
     }

     public StsQnion(double[] axis, double radians)
     {
         quat = new double[4];
         setValue(axis, radians);
     }

     public StsQnion(float f[])
     {
         quat = new double[4];
         for (int n = 0; n < 4; n++)
             quat[n] = f[n];
     }

     private StsQnion setValue(StsPoint axis, float radians)
     {
         StsPoint q = new StsPoint(3);
         q.copyFrom(axis);
         q.normalize();
         q.multiply((float) Math.sin(radians / 2.0f));

         quat[0] = q.v[0];
         quat[1] = q.v[1];
         quat[2] = q.v[2];

         quat[3] = (float) Math.cos(radians / 2.0f);

         return this;
     }

     private StsQnion setValue(double[] axis, double radians)
     {
         StsMath.normalize(axis);
         double sin = Math.sin(radians / 2.0);
         for (int n = 0; n < 3; n++)
             quat[n] = axis[n] * sin;
         quat[3] = (float) Math.cos(radians / 2.0f);
         return this;
     }

     public float getAxisValue(StsPoint axis)
     {
         float len;
         StsPoint q = new StsPoint(quat[0], quat[1], quat[2]);

         if ((len = q.length()) > 0.00001f)
         {
             q.multiply(1.0f / len);
             axis.copyFrom(q);
             return (float) (2.0 * Math.acos(quat[3]));
         }

         else
         {
             axis.setXYZ(0.0f, 0.0f, 1.0f);
             return 0.0f;
         }
     }

     /**
      * return the value of the quaternion the way VRML likes it. i.e. an array
      * of 4 floats representing the 3 vector values and the angle.
      */
     public float[] getOrientation()
     {
         float o[] = new float[4];
         StsPoint axis = new StsPoint(3);
         o[3] = getAxisValue(axis);
         o[0] = axis.v[0];
         o[1] = axis.v[1];
         o[2] = axis.v[2];
         return o;
     }

     // q2*this ??
     public StsQnion times(StsQnion q2)
     {
         StsQnion q = new StsQnion(q2.quat[3] * quat[0] + q2.quat[0] * quat[3] +
             q2.quat[1] * quat[2] - q2.quat[2] * quat[1],

             q2.quat[3] * quat[1] + q2.quat[1] * quat[3] +
                 q2.quat[2] * quat[0] - q2.quat[0] * quat[2],

             q2.quat[3] * quat[2] + q2.quat[2] * quat[3] +
                 q2.quat[0] * quat[1] - q2.quat[1] * quat[0],

             q2.quat[3] * quat[3] - q2.quat[0] * quat[0] -
                 q2.quat[1] * quat[1] - q2.quat[2] * quat[2]);
         //        q.normalize();
         return (q);
     }

     // Implementation details

     public void normalize()
     {
         float dist = (float) (1.0 / magnitude());

         quat[0] *= dist;
         quat[1] *= dist;
         quat[2] *= dist;
         quat[3] *= dist;
     }

     public float magnitude()
     {
         return (float) (Math.sqrt(quat[0] * quat[0] +
             quat[1] * quat[1] +
             quat[2] * quat[2] +
             quat[3] * quat[3]));
     }

     public StsQnion inverse()
     {
         //        normalize();
         return (conjugate());
     }

     public StsQnion conjugate()
     {
         return (new StsQnion(-quat[0], -quat[1], -quat[2], quat[3]));
     }

     public StsPoint rotateVec(StsPoint v)
     {
         StsQnion q = new StsQnion();

         q.quat[0] = v.v[0];
         q.quat[1] = v.v[1];
         q.quat[2] = v.v[2];
         q.quat[3] = 0F;
         //       q.normalize();
         q = inverse().times(q.times(this));
         StsPoint v2 = new StsPoint(q.quat[0], q.quat[1], q.quat[2]);
         return (v2);
     }

     public StsPoint leftRotateVec(StsPoint v)
     {
         StsQnion q = new StsQnion();

         q.quat[0] = v.v[0];
         q.quat[1] = v.v[1];
         q.quat[2] = v.v[2];
         q.quat[3] = 0F;
         //       q.normalize();
         q = times(q.times(inverse()));
         StsPoint v2 = new StsPoint(q.quat[0], q.quat[1], q.quat[2]);
         return (v2);
     }

     public double[] leftRotateVec(double[] v)
     {
         StsQnion q = new StsQnion();

         q.quat[0] = v[0];
         q.quat[1] = v[1];
         q.quat[2] = v[2];
         q.quat[3] = 0F;
         //       q.normalize();
         q = times(q.times(inverse()));
         return q.quat;
     }

     /**
      * Given a point, rotate it thru an angle about an axis thru the origin
      *
      * @return rotated point
      */
     static public StsPoint rotatePointAroundAxis(StsPoint point, StsPoint axis, float angleRad)
     {
         StsQnion q = new StsQnion(axis, angleRad);
         return q.leftRotateVec(point);
     }

     public String toString()
     {
         return ("[" + quat[3] + ",(" + quat[0] + "," + quat[1] + "," + quat[2]
             + ")]");
     }

     static public void main(String[] args)
     {
         float angle = (float) StsMath.RADperDEG * 45;
         StsPoint axis = new StsPoint(0.0f, 1.0f, 0.0f);
         StsPoint point = new StsPoint(10.0f, 0.0f, 0.0f);
         System.out.println(" axis: " + axis.toString() + " angle: " + angle + " result should be right: 7.07, 0, -7.07 or left: 7.07, 0, 7.07 ");
         debugTest(axis, angle, point);

         axis = new StsPoint(1.0f, 0.0f, 0.0f);
         point = new StsPoint(0.0f, 10.0f, 0.0f);
         System.out.println(" axis: " + axis.toString() + " angle: " + angle + " result should be right: 0, 7.07, 7.07 or left: 7.07, 0, -7.07 ");
         debugTest(axis, angle, point);

         axis = new StsPoint(0.0f, 0.0f, 1.0f);
         point = new StsPoint(10.0f, 0.0f, 0.0f);
         System.out.println(" axis: " + axis.toString() + " angle: " + angle + " result should be right: 7.07, 7.07, 0 or left:  7.07, -7.07, 0");
         debugTest(axis, angle, point);
     }

     static private void debugTest(StsPoint axis, float angle, StsPoint p)
     {
         StsQnion q = new StsQnion(axis, angle);
         StsPoint r = q.rotateVec(p);
         System.out.println("rite-hand rotation");
         r.print();
         r = q.leftRotateVec(p);
         System.out.println("left-hand rotation");
         r.print();
     }
 }
