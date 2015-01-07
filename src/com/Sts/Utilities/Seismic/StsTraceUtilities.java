package com.Sts.Utilities.Seismic;

import com.Sts.Utilities.*;
 import com.Sts.DBTypes.*;
 import com.Sts.UI.*;
 import com.Sts.Types.StsRotatedGridBoundingBox;

 import javax.media.opengl.*;
 import java.nio.*;

/**
 * Created by IntelliJ IDEA.
  * User: Tom Lasseter
  * Date: May 15, 2009
  * Time: 6:42:54 AM
  * To change this template use File | Settings | File Templates.
  */
 public class StsTraceUtilities
 {
     /** fractional max error allowed in curve fitting wavelets */
     static double maxError = 0.01;

     public final static byte POINT_ORIGINAL = 0;

     public final static byte POINT_PLUS_ZERO_CROSSING = 1;
     public final static byte POINT_MAXIMUM = 2;
     public final static byte POINT_MINUS_ZERO_CROSSING = 3;
     public final static byte POINT_MINIMUM = 4;

     public final static byte POINT_FALSE_MAXIMUM = 5;
     public final static byte POINT_FALSE_MINIMUM = 6;

     public final static byte POINT_INTERPOLATED = 7;
     public final static byte POINT_FLAT_ZERO = 8;

     public final static byte POINT_ANY = 9;

     public final static byte[] pointTypesBefore = { POINT_ANY, POINT_MINIMUM, POINT_PLUS_ZERO_CROSSING, POINT_MAXIMUM, POINT_MINUS_ZERO_CROSSING, POINT_ANY, POINT_ANY, POINT_ANY, POINT_ANY, POINT_ANY };
     public final static byte[] pointTypesAfter = { POINT_ANY, POINT_MAXIMUM, POINT_MINUS_ZERO_CROSSING, POINT_MINIMUM, POINT_PLUS_ZERO_CROSSING, POINT_ANY, POINT_ANY, POINT_ANY, POINT_ANY, POINT_ANY };

     public final static String[] typeStrings = { "Original", "Zero Plus", "Max", "Zero Minus", "Min", "False Max", "False Min", "Interpolated", "Flat Zero", "Any" };
     public final static StsColor[] typeColors = { StsColor.GRAY, StsColor.MAGENTA, StsColor.RED, StsColor.CYAN, StsColor.BLUE, StsColor.GREEN, StsColor.GREEN, StsColor.PINK, StsColor.AQUAMARINE, StsColor.YELLOW };
     public static double[][] tracePoints;
     public static byte[] tracePointTypes;
     public static final StsColor outsideMuteColor = StsColor.GRAY;

     public static final boolean debug = false;

     /**
      * create a vector points[nPoints][3]: time, amplitude, and slope for each point.
      * Insert the first point; if it is zero, insert another zero point just before the first non-zero point
      * but only if it is greater than the time value of the third point.  So if the first point is zero, and the second point
      * is not, we insert all points.  If the first n points are zero, we insert a zero at index 0 and at index n-1, and then
      * all points after that.
      *
      * @param floatBuffer contains data for this trace; must be positioned to start of trace if offset is not zero
      * @param zMin        time or depth of first sample
      * @param zInc        time or depth increment between samples
      * @param nSamples    number of samples in trace
      * @return
      */
     static public double[][] getDataPoints(FloatBuffer floatBuffer, float zMin, float zInc, int nSamples)
     {
         float[] floats = new float[nSamples];
         floatBuffer.get(floats);
         return getDataPoints(floats, 0, zMin, zInc, nSamples);
     }

     static public double[][] getDataPoints(float[] floats, float zMin, float zInc, int nPoints)
     {
         return getDataPoints(floats, 0, zMin, zInc, nPoints);
     }

     static public double[][] getDataPoints(float[] floats, int sliceMin, float sliceZMin, float zInc, int nSlices)
     {
         int n = 0, s = 0;
         try
         {
             double[][] points = new double[nSlices][3];
             float t = sliceZMin;
             for (n = 0, s = sliceMin; n < nSlices; n++, s++, t += zInc)
             {
                 points[n][0] = t;
                 double v = floats[s];
                 if (Math.abs(v) < StsParameters.roundOff) v = 0;
                 points[n][1] = v;
             }
             // compute slopes (0 curvature as end conditions)
             for (n = 1; n < nSlices - 1; n++)
                 points[n][2] = (points[n + 1][1] - points[n - 1][1]) / 2;

             points[0][2] = 1.5 * (points[1][1] - points[0][1]) - points[1][2] / 2;
             points[nSlices - 1][2] = 1.5 * (points[nSlices - 1][1] - points[nSlices - 2][1]) - points[nSlices - 2][2] / 2;
             return points;
         }
         catch (Exception e)
         {
             StsException.outputWarningException(StsTraceUtilities.class, "getWigglePoints", " sliceMin " + sliceMin + " slice " + s + " index " + n, e);
             return null;
         }
     }
     /*
         static public double[][] computeInterpolatedPoints(FloatBuffer floatBuffer, float zMin, float zInc, int nSamples)
         {
             return computeInterpolatedPoints(floatBuffer, 0, zMin, zInc, nSamples, nSamples, 0.0f, 0.0f);
         }

         static public double[][] computeInterpolatedPoints(FloatBuffer floatBuffer, int offset, float zMin, float zInc, int nSamples)
         {
             return computeInterpolatedPoints(floatBuffer, offset, zMin, zInc, nSamples, nSamples, 0.0f, 0.0f);
         }

         static public double[][] computeInterpolatedPoints(FloatBuffer floatBuffer, int inputSliceMin, float zMin, float zInc, int nInputSlices, int nOutputSlices, float valueScale, float valueOffset)
         {
             float[] floats = new float[nInputSlices];
             floatBuffer.get(floats);
             return computeInterpolatedPoints(floats, inputSliceMin, zMin, zInc, nOutputSlices, valueScale, valueOffset);
         }

         static public void displayInterpolatedPoints(GL gl, FloatBuffer floatBuffer, float y, int inputSliceMin, float zMin, float zInc, int nInputSlices, int nOutputSlices, float valueScale, float valueOffset, StsWiggleDisplayProperties wiggleProperties, double[] muteRange)
         {
             float[] floats = new float[nInputSlices];
             floatBuffer.get(floats);
             displayInterpolatedPoints(gl, floats, y, inputSliceMin, zMin, zInc, nOutputSlices, valueScale, valueOffset, wiggleProperties, muteRange);
         }
     */

     static public void displayInterpolatedPoints(GL gl, float[] traceValues, float y, int displayMin, float displayZMin, float zInc, int nDisplaySlices, float valueScale, float valueOffset, StsWiggleDisplayProperties wiggleProperties, double[] muteRange, int displayInc)
     {
         float[] scaledValues = new float[nDisplaySlices];
         if (valueOffset != 0.0f)
         {
             for (int n = 0, i = displayMin; n < nDisplaySlices; n++, i++)
                 scaledValues[n] = valueOffset + valueScale * traceValues[i];
         }
         else
             for (int n = 0, i = displayMin; n < nDisplaySlices; n++, i++)
                 scaledValues[n] = valueScale * traceValues[i];
         if (wiggleProperties.hasFill())
         {
             drawFilledWiggleTraces(gl, scaledValues, 0, nDisplaySlices, y, displayZMin, zInc, wiggleProperties, muteRange, displayInc);
             drawFilledWiggleTracesLine(gl, scaledValues, 0, nDisplaySlices, y, displayZMin, zInc, wiggleProperties, muteRange, displayInc);
         }
         if (wiggleProperties.getWiggleDrawLine())
             drawWiggleTraces(gl, scaledValues, 0, nDisplaySlices, y, displayZMin, zInc, wiggleProperties, displayInc);
         if (wiggleProperties.getWiggleDrawPoints())
             drawWigglePoints(gl, scaledValues, 0, nDisplaySlices, y, displayZMin, zInc, displayInc); // uses tracePoints and tracePointTypes in StsTraceUtilities
     }

     static public void displayInterpolatedPoints(GL gl, float[] traceValues, float y, int displayMin, float displayZMin, float zInc, int nDisplaySlices, StsWiggleDisplayProperties wiggleProperties, double[] muteRange, int displayInc)
     {
         if (wiggleProperties.hasFill())
         {
             drawFilledWiggleTraces(gl, traceValues, displayMin, nDisplaySlices, y, displayZMin, zInc, wiggleProperties, muteRange, displayInc);
             drawFilledWiggleTracesLine(gl, traceValues, displayMin, nDisplaySlices, y, displayZMin, zInc, wiggleProperties, muteRange, displayInc);
         }
         if (wiggleProperties.getWiggleDrawLine())
             drawWiggleTraces(gl, traceValues, displayMin, nDisplaySlices, y, displayZMin, zInc, wiggleProperties, displayInc);
         if (wiggleProperties.getWiggleDrawPoints())
             drawWigglePoints(gl, traceValues, displayMin, nDisplaySlices, y, displayZMin, zInc, displayInc); // uses tracePoints and tracePointTypes in StsTraceUtilities
     }

     /*
         static public double[][] computeInterpolatedPoints(float[] floats, int inputSliceMin, float zMin, float zInc, int nOutputSlices, float valueScale, float valueOffset)
         {
             int nIntervals = StsTraceUtilities.computeInterpolationInterval(zInc, 5);
             double[][] dataPoints = StsTraceUtilities.getDataPoints(floats, inputSliceMin, zMin, zInc, nOutputSlices);
             tracePoints = computeInterpolatedDataPoints(dataPoints, nIntervals, zInc, valueScale, valueOffset);
             if(tracePoints == null) return null;
             defineTracePointTypes(tracePoints);
             return tracePoints;
         }
     */
     static public float[] computeCubicInterpolatedPoints(float[] inputPoints, int nIntervals)
     {
         int nInputPoints = inputPoints.length;
         int nPoints = (nInputPoints - 1) * nIntervals + 1;
         try
         {
             StsCubicForwardDifference cubicFd = new StsCubicForwardDifference(nIntervals);
             float[] points = new float[nPoints];
             double z0;
             double z1 = inputPoints[0];
             double z2 = inputPoints[1];
             double z3 = inputPoints[2];
             int out = 0;
             double slope1 = z2 - z1;
             double slope2 = 0.5 * (z3 - z1);
             cubicFd.hermiteInitialize(z1, z2, slope1, slope2);
             double[] interpolatedPoints = cubicFd.evaluate();
             for (int j = 0; j < nIntervals; j++)
                 points[out++] = (float) interpolatedPoints[j];
             for (int i = 3; i < nInputPoints; i++)
             {
                 z0 = z1;
                 z1 = z2;
                 z2 = z3;
                 z3 = inputPoints[i];
                 cubicFd.splineInitialize(z0, z1, z2, z3);
                 interpolatedPoints = cubicFd.evaluate();
                 for (int j = 0; j < nIntervals; j++)
                     points[out++] = (float) interpolatedPoints[j];
             }
             z0 = z1;
             z1 = z2;
             z2 = z3;
             slope1 = slope2;
             slope2 = z3 - z1;
             cubicFd.hermiteInitialize(z1, z2, slope1, slope2);
             interpolatedPoints = cubicFd.evaluate();
             for (int j = 0; j < nIntervals; j++)
                 points[out++] = (float) interpolatedPoints[j];
             points[out++] = (float) z2;
             return points;
         }
         catch (Exception e)
         {
             StsException.outputWarningException(StsTraceUtilities.class, "computeInterpolatedPoints", e);
             return null;
         }
     }

     static public float[] computeCubicInterpolatedPoints(double[] inputPoints, int nIntervals)
     {
         int nInputPoints = inputPoints.length;
         int nPoints = (nInputPoints - 1) * nIntervals + 1;
         try
         {
             StsCubicForwardDifference cubicFd = new StsCubicForwardDifference(nIntervals);
             float[] points = new float[nPoints];
             double z0;
             double z1 = inputPoints[0];
             double z2 = inputPoints[1];
             double z3 = inputPoints[2];
             int out = 0;
             double slope1 = z2 - z1;
             double slope2 = 0.5 * (z3 - z1);
             cubicFd.hermiteInitialize(z1, z2, slope1, slope2);
             double[] interpolatedPoints = cubicFd.evaluate();
             for (int j = 0; j < nIntervals; j++)
                 points[out++] = (float) interpolatedPoints[j];
             for (int i = 3; i < nInputPoints; i++)
             {
                 z0 = z1;
                 z1 = z2;
                 z2 = z3;
                 z3 = inputPoints[i];
                 cubicFd.splineInitialize(z0, z1, z2, z3);
                 interpolatedPoints = cubicFd.evaluate();
                 for (int j = 0; j < nIntervals; j++)
                     points[out++] = (float) interpolatedPoints[j];
             }
             z0 = z1;
             z1 = z2;
             z2 = z3;
             slope1 = slope2;
             slope2 = z3 - z1;
             cubicFd.hermiteInitialize(z1, z2, slope1, slope2);
             interpolatedPoints = cubicFd.evaluate();
             for (int j = 0; j < nIntervals; j++)
                 points[out++] = (float) interpolatedPoints[j];
             points[out++] = (float) z2;
             return points;
         }
         catch (Exception e)
         {
             StsException.outputWarningException(StsTraceUtilities.class, "computeInterpolatedPoints", e);
             return null;
         }
     }

     static public float[] computeLinearInterpolatedPoints(double[] inputPoints, int nIntervals)
     {
         int nInputPoints = inputPoints.length;
         int nPoints = (nInputPoints - 1) * nIntervals + 1;
         try
         {
             StsCubicForwardDifference cubicFd = new StsCubicForwardDifference(nIntervals);
             float[] points = new float[nPoints];
             double z0;
             int out = 0;
             double z1 = inputPoints[0];
             double z = 0.0;
             for (int i = 1; i < nInputPoints; i++)
             {
                 z0 = z1;
                 z1 = inputPoints[i];
                 z = z0;
                 double dz = (z1 - z0)/nIntervals;
                 for (int j = 0; j < nIntervals; j++)
                 {
                     points[out++] = (float)z;
                     z += dz;
                 }
             }
             points[out++] = (float) z;
             return points;
         }
         catch (Exception e)
         {
             StsException.outputWarningException(StsTraceUtilities.class, "computeInterpolatedPoints", e);
             return null;
         }
     }

     static public byte[] defineTracePointTypes(float[] pointValues)
     {
         int nPoints = pointValues.length;
         byte[] tracePointTypes = new byte[nPoints];
         // for first point, assign as Max or min, or flat-zero
         float v = pointValues[0];
         if (v > 0)
             tracePointTypes[0] = POINT_FALSE_MAXIMUM;
         else if (v < 0)
             tracePointTypes[0] = POINT_FALSE_MINIMUM;
         else // v == 0
             tracePointTypes[0] = POINT_FLAT_ZERO;
         // assign maxima and minima
         float vp = pointValues[1];
         for (int n = 1; n < nPoints - 1; n++)
         {
             float vm = v;
             v = vp;
             vp = pointValues[n + 1];

             if (v < vm && v <= vp && v < 0)
                 tracePointTypes[n] = POINT_MINIMUM;
             else if (v > vm && v >= vp && v > 0)
                 tracePointTypes[n] = POINT_MAXIMUM;
             else
             {
                 if (vm <= 0.0 && v > 0.0)
                 {
                     if (-vm < v)
                         tracePointTypes[n - 1] = POINT_PLUS_ZERO_CROSSING;
                     else
                         tracePointTypes[n] = POINT_PLUS_ZERO_CROSSING;
                 }
                 else if (vm >= 0.0 && v < 0.0)
                 {
                     if (vm < -v)
                         tracePointTypes[n - 1] = POINT_MINUS_ZERO_CROSSING;
                     else
                         tracePointTypes[n] = POINT_MINUS_ZERO_CROSSING;
                 }
             }
         }
         v = vp;
         if (v > 0)
             tracePointTypes[nPoints - 1] = POINT_FALSE_MAXIMUM;
         else if (v < 0)
             tracePointTypes[nPoints - 1] = POINT_FALSE_MINIMUM;
         else
             tracePointTypes[nPoints - 1] = POINT_FLAT_ZERO;

         return tracePointTypes;
     }

     static public boolean computeInstantAmpAndPhase(float[] gatherData, byte[] gatherPointTypes, float[] amplitudes, float[] phases)
     {
         int nPoints = gatherData.length;
         int[] eventIndexes = new int[nPoints];
         int nEvents = 0;
         for (int n = 0; n < nPoints; n++)
         {
             byte newPointType = gatherPointTypes[n];
             if (isPointTypeEvent(newPointType))
                 eventIndexes[nEvents++] = n;
         }
         for (int n = 0; n < nEvents; n++)
         {
             int i = eventIndexes[n];
             byte pointType = gatherPointTypes[i];
             switch (pointType)
             {
                 case POINT_MAXIMUM:
                     phases[i] = 0.0f;
                     break;
                 case POINT_MINIMUM:
                     phases[i] = 180;
                     break;
                 case POINT_PLUS_ZERO_CROSSING:
                     phases[i] = -90;
                     break;
                 case POINT_MINUS_ZERO_CROSSING:
                     phases[i] = 90;
                     break;
             }
         }
         int prevIndex = eventIndexes[0];
         for (int n = 1; n < nEvents; n++)
         {
             int i = eventIndexes[n];
             byte pointType = gatherPointTypes[i];
             switch (pointType)
             {
                 case POINT_MAXIMUM:
                     amplitudes[i] = gatherData[i];
                     interpolateAmplitudes(prevIndex, i, amplitudes);
                     prevIndex = i;
                     break;
                 case POINT_MINIMUM:
                     amplitudes[i] = -gatherData[i];
                     interpolateAmplitudes(prevIndex, i, amplitudes);
                     prevIndex = i;
                     break;
             }
         }
         int nextIndex = eventIndexes[0];
         for (int n = 1; n < nEvents; n++)
         {
             int index = nextIndex;
             nextIndex = eventIndexes[n];
             interpolatePhases(index, nextIndex, phases);
         }
         return true;
     }

     static public void interpolateAmplitudes(int index, int nextIndex, float[] amplitudes)
     {
         double df = 1.0 / (nextIndex - index);
         double dAmp = df * (amplitudes[nextIndex] - amplitudes[index]);
         for (int i = index + 1; i < nextIndex; i++)
             amplitudes[i] = (float) (amplitudes[i - 1] + dAmp);
     }

     static public void interpolatePhases(int index, int nextIndex, float[] phases)
     {
         double df = 1.0 / (nextIndex - index);
         double dPhase = phases[nextIndex] - phases[index];
         if (dPhase < 0.0)
             dPhase += 360;
         dPhase *= df;
         for (int i = index + 1; i < nextIndex; i++)
         {
             phases[i] = (float) (phases[i - 1] + dPhase);
             if (phases[i] > 180)
                 phases[i] -= 360;
         }
     }

     static public final boolean isPointTypeEvent(byte type)
     {
         return type >= POINT_MINIMUM && type <= POINT_MINUS_ZERO_CROSSING;
     }

     /**
      * Determining min, max, plus_zero_crossing, and minus_zero_crossing
      * Given a point value v, and the previous point value vm, and the next point value vp:
      * if v < vm && v <= vp : v is a minimum if v < 0
      * if v > vm && v >= vp : v is a maximum if v > 0
      * Having assigned minima and maxima, we can go back thru and assign zero crossings
      * Given a point v and a previous point vm:
      * If vm <= 0 and v > 0: if -vm < v && vm isn't a min or a max: vm is a plus-zero-crossing; otherwise v is a plus zero-crossing
      * If vm >= 0 and v < 0: if vm < -v && vm isn't a min or a max: vm is a minus-zero-crossing; otherwise v is a minus zero-crossing
      * /*
      *
      * @return
      */
     static private void assignMaxAndMin(int n, float[] pointValues, byte[] tracePointTypes)
     {
         double vm = pointValues[n - 1];
         double v = pointValues[n];
         double vp = pointValues[n + 1];

         if (v < vm && v <= vp && v < 0)
             tracePointTypes[n] = POINT_MINIMUM;
         else if (v > vm && v >= vp && v > 0)
             tracePointTypes[n] = POINT_MAXIMUM;
         else
             tracePointTypes[n] = POINT_INTERPOLATED;
     }

     /**
      * Two successive points are either side of zero axis; assign one as plus or minus zero-crossing
      * depending on which is closer and whether increasing or decreasing.
      * Don't assign if already assigned as a max or min.
      *
      * @param n
      * @param interpolatedPoints
      * @param tracePointTypes
      */
     static private void assignZeroCrossings(int n, double[][] interpolatedPoints, byte[] tracePointTypes)
     {
         double vm = interpolatedPoints[n - 1][1];
         double v = interpolatedPoints[n][1];
         if (vm <= 0.0 && v > 0.0)
         {
             if (-vm < v && tracePointTypes[n - 1] == POINT_INTERPOLATED)
                 tracePointTypes[n - 1] = POINT_PLUS_ZERO_CROSSING;
             else if (tracePointTypes[n] == POINT_INTERPOLATED)
                 tracePointTypes[n] = POINT_PLUS_ZERO_CROSSING;
         }
         else if (vm >= 0.0 && v < 0.0)
         {
             if (vm < -v && tracePointTypes[n - 1] == POINT_INTERPOLATED)
                 tracePointTypes[n - 1] = POINT_MINUS_ZERO_CROSSING;
             else if (tracePointTypes[n] == POINT_INTERPOLATED)
                 tracePointTypes[n] = POINT_MINUS_ZERO_CROSSING;
         }
     }

     static public boolean isMaxMinOrZero(byte pointType)
     {
         return pointType >= POINT_PLUS_ZERO_CROSSING && pointType <= POINT_MINIMUM;
     }

     static private void adjustPointTypes()
     {
         byte pointType = POINT_INTERPOLATED;
         byte nextPointType = POINT_INTERPOLATED;
         int nPoint, nNextPoint;
         int nPoints = tracePointTypes.length;
         int n = 0;
         for (n = 0; n < nPoints; n++)
         {
             pointType = tracePointTypes[n];
             if (pointType == POINT_INTERPOLATED || pointType == POINT_ORIGINAL) continue;
             nPoint = n;
             break;
         }
         nPoint = n;
         n++;
         for (; n < nPoints; n++)
         {
             nextPointType = tracePointTypes[n];
             if (nextPointType == POINT_INTERPOLATED || nextPointType == POINT_ORIGINAL)
                 nNextPoint = n;
             break;
         }
         nNextPoint = n;
         n++;
         for (; n < nPoints; n++)
         {
             byte type = tracePointTypes[n];
             if (type != POINT_INTERPOLATED && type != POINT_ORIGINAL)
             {
                 byte prevPointType = pointType;
                 int nPrevPoint = nPoint;
                 pointType = nextPointType;
                 nPoint = nNextPoint;
                 nextPointType = type;
                 nNextPoint = n;
                 if (pointType == POINT_MINIMUM)
                 {
                     if (prevPointType != POINT_MINUS_ZERO_CROSSING && nextPointType != POINT_PLUS_ZERO_CROSSING)
                     {
                         pointType = POINT_INTERPOLATED;
                         tracePointTypes[nPoint] = POINT_INTERPOLATED;
                     }
                 }
                 else if (pointType == POINT_MAXIMUM)
                 {
                     if (prevPointType != POINT_PLUS_ZERO_CROSSING && nextPointType != POINT_MINUS_ZERO_CROSSING)
                     {
                         pointType = POINT_INTERPOLATED;
                         tracePointTypes[nPoint] = POINT_INTERPOLATED;
                     }
                 }
                 else if (pointType == POINT_PLUS_ZERO_CROSSING)
                 {
                     if (prevPointType != POINT_MINIMUM && nextPointType != POINT_MAXIMUM)
                     {
                         pointType = POINT_INTERPOLATED;
                         tracePointTypes[nPoint] = POINT_INTERPOLATED;
                     }
                     else
                         tracePoints[nPoint][1] = 0.0;
                 }
                 else if (pointType == POINT_MINUS_ZERO_CROSSING)
                 {
                     if (prevPointType != POINT_MAXIMUM && nextPointType != POINT_MINIMUM)
                     {
                         pointType = POINT_INTERPOLATED;
                         tracePointTypes[nPoint] = POINT_INTERPOLATED;
                     }
                     else
                         tracePoints[nPoint][1] = 0.0;
                 }
             }
         }
     }

     static public double interpolateCubic(double[] point0, double[] point1, double z, double zInc)
     {
         double z0 = point0[0];
         double v0 = point0[1];
         double s0 = point0[2];
         double z1 = point1[0];
         double v1 = point1[1];
         double s1 = point1[2];
         double f = (z - z0) / (z1 - z0);
         double slopeMult = (z1 - z0) / zInc;
         return StsMath.hermiteCubic(v0, v1, slopeMult * s0, slopeMult * s1, f);
     }

     static public double[] interpolateCubicValueAndSlope(double[] point0, double[] point1, double z, double zInc)
     {
         double z0 = point0[0];
         double v0 = point0[1];
         double s0 = point0[2];
         double z1 = point1[0];
         double v1 = point1[1];
         double s1 = point1[2];
         double f = (z - z0) / (z1 - z0);
         double slopeMult = (z1 - z0) / zInc;
         return StsMath.hermiteValueAndSlope(z, v0, v1, slopeMult * s0, slopeMult * s1, f);
     }

     static public void drawWiggleTraces(GL gl, float[] values, int nValueMin, int nValues, float x0, float displayZMin, float zInc, StsWiggleDisplayProperties wiggleProperties, int displayInc)
     {
         try
         {
             if (values == null || values.length < 3) return;
             StsColor lineColor = wiggleProperties.getLineColor();
             gl.glDisable(GL.GL_LIGHTING);
             gl.glLineWidth(1.0f);
             lineColor.setGLColor(gl);
             gl.glBegin(GL.GL_LINE_STRIP);
             float z = displayZMin;
             for (int i = 0, n = nValueMin; i < nValues; i++, n += displayInc, z += zInc)
                 gl.glVertex2f(values[n] + x0, z);
         }
         catch (Exception e)
         {
             StsException.outputWarningException(StsTraceUtilities.class, "drawWiggleTraces", e);
         }
         finally
         {
             gl.glEnd();
             gl.glEnable(GL.GL_LIGHTING);
         }
     }

     static public void drawWiggleTraces(GL gl, float[] values, float[] zValues, float x0, StsWiggleDisplayProperties wiggleProperties)
     {
         try
         {
             if (values == null || values.length < 3) return;
             int nValues = values.length;
             StsColor lineColor = wiggleProperties.getLineColor();
             gl.glDisable(GL.GL_LIGHTING);
             gl.glLineWidth(1.0f);
             lineColor.setGLColor(gl);
             gl.glBegin(GL.GL_LINE_STRIP);
             for (int i = 0; i < nValues; i++)
                 gl.glVertex2f(values[i] + x0, zValues[i]);
         }
         catch (Exception e)
         {
             StsException.outputWarningException(StsTraceUtilities.class, "drawWiggleTraces", e);
         }
         finally
         {
             gl.glEnd();
             gl.glEnable(GL.GL_LIGHTING);
         }
     }
     static public void drawWigglePoints(GL gl, float[] values, int nValueMin, int nValues, float x0, float displayZMin, float zInc, int displayInc)
     {
         if (nValues < 3) return;
         gl.glDisable(GL.GL_LIGHTING);
         gl.glPointSize(4);
         gl.glBegin(GL.GL_POINTS);
         float z = displayZMin;
         for (int i = 0, n = nValueMin; i < nValues; i++, n += displayInc, z += zInc)
         {
             byte pointType = getPointType(values, n);
             if (pointType <= POINT_INTERPOLATED) continue;
             getPointTypeColor(pointType).setGLColor(gl);
             gl.glVertex2f(values[n] + x0, z);
         }
         gl.glEnd();
         gl.glEnable(GL.GL_LIGHTING);
     }

     static public byte getPointType(float[] values, int n)
     {
         if (n == 0)
             return getEndPointType(values[0]);
         int nValues = values.length;
         if (n == nValues - 1)
             return getEndPointType(values[nValues - 1]);
         else
             return getPointType(values[n - 1], values[n], values[n + 1]);
     }

     static private byte getEndPointType(float value)
     {
         if (value > 0)
             return POINT_FALSE_MAXIMUM;
         else if (value < 0)
             return POINT_FALSE_MINIMUM;
         else // value == 0
             return POINT_FLAT_ZERO;
     }

     static private byte getPointType(float vm, float v, float vp)
     {
         if (v < vm && v <= vp && v < 0)
             return POINT_MINIMUM;
         else if (v > vm && v >= vp && v > 0)
             return POINT_MAXIMUM;
         else
         {
             if (vm <= 0.0 && v > 0.0)
             {
                 if (-vm < v)
                     return POINT_PLUS_ZERO_CROSSING;
                 else
                     return POINT_PLUS_ZERO_CROSSING;
             }
             else if (vm >= 0.0 && v < 0.0)
             {
                 if (vm < -v)
                     return POINT_MINUS_ZERO_CROSSING;
                 else
                     return POINT_MINUS_ZERO_CROSSING;
             }
             else
                 return POINT_INTERPOLATED;
         }
     }

     static public void drawFilledWiggleTraces(GL gl, float[] values, int nValueMin, int nValues, float x0, float displayZMin, float zInc, StsWiggleDisplayProperties wiggleProperties, int displayInc)
     {
         if (values == null || values.length < 3) return;
         try
         {
             gl.glDisable(GL.GL_LIGHTING);
             StsColor plusColor = wiggleProperties.getWigglePlusColor();
             StsColor minusColor = wiggleProperties.getWiggleMinusColor();
             float amp1 = values[nValueMin];
             float t1 = displayZMin;
             boolean plus1 = amp1 >= 0;
             if (plus1)
                 plusColor.setGLColor(gl);
             else
                 minusColor.setGLColor(gl);

             gl.glBegin(GL.GL_QUAD_STRIP);
             gl.glVertex2d(x0, t1);
             gl.glVertex2d(x0 + amp1, t1);
             for (int i = 1, n = nValueMin + displayInc; i < nValues; i++, n += displayInc)
             {
                 float t0 = t1;
                 t1 += zInc;
                 float amp0 = amp1;
                 amp1 = values[n];
                 boolean plus0 = plus1;
                 plus1 = amp1 >= 0;

                 if (plus0 && plus1 || !plus0 && !plus1) // line does not cross zero, so continue drawing fill color
                 {
                     gl.glVertex2f(x0, t1);
                     gl.glVertex2f(x0 + amp1, t1);
                 }
                 else // line crosses at tz; draw old color to tz and then new color to t1
                 {
                     float tz = t0 + amp0 * (t1 - t0) / (amp0 - amp1);
                     gl.glVertex2f(x0, tz);
                     gl.glVertex2f(x0, tz);
                     if (plus1)
                         plusColor.setGLColor(gl);
                     else
                         minusColor.setGLColor(gl);
                     gl.glVertex2f(x0, tz);
                     gl.glVertex2f(x0, tz);
                     gl.glVertex2f(x0, t1);
                     gl.glVertex2f(x0 + amp1, t1);
                 }
             }
             gl.glEnd();
             gl.glEnable(GL.GL_LIGHTING);
         }
         catch (Exception e)
         {
             StsException.outputException("StsPreStackLine.drawFilledWiggleTraces() failed.", e, StsException.WARNING);
         }
         finally
         {
             gl.glEnd();
             gl.glEnable(GL.GL_LIGHTING);
         }
     }

     static public void drawFilledWiggleTraces(GL gl, float[] values, int nValueMin, int nValues, float x0, float displayZMin, float zInc, StsWiggleDisplayProperties wiggleProperties, double[] muteRange, int displayInc)
     {
         if (values == null || values.length < 3) return;
         try
         {
             gl.glDisable(GL.GL_LIGHTING);
             StsColor plusColor = wiggleProperties.getWigglePlusColor();
             StsColor minusColor = wiggleProperties.getWiggleMinusColor();
             float amp1 = values[nValueMin];
             float t1 = displayZMin;
             boolean plus1 = amp1 >= 0;
             float topMute = (float) muteRange[0];
             float botMute = (float) muteRange[1];
             boolean muted = (t1 < topMute || t1 > botMute);
             if (muted) // start color is mute
                 outsideMuteColor.setGLColor(gl);
             else // start color is plus or minus color
             {
                 if (plus1)
                     plusColor.setGLColor(gl);
                 else
                     minusColor.setGLColor(gl);
             }
             gl.glBegin(GL.GL_QUAD_STRIP);
             gl.glVertex2d(x0, t1);
             gl.glVertex2d(x0 + amp1, t1);
             for (int i = 1, n = nValueMin + displayInc; i < nValues; i++, n += displayInc)
             {
                 float t0 = t1;
                 t1 += zInc;
                 float amp0 = amp1;
                 amp1 = values[n];
                 boolean plus0 = plus1;
                 plus1 = amp1 >= 0;

                 if (muted) // previous point is muted
                 {
                     muted = (t1 < topMute || t1 > botMute);
                     if (muted) // still muted
                     {
                         if (plus0 != plus1) // line does not cross zero, so continue drawing muted color to t1
                         {
                             float tz = t0 + amp0 * (t1 - t0) / (amp0 - amp1);
                             if (plus1)
                                 outsideMuteColor.setGLColor(gl);
                             else
                                 minusColor.setGLColor(gl); //don't fill negative wiggles with mute color!! SWC 10/20/09
                             gl.glVertex2f(x0, tz);
                             gl.glVertex2f(x0, tz);
                         }
                         gl.glVertex2d(x0, t1);
                         gl.glVertex2d(x0 + amp1, t1);
                     }
                     else // have crossed into unmuted interval; only occurs if previous point is muted and current point is not
                     {
                         float ampMute = amp0 + (amp1 - amp0) * (topMute - t0) / (t1 - t0);
                         if (plus0 && plus1 || !plus0 && !plus1) // line does not cross zero, so just draw unmuted color to t1
                         {
                             // complete muted color for first part of interval
                             gl.glVertex2d(x0, topMute);
                             gl.glVertex2d(x0 + ampMute, topMute);

                             // complete unmuted plus or minus color for rest of interval
                             if (plus1)
                                 plusColor.setGLColor(gl);
                             else
                                 minusColor.setGLColor(gl);

                             gl.glVertex2f(x0, topMute);
                             gl.glVertex2f(x0 + ampMute, topMute);
                             gl.glVertex2f(x0, t1);
                             gl.glVertex2f(x0 + amp1, t1);
                         }
                         else // line crosses at tz; so complete mute color to tz if tz > topMute
                         {
                             float tz = t0 + amp0 * (t1 - t0) / (amp0 - amp1);

                             if (topMute < tz)
                             {
                                 // complete muted color for first part of interval
                                 gl.glVertex2f(x0, topMute);
                                 gl.glVertex2f(x0 + ampMute, topMute);
                                 // switch to old color and complete to tz
                                 if (plus0)
                                     plusColor.setGLColor(gl);
                                 else
                                     minusColor.setGLColor(gl);
                                 gl.glVertex2f(x0, topMute);
                                 gl.glVertex2f(x0 + ampMute, topMute);
                                 gl.glVertex2f(x0, tz);
                                 gl.glVertex2f(x0, tz);
                                 // switch to new color and complete from tz to t1
                                 if (plus1)
                                     plusColor.setGLColor(gl);
                                 else
                                     minusColor.setGLColor(gl);
                                 gl.glVertex2f(x0, tz);
                                 gl.glVertex2f(x0, tz);
                                 gl.glVertex2f(x0, t1);
                                 gl.glVertex2f(x0 + amp1, t1);
                             }
                             else // topMute is > tz; so complete mute color to tz and then to topMute and new color beyond
                             {
                                 gl.glVertex2f(x0, tz);
                                 gl.glVertex2f(x0, tz);
                                 gl.glVertex2f(x0, topMute);
                                 gl.glVertex2f(x0 + ampMute, topMute);

                                 if (plus1)
                                     plusColor.setGLColor(gl);
                                 else
                                     minusColor.setGLColor(gl);

                                 gl.glVertex2f(x0, topMute);
                                 gl.glVertex2f(x0 + ampMute, topMute);
                                 gl.glVertex2f(x0, t1);
                                 gl.glVertex2f(x0 + amp1, t1);
                             }
                         }
                     }
                 }
                 else // previous point not muted
                 {
                     muted = (t1 > botMute);
                     if (!muted) // still not muted
                     {
                         if (plus0 && plus1 || !plus0 && !plus1) // line does not cross zero, so continue drawing fill color
                         {
                             gl.glVertex2f(x0, t1);
                             gl.glVertex2f(x0 + amp1, t1);
                         }
                         else // line crosses at tz; draw old color to tz and then new color to t1
                         {
                             float tz = t0 + amp0 * (t1 - t0) / (amp0 - amp1);
                             gl.glVertex2f(x0, tz);
                             gl.glVertex2f(x0, tz);
                             if (plus1)
                                 plusColor.setGLColor(gl);
                             else
                                 minusColor.setGLColor(gl);
                             gl.glVertex2f(x0, tz);
                             gl.glVertex2f(x0, tz);
                             gl.glVertex2f(x0, t1);
                             gl.glVertex2f(x0 + amp1, t1);
                         }
                     }
                     else // have crossed into muted interval; occurs only as we cross the botMute time
                     {
                         float ampMute = amp0 + (amp1 - amp0) * (botMute - t0) / (t1 - t0);
                         if (plus0 && plus1 || !plus0 && !plus1) // line does not cross zero, so just draw unmuted color to botMute and
                         {
                             // complete unmuted color for first part of interval
                             gl.glVertex2f(x0, botMute);
                             gl.glVertex2f(x0 + ampMute, botMute);

                             // switch to mute color for rest of interval
                             outsideMuteColor.setGLColor(gl);

                             gl.glVertex2f(x0, botMute);
                             gl.glVertex2f(x0 + ampMute, botMute);
                             gl.glVertex2f(x0, t1);
                             gl.glVertex2f(x0 + amp1, t1);
                         }
                         else // line crosses at tz; so complete color to tz if tz < botMute
                         {
                             float tz = t0 + amp0 * (t1 - t0) / (amp0 - amp1);

                             if (tz < botMute)
                             {
                                 // complete color for first part of interval
                                 gl.glVertex2f(x0, tz);
                                 gl.glVertex2f(x0, tz);
                                 gl.glVertex2f(x0, botMute);
                                 gl.glVertex2f(x0 + ampMute, botMute);

                                 // switch to mute color for rest of interval
                                 outsideMuteColor.setGLColor(gl);

                                 gl.glVertex2f(x0, botMute);
                                 gl.glVertex2f(x0 + ampMute, botMute);
                                 gl.glVertex2f(x0, t1);
                                 gl.glVertex2f(x0 + amp1, t1);
                             }
                             else  // tz > botMute, so draw old color to botMute, mute color to tz and then to t1
                             {
                                 gl.glVertex2f(x0, botMute);
                                 gl.glVertex2f(x0 + ampMute, botMute);

                                 // switch to mute color for rest of interval
                                 outsideMuteColor.setGLColor(gl);

                                 gl.glVertex2f(x0, botMute);
                                 gl.glVertex2f(x0 + ampMute, botMute);
                                 gl.glVertex2f(x0, t1);
                                 gl.glVertex2f(x0 + amp1, t1);
                             }
                         }
                     }
                 }
             }
             gl.glEnd();
             gl.glEnable(GL.GL_LIGHTING);
         }
         catch (Exception e)
         {
             StsException.outputException("StsPreStackLine.drawFilledWiggleTraces() failed.", e, StsException.WARNING);
         }
         finally
         {
             gl.glEnd();
             gl.glEnable(GL.GL_LIGHTING);
         }
     }

     static public void drawFilledWiggleTraces(GL gl, float[] values, float[] zValues, float x0, StsWiggleDisplayProperties wiggleProperties)
     {
         if (values == null || values.length < 3) return;
         try
         {
             int nValues = values.length;
             gl.glDisable(GL.GL_LIGHTING);
             StsColor plusColor = wiggleProperties.getWigglePlusColor();
             StsColor minusColor = wiggleProperties.getWiggleMinusColor();
             float amp1 = values[0];
             float t1 = zValues[0];
             boolean plus1 = amp1 >= 0;
             if (plus1)
                 plusColor.setGLColor(gl);
             else
                 minusColor.setGLColor(gl);

             gl.glBegin(GL.GL_QUAD_STRIP);
             gl.glVertex2d(x0, t1);
             gl.glVertex2d(x0 + amp1, t1);
             for (int i = 1; i < nValues; i++)
             {
                 float t0 = t1;
                 t1 = zValues[i];
                 float amp0 = amp1;
                 amp1 = values[i];
                 boolean plus0 = plus1;
                 plus1 = amp1 >= 0;

                 if (plus0 && plus1 || !plus0 && !plus1) // line does not cross zero, so continue drawing fill color
                 {
                     gl.glVertex2f(x0, t1);
                     gl.glVertex2f(x0 + amp1, t1);
                 }
                 else // line crosses at tz; draw old color to tz and then new color to t1
                 {
                     float tz = t0 + amp0 * (t1 - t0) / (amp0 - amp1);
                     gl.glVertex2f(x0, tz);
                     gl.glVertex2f(x0, tz);
                     if (plus1)
                         plusColor.setGLColor(gl);
                     else
                         minusColor.setGLColor(gl);
                     gl.glVertex2f(x0, tz);
                     gl.glVertex2f(x0, tz);
                     gl.glVertex2f(x0, t1);
                     gl.glVertex2f(x0 + amp1, t1);
                 }
             }
             gl.glEnd();
             gl.glEnable(GL.GL_LIGHTING);
         }
         catch (Exception e)
         {
             StsException.outputException("StsPreStackLine.drawFilledWiggleTraces() failed.", e, StsException.WARNING);
         }
         finally
         {
             gl.glEnd();
             gl.glEnable(GL.GL_LIGHTING);
         }
     }

     static public void drawFilledWiggleTracesLine(GL gl, float[] values, int nValueMin, int nValues, float x0, float displayZMin, float zInc, StsWiggleDisplayProperties wiggleProperties, int displayInc)
     {
         if (values == null || values.length < 3) return;

         try
         {
             if (!wiggleProperties.hasFill()) return;
             gl.glDisable(GL.GL_LIGHTING);
             StsColor plusColor = wiggleProperties.getWigglePlusColor();
             StsColor minusColor = wiggleProperties.getWiggleMinusColor();
             gl.glLineWidth(1.0f);
             float amp1 = values[nValueMin];
             float t1 = displayZMin;
             boolean plus1 = amp1 >= 0;
             if (plus1)
                 plusColor.setGLColor(gl);
             else
                 minusColor.setGLColor(gl);

             gl.glBegin(GL.GL_LINE_STRIP);
             gl.glVertex2d(x0, t1);
             gl.glVertex2d(x0 + amp1, t1);
             for (int i = 1, n = nValueMin + displayInc; i < nValues; i++, n += displayInc)
             {
                 float t0 = t1;
                 t1 += zInc;
                 float amp0 = amp1;
                 amp1 = values[n];
                 boolean plus0 = plus1;
                 plus1 = amp1 >= 0;

                 if (plus0 && plus1 || !plus0 && !plus1) // line does not cross zero, so continue drawing fill color
                 {
                     gl.glVertex2d(x0, t1);
                     //                            gl.glVertex2d(x0 + amp1, t1);
                 }
                 else // line crosses at tz; draw old color to tz and then new color to t1
                 {
                     double tz = t0 + amp0 * (t1 - t0) / (amp0 - amp1);
                     gl.glVertex2d(x0, tz);
                     //                            gl.glVertex2d(x0, tz);
                     if (plus1)
                         plusColor.setGLColor(gl);
                     else
                         minusColor.setGLColor(gl);
                     gl.glVertex2d(x0, tz);
                     //                            gl.glVertex2d(x0, tz);
                     gl.glVertex2d(x0, t1);
                     //                            gl.glVertex2d(x0 + amp1, t1);
                 }
             }
             gl.glEnd();
             gl.glEnable(GL.GL_LIGHTING);
         }
         catch (Exception e)
         {
             StsException.outputException("StsPreStackLine.drawFilledWiggleTraces() failed.", e, StsException.WARNING);
         }
         finally
         {
             gl.glEnd();
             gl.glEnable(GL.GL_LIGHTING);
         }
     }

     static public void drawFilledWiggleTracesLine(GL gl, float[] values, int nValueMin, int nValues, float x0, float displayZMin, float zInc, StsWiggleDisplayProperties wiggleProperties, double[] muteRange, int displayInc)
     {
         if (values == null || values.length < 3) return;

         try
         {
             if (!wiggleProperties.hasFill()) return;
             gl.glDisable(GL.GL_LIGHTING);
             StsColor plusColor = wiggleProperties.getWigglePlusColor();
             StsColor minusColor = wiggleProperties.getWiggleMinusColor();
             gl.glLineWidth(1.0f);
             float amp1 = values[nValueMin];
             float t1 = displayZMin;
             boolean plus1 = amp1 >= 0;
             double topMute = muteRange[0];
             double botMute = muteRange[1];
             boolean muted = (t1 < topMute || t1 > botMute);
             if (muted) // start color is mute
                 outsideMuteColor.setGLColor(gl);
             else // start color is plus or minus color
             {
                 if (plus1)
                     plusColor.setGLColor(gl);
                 else
                     minusColor.setGLColor(gl);
             }
             gl.glBegin(GL.GL_LINE_STRIP);
             gl.glVertex2d(x0, t1);
             gl.glVertex2d(x0 + amp1, t1);
             for (int i = 1, n = nValueMin + displayInc; i < nValues; i++, n += displayInc)
             {
                 float t0 = t1;
                 t1 += zInc;
                 float amp0 = amp1;
                 amp1 = values[n];
                 boolean plus0 = plus1;
                 plus1 = amp1 >= 0;

                 if (muted) // previous point is muted
                 {
                     muted = (t1 < topMute || t1 > botMute);
                     if (muted) // still muted
                     {
                         if (plus0 != plus1) // line does not cross zero, so continue drawing mmuted color to t1
                         {
                             float tz = t0 + amp0 * (t1 - t0) / (amp0 - amp1);
                             //                            gl.glVertex2d(x0, tz);
                             gl.glVertex2f(x0, tz);
                         }
                         gl.glVertex2d(x0, t1);
                         //                        gl.glVertex2d(x0 + amp1, t1);
                     }
                     else // have crossed into unmuted interval; only occurs if previous point is muted and current point is not
                     {
                         double ampMute = amp0 + (amp1 - amp0) * (topMute - t0) / (t1 - t0);
                         if (plus0 && plus1 || !plus0 && !plus1) // line does not cross zero, so just draw unmmuted color to t1
                         {
                             // complete muted color for first part of interval
                             gl.glVertex2d(x0, topMute);
                             //                            gl.glVertex2d(x0 + ampMute, topMute);

                             // complete unmuted plus or minus color for rest of interval
                             if (plus1)
                                 plusColor.setGLColor(gl);
                             else
                                 minusColor.setGLColor(gl);

                             gl.glVertex2d(x0, topMute);
                             //                            gl.glVertex2d(x0 + ampMute, topMute);
                             gl.glVertex2d(x0, t1);
                             //                            gl.glVertex2d(x0 + amp1, t1);
                         }
                         else // line crosses at tz; so complete mute color to tz if tz > topMute
                         {
                             double tz = t0 + amp0 * (t1 - t0) / (amp0 - amp1);

                             if (topMute < tz)
                             {
                                 // complete muted color for first part of interval
                                 gl.glVertex2d(x0, topMute);
                                 //                                gl.glVertex2d(x0 + ampMute, topMute);
                                 // switch to old color and complete to tz
                                 if (plus0)
                                     plusColor.setGLColor(gl);
                                 else
                                     minusColor.setGLColor(gl);
                                 gl.glVertex2d(x0, topMute);
                                 //                                gl.glVertex2d(x0 + ampMute, topMute);
                                 gl.glVertex2d(x0, tz);
                                 //                                gl.glVertex2d(x0, tz);
                                 // switch to new color and complete from tz to t1
                                 if (plus1)
                                     plusColor.setGLColor(gl);
                                 else
                                     minusColor.setGLColor(gl);
                                 gl.glVertex2d(x0, tz);
                                 //                                gl.glVertex2d(x0, tz);
                                 gl.glVertex2d(x0, t1);
                                 //                                gl.glVertex2d(x0 + amp1, t1);
                             }
                             else // topMute is > tz; so complete mute color to tz and then to topMute and new color beyond
                             {
                                 gl.glVertex2d(x0, tz);
                                 //                                gl.glVertex2d(x0, tz);
                                 gl.glVertex2d(x0, topMute);
                                 //                                gl.glVertex2d(x0 + ampMute, topMute);

                                 if (plus1)
                                     plusColor.setGLColor(gl);
                                 else
                                     minusColor.setGLColor(gl);

                                 gl.glVertex2d(x0, topMute);
                                 //                                gl.glVertex2d(x0 + ampMute, topMute);
                                 gl.glVertex2d(x0, t1);
                                 //                                gl.glVertex2d(x0 + amp1, t1);
                             }
                         }
                     }
                 }
                 else // previous point not muted
                 {
                     muted = (t1 > botMute);
                     if (!muted) // still not muted
                     {
                         if (plus0 && plus1 || !plus0 && !plus1) // line does not cross zero, so continue drawing fill color
                         {
                             gl.glVertex2d(x0, t1);
                             //                            gl.glVertex2d(x0 + amp1, t1);
                         }
                         else // line crosses at tz; draw old color to tz and then new color to t1
                         {
                             double tz = t0 + amp0 * (t1 - t0) / (amp0 - amp1);
                             gl.glVertex2d(x0, tz);
                             //                            gl.glVertex2d(x0, tz);
                             if (plus1)
                                 plusColor.setGLColor(gl);
                             else
                                 minusColor.setGLColor(gl);
                             gl.glVertex2d(x0, tz);
                             //                            gl.glVertex2d(x0, tz);
                             gl.glVertex2d(x0, t1);
                             //                            gl.glVertex2d(x0 + amp1, t1);
                         }
                     }
                     else // have crossed into muted interval; occurs only as we cross the botMute time
                     {
                         double ampMute = amp0 + (amp1 - amp0) * (botMute - t0) / (t1 - t0);
                         if (plus0 && plus1 || !plus0 && !plus1) // line does not cross zero, so just draw unmuted color to botMute and
                         {
                             // complete unmuted color for first part of interval
                             gl.glVertex2d(x0, botMute);
                             //                            gl.glVertex2d(x0 + ampMute, botMute);

                             // switch to mute color for rest of interval
                             outsideMuteColor.setGLColor(gl);

                             gl.glVertex2d(x0, botMute);
                             //                           gl.glVertex2d(x0 + ampMute, botMute);
                             gl.glVertex2d(x0, t1);
                             //                           gl.glVertex2d(x0 + amp1, t1);
                         }
                         else // line crosses at tz; so complete color to tz if tz < botMute
                         {
                             double tz = t0 + amp0 * (t1 - t0) / (amp0 - amp1);

                             if (tz < botMute)
                             {
                                 // complete color for first part of interval
                                 gl.glVertex2d(x0, tz);
                                 //                                gl.glVertex2d(x0, tz);
                                 gl.glVertex2d(x0, botMute);
                                 //                                gl.glVertex2d(x0 + ampMute, botMute);

                                 // switch to mute color for rest of interval
                                 outsideMuteColor.setGLColor(gl);

                                 gl.glVertex2d(x0, botMute);
                                 //                                gl.glVertex2d(x0 + ampMute, botMute);
                                 gl.glVertex2d(x0, t1);
                                 //                                gl.glVertex2d(x0 + amp1, t1);
                             }
                             else  // tz > botMute, so draw old color to botMute, mute color to tz and then to t1
                             {
                                 gl.glVertex2d(x0, botMute);
                                 //                                gl.glVertex2d(x0 + ampMute, botMute);

                                 // switch to mute color for rest of interval
                                 outsideMuteColor.setGLColor(gl);

                                 gl.glVertex2d(x0, botMute);
                                 //                                gl.glVertex2d(x0 + ampMute, botMute);
                                 gl.glVertex2d(x0, t1);
                                 //                                gl.glVertex2d(x0 + amp1, t1);
                             }
                         }
                     }
                 }
             }
             gl.glEnd();
             gl.glEnable(GL.GL_LIGHTING);
         }
         catch (Exception e)
         {
             StsException.outputException("StsPreStackLine.drawFilledWiggleTraces() failed.", e, StsException.WARNING);
         }
         finally
         {
             gl.glEnd();
             gl.glEnable(GL.GL_LIGHTING);
         }
     }

     static public int[] getGatherDataRange(float[][] axisRanges, double[][] traceData)
     {
         return new int[]{0, traceData.length - 1};
     }

     static public int[] getGatherDataRange(float[][] axisRanges, double[][] traceData, float minVelocity, boolean isFlattened)
     {
         return new int[]{0, traceData.length - 1};
         /*
             int nPoints = traceData.length - 1;
             double z1 = traceData[0][0];
             double z2 = traceData[nPoints - 1][0];
             double maxOffset = 30000;
             if(isFlattened)
                 z2 = Math.sqrt(z2 * z2 + maxOffset * maxOffset / (minVelocity * minVelocity));

             double dzGuess = (z2 - z1) / (nPoints - 1);
             float rangeZMin = axisRanges[1][1];
             float rangeZMax = axisRanges[1][0];

             int min = StsMath.floor((rangeZMin - z1) / dzGuess);
             min = StsMath.minMax(min, 0, nPoints - 1);
             double z = traceData[min][0];
             if(z < rangeZMin)
             {
                 while(z < rangeZMin && min < nPoints - 1)
                     z = traceData[++min][0];
                 min--;
             }
             else if(z > rangeZMin)
             {
                 while(z > rangeZMin && min > 0)
                     z = traceData[--min][0];
             }

             int max = StsMath.ceiling((rangeZMax - z1) / dzGuess);
             max = StsMath.minMax(max, 0, nPoints - 1);

             z = traceData[max][0];
             if(z < rangeZMax)
             {
                 while(z < rangeZMax && max < nPoints - 1)
                     z = traceData[++max][0];
             }
             else if(z > rangeZMax)
             {
                 while(z > rangeZMax && max > 0)
                     z = traceData[--max][0];
                 max++;
             }
             return new int[]{min, max};
         */
     }

     static public int computeInterpolationInterval(float zInc, int approximateNIntervals)
     {
         double[] scale = StsMath.niceScale(0.0, zInc, approximateNIntervals, true);
         return (int) Math.round(zInc / scale[2]);
     }

     static public double[][] getDisplayedData(double[][] points, int[] displayedDataRange, float valueScale, float valueOffset)
     {
         int min = displayedDataRange[0];
         int max = displayedDataRange[1];
         int nDisplayedPoints = max - min + 1;
         if (nDisplayedPoints < 0)
         {
             StsException.systemError(StsTraceUtilities.class, "getDisplayedData");
             return null;
         }
         double[][] displayedPoints = new double[nDisplayedPoints][2];
         for (int n = min, i = 0; n <= max; n++, i++)
         {
             displayedPoints[i][0] = points[n][0];
             displayedPoints[i][1] = points[n][1] * valueScale + valueOffset;
         }
         return displayedPoints;
     }

     static public float[] getDisplayZRange(float[][] axisRanges, float zInc)
     {
         float displayZMin = axisRanges[1][1] - 3 * zInc;
         float displayZMax = axisRanges[1][0] + 3 * zInc;
         return new float[]{displayZMin, displayZMax};
     }

     static public int[] getDisplayRange(StsRotatedGridBoundingBox displayBoundingBox, float[] displayZRange)
     {

         int displayMin = displayBoundingBox.getNearestBoundedSliceCoor(displayZRange[0]);
         int displayMax = displayBoundingBox.getNearestBoundedSliceCoor(displayZRange[1]);
         return new int[]{displayMin, displayMax};
     }

     static public double

			 [][] computeInterpolatedPoints(double[][] dataPoints, int nIntervals, float zInc)
     {
         return computeInterpolatedDataPoints(dataPoints, nIntervals, zInc, 0.0, 0.0);
     }

     static public double[][] computeInterpolatedDataPoints(double[][] dataPoints, int nIntervals, float zInc, double valueScale, double valueOffset)
     {
         int n = 0;

         try
         {
             if (dataPoints == null) return null;
             int nDataPoints = dataPoints.length;
             if (nDataPoints == 0) return dataPoints;

             int nTotalPoints = (nDataPoints - 1) * nIntervals + 1;
             tracePoints = new double[nTotalPoints][2];
             tracePointTypes = new byte[nTotalPoints];
             double t1 = dataPoints[0][0];
             double v1 = dataPoints[0][1];
             int p = 0;
             double df = 1.0 / (nIntervals);
             StsCubicForwardDifference cubicFd = new StsCubicForwardDifference(nIntervals);
             for (n = 0; n < nDataPoints - 1; n++)
             {
                 double t0 = t1;
                 t1 = dataPoints[n + 1][0];
                 double v0 = v1;
                 v1 = dataPoints[n + 1][1];
                 double dt = df * (t1 - t0);
                 double slopeScale = (t1 - t0) / zInc;

                 double f = df;
                 cubicFd.hermiteInitialize(dataPoints[n][1], dataPoints[n + 1][1], slopeScale * dataPoints[n][2], slopeScale * dataPoints[n + 1][2]);
                 double[] interpolatedPoints = cubicFd.evaluate();
                 double t = t0;
                 for (int i = 0; i < nIntervals; i++, p++)
                 {
                     tracePoints[p][0] = t;
                     t += dt;
                     if (valueScale != 0.0)
                         tracePoints[p][1] = interpolatedPoints[i] * valueScale + valueOffset;
                     else
                         tracePoints[p][1] = interpolatedPoints[i];

                 }
             }
             tracePoints[p] = dataPoints[nDataPoints - 1];
             p = 0;
             for (n = 0; n < nDataPoints - 1; n++)
             {
                 tracePointTypes[p++] = POINT_ORIGINAL;
                 for (int i = 0; i < nIntervals - 1; i++, p++)
                     tracePointTypes[p] = POINT_INTERPOLATED;
             }
             tracePoints[p] = dataPoints[nDataPoints - 1];
             return tracePoints;
         }
         catch (Exception e)
         {
             StsException.outputWarningException(StsTraceUtilities.class, "computeInterpolatedDataPoints", "Failed at point index: " + n, e);
             return null;
         }
     }

     static public boolean isZeroCrossingOrMinOrMax(byte pointType)
     {
         return pointType >= POINT_MAXIMUM && pointType <= POINT_MINUS_ZERO_CROSSING;
     }

     static private byte[] scratchByteData;

     static public float[] computeWigglePoints(ByteBuffer byteBuffer, float horizScale, int displayMin, int nSlices, int nDisplaySlices)
     {
         float[] points = new float[nDisplaySlices];
         if (scratchByteData == null || scratchByteData.length < nSlices)
             scratchByteData = new byte[nSlices];
         byteBuffer.get(scratchByteData, 0, nSlices);
         for (int n = 0, i = displayMin; n < nDisplaySlices; n++, i++)
             if (scratchByteData[i] != -1)
                 points[n] = horizScale * StsMath.unsignedByteToSignedInt(scratchByteData[i]);
         return points;
     }

     static public StsColor getPointTypeColor(byte pointType) { return typeColors[pointType]; }

     /**
      * linearly interpolates to get data value.
      * assumes data is regularly spaced in time (zInc is constant)
      *
      * @param traceData
      * @param t
      * @param zInc
      * @return
      */
     public static double getInterpolatedValue(float[] traceData, double t, double zInc)
     {
         int index = (int) Math.floor(t / zInc);
         if (index < traceData.length - 1)
         {
             double r = t - index * zInc;
             double v0 = traceData[index];
             double v1 = traceData[index + 1];
             return v0 + r * (v1 - v0) / zInc;
         }
         if (index >= traceData.length - 1) return traceData[traceData.length - 1];
         return 0;
     }

     public static void stack(float[] traces, int nTraces, float[] stackTrace)
     {
         stack(traces, 0, nTraces - 1, stackTrace);
     }


     /**
      * Interpolates value between samples along a trace.
      *
      * @param traces
      * @param t        time to be sampled (in samples)
      * @param trace
      * @param nSamples
      * @return
      */
     public static double getInterpolatedValue(float[] traces, double t, int trace, int nSamples) //, double zInc)
     {
         int index = (int) (t);///zInc); //Math.floor(tx/zInc);
         if (index < nSamples - 1)
         {
             double r = t - index;//*zInc;
             double v0 = traces[trace * nSamples + index];
             double v1 = traces[trace * nSamples + index + 1];
             return v0 + r * (v1 - v0);///zInc;
         }
         return traces[trace * nSamples + nSamples - 1];
     }

     public static void stack(float[] inputTraces, int iTrace0, int iTraceN, float[] outputTrace)
     {
         if (inputTraces == null || outputTrace == null) return;
         if (inputTraces.length == 0 || inputTraces.length == 0 || outputTrace.length == 0) return;
         int nSamples = outputTrace.length;
         int nTraces = iTraceN - iTrace0 + 1;
         if (nTraces * nSamples > inputTraces.length)
         {
             StsMessage.printMessage("StsTraceUtilities.stack(): nTraces*nSamples >= traces.length");
             return;
         }
         double sum = 0;
         double val = 0;
         int nVals = 0;
         for (int i = 0; i < nSamples; i++)
         {
             sum = 0;
             nVals = 0;
             for (int j = iTrace0; j < iTraceN; j++)
             {
                 val = inputTraces[j * nSamples + i];
                 if (val == 0) continue;
                 sum += val;
                 nVals += 1;
             }
             if (nVals > 0)
                 outputTrace[i] = (float) (sum / nVals);
             else
                 outputTrace[i] = 0;
         }

     }

     public static void stack(float[][] inputTraces, int iTrace0, int iTraceN, float[] outputTrace)
     {
         if (inputTraces == null || outputTrace == null) return;
         if (inputTraces.length < iTraceN || outputTrace.length == 0) return;
         int nSamples = outputTrace.length;
         int nTraces = iTraceN - iTrace0 + 1;
         if (nTraces > inputTraces.length)
         {
             StsMessage.printMessage("StsTraceUtilities.stack(): nTraces > inputTraces.length");
             return;
         }
         if (inputTraces[iTrace0] != null && nSamples > inputTraces[iTrace0].length)
         {
             StsMessage.printMessage("StsTraceUtilities.stack(): nSamples > inputTraces[0].length");
             return;
         }
         double sum = 0;
         double val = 0;
         int nVals = 0;
         for (int i = 0; i < nSamples; i++)
         {
             sum = 0;
             nVals = 0;
             for (int j = iTrace0; j < iTraceN; j++)
             {
                 if (inputTraces[j] == null) continue;
                 val = inputTraces[j][i];
                 if (val == 0) continue;
                 sum += val;
                 nVals += 1;
             }
             if (nVals > 0)
                 outputTrace[i] = (float) (sum / nVals);
             else
                 outputTrace[i] = 0;
         }

     }
 }
