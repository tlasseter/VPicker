package com.Sts.Types.Seismic;

import com.Sts.Utilities.Seismic.*;
import com.Sts.Utilities.*;

/**
 * Created by IntelliJ IDEA.
  * User: Tom Lasseter
  * Date: Nov 18, 2009
  * Time: 2:40:42 PM
  * To change this template use File | Settings | File Templates.
  */
 public class StsTrace
 {
     /** data for traces */
     public float[] data;
     public byte[] pointTypes;
     public float[] instantAmplitude;
     public float[] instantPhase;
     public float[] hilbert;
     static final boolean debug = false;

     public StsTrace(float[] data)
     {
         this.data = data;
     }

     public StsTrace(float[] gatherData, int offset, int nValues)
     {
         if (data == null || data.length != nValues)
             data = new float[nValues];
         System.arraycopy(gatherData, offset, data, 0, nValues);
     }

     public void interpolateData(int nInterpolationIntervals, float zInc)
     {
         data = StsTraceUtilities.computeCubicInterpolatedPoints(data, nInterpolationIntervals);
     }

     public void computePointTypes()
     {
         if (pointTypes != null) return;
         pointTypes = StsTraceUtilities.defineTracePointTypes(data);
     }

     public void scaleData(float scaleMultiplier)
     {
         StsMath.scale(data, scaleMultiplier);
     }

     public synchronized boolean computeInstantaneousAmpAndPhaseGatherData()
     {
         try
         {
             computePointTypes();
             int nSamples = data.length;
             instantAmplitude = new float[nSamples];
             instantPhase = new float[nSamples];
             if (debug) StsException.systemDebug(this, "computeGatherData", "Computing gatherData for gather view");
             if(!StsTraceUtilities.computeInstantAmpAndPhase(data, pointTypes, instantAmplitude, instantPhase)) return false;
             computeComplexComponents();
             return true;
         }
         catch (Exception e)
         {
             StsException.outputWarningException(this, "computeInstantaneousAmpAndPhaseGatherData", e);
             return false;
         }
     }

     private void computeComplexComponents()
     {
        int nSamples = data.length;
        hilbert = new float[nSamples];
        for(int n = 0; n < nSamples; n++)
             hilbert[n] = (float)StsMath.sind(instantPhase[n])*instantAmplitude[n];
     }
     
     public synchronized boolean computeHilbertAmpAndPhase()
     {
         try
         {
             computeHilbertTransform();
             instantAmplitude = StsHilbertTransform.computeAmplitude(hilbert, data);
             instantPhase = StsHilbertTransform.computePhase(hilbert, data);
             if (debug) StsException.systemDebug(this, "computeGatherData", "Computing gatherData for gather view");
             return true;
         }
         catch (Exception e)
         {
             StsException.outputWarningException(this, "computeInstantaneousAmpAndPhaseGatherData", e);
             return false;
         }
     }

     public synchronized boolean computeHilbertTransform()
     {
         StsHilbertTransform.staticInitialize();
         return calculateHilbertTransform();
     }

     public synchronized boolean computeHilbertTransform(int halfWindowSize)
     {
         StsHilbertTransform.staticInitialize(halfWindowSize);
         return calculateHilbertTransform();
     }

     private boolean calculateHilbertTransform()
     {
         try
         {
             hilbert = StsHilbertTransform.computeHilbert(data);
             return true;
         }
         catch (Exception e)
         {
             StsException.outputWarningException(this, "computeInstantaneousAmpAndPhaseGatherData", e);
             return false;
         }
     }

     public double computeRmsAmplitude()
     {
         return StsMath.rmsIgnoreZero(data);

     }

     public void normalizeAmplitude(double amplitude)
     {
         StsMath.normalizeAmplitude(data, amplitude);
     }

     public double getValue(double f, int nDataSamples)
     {
         if (f < 0.0 || f > nDataSamples) return StsParameters.doubleNullValue;
         double value;
         int index = StsMath.floor(f);
         if (index == nDataSamples - 1)
         {
             value = data[index];
         }
         else
         {
             f = f - index;
             float v0 = data[index];
             float v1 = data[index + 1];
             value = v0 + f * (v1 - v0);
         }
         return value;
     }

     public double[] getComplexComponents(double f, int nDataSamples)
     {
         if (f < 0.0 || f > nDataSamples) return new double[] { 0.0, StsParameters.doubleNullValue };
         int index = StsMath.floor(f);
         if (index >= nDataSamples - 1)
             return getComplexComponents(index);
         else
         {
             f = f - index;
             return getComplexComponents(index, f);
         }
     }

     private double[] getComplexComponents(int index)
     {
            return new double[]{data[index], hilbert[index]};
     }

     private double[] getComplexComponents(int index, double f)
     {
         double data0 = data[index];
         double data1 = data[index + 1];
         double d = data0 + f * (data1 - data0);
         double h0 = hilbert[index];
         double h1 = hilbert[index + 1];
         double h = h0 + f * (h1 - h0);
         return new double[]{d, h};
     }

     public double[] getPhaseAmp(double f, int nDataSamples)
     {
         if (f < 0.0 || f > nDataSamples) return new double[] { 0.0, StsParameters.doubleNullValue };
         int index = StsMath.floor(f);
         if (index >= nDataSamples - 1)
             return getPhaseAmp(index);
         else
         {
             f = f - index;
             return getPhaseAmp(index, f);
         }
     }

     private double[] getPhaseAmp(int index)
     {
         return new double[]{instantPhase[index], instantAmplitude[index]};
     }

     private double[] getPhaseAmp(int index, double f)
     {
         double amp0 = instantAmplitude[index];
         double amp1 = instantAmplitude[index + 1];
         double amplitude = amp0 + f * (amp1 - amp0);
         double phase0 = instantPhase[index];
         double phase1 = instantPhase[index + 1];
         double phase = phase0 + f * (phase1 - phase0);
         return new double[]{phase, amplitude};
     }
 }