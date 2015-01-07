package com.Sts.Types.PreStack;

import com.Sts.Types.Seismic.*;
import com.Sts.Types.Seismic.StsTrace;
import com.Sts.Types.Seismic.StsTraceSet;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

 import java.util.*;


/**
 * Created by IntelliJ IDEA.
  * User: Tom Lasseter
  * Date: Dec 19, 2009
  * Time: 8:46:05 AM
  * To change this template use File | Settings | File Templates.
  */
 public class StsGatherTraceSet extends StsTraceSet
 {
     StsGather gather;
     int minOffsetTraceIndex;
     StsGatherTrace[] absIndexTraces;

     public StsGatherTraceSet(StsGather gather)
     {
         super(gather.lineSet.currentLine);
         traces = gather.constructGatherTraces();
         this.gather = gather;
         nTraces = gather.nGatherTraces;
         setMinimumOffsetTrace();
     }

     public StsGatherTraceSet(StsGather gather, int nApproxInterpIntervals)
     {
         super(gather.lineSet.currentLine, nApproxInterpIntervals);
         traces = gather.constructGatherTraces();
         this.gather = gather;
         nTraces = gather.nGatherTraces;
         setMinimumOffsetTrace();
     }

     public StsGatherTraceSet(StsRotatedGridBoundingBox boundingBox, int nApproxInterpIntervals)
     {
         super(boundingBox, nApproxInterpIntervals);
     }

     StsGatherTrace[] getGatherTraces() { return (StsGatherTrace[]) traces; }

     public void setMinimumOffsetTrace()
     {
         double minOffsetSq = StsParameters.largeDouble;
         for (int n = 0; n < traces.length; n++)
         {
             StsGatherTrace trace = (StsGatherTrace) traces[n];
             if (trace.xsq < minOffsetSq)
             {
                 minOffsetTraceIndex = n;
                 minOffsetSq = trace.xsq;
             }
         }
     }

     public void constructAbsIndexTraces()
     {
         absIndexTraces = (StsGatherTrace[]) StsMath.arraycopy(traces);
         Arrays.sort(absIndexTraces, StsGatherTrace.getAbsOffsetComparator());
     }

     public double getTOffsetValue(double tOffset, StsTrace trace)
     {
         double f;
         int nDataSamples;
         if (isInterpolated)
         {
             f = (tOffset - zMin) / zIncInterpolated;
             nDataSamples = nInterpolatedSamples;
         }
         else
         {
             f = (tOffset - zMin) / zInc;
             nDataSamples = nSamples;
         }
         return trace.getValue(f, nDataSamples);
     }

     public double[] getPhaseAmp(double tOffset, StsTrace trace)
     {
         double f;
         int nDataSamples;
         if(tOffset < zMin || tOffset > zMax) return null;
         if (isInterpolated)
         {
             f = (tOffset - zMin) / zIncInterpolated;
             nDataSamples = nInterpolatedSamples;
         }
         else
         {
             f = (tOffset - zMin) / zInc;
             nDataSamples = nSamples;
         }
         return trace.getPhaseAmp(f, nDataSamples);
     }

     public double[] getComplexComponents(double tOffset, StsTrace trace)
     {
         double f;
         int nDataSamples;
         if(tOffset < zMin || tOffset > zMax) return null;
         if (isInterpolated)
         {
             f = (tOffset - zMin) / zIncInterpolated;
             nDataSamples = nInterpolatedSamples;
         }
         else
         {
             f = (tOffset - zMin) / zInc;
             nDataSamples = nSamples;
         }
         return trace.getComplexComponents(f, nDataSamples);
     }
 }
