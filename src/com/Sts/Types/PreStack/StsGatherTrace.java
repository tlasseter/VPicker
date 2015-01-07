package com.Sts.Types.PreStack;

import com.Sts.Types.Seismic.*;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Dec 22, 2009
 * Time: 2:50:47 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsGatherTrace extends StsTrace
{
    /** index in offset order (max minus to max plus) */
    public int index;
    /** offset */
    public double x;
    /** offset squared */
    public double xsq;

    static public Comparator traceComparator;
    static public Comparator traceAbsOffsetComparator;
    static
    {
        traceComparator = new TraceComparator();
        traceAbsOffsetComparator = new TraceAbsOffsetComparator();
    }

    public StsGatherTrace(float[] data)
    {
        super(data);
    }

    public StsGatherTrace(float[] gatherData, int dataOffset, int nValues, double traceOffset)
    {
        super(gatherData, dataOffset, nValues);
        x = traceOffset;
        xsq = x*x;
    }

    static public StsGatherTrace[] constructGatherTraces(StsGather gather)
    {
        int nSlices = gather.line.nSlices;
        int nTraces = gather.nGatherTraces;
        StsGatherTrace[] traces = new StsGatherTrace[nTraces];
        float[] dataF = gather.getGatherDataf();
        double[] traceOffsets = gather.getGatherOffsets();
        if (dataF != null)
        {
            for(int t = 0, i = 0; t < nTraces; t++, i += nSlices)
                traces[t] = new StsGatherTrace(dataF, i, nSlices, traceOffsets[t]);
        }
        return traces;
    }

    static public Comparator getComparator()
    {
        return traceComparator;
    }

    static class TraceComparator implements Comparator
    {
        TraceComparator()
        {

        }

        public int compare(Object t1, Object t2)
        {
            double x1 = ((StsGatherTrace)t1).x;
            double x2 = ((StsGatherTrace)t2).x;
            if(x1 > x2) return 1;
            else if(x1 < x2) return -1;
            else return 0;
        }
    }

    static public Comparator getAbsOffsetComparator()
    {
        return traceAbsOffsetComparator;
    }

    static class TraceAbsOffsetComparator implements Comparator
    {
        TraceAbsOffsetComparator()
        {

        }

        public int compare(Object t1, Object t2)
        {
            double x1 = Math.abs(((StsGatherTrace)t1).x);
            double x2 = Math.abs(((StsGatherTrace)t2).x);
            if(x1 > x2) return 1;
            else if(x1 < x2) return -1;
            else return 0;
        }
    }
}
