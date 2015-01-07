package com.Sts.Types.PreStack;

import com.Sts.Types.Seismic.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Dec 19, 2009
 * Time: 8:46:05 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsSuperGatherTraceSet extends StsGatherTraceSet
{
    StsSuperGather superGather;
    StsGatherTraceSet[] gatherTraceSets;

    public StsSuperGatherTraceSet(StsSuperGather superGather, int nApproxInterpIntervals)
    {
        super(superGather.lineSet.currentLine, nApproxInterpIntervals);
        this.superGather = superGather;
        this.gather = superGather.centerGather;
        StsGather[] gathers = superGather.gathers;
        if (gathers == null || gathers.length == 0) return;
        int nGathers = gathers.length;
        gatherTraceSets = new StsGatherTraceSet[nGathers];
        for(int n = 0; n < nGathers; n++)
        {
            gatherTraceSets[n] = new StsGatherTraceSet(gathers[n], nApproxInterpIntervals);
            nTraces += gathers[n].nGatherTraces;
            traces = (StsGatherTrace[])StsMath.arrayAddArray(traces, gatherTraceSets[n].getGatherTraces());
        }
        Arrays.sort(traces, StsGatherTrace.getComparator());
        setMinimumOffsetTrace();
    }

    public StsGatherTrace[] getGatherTraces(int offsetAxisType)
    {
        if(offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_ABS_INDEX)
        {
            if(absIndexTraces == null) constructAbsIndexTraces();
            return absIndexTraces;
        }
        else
            return (StsGatherTrace[])traces;
    }

    public void rmsNormalizeAmplitudes()
    {
        double rmsAmplitude = computeRmsAmplitude();
        normalizeAmplitude(rmsAmplitude);
    }
}
