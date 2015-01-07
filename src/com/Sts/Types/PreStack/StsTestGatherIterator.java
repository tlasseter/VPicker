package com.Sts.Types.PreStack;

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.Types.*;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: May 23, 2008
 * Time: 6:09:32 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsTestGatherIterator extends StsAbstractGatherIterator implements Iterator
{
    StsPreStackLineSet lineSet;
    StsPreStackLine[] lines;
    int nLines;
    int nLine;
    int nLineGathers;
    StsPreStackLine currentLine;
    int minGatherIndex;
    int maxGatherIndex;
    int nGather;
    StsGather gather;

    StsMappedDoubleBuffer lineAttributesArrayBuffer;

    public StsTestGatherIterator(StsPreStackLineSet testLineSet)
    {
        super(testLineSet);
    }

    public void initialize()
    {
        lineSet = (StsPreStackLineSet)dataLineSet;
        lines = lineSet.lines;
        nLines = lines.length;
        gather = new StsGather(lineSet);
        setNextLine();
    }

    private boolean setNextLine()
    {
        while(nLine < nLines-1)
        {
            nLine++;
            currentLine = lines[nLine];
            if(currentLine == null) continue;
            minGatherIndex = currentLine.minGatherIndex;
            maxGatherIndex = currentLine.maxGatherIndex;
            nGather = minGatherIndex;
            if(maxGatherIndex < minGatherIndex) continue;
            int nLineTraces = currentLine.getNTraces();
            lineAttributesArrayBuffer = currentLine.getAllAttributesArrayBuffer(nLineTraces);
            return true;
        }
        return false;
    }

    public boolean hasNext()
    {
        if(nGather > maxGatherIndex)
            if(!setNextLine()) return false;
        return setNextGather();
    }

    private boolean setNextGather()
    {
        while(nGather <= maxGatherIndex)
        {
            gather.initializeLineGather(nLine, nGather);
            nGather++;
            if(gather.nGatherTraces > 0)
            {
                gatherData.setTraceData(gather.getGatherDataf());
                double[] traceOrderedGatherAttributes = gather.getTraceOrderedAttributes(lineAttributesArrayBuffer);
                gatherData.setTraceOrderedAttributes(traceOrderedGatherAttributes);
                return true;
            }
        }
        return false;
    }
}