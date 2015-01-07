package com.Sts.DBTypes;

import com.Sts.Types.*;
import com.Sts.Utilities.*;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Mar 15, 2010
 * Time: 12:54:25 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsBlockCellColumn
{
    byte type;
    LayerInterval[] layerIntervals;
    ArrayList gridPolygons = new ArrayList();

    static public final byte EMPTY = StsParameters.CELL_EMPTY;
    static public final byte FULL = StsParameters.CELL_FULL;
    static public final byte EDGE = StsParameters.CELL_EDGE;

    void addLayerRange(StsPolygon gridPolygon, int[] layerRange)
    {
        addLayerInterval(new LayerInterval(layerRange));
        gridPolygons.add(gridPolygon);
    }

    void addLayerInterval(LayerInterval layerInterval)
    {
        if(layerIntervals == null)
        {
            layerIntervals = new LayerInterval[] { layerInterval };
            return;
        }
        for(int n = 0; n < layerIntervals.length; n++)
        {
            byte position = layerInterval.getPosition(layerIntervals[n]);
            if(position == BEFORE)
            {
                layerIntervals = (LayerInterval[]) StsMath.arrayInsertElementBefore(layerIntervals, layerInterval, n);
                return;
            }
            else if(position == OVERLAP)
            {
                layerIntervals[n].combineLayerInterval(layerInterval);
                return;
            }
        }
        layerIntervals = (LayerInterval[]) StsMath.arrayAddElement(layerIntervals, layerInterval);
    }

    static final byte BEFORE = -1;
    static final byte EQUAL = 0;
    static final byte AFTER = 1;
    static final byte OVERLAP = 2;

    class LayerInterval
    {
        int topLayerNumber;
        int botLayerNumber;

        LayerInterval(int[] layerRange)
        {
            topLayerNumber = layerRange[0];
            botLayerNumber = layerRange[1];
        }

        byte getPosition(LayerInterval otherLayerInterval)
        {
            if(topLayerNumber == otherLayerInterval.topLayerNumber)
            {
                if(botLayerNumber >= otherLayerInterval.botLayerNumber)
                    return EQUAL;
                else
                    return OVERLAP;
            }
            else if(topLayerNumber < otherLayerInterval.topLayerNumber)
            {
                if(botLayerNumber >= otherLayerInterval.botLayerNumber)
                    return EQUAL;
                if(topLayerNumber >= otherLayerInterval.botLayerNumber-1)
                    return OVERLAP;
                else
                    return BEFORE;
            }
            else // topLayerNumber > otherLayerInterval.topLayerNumber
            {
               if(topLayerNumber > otherLayerInterval.botLayerNumber+1)
                    return AFTER;
                else
                    return OVERLAP;
            }
        }

        void combineLayerInterval(LayerInterval layerInterval)
        {
            topLayerNumber = Math.min(topLayerNumber, layerInterval.topLayerNumber);
            botLayerNumber = Math.max(botLayerNumber, layerInterval.botLayerNumber);
        }
    }
}
