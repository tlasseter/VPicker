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
public class StsBlockCell
{
    byte type;
    ArrayList gridPolygons = new ArrayList();

    static public final byte EMPTY = StsParameters.CELL_EMPTY;
    static public final byte FULL = StsParameters.CELL_FULL;
    static public final byte EDGE = StsParameters.CELL_EDGE;

    void addPolygon(StsPolygon layerPolygon)
    {
        gridPolygons.add(layerPolygon);
        type = EDGE;
    }
}