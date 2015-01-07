
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

/** StsFaultGap is a temporary data structure which holds information about
 *  the gap geometry along a fault */

public class StsFaultGap
{
    public int rowOrCol;        /** point is on a ROW or COL                                         */
    public int rightDirection;  /** The direction away from the right side of the fault              */
    public float rowF, colF;    /** Float coordinates of the gap centerLine point                    */
    public StsSection section;  /** The section this point is on                                     */
    public float sectionColF;   /** The distance in positive direction along the section (0.0 to nCols-1 */
    public int minIndex;        /** Index of max range in minus direction                            */
    public int maxIndex;        /** Index of max range in plus direction                             */
    public float minIndexF;     /** Index of next edge or grid side in minus direction               */
    public float maxIndexF;     /** Index of next edge or grid side in plus direction                */

    // Convenience copies of flags
    static final int ROW = StsParameters.ROW;
    static final int COL = StsParameters.COL;
    static final int MINUS = StsParameters.MINUS;
    static final int PLUS = StsParameters.PLUS;

// constructTraceAnalyzer

    public StsFaultGap(StsEdgeLoopRadialGridLink link, StsSection section, int rowOrCol)
    {
        this.section = section;
        this.rowOrCol = rowOrCol;
        this.rightDirection = link.getInsideDirection(rowOrCol);
        this.rowF = link.getRowF();
        this.colF = link.getColF();
        this.sectionColF = link.getPoint().getColF(section);
    }
}



