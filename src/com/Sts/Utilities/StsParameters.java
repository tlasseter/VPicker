//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Utilities;

import com.Sts.DBTypes.*;

public class StsParameters
{
	/* Changed to max values of types to allow wider range of numbers needed for sensors where 1e-20 is common
    public static final float roundOff = 1.e-5f;
    public static final float doubleRoundOff = 1.e-10f;
    public static final float nullValue = 1.e30f;
    public static final float largeFloat = 1.e+20f;
    public static final float veryLargeFloat = 1.e+30f;
    public static final long largeLong = Long.MAX_VALUE;
    public static final long nullLongValue = -Long.MAX_VALUE;
    public static final double largeDouble = 1.e+10;
    public static final float smallFloat = 1.e-10f;
	public static final int largeInt = 100000;
    public static final int nullInteger = -99999;
    public static final double doubleNullValue = 1.e30;	
	*/
    public static final float roundOff = 1.e-5f;
    public static final float doubleRoundOff = 1.e-10f;
    public static final float nullValue = Float.MAX_VALUE;
    public static final float largeFloat = Float.MAX_VALUE;
    public static final float veryLargeFloat = Float.MAX_VALUE;
    public static final long largeLong = Long.MAX_VALUE;
    public static final long nullLongValue = -Long.MAX_VALUE;
    public static final double largeDouble = Double.MAX_VALUE;
    public static final float smallFloat = Float.MIN_VALUE;
	public static final int largeInt = Integer.MAX_VALUE;
    public static final int nullInteger = Integer.MIN_VALUE;
    public static final double doubleNullValue = Double.MAX_VALUE;
    
    /**/
    public static final double sqrt2 = 1.41421356237309504880;
    public static final float sqrt2f = 1.41421356f;
    public static final double halfSqrt2 = 0.70710678118654752440;
    public static final float halfSqrt2f = 0.70710678f;

    /** Used primarily by StsList to indicate item not found in list. */
    public static final int NO_MATCH = -1;

    /** Model status flags set by StsModel, set by StsBuildSections, etc */
    static public final byte NONE = 0;
    static public final int DATA = 1;
    static public final int FRAME = 2;
    static public final int MODEL = 3;

    /** Grid directions used by StsSectionPoint, StsSurfaceVertex */
    static public final int MINUS = -1;
    static public final int PLUS = 1;
    static public final int PLUS_AND_MINUS = 2;

    static public final int OK = 1;
    static public final int NOT_OK = -1;

    /** Next or prev used by  by StsHorizonPointLinkedGrid */
    static public final int PREV = -1;
    static public final int NEXT = 1;

    static public final String NONE_STRING = "none";

    /** unsigned bytes range in value from 0 to 254; 255 if a nullByte
     *  Because Java doesn't support unsigned bytes, these are equivalent to:
     *  unsigned 0 is signed 0; unsigned 254 is signed -2; unsigned 255 is signed -1
     *  So we have the definitions below:
     */
    static public final byte nullByte = -1;
    static public final byte nullMinByte = 0;
    static public final byte nullMaxByte = -2;

    /** Side flags used by StsSection and others */
    public static final int RIGHT = 1;
    public static final int LEFT = -1;

    public static final String sideLabel(int side)
    {
        switch (side)
        {
            case NONE:
                return new String("side: NONE");
            case RIGHT:
                return new String("side: RIGHT");
            case LEFT:
                return new String("side: LEFT");
            default:
                return new String("side: UNDEFINED");
        }
    }

    /** Well types: used by StsWell, StsWellLine, StsSection */
    public static final byte WELL = 0;
    public static final byte PSEUDO = 1;
    public static final byte FAULT = 2;
    public static final byte HORIZ = 3;
    public static final byte PSEUDO_HORIZ = 4;
    public static final byte BOUNDARY = 5;
    public static final byte FRACTURE = 6;
    public static final byte AUXILIARY = 7;
    public static final byte REFERENCE = 8;

     public static final byte RIBBON = 9;
     // A curved section is two wells connected by curved sectionEdges
     public static final byte CURVED = 10;

     public static final String WELL_STRING = "Well";
     public static final String PSEUDO_STRING = "Pseudo";
     public static final String FAULT_STRING = "Fault";
     public static final String HORIZ_STRING = "Horizontal";
     public static final String PSEUDO_HORIZ_STRING = "Pseudo-horizontal";
     public static final String BOUNDARY_STRING = "Boundary";
     public static final String FRACTURE_STRING = "Fracture";
     public static final String AUXILIARY_STRING = "Auxiliary";
     public static final String REFERENCE_STRING = "Reference";
     // A ribbon section has one or two non-vertical wells connected by straight sectionEdges
     public static final String RIBBON_STRING = "Ribbon";
     // A curved section is two wells connected by curved sectionEdges
     public static final String CURVED_STRING = "Curved";

     public static final String[] typeStrings = new String[]{WELL_STRING, PSEUDO_STRING, FAULT_STRING, HORIZ_STRING, PSEUDO_HORIZ_STRING,
         BOUNDARY_STRING, FRACTURE_STRING, AUXILIARY_STRING, REFERENCE_STRING, RIBBON_STRING, CURVED_STRING};

     public static final String typeToString(int type)
     {
         return typeStrings[type];
     }

     /** Panel orientation */
     public static final int HORIZONTAL = 1;
     public static final int VERTICAL = 2;

    /** Grid point types: used in StsGrid, StsGridPoint */
     static public final byte SURF_NONE = 0;
     /** good point */
     static public final byte SURF_PNT = 1;
     /** nulled by gapping routine */
     static public final byte SURF_GAP_SET = 2;
     /** originally null */
     static public final byte SURF_GAP = 3;
     /** originally null and is outside boundary */
     static public final byte SURF_BOUNDARY = 4;
     /** NULL_GAP type filled in. */
     static public final byte SURF_GAP_FILLED = 5;
     /** Needs to be filled by interpolation */
     static public final byte SURF_GAP_NOT_FILLED = 7;
     /** uninitialized */
     static public final String SURF_NONE_STRING = "SURF_NONE";
     /** good point */
     static public final String SURF_PNT_STRING = "SURF_PNT";
     /** nulled by gapping routine */
     static public final String SURF_GAP_SET_STRING = "SURF_GAP_SET";
     /** originally null */
     static public final String SURF_GAP_STRING = "SURF_GAP";
     /** originally null and is outside boundary */
     static public final String SURF_BOUNDARY_STRING = "SURF_BOUNDARY";
     /** NULL_GAP type filled in. */
     static public final String SURF_GAP_FILLED_STRING = "SURF_GAP_FILLED";
     /** Needs to be filled by interpolation */
     static public final String SURF_GAP_NOT_FILLED_STRING = "SURF_GAP_NOT_FILLED";

     static public final String[] surfacePointTypeNames = new String[]{SURF_NONE_STRING, SURF_PNT_STRING, SURF_GAP_SET_STRING, SURF_GAP_STRING,
         SURF_BOUNDARY_STRING, SURF_GAP_FILLED_STRING, SURF_GAP_NOT_FILLED_STRING};

     static public final StsColor[] surfacePointTypeColors = new StsColor[]{StsColor.GREY, StsColor.BLUE, StsColor.PURPLE, StsColor.CYAN,
         StsColor.RED, StsColor.ORANGE, StsColor.MAGENTA};

     public static final String getSurfacePointTypeName(byte nullType)
     {
         if (nullType <= 0 || nullType > SURF_GAP_NOT_FILLED)
             return "Unidentified: " + nullType;
         else
             return surfacePointTypeNames[nullType];
     }

     public static final StsColor getSurfacePointTypeColor(byte nullType)
     {
         if (nullType <= 0 || nullType > SURF_GAP_NOT_FILLED)
             return StsColor.BLACK;
         else
             return surfacePointTypeColors[nullType];
     }

    // pointsType flags used by StsBlockGrid
    static public final byte GAP_NULL = 0; // Not initialized
    static public final byte GAP_GRID = 1; // gridPoint for interpolating gapPoints
    static public final byte GAP_SURF_GRID = 2; // gridPoint from original surface: not used now
    static public final byte GAP_NOT_FILLED = 3; // gapPoint we need to fill
    static public final byte GAP_FILLED = 4; // filled by least-sqs interpolated value
    static public final byte GAP_CANT_FILL = 5; // couldn't compute least-sq value
    static public final byte GAP_CUT = 6; // gapPoint, but has been cutoff by fault; not used
    static public final byte GAP_EXTRAP = 9; // extrapolated by least-sqs interpolated value
    static public final byte GAP_NONE = 10; // passed when type not to be changed
    static public final byte GAP_FILL_HOLE = 11; // back-filled hole (extrapolated): not used

     static public final String GAP_NULL_STRING = "GAP_NULL"; // Not initialized
     static public final String GAP_GRID_STRING = "GAP_GRID"; // gridPoint for interpolating gapPoints
     static public final String GAP_SURF_GRID_STRING = "GAP_SURF_GRID"; // gridPoint from original surface: not used now
     static public final String GAP_NOT_FILLED_STRING = "GAP_NOT_FILLED"; // gapPoint we need to fill
     static public final String GAP_FILLED_STRING = "GAP_FILLED"; // filled by least-sqs interpolated value
     static public final String GAP_CANT_FILL_STRING = "GAP_CANT_FILL"; // couldn't compute least-sq value
     static public final String GAP_CUT_STRING = "GAP_CUT"; // gapPoint, but has been cutoff by fault; not used
     static public final String GAP_EXTRAP_STRING = "GAP_EXTRAP"; // extrapolated by least-sqs interpolated value
     static public final String GAP_NONE_STRING = "GAP_NONE"; // passed when type not to be changed
     static public final String GAP_FILL_HOLE_STRING = "GAP_FILL_HOLE"; // back-filled hole (extrapolated): not used

     static public final String[] gapPointTypeNames = new String[]{GAP_NULL_STRING, GAP_GRID_STRING, GAP_SURF_GRID_STRING,
         GAP_NOT_FILLED_STRING, GAP_FILLED_STRING, GAP_CANT_FILL_STRING, GAP_CUT_STRING, GAP_NONE_STRING, GAP_NONE_STRING, GAP_EXTRAP_STRING, GAP_NONE_STRING, GAP_FILL_HOLE_STRING};

     static public final StsColor[] gapPointTypeColors = new StsColor[]{StsColor.GREY, StsColor.BLUE, StsColor.GREY,
         StsColor.PURPLE, StsColor.CYAN, StsColor.BROWN, StsColor.MAGENTA, StsColor.RED, StsColor.GREY, StsColor.ORANGE};

     public static final String getGapTypeName(byte nullType)
     {
         if (nullType < 0 || nullType > GAP_FILL_HOLE)
             return "Unidentified: " + nullType;
         else
             return gapPointTypeNames[nullType];
     }

     public static final StsColor getGapPointTypeColor(byte nullType)
     {
         if (nullType <= 0 || nullType > GAP_FILL_HOLE)
             return StsColor.BLACK;
         else
             return gapPointTypeColors[nullType];
     }


    public static final byte CELL_EMPTY = 0;
    public static final byte CELL_FULL = 1;
    public static final byte CELL_EDGE = 2;

    public static final String cellTypeString(byte cellType)
    {
        switch (cellType)
        {
            case CELL_EMPTY:
                return new String("CELL_EMPTY");
            case CELL_FULL:
                return new String("CELL_FILL");
            case CELL_EDGE:
                return new String("CELL_EDGE");
            default:
                return new String("Unknown cellType ");
        }
    }

    /** Grid orientation: used by StsSurfaceVertex */
    static public final int ROW = 1;
    static public final int COL = 2;
    static public final int ROWCOL = 3;

    static public final String ROW_STRING = "ROW";
    static public final String COL_STRING = "COL";
    static public final String ROWCOL_STRING = "ROW&COL";
    static public final String[] rowOrColStrings = new String[] { "None", ROW_STRING, COL_STRING, ROWCOL_STRING };

    public static final String getRowOrColString(int rowOrCol)
    {
        return rowOrColStrings[rowOrCol];
    }

    public static final String rowCol(int rowOrCol)
    {
        switch (rowOrCol)
        {
            case ROW:
                return new String("Row: ");
            case COL:
                return new String("Col: ");
            case ROWCOL:
                 return new String("RowCol: ");
            default:
                return new String("No rowCol ");
        }
    }

    public static final int getCrossingRowOrCol(int rowOrCol)
    {
        switch (rowOrCol)
        {
            case ROW:
                return COL;
            case COL:
                return ROW;
            default:
                return NONE;
        }
    }

    /** TStrip flags: used by StsGrid and StsLineLinkGrid */
    static public final byte STRIP_INVALID = 0;
    static public final byte STRIP_BOTH = 1;
    static public final byte STRIP_BOT = 2;
    static public final byte STRIP_TOP = 3;
    static public final byte STRIP_BOT_TOP = 4;

    /** Position flags used by StsSection.compareWellPositions */
    public static final int BEFORE = -1;
    public static final int AFTER = 1;
    public static final int UNKNOWN = 0;

    /** Comparison flags used in comparing StsObjects; see StsLink in LinkLists */
    static public final int LESS = -1;
    static public final int EQUALS = 0;
    static public final int GREATER = 1;
    static public final int UNDEFINED = -99;

    /** Type flags used by StsLineLinkedGrid */
    static public final byte GAP_LINE = 1; /** These points are in the gap  */
    static public final byte GRID_LINE = 2; /** These points are on the grid */
    static public final byte BOUNDARY_LINE = 3; /** These points are on the grid/gap boundary */

    /** Type flags used by StsGridLine
     *  flags must be negative since >= 0 indicates an edgeGroupIndex
     */
    static public final int GRIDLINE_GRID = -1;
    static public final int GRIDLINE_BOUNDARY = -2;
    static public final int GRIDLINE_NONE = -99;

    /** Skip flags used by StsSurfaceEdge.computeEdgePoints */
    static public final int FIRST = 1;
    static public final int LAST = 2;

    /** Loop type flags used by StsEdgeGroup */
    static public final byte LOOP = 1;
    static public final byte OPEN_LOOP = 2;

    /** Edge types used by StsGridSectionPoint */
    static public final byte X_LINE = 1;
    static public final byte Y_LINE = 2;
    static public final byte WELL_LINE = 3;
    static public final byte HORIZON_LINE = 4;
    static public final byte SUBHORIZON_LINE = 5;

    /** Direction types used by StsSection and pointLinkedGrid polygon construction */
    static public final int INSIDE = 1;
    static public final int OUTSIDE = -1;

    /** StsLinkedGrid types used by StsGridLink, StsLinkedGrid, and StsGridSectionPoint */
    static public final byte SECTION = 1;
    static public final byte HORIZON = 2;
    static public final byte SECTION_VERTICAL = 3;

    static public final byte LINK_NULL = 0;
    static public final byte LINK = 1;
    static public final byte LINK_SPLIT = 2;

    /** StsGridSectionPoint linked points direction indices used by StsLinkedGrid */
    static public final byte ROW_PLUS = 0;
    static public final byte COL_PLUS = 1;
    static public final byte ROW_MINUS = 2;
    static public final byte COL_MINUS = 3;
    static public final byte PREV_POINT = 4;
    static public final byte NEXT_POINT = 5;
    static public final byte GRID_ROW_PLUS = 6;
    static public final byte GRID_COL_PLUS = 7;
    static public final byte GRID_ROW_MINUS = 8;
    static public final byte GRID_COL_MINUS = 9;
    static public final byte NO_LINK = 10;
    static public final byte MULTI_INSIDE_LINKS = 11;
    
     /** StsGridSectionPoint and StsGridLink direction index strings */
     static public final String ROW_PLUS_STRING = "Row Plus";
     static public final String COL_PLUS_STRING = "Col Plus";
     static public final String ROW_MINUS_STRING = "Row Minus";
     static public final String COL_MINUS_STRING = "Row Plus";
     static public final String PREV_POINT_STRING = "Prev Point";
     static public final String NEXT_POINT_STRING = "Next Point";
     static public final String GRID_ROW_PLUS_STRING = "Grid Row Plus";
     static public final String GRID_COL_PLUS_STRING = "Grid Col Plus";
     static public final String GRID_ROW_MINUS_STRING = "Grid Row Minus";
     static public final String GRID_COL_MINUS_STRING = "Grid Col Minus";
     static public final String NO_LINK_STRING = "No link";
     static public final String MULTI_LINK_STRING = "Multi links";
     static public final String[] gridLinkNames = new String[] { ROW_PLUS_STRING, COL_PLUS_STRING, ROW_MINUS_STRING, COL_MINUS_STRING,
                                  PREV_POINT_STRING, NEXT_POINT_STRING, GRID_ROW_PLUS_STRING, GRID_COL_PLUS_STRING, GRID_ROW_MINUS_STRING,
                                  GRID_COL_MINUS_STRING, NO_LINK_STRING, MULTI_LINK_STRING };


     public static final String getGridLinkName(int linkType)
     {
         if (linkType < 0 || linkType > MULTI_INSIDE_LINKS)
             return "Unidentified: " + linkType;
         else
             return gridLinkNames[linkType];
     }
    /** WellZone types used by StsWell and StsWellZone */
    static public final int STRAT = 0;
    static public final int LITH = 1;

    /** relationship of edge to section used in StsSection */
    static public final byte MATCH = 0;
    static public final byte MATCH_NOT = 1;
    static public final byte MATCH_REVERSED = 2;
    static public final byte MATCH_UNKNOWN = 3;

    static public final String EOL = System.getProperty("line.separator");

    static public final byte NORTH = 0;
    static public final byte EAST = 1;
    static public final byte SOUTH = 2;
    static public final byte WEST = 3;

    static public final byte TD_DEPTH = 0;
    static public final byte TD_TIME = 1;
    static public final byte TD_TIME_DEPTH = 2;
    static public final byte TD_APPROX_DEPTH = 3;
    static public final byte TD_APPROX_DEPTH_AND_DEPTH = 4;
    static public final byte TD_NONE = 5;

    static public final String TD_DEPTH_STRING = "Depth";
    static public final String TD_TIME_STRING = "Time";
    static public final String TD_TIME_DEPTH_STRING = "Time & Depth";
     static public final String TD_APPROX_DEPTH_STRING = "Seismic Depth";
     static public final String TD_APPROX_DEPTH_AND_DEPTH_STRING = "Seismic depth & Depth";
    static public final String TD_NONE_STRING = "None";
    static public final String[] TD_ALL_STRINGS = new String[]
    {
        TD_DEPTH_STRING, TD_TIME_STRING, TD_TIME_DEPTH_STRING, TD_APPROX_DEPTH_STRING, TD_APPROX_DEPTH_AND_DEPTH_STRING, TD_NONE_STRING
    };
    static public final String[] TD_STRINGS = new String[]
    {
        TD_DEPTH_STRING, TD_TIME_STRING
    };

    static public final String[] TD_BOTH_STRINGS = new String[] { TD_DEPTH_STRING, TD_TIME_STRING, TD_TIME_DEPTH_STRING };
    static public final String[] TD_SELECT_STRINGS = new String[]{TD_DEPTH_STRING, TD_TIME_STRING, TD_APPROX_DEPTH_STRING};

    static public final byte TD_CONV_NONE = 00;
    static public final byte TD_CONV_DEPTH_TIME = 01;
    static public final byte TD_CONV_TIME_DEPTH = 10;
    static public final byte TD_CONV_APPROX_DEPTH_TO_DEPTH = 30;
    static public final byte TD_CONV_DEPTH_TO_APPROX_DEPTH = 03;

    static public final String TD_CONV_NONE_STRING = "None";
    static public final String TD_CONV_DEPTH_TIME_STRING = "Depth to Time";
    static public final String TD_CONV_TIME_DEPTH_STRING = "Time to Depth";
    static public final String TD_CONV_APPROX_DEPTH_TO_DEPTH_STRING = "Seismic Depth to Depth";
    static public final String TD_CONV_DEPTH_TO_APPROX_DEPTH_STRING = "Depth to Seismic Depth";

    static public byte computeDomainConversionFlag(byte fromDomain, byte toDomain)
    {
        if(fromDomain == toDomain) return (byte)0;
        else return (byte)(fromDomain*10 + toDomain);
    }

    static public final byte SAMPLE_TYPE_NONE = 0;
    static public final byte SAMPLE_TYPE_VEL_AVG = 1;
    static public final byte SAMPLE_TYPE_VEL_INTERVAL = 2;
    static public final byte SAMPLE_TYPE_VEL_RMS = 3;
    static public final byte SAMPLE_TYPE_VEL_INSTANT = 4;

    static public final String V_NONE_STRING = "Amplitude";
    static public final String V_AVG_STRING = "Average Vel";
    static public final String V_INTERVAL_STRING = "Interval Vel";
    static public final String V_RMS_STRING = "RMS Vel";
    static public final String V_INSTANT_STRING = "Instant Vel";

    static public final String[] SAMPLE_TYPE_STRINGS = new String[]
    {
        V_NONE_STRING, V_AVG_STRING, V_INTERVAL_STRING, V_RMS_STRING, V_INSTANT_STRING
    };

    static public final boolean volumeTypeIsVelocity(byte volumeType) { return volumeType > SAMPLE_TYPE_NONE; }

    static public String getSampleTypeString(byte sampleType) { return SAMPLE_TYPE_STRINGS[sampleType]; }

    static public final String[] VEL_STRINGS = new String[]
    {
        V_AVG_STRING, V_INTERVAL_STRING, V_RMS_STRING, V_INSTANT_STRING
    };

    static public final byte TIME_SECOND = 0;
    static public final String T_SECOND_STRING = "Seconds";
    static public final float TIME_SECOND_SCALE = 1.0f;
    static public final byte TIME_MSECOND = 1;
    static public final String T_MSECOND_STRING = "MilliSeconds";
    static public final float TIME_MSECOND_SCALE = 1000.0f;
    static public final byte TIME_USECOND = 2;
    static public final String T_USECOND_STRING = "MicroSeconds";
    static public final float TIME_USECOND_SCALE = 1000000.0f;
    static public final byte TIME_NONE = 3;
    static public final String T_NONE_STRING = "None";
    static public final float TIME_NONE_SCALE = nullValue;
    static public final String[] TIME_STRINGS = new String[]
    {
        T_SECOND_STRING, T_MSECOND_STRING, T_USECOND_STRING, T_NONE_STRING
    };
    static public final float[] TIME_SCALES = new float[]
    {
        TIME_SECOND_SCALE, TIME_MSECOND_SCALE, TIME_USECOND_SCALE, TIME_NONE_SCALE
    };

    static public final byte DIST_METER = 0;
    static public final float DIST_METER_SCALE = 1.0f;
    static public final String D_METER_STRING = "Meters";
    static public final byte DIST_FEET = 1;
    static public final float DIST_FEET_SCALE = 3.28084f;
    static public final String D_FEET_STRING = "Feet";
    static public final byte DIST_INCH = 2;
    static public final float DIST_INCH_SCALE = 39.37008f;
    static public final String D_INCH_STRING = "Inches";
    static public final byte DIST_DMETER = 3;
    static public final float DIST_DMETER_SCALE = 10.0f;
    static public final String D_DMETER_STRING = "Deci-Meters";
    static public final byte DIST_NONE = 4;
    static public final float DIST_NONE_SCALE = nullValue;
    static public final String D_NONE_STRING = "None";
    static public final String[] DIST_STRINGS = new String[]
    {
        D_METER_STRING, D_FEET_STRING, D_INCH_STRING, D_DMETER_STRING, D_NONE_STRING
    };
    static public final float[] DIST_SCALES = new float[]
    {
        DIST_METER_SCALE, DIST_FEET_SCALE, DIST_INCH_SCALE, DIST_DMETER_SCALE, DIST_NONE_SCALE
    };

    static public final byte TEMP_F = 0;
    static public final byte TEMP_C = 1;
    static public final String TEMP_F_STRING = "F";
    static public final String TEMP_C_STRING = "C";
    static public final String[] TEMP_STRINGS = { TEMP_F_STRING, TEMP_C_STRING };

    static public final String[] SEISMIC_VERTICAL_UNITS = new String[] { T_SECOND_STRING, T_MSECOND_STRING, D_FEET_STRING, D_METER_STRING };

    static public final String VEL_FT_PER_MSEC = "ft/msec";
    static public final String VEL_FT_PER_SEC = "ft/sec";
    static public final String VEL_M_PER_MSEC = "m/msec";
    static public final String VEL_M_PER_SEC = "m/sec";
    static public final String VEL_UNITS_NONE = "none";
    static public final String[] VEL_UNITS = new String[]
    {
        VEL_FT_PER_MSEC, VEL_FT_PER_SEC, VEL_M_PER_MSEC, VEL_M_PER_SEC, VEL_UNITS_NONE
    };

    static public final String ONE_WAY_VELOCITY = "one-way";
    static public final String TWO_WAY_VELOCITY = "two-way";
    static public final String[] ONE_OR_TWO_WAY = new String[] { ONE_WAY_VELOCITY, TWO_WAY_VELOCITY };

    static public byte getZDomainFromString(String zDomainString)
    {
        for (int n = 0; n < TD_ALL_STRINGS.length; n++)
        {
            if (zDomainString.equals(TD_ALL_STRINGS[n]))
            {
                return (byte) n;
            }
        }
        return TD_NONE;
    }

    static public String[] getSupportedDomainStrings(byte zDomain)
    {
        switch(zDomain)
        {
            case TD_NONE:
                return new String[] { TD_NONE_STRING };
            case TD_DEPTH:
                return new String[] { TD_DEPTH_STRING };
            case TD_TIME:
                return new String[] { TD_TIME_STRING };
            case TD_APPROX_DEPTH:
                return new String[] { TD_APPROX_DEPTH_STRING };
            case TD_TIME_DEPTH:
                return new String[] { TD_DEPTH_STRING, TD_TIME_STRING };
            case TD_APPROX_DEPTH_AND_DEPTH:
                return new String[] { TD_APPROX_DEPTH_STRING, TD_DEPTH_STRING };
            default:
                return new String[] { TD_NONE_STRING };
        }
    }

    static public byte getDistanceUnitsFromString(String dUnitString)
    {
        for (int n = 0; n < DIST_STRINGS.length; n++)
        {
            if (dUnitString.equals(DIST_STRINGS[n]))
            {
                return (byte) n;
            }
        }
        return DIST_NONE;
    }
    static public String getDistanceUnitString(byte type)
    {
        return DIST_STRINGS[type];
    }

    static public byte getTimeUnitsFromString(String tUnitString)
    {
        for (int n = 0; n < TIME_STRINGS.length; n++)
        {
            if (tUnitString.equals(TIME_STRINGS[n]))
            {
                return (byte) n;
            }
        }
        return TIME_NONE;
    }
    static public String getTimeUnitString(byte type)
    {
        return TIME_STRINGS[type];
    }

    static public String getTDString(byte type)
    {
        return TD_BOTH_STRINGS[type];
    }

    static public String estimateVelocityUnits(float dataItem)
    {
        if(dataItem < 1.0f)
         {
             return VEL_M_PER_MSEC;
         }
         else if(dataItem < 2.0f)
             return VEL_M_PER_MSEC;
         else if(dataItem < 10.0f)
             return VEL_FT_PER_MSEC;
        else if(dataItem < 2000.0f)
             return VEL_M_PER_SEC;
         else
             return VEL_FT_PER_SEC;
    }

    static public byte getVelocityTypeFromString(String vTypeString)
    {
        if(vTypeString == null) return SAMPLE_TYPE_NONE;
        for (int n = 0; n < VEL_STRINGS.length; n++)
        {
            if (VEL_STRINGS[n].equals(vTypeString))
            {
                return (byte) n;
            }
        }
        return SAMPLE_TYPE_NONE;
    }
    static public String getVelocityTypeString(byte type)
    {
        return VEL_STRINGS[type];
    }

    static public byte getStringMatchByteIndex(String[] strings, String string)
    {
        for(byte n = 0; n < strings.length; n++)
            if(strings[n] == string) return n;
        return -1;
    }

    static public int getStringMatchIndex(String[] strings, String string)
    {
        for(int n = 0; n < strings.length; n++)
            if(strings[n] == string) return n;
        return -1;
    }

    /** Types for seismic volumes */
	static final public byte VELOCITY = 1;
	static final public byte HILBERT_AMP = 2;
	static final public byte HILBERT_FREQ = 3;
	static final public byte HILBERT_PHASE = 4;
	static final public byte HILBERT_TRANS = 5;
	static final public byte ANALOGUE = 6;
	static final public byte STACKED = 7;
	static final public byte SEMBLANCE = 8;

	public static final int XDIR = 0;
	/** static variable YDIR = 1 - Y oriented cursor section */
	public static final int YDIR = 1;
	/** static variable ZDIR = 2 - Z oriented cursor section */
	public static final int ZDIR = 2;
	/** Default orientation labels */
    public static final String[] coorLabels = new String[] { "X", "Y", "Z" };

	public static float convertVelocity(float velIn, String distUnitsInput, String timeUnitsInput, String distUnitsOutput, String timeUnitsOutput)
	{
		byte distUnitsInputIndex = getDistanceUnitsFromString(distUnitsInput);
		byte distUnitsOutputIndex = getDistanceUnitsFromString(distUnitsOutput);
		byte timeUnitsInputIndex = getTimeUnitsFromString(timeUnitsInput);
		byte timeUnitsOutputIndex = getTimeUnitsFromString(timeUnitsOutput);
		float velOut = convertVelocity(velIn, distUnitsInputIndex, distUnitsOutputIndex, timeUnitsInputIndex, timeUnitsOutputIndex);
		System.out.println("input velocity " + velIn + " " + distUnitsInput + "/" + timeUnitsInput + ". output velocity " + velOut + " " + distUnitsOutput + "/" + timeUnitsOutput);
		return velOut;
	}

	public static boolean velocityUnitsChanged(byte distUnitsInputIndex, byte timeUnitsInputIndex, byte distUnitsOutputIndex, byte timeUnitsOutputIndex)
	{
		return distUnitsInputIndex != distUnitsOutputIndex || timeUnitsInputIndex != timeUnitsOutputIndex;
	}

	public static float convertVelocity(float velIn, byte distUnitsInputIndex, byte timeUnitsInputIndex, byte distUnitsOutputIndex, byte timeUnitsOutputIndex)
	{
		float distUnitsInputScale = DIST_SCALES[distUnitsInputIndex];
		float distUnitsOutputScale = DIST_SCALES[distUnitsOutputIndex];
		float timeUnitsInputScale = TIME_SCALES[timeUnitsInputIndex];
		float timeUnitsOutputScale = TIME_SCALES[timeUnitsOutputIndex];
		return velIn*distUnitsOutputScale*timeUnitsInputScale/(distUnitsInputScale*timeUnitsOutputScale);
	}

	public static boolean timeUnitsChanged(byte timeUnitsInputIndex, byte timeUnitsOutputIndex)
	{
		return timeUnitsInputIndex != timeUnitsOutputIndex;
	}

	public static float convertTime(float timeIn, byte timeUnitsInputIndex, byte timeUnitsOutputIndex)
	{
		float timeUnitsInputScale = TIME_SCALES[timeUnitsInputIndex];
		float timeUnitsOutputScale = TIME_SCALES[timeUnitsOutputIndex];
		return timeIn*timeUnitsOutputScale/timeUnitsInputScale;
	}

	public static boolean distanceUnitsChanged(byte depthUnitsInputIndex, byte depthUnitsOutputIndex)
	{
		return depthUnitsInputIndex != depthUnitsOutputIndex;
	}

	public static float convertDistance(float depthIn, byte depthUnitsInputIndex, byte depthUnitsOutputIndex)
	{
		float depthUnitsInputScale = DIST_SCALES[depthUnitsInputIndex];
		float depthUnitsOutputScale = DIST_SCALES[depthUnitsOutputIndex];
		return depthIn*depthUnitsOutputScale/depthUnitsInputScale;
	}
	static public String getVelocityString(byte distUnitsIndex, byte timeUnitsIndex)
	{
		return new String(DIST_STRINGS[distUnitsIndex] + "/" +  TIME_STRINGS[timeUnitsIndex]);
	}

	static public String getTimeString(byte timeUnitsIndex)
	{
		return TIME_STRINGS[timeUnitsIndex];
	}

	static public String getDepthString(byte depthUnitsIndex)
	{
		return DIST_STRINGS[depthUnitsIndex];
	}

     static public final String TVD_STRING = "TVD";
     static public final String SUBSEA_STRING = "SubSea";
     static public final String[] TVD_SUBSEA_STRINGS = new String[]{TVD_STRING, SUBSEA_STRING};
	static public void main(String[] args)
	{
		float vi = 10.0f; // m/sec
		convertVelocity(vi, StsParameters.D_METER_STRING, T_SECOND_STRING, D_FEET_STRING, T_SECOND_STRING);
		convertVelocity(vi, StsParameters.D_METER_STRING, T_SECOND_STRING, D_FEET_STRING, T_MSECOND_STRING);
	}
}
