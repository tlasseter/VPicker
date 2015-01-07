package com.Sts.Types.PreStack;

import com.Sts.Utilities.*;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: May 23, 2008
 * Time: 6:23:10 AM
 * To change this template use File | Settings | File Templates.
 */
public interface StsDataLineSetFace
{
    public int getNSamples(int row);
    public Iterator getGatherIterator();
    public String[] getAttributeNames();
    public int getNTraces();
    public String getStemname();
    public String getDescription();
    public byte getTimeUnits(); // returns seconds, milliseconds, or microseconds
    public byte getDepthUnits(); // returns feet or meters
    public byte getHorizontalUnits(); // returns feet or meters
    public float getAngle(); //angle counterclockwise from due East of first inline
    public float getDataMax(); //maximum sample value
    public float getDataMin(); //minimum sample value
    public boolean getIsNMOed(); //has NMO been applied?
    public boolean getIsXLineCCW(); //is the first crossline counterclockwise to the first inline?
    public int getNAttributes(); //number of attributes
    public int getNSlices(); //number of samples per trace
    public float getXInc(); //bin increment in crossline direction
    public float getYInc(); //bin increment in inline direction
    public byte getZDomain(); //are z units time or depth?
    public float getZInc(); //sample increment (milliseconds)
    public float getZMax(); //time of last sample (milliseconds)
    public float getZMin(); //time of first sample (milliseconds)
    public void initializeDataSet(); //do any prepwork necessary to get all of this stuff
    public float getColNumInc();  //increment in crossline numbers
    public float getRowNumInc();  //increment in inline numbers
    public String getSegyFilename();  //name of SEG-Y input file (or other useful information about source of data)

    static public final byte TIME_SECOND = StsParameters.TIME_SECOND;
    static public final byte TIME_MSECOND = StsParameters.TIME_MSECOND;
    static public final byte TIME_USECOND = StsParameters.TIME_USECOND;

    static public final byte DIST_METER = StsParameters.DIST_METER;
    static public final byte DIST_FEET = StsParameters.DIST_FEET;
    
    static public final byte TIME = StsParameters.TD_TIME;
    static public final byte DEPTH = StsParameters.TD_DEPTH;
}
