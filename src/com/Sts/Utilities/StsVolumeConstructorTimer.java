package com.Sts.Utilities;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */
public class StsVolumeConstructorTimer extends StsSumTimer
{
    static String TIMER_PROCESS_BLOCK_OUTPUT = "block input";
    static String TIMER_MAP_BLOCK = "map block";
    static String TIMER_CLEAR_INPUT_BLOCK = "clear input block";
    static String TIMER_CLEAR_OUTPUT_BLOCK = "clear output block";
    static String TIMER_FLOAT_ROW_OUTPUT = "float row output";
    static String TIMER_BYTE_ROW_OUTPUT = "row byte output";
	static String TIMER_BYTE_COL_OUTPUT = "col byte output";
	static String TIMER_SLICE_ROW_OUTPUT = "slice byte output";
    static String TIMER_CALC_VEL = "calculate vel";
	static String TIMER_STACK = "stack";

    public static StsTimer getBlockInputTimer;
    public static StsTimer getMapBlockTimer;
    public static StsTimer getClearInputBlockTimer;
    public static StsTimer getClearOutputBlockTimer;
    public static StsTimer getFloatRowTimer;
    public static StsTimer getByteRowTimer;
	public static StsTimer getByteColTimer;
	public static StsTimer getByteSliceTimer;
    public static StsTimer overallTimer;
    public static StsTimer calcVelTimer;
	public static StsTimer stackTimer;

    static
	{
        getBlockInputTimer = StsSumTimer.addTimer(TIMER_PROCESS_BLOCK_OUTPUT);
        getMapBlockTimer = StsSumTimer.addTimer(TIMER_MAP_BLOCK);
        getClearInputBlockTimer = StsSumTimer.addTimer(TIMER_CLEAR_INPUT_BLOCK);
        getClearOutputBlockTimer = StsSumTimer.addTimer(TIMER_CLEAR_OUTPUT_BLOCK);
        getFloatRowTimer = StsSumTimer.addTimer(TIMER_FLOAT_ROW_OUTPUT);
        getByteRowTimer = StsSumTimer.addTimer(TIMER_BYTE_ROW_OUTPUT);
		getByteColTimer = StsSumTimer.addTimer(TIMER_BYTE_COL_OUTPUT);
		getByteSliceTimer = StsSumTimer.addTimer(TIMER_SLICE_ROW_OUTPUT);
		calcVelTimer = StsSumTimer.addTimer(TIMER_CALC_VEL);
		stackTimer = StsSumTimer.addTimer(TIMER_STACK);
        overallTimer = new StsTimer("overall");
    }

	static public void clear()
	{
		StsSumTimer.clear();
		overallTimer.clear();
	}
}