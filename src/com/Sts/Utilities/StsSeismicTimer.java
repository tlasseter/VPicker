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
public class StsSeismicTimer extends StsSumTimer
{
    static String TIMER_MAP_BLOCK = "map block";
    static String TIMER_OUTPUT_BLOCK = "get output block";
	static String TIMER_INPUT_BLOCK = "get input block";
	static String TIMER_READ_INPUT_BLOCK = "read in input block";
	static String TIMER_NORMALIZE = "normalize/filter/AGC";
	static String TIMER_CALC_VEL = "calculate vel";
	static String TIMER_STACK = "stack";
    public static StsTimer getMapBlockTimer;
    public static StsTimer getOutputBlockTimer;
	public static StsTimer getInputBlockTimer;
	public static StsTimer readInputBlockTimer;
	public static StsTimer normalizeTimer;
	public static StsTimer calcVelTimer;
	public static StsTimer stackTimer;
	public static StsTimer overallTimer;

	static
	{
        getMapBlockTimer = StsSumTimer.addTimer(TIMER_MAP_BLOCK);
        getOutputBlockTimer = StsSumTimer.addTimer(TIMER_OUTPUT_BLOCK);
		getInputBlockTimer = StsSumTimer.addTimer(TIMER_INPUT_BLOCK);
		readInputBlockTimer = StsSumTimer.addTimer(TIMER_READ_INPUT_BLOCK);
		normalizeTimer = StsSumTimer.addTimer(TIMER_NORMALIZE);
		calcVelTimer = StsSumTimer.addTimer(TIMER_CALC_VEL);
		stackTimer = StsSumTimer.addTimer(TIMER_STACK);
		overallTimer = new StsTimer("overall");
	}

	static public void clear()
	{
		StsSumTimer.clear();
		overallTimer.clear();
	}

	public static void main(String[] args)
	{
		try
		{
			com.Sts.MVC.Main.setVersion();
			overallTimer.start();
			for(int i = 0; i < 10; i++)
			{
				getOutputBlockTimer.start();
				Thread.currentThread().sleep(100);
				getOutputBlockTimer.stopAccumulateIncrementCount();
				getInputBlockTimer.start();
				Thread.currentThread().sleep(200);
				getInputBlockTimer.stopAccumulateIncrementCount();
			}
			printTimers("test");
			overallTimer.stopPrint("overall time");
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}
}
