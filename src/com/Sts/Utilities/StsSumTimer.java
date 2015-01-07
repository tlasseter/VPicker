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
public class StsSumTimer
{
	static StsTimer[] timers;
	static int nTimers = 0;
	static int incSize = 10;

	static
	{
	    timers = new StsTimer[incSize];
	}

	static public StsTimer addTimer(String name)
	{
		StsTimer timer = new StsTimer(name);
        addTimer(timer);
        return timer;
    }

    static public void addTimer(StsTimer timer)
    {
		if(nTimers+1 > timers.length)
		{
			StsTimer[] newTimers = new StsTimer[timers.length + 10];
			System.arraycopy(timers, 0, newTimers, 0, timers.length);
			timers = newTimers;
		}
		timers[nTimers++] = timer;
	}

	static public void startTimer(String name)
	{
		StsTimer timer = getTimer(name);
		if(timer == null)
		{
			StsException.systemError("StsSumTimer.startTimer() failed.  Couldn't find timer " + name);
			return;
		}
		timer.start();
	}

	static public void stopTimer(String name)
	{
		StsTimer timer = getTimer(name);
		if(timer == null)
		{
			StsException.systemError("StsSumTimer.startTimer() failed.  Couldn't find timer " + name);
			return;
		}
		timer.stopAccumulateIncrementCount();
	}

	static public StsTimer getTimer(String name)
	{
		for(int n = 0; n < nTimers; n++)
		{
			StsTimer timer = timers[n];
			if(timer.name == name) return timer;
		}
		StsException.systemError("StsSumTimer.getTimer() failed. Couldn't find name " + name);
		return null;
	}

	static public void printTimers(String name)
	{
		double sum = 0.0;
		System.out.println("Summary for timers " + name);
		for(int n = 0; n < nTimers; n++)
		{
			timers[n].printElapsedTime();
			sum += timers[n].getElapsedTime();
		}
        String sumTimeString = StsTimer.getTimeString(sum);
	    System.out.println("Check sum timers: " + sumTimeString);
	}

	static public void clear()
	{
		for(int n = 0; n < nTimers; n++)
			timers[n].clear();
	}

    static public void printIntervalElapsedTimes(boolean printIntervals)
    {
        for(int n = 0; n < nTimers; n++)
            timers[n].printIntervalElapsedTimes = printIntervals;
    }
/*
	public static void main(String[] args)
	{
		com.Sts.MVC.Main.setVersion();
		StsSumTimer.classInitialize("Sum test");
		String timerOneString = "timer one";
		String timerTwoString = "timer two";
		StsSumTimer.addTimer(timerOneString);
		StsSumTimer.addTimer(timerTwoString);
		try
		{
			for(int i = 0; i < 10; i++)
			{
				StsSumTimer.startTimer(timerOneString);
				Thread.currentThread().sleep(100);
				StsSumTimer.stopTimer(timerOneString);
				StsSumTimer.startTimer(timerTwoString);
				Thread.currentThread().sleep(200);
				StsSumTimer.stopTimer(timerTwoString);
		    }
			StsSumTimer.printTimers();
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}
*/
}
