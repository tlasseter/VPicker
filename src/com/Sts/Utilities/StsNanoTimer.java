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
public class StsNanoTimer extends StsTimer
{
	static public boolean gaveErrorMessage = false;
	StsNanoTimer()
	{
	}

	StsNanoTimer(String name)
	{
		this.name = name;
		System.out.println("Created nanoTimer " + name);
	}

	public void start()
	{

		if(gaveErrorMessage)
		{
			System.err.println("Nanotimer called; Java versions is 1.5 so uncomment lines in StsNanoTimer start and stop methods.");
			gaveErrorMessage = true;
		}

//		startTime = System.nanoTime();
	}

	public double stop()
	{
		stopTime = System.currentTimeMillis();
//		stopTime = System.nanoTime();
		elapsedTime += 1.e-6*(stopTime - startTime);
		count++;
		return elapsedTime;
    }

	public static void main(String[] args)
	{
		StsNanoTimer stsnanotimer = new StsNanoTimer();
	}
}
