package com.Sts.UI.Progress;

import com.Sts.Utilities.*;

/** StsProgressBar can be run in two modes:  a simple progress from minValue to maxValue or
 *  a series of intervals with optional subIntervals.  The intervals are of equal length,
 *  but the subIntervals are defined as fractions of the interval.  For each subInterval,
 *  you must specify the subIntervalCountTotal which is the number of steps the application
 *  will run in completing the subInterval.  The progressBar is given the currentCount and
 *  computes where the progressBar should be based on the currentInterval, the subInterval,
 *  and the countTotal for the subInterval.
 */
public class StsIntervalProgressBar extends StsProgressBar
{
    protected int nIntervals = 1;
    protected float intervalMin = 0;
    protected float intervalMax = 100;
    protected int intervalSize = 100;
    protected float countPerIncrement = 1;
    protected int nInterval = 0;
    protected float intervalCount = 0;

    private StsIntervalProgressBar(boolean addInsets, boolean addCancelButton)
    {
        super(addInsets, addCancelButton);
    }

    static public StsIntervalProgressBar constructor(boolean addInsets, boolean addCancelButton)
    {

        try
        {
            return new StsIntervalProgressBar(addInsets, addCancelButton);
        }
        catch(Exception e)
        {
            return null;
        }
    }

    static public StsIntervalProgressBar constructor()
     {
         return constructor(false, false);
     }

    static public StsIntervalProgressBar constructorWithInsets()
    {
        return constructor(true, false);
    }

    static public StsIntervalProgressBar constructorWithCancel()
    {
        return constructor(false, true);
    }

    static public StsIntervalProgressBar constructorWithInsetsAndCancel()
    {
        return constructor(true, true);
    }

    public void initializeIntervals(int nIntervals, int intervalSize)
    {
        this.nIntervals = nIntervals;
        this.intervalSize = intervalSize;
        nInterval = 0;
        initialize(nIntervals*intervalSize);
    }

    public void initializeIntervals(int nIntervals)
    {
        initializeIntervals(nIntervals, 100);
    }

    public void incrementInterval()
    {
        setInterval(nInterval + 1);
    }

    public void setInterval(int nInterval, long intervalTotalCount)
    {
        setInterval(nInterval);
        setIntervalCount(intervalTotalCount);
    }

    public void setInterval(int n)
    {
        nInterval = Math.min(n, nIntervals-1);
        intervalMin = nInterval*intervalSize;
        intervalMax = intervalMin + intervalSize;
        setValue(intervalSize*nInterval);
        intervalCount = intervalMin;
    }

    public void setIntervalCount(long intervalTotalCount)
    {
        countPerIncrement = (intervalMax - intervalMin)/intervalTotalCount;
    }

    public void setSubInterval(float startFraction, float endFraction, long subIntervalTotalCount)
    {
        setSubInterval(startFraction, endFraction);
        setSubIntervalTotalCount(subIntervalTotalCount);
    }

    public void setSubInterval(float startFraction, float endFraction)
    {
        intervalMin = (nInterval + startFraction)*intervalSize;
        intervalMax = (nInterval + endFraction)*intervalSize;
        intervalCount = intervalMin;
    }

    public void setSubIntervalTotalCount(long subIntervalTotalCount)
    {
        countPerIncrement = (intervalMax - intervalMin)/subIntervalTotalCount;
    }

    public void incrementCount()
    {
        intervalCount += countPerIncrement;
        int newValue = Math.round(intervalCount);
        if(newValue > progressBarValue) setValue(newValue);
    }

    public void setCount(long count)
    {
        intervalCount = intervalMin + count/countPerIncrement;
        int newValue = (int)Math.min(intervalMax, Math.round(intervalCount));
        if(newValue > progressBarValue) setValue(newValue);
    }

    public void resetProgressBar()
    {
        super.resetProgressBar();
        intervalCount = intervalMin;
    }

    static public void main(String[] args)
	{
        final StsProgressPanel progressBar = StsProgressPanel.constructorWithCancelButton(5, 50);
        final IntervalProgressBarTest test =  new IntervalProgressBarTest(progressBar);
        Runnable winThread = new Runnable()
        {
            public void run()
            {
                StsToolkit.createDialog(progressBar, false);
            }
        };
        StsToolkit.runWaitOnEventThread(winThread);
        Thread testThread = new Thread(test);
        testThread.start();
    }
}

class IntervalProgressBarTest implements Runnable
{
	StsProgressPanel progressPanel;
    int max = 100;
    boolean interrupt = false;

    IntervalProgressBarTest(StsProgressPanel progressBar)
	{
        this.progressPanel = progressBar;
    }

    public void setMax(int max) { this.max = max; }
    public int getMax() { return max; }

    public void run()
	{
		Runnable progressTest = new Runnable()
		{
			public void run()
			{
				runSomething();
			}
		};

		Thread progressTestThread = new Thread(progressTest);
		progressTestThread.start();
	}

    public void cancel()
    {
        interrupt = true;
    }

    void runSomething()
	{
        int nIntervals, nInterval, intervalCount;
        try
		{
            // this is a test of two intervals with no subintervals
             progressPanel.appendLine("ProgressBarInterval Test: 2 intervals with no subIntervals. ");
             nIntervals = 2;
             nInterval = 0;
             progressPanel.initializeIntervals(nIntervals);
             intervalCount = 50;
             progressPanel.setIntervalCount(intervalCount);
             progressPanel.appendLine("Doing interval " + progressPanel.progressBar.nInterval + " with count " + intervalCount);
             for(int i = 0; i < intervalCount; i++)
             {
                if(progressPanel.isCanceled()) return;
                progressPanel.incrementCount();
                StsToolkit.sleep(20);
             }
             StsToolkit.sleep(200);
             intervalCount = 10;
            nInterval++;
            progressPanel.incrementInterval();
            progressPanel.setIntervalCount(intervalCount);
            progressPanel.appendLine("Doing interval " + nInterval + " with count " + intervalCount);
            for(int i = 0; i < intervalCount; i++)
            {
               if(progressPanel.isCanceled()) return;
               progressPanel.incrementCount();
               StsToolkit.sleep(20);
            }
            StsToolkit.sleep(2000);
        // ----------------------------------
            progressPanel.appendLine("ProgressBarInterval Test: 2 intervals with subIntervals. ");
             nIntervals = 2;
             nInterval = 0;
             progressPanel.initializeIntervals(nIntervals, 100);
             intervalCount = 50;
             progressPanel.setSubInterval(0.0f, 0.4f, intervalCount);
             progressPanel.appendLine("Doing interval " + progressPanel.progressBar.nInterval + " with count " + intervalCount);
             for(int i = 0; i < intervalCount; i++)
             {
                if(progressPanel.isCanceled()) return;
                progressPanel.incrementCount();
                StsToolkit.sleep(20);
             }
             StsToolkit.sleep(200);
             intervalCount = 10;
             progressPanel.setSubInterval(0.4f, 1.0f, intervalCount);
             progressPanel.appendLine("Doing interval " + progressPanel.progressBar.nInterval + " with count " + intervalCount);
             for(int i = 0; i < intervalCount; i++)
             {
                if(progressPanel.isCanceled()) return;
                progressPanel.incrementCount();
                StsToolkit.sleep(20);
             }

            progressPanel.incrementInterval();
            intervalCount = 150;
            progressPanel.setSubInterval(0.0f, 0.2f, intervalCount);
            progressPanel.appendLine("Doing interval " + progressPanel.progressBar.nInterval + " with count " + intervalCount);
            for(int i = 0; i < intervalCount; i++)
            {
               if(progressPanel.isCanceled()) return;
               progressPanel.incrementCount();
               StsToolkit.sleep(20);
            }
            StsToolkit.sleep(200);
            intervalCount = 75;
            progressPanel.setSubInterval(0.2f, 0.4f, intervalCount);
            progressPanel.appendLine("Doing interval " + progressPanel.progressBar.nInterval + " with count " + intervalCount);
            for(int i = 0; i < intervalCount; i++)
            {
               if(progressPanel.isCanceled()) return;
               progressPanel.incrementCount();
               StsToolkit.sleep(20);
            }
            StsToolkit.sleep(200);
            intervalCount = 100;
            progressPanel.setSubInterval(0.4f, 1.0f, intervalCount);
            progressPanel.appendLine("Doing interval " + progressPanel.progressBar.nInterval + " with count " + intervalCount);
            for(int i = 0; i < intervalCount; i++)
            {
               if(progressPanel.isCanceled()) return;
               progressPanel.incrementCount();
               StsToolkit.sleep(20);
            }
        }
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}
}
