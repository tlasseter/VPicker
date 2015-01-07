package com.Sts.UI.Progress;

import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import java.awt.*;

/**
 * <p>Title: </p>
 * <p/>
 * <p>Description: </p>
 * <p/>
 * <p>Copyright: Copyright (c) 2006</p>
 * <p/>
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public class StsProgressPanel extends StsJPanel
{
    boolean addCancelButton = false;
    int nRows = 0;
    int nCols = 0;

    public StsProgressTextPanel textPanel;
    public StsIntervalProgressBar progressBar;

    private StsProgressPanel(boolean addCancelButton)
    {
        super(true);
        this.addCancelButton = addCancelButton;
        constructPanel();
    }

    private StsProgressPanel(boolean addCancelButton, int nRows, int nCols)
    {
        super();
        this.addCancelButton = addCancelButton;
        this.nRows = nRows;
        this.nCols = nCols;
        constructPanel();
    }

    static public StsProgressPanel constructor()
    {
        return new StsProgressPanel(false);
    }

    static public StsProgressPanel constructor(boolean addCancel, int nRows, int nCols)
    {
        return new StsProgressPanel(addCancel, nRows, nCols);
    }

    static public StsProgressPanel constructor(int nRows, int nCols)
    {
        return new StsProgressPanel(false, nRows, nCols);
    }

    static public StsProgressPanel constructorWithCancelButton(int nRows, int nCols)
    {
        return new StsProgressPanel(true, nRows, nCols);
    }

    static public StsProgressPanel constructorWithCancelButton()
    {
        return new StsProgressPanel(true);
    }

    public void addCancelButton() { addCancelButton = true; }

    public void initialize(int maxValue)
    {
        setMaximum(maxValue);
        setValue(0);
    }

    public void setMaximum(int maxValue)
    {
        progressBar.setMaximum(maxValue);
    }

    public void setValue(int value)
    {
        progressBar.setValue(value);
    }

    public void finished()
    {
        progressBar.finished();
    }

    public void constructPanel()
    {
        if(nRows > 0)
        {
            textPanel = StsProgressTextPanel.constructor(nRows, nCols);
            int width = Math.max(50, 40+nCols*8);
            int height = 60 + nRows*20;
            textPanel.setPreferredSize(width, height);
            gbc.fill = GridBagConstraints.BOTH;
            gbc.weighty = 1.0;
            add(textPanel);
        }

        if(addCancelButton)
            progressBar = StsIntervalProgressBar.constructorWithCancel();
        else
            progressBar = StsIntervalProgressBar.constructor();
        
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.weighty = 0.0;
        add(progressBar);
    }

    public void setValue(double value)
    {
        progressBar.setValue((int)(value * 100));
    }

    public void appendLine(String line)
    {
        if(textPanel == null) return;
        textPanel.appendLine(line);
    }

    public void appendErrorLine(String line)
    {
        if(textPanel == null) return;
        textPanel.appendErrorLine(line);
    }

    public void setDescriptionAndLevel(String progressDescription, int level)
    {
        progressBar.setDescriptionAndLevel(progressDescription, level);
    }

    public void setDescription(String progressDescription)
    {
        progressBar.setDescription(progressDescription);
    }

    public void incrementProgress()
    {
        progressBar.incrementProgress();
    }

    public void setLevel(int level)
    {
        progressBar.setLevel(level);
    }

    public void resetProgressBar()
    {
        progressBar.resetProgressBar();
    }

    public void initializeIntervals(int nIntervals, int intervalSize)
    {
        progressBar.initializeIntervals(nIntervals, intervalSize);
    }

    public void initializeIntervals(int nIntervals)
    {
        progressBar.initializeIntervals(nIntervals);
    }

    public void incrementInterval()
    {
        progressBar.incrementInterval();
    }

    public void setInterval(int nInterval)
    {
        progressBar.setInterval(nInterval);
    }

    public void setIntervalCount(long intervalTotalCount)
    {
        progressBar.setIntervalCount(intervalTotalCount);
    }

    public void setSubInterval(float startFraction, float endFraction, long subIntervalTotalCount)
    {
        progressBar.setSubInterval(startFraction, endFraction, subIntervalTotalCount);
    }

    public void setSubInterval(float startFraction, float endFraction)
    {
        progressBar.setSubInterval(startFraction, endFraction);
    }

    public void setSubIntervalTotalCount(long subIntervalTotalCount)
    {
        progressBar.setSubIntervalTotalCount(subIntervalTotalCount);
    }

    public void incrementCount()
    {
        progressBar.incrementCount();
    }

    public void setCount(long count)
    {
        progressBar.setCount(count);
    }

    public boolean isCanceled() { return progressBar.canceled; }

    public void cancel()
    {
        progressBar.cancel();
    }

    public void clearProcess()
    {
        progressBar.clearProcess();
    }

    static public void main(String[] args)
	{
        Runnable winThread = new Runnable()
        {
            public void run()
            {
                StsProgressPanel panel = constructorWithCancelButton(5, 20);
                TestProgressStatusPanel test = new TestProgressStatusPanel(panel);
                StsToolkit.createDialog(panel, false);
                test.run();
            }
        };
        StsToolkit.runWaitOnEventThread(winThread); 
    }
}

class TestProgressStatusPanel implements Runnable
{
    StsProgressPanel progressPanel;
    int max = 100;

    TestProgressStatusPanel(StsProgressPanel panel)
    {
        progressPanel = panel;
    }

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

    public void runSomething()
    {
        progressPanel.setMaximum(max);
        for(int n = 0; n < max; n++)
        {
            if(progressPanel.progressBar.canceled) return;
            progressPanel.setValue(n);
            progressPanel.appendLine("doing " + n);
            StsToolkit.sleep(100);
            if(n == 50) progressPanel.setDescriptionAndLevel("Warning", StsProgressBar.WARNING);
        }
    }
}
