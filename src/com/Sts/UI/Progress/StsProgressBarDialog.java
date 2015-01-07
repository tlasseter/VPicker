package com.Sts.UI.Progress;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version beta 1.0
 */

import javax.swing.*;
import java.awt.*;

public class StsProgressBarDialog extends JDialog
{
	public StsProgressPanel progressPanel;

    private StsProgressBarDialog(Frame frame, String title, boolean modal, boolean addCancel, int nRows, int nCols, int width, int height)
    {
        super(frame, title, modal);
        if (width != 0 && height != 0) setSize(width, height);
        progressPanel = StsProgressPanel.constructor(addCancel, nRows, nCols);
        getContentPane().add(progressPanel);
        setLocationRelativeTo(frame);
        pack();
		setVisible(true);
    }

    static public StsProgressBarDialog constructor(Frame frame, String title)
    {
        return new StsProgressBarDialog(frame, title, false, false, 0, 0, 200, 40);
    }

    static public StsProgressBarDialog constructor(Frame frame, String title, boolean modal)
    {
        return new StsProgressBarDialog(frame, title, modal, false, 0, 0, 200, 40);
    }

    static public StsProgressBarDialog constructor(Frame frame, String title, boolean modal, int nRows, int nCols, int width, int height)
    {
        return new StsProgressBarDialog(frame, title, modal, false, nRows, nCols, width, height);
    }

    static public StsProgressBarDialog constructorWithCancel(Frame frame, String title, boolean modal, int nRows, int nCols, int width, int height)
    {
        return new StsProgressBarDialog(frame, title, modal, true, nRows, nCols, width, height);
    }

	public void setLabelText(String text)
	{
        progressPanel.setDescription(text);
	}

	public void setProgressMax(int max)
	{
		progressPanel.setMaximum(max);
	}

	public void setProgress(int progress)
	{
		progressPanel.setValue(progress);
	}

    public boolean isCanceled()
    {
        return progressPanel.isCanceled();
    }

    public void finished()
	{
		progressPanel.finished();
        dispose();
	}

	public static void main(String[] args)
	{
		TestProgress testProgress = new TestProgress();
		testProgress.run();
	}
}

class TestProgress
{
	StsProgressBarDialog progressBarDialog;

	TestProgress()
	{
		progressBarDialog = StsProgressBarDialog.constructor(null, "Test progress", false, 5, 20, 500, 40);
	}

	void run()
	{
		Runnable testProgress = new Runnable()
		{
			public void run()
			{
				runSomething();
			}
		};

		Thread testProgressThread = new Thread(testProgress);
		testProgressThread.start();
	}
	void runSomething()
	{
		try
		{
			progressBarDialog.progressPanel.setMaximum(10);
			progressBarDialog.progressPanel.appendLine("Start message.");
			for(int n = 1; n <= 10; n++)
			{
				Thread.sleep(1000);
				progressBarDialog.progressPanel.setValue(n);
				progressBarDialog.progressPanel.appendLine("Step " + n);
				if(n <= 5)
					progressBarDialog.progressPanel.setDescriptionAndLevel("Step " + n, StsProgressBar.WARNING);
				else
					progressBarDialog.progressPanel.setDescriptionAndLevel("Step " + n, StsProgressBar.WARNING);
			}
			progressBarDialog.progressPanel.finished();
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}
}
