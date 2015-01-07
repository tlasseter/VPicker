package com.Sts.UI.Progress;

import com.Sts.DBTypes.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
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

/**
 * StsProgressBar can be run in two modes:  a simple progress from minValue to maxValue or
 * a series of intervals with optional subIntervals.  The intervals are of equal length,
 * but the subIntervals are defined as fractions of the interval.  For each subInterval,
 * you must specify the subIntervalCountTotal which is the number of steps the application
 * will run in completing the subInterval.  The progressBar is given the currentCount and
 * computes where the progressBar should be based on the currentInterval, the subInterval,
 * and the countTotal for the subInterval.
 */

//TODO for Runnables, use StsToolkit methods which check for thread (like StsToolkit.runLaterOnEventThread(Runnable)
public class StsProgressBar extends StsJPanel
{
    String description;
    public boolean canceled = false;

    protected int progressBarMaxValue = 100;
    protected int progressBarMinValue = 0;
    protected int progressBarValue = 0;

    JProgressBar progressBar = new JProgressBar(0, 100);
    boolean addCancelButton = false;
    StsButton cancelButton;

    public static final int INFO = 0;
    public static final int WARNING = 1;
    public static final int FATAL = 2;
    public static final int ERROR = 3;

    protected StsProgressBar(boolean addInsets, boolean addCancelButton)
    {
        super(addInsets);
        this.addCancelButton = addCancelButton;
        constructPanel();
    }

    static public StsProgressBar constructor(boolean addInsets, boolean addCancelButton)
    {
        try
        {
            return new StsProgressBar(addInsets, addCancelButton);
        }
        catch(Exception e)
        {
            return null;
        }
    }

    static public StsProgressBar constructor()
    {
        return constructor(false, false);
    }

    static public StsProgressBar constructorWithInsets()
    {
        return constructor(true, false);
    }

    static public StsProgressBar constructorWithCancel()
    {
        return constructor(false, true);
    }

    static public StsProgressBar constructorWithInsetsAndCancel()
    {
        return constructor(true, true);
    }

    protected void constructPanel()
    {
        progressBar.setBackground(Color.WHITE);
        progressBar.setStringPainted(true);
        progressBar.setFont(new java.awt.Font("Serif", 1, 12));
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.anchor = GridBagConstraints.CENTER;
        addToRow(progressBar);
        gbc.fill = GridBagConstraints.NONE;
        gbc.weightx = 0;
        if(addCancelButton)
        {
            cancelButton = new StsButton("Cancel", "Hit button to terminate process.", this, "cancel");
            addEndRow(cancelButton);
        }
    }

    public void setMaximum(int maxValue)
    {
        progressBarMaxValue = maxValue;
        StsToolkit.runLaterOnEventThread(new Runnable() { public void run() { progressBar.setMaximum(progressBarMaxValue); }});
    }

    public void setMinimum(int minValue)
    {
        progressBarMinValue = minValue;
        Runnable setMinimumRunnable = new Runnable()
        {
            public void run()
            {
                int percent = (int)(100.0f * (progressBarValue - progressBarMinValue) / (progressBarMaxValue - progressBarMinValue));
                progressBar.setMaximum(100);
                progressBar.setValue(percent);
            }
        };
        StsToolkit.runLaterOnEventThread(setMinimumRunnable);
    }

    public void setValue(int value)
    {
        progressBarValue = value;
        Runnable setValueRunnable = new Runnable()
        {
            public void run()
            {
                progressBar.setString(getString());
                progressBar.setValue(progressBarValue);
            }
        };
        StsToolkit.runLaterOnEventThread(setValueRunnable);
    }

    public void setValueImmediate(int value)
    {
        progressBar.setValue(value);
        paintImmediately();
    }

    public void setValueAndStringImmediate(int value, String string)
    {
        progressBar.setValue(value);
        progressBar.setString(string);
        paintImmediately();
    }

    public void paintImmediately()
    {
        progressBar.paintImmediately(0, 0, progressBar.getWidth(), progressBar.getHeight());
    }

    public void setLevel(int level_)
    {
        final int level = level_;
        Runnable runnable = new Runnable()
        {
            public void run()
            {
                switch(level)
                {
                    case ERROR:
                        progressBar.setForeground(Color.RED);
                        break;
                    case WARNING:
                        progressBar.setForeground(Color.ORANGE);
                        break;
                    case INFO:
                        progressBar.setForeground(Color.GREEN);
                        break;
                    default:
                        progressBar.setForeground(Color.BLUE);
                }
            }
        };
        StsToolkit.runLaterOnEventThread(runnable);
    }

    public void setDescriptionAndLevel(String description, int level_)
    {
        this.description = description;
        final int level = level_;

        Runnable setDescriptionRunnable = new Runnable()
        {
            public void run()
            {
                switch(level)
                {
                    case ERROR:
                        progressBar.setForeground(StsColor.RED.getColor());
                        break;
                    case WARNING:
                        progressBar.setForeground(StsColor.DARKORANGE.getColor());
                        break;
                    case INFO:
                        progressBar.setForeground(StsColor.GREEN.getColor());
                        break;
                    default:
                        progressBar.setForeground(StsColor.BLUE.getColor());
                }
                progressBar.setString(getString());
            }
        };
        StsToolkit.runLaterOnEventThread(setDescriptionRunnable);
    }

    public void setDescription(String description)
    {
        this.description = description;

        Runnable setDescriptionRunnable = new Runnable()
        {
            public void run()
            {
                progressBar.setString(getString());
            }
        };
        StsToolkit.runLaterOnEventThread(setDescriptionRunnable);
    }

    public void setStringImmediate(String string)
    {
        progressBar.setString(string);
        paintImmediately();
    }

    public void setDescriptionandValue(String description, int value)
    {
        progressBarValue = value;
        this.description = description;

        Runnable setDescriptionRunnable = new Runnable()
        {
            public void run()
            {
                progressBar.setString(getString());
                progressBar.setValue(progressBarValue);
            }
        };
        StsToolkit.runLaterOnEventThread(setDescriptionRunnable);
    }

    private String getString()
    {
        int percent = (int)(100.0f * progressBarValue / progressBarMaxValue);
        if(description == null)
            return new String(percent + "%");
        else
            return new String(percent + "% " + description);
    }

    public void initialize(int maxValue)
    {
        setMaximum(maxValue);
        resetProgressBar();
    }

    public void initializeImmediate(int maxValue)
    {
        progressBar.setMaximum(maxValue);
        setValueImmediate(0);
    }

    public void resetProgressBar()
    {
        setValue(0);
        setDescriptionAndLevel("", INFO);
        canceled = false;
    }

    public void finished()
    {
        setValue(progressBarMaxValue);
    }

    public void incrementProgress()
    {
        progressBarValue++;
        setValue(progressBarValue);
    }

    public void cancel()
    {
        canceled = true;
        setDescriptionAndLevel("Interrupted", WARNING);
    }

    public void clearProcess()
    {
        canceled = true;
        StsToolkit.sleep(200);
    }

    public boolean isCanceled() { return canceled; }

    public void appendLine(String line) { }

    static public void main(String[] args)
    {
        final StsProgressBar stsProgressBar = StsProgressBar.constructorWithCancel();
        Runnable winThread = new Runnable()
        {
            public void run()
            {
                StsToolkit.createDialog(stsProgressBar, false);
            }
        };
        StsToolkit.runWaitOnEventThread(winThread);
        ProgressBarTest test = new ProgressBarTest(stsProgressBar);
        Thread testThread = new Thread(test);
        testThread.start();
    }
}

class ProgressBarTest implements Runnable
{
    StsProgressBar progressBar;
    int max = 100;
    int[] intervals = new int[]{0, 20, 50, 60, 100};
    boolean interrupt = false;

    ProgressBarTest(StsProgressBar progressBar)
    {
        this.progressBar = progressBar;
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
        try
        {
            System.out.println("Doing cancel test.");
            System.out.println("First progress run: canceling at 25%");
            for(int n = 1; n <= max; n++)
            {
                if(progressBar.canceled) break;

                progressBar.setValue(n);
                if(n > max / 4)
                    progressBar.canceled = true;
                StsToolkit.sleep(100);
            }
            StsToolkit.sleep(1000);
            System.out.println("Second progress run: testing various levels/colors");
            progressBar.resetProgressBar();
            for(int n = 1; n <= max; n++)
            {
                if(progressBar.canceled) return;

                progressBar.setValue(n);
                if(n <= max / 4)
                    progressBar.setDescriptionAndLevel("Step " + n, StsProgressBar.INFO);
                else if(n < max / 2)
                    progressBar.setDescriptionAndLevel("Step " + n, StsProgressBar.ERROR);
                else if(n < 3 * max / 4)
                    progressBar.setDescriptionAndLevel("Step " + n, StsProgressBar.WARNING);
                StsToolkit.sleep(100);
            }
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }
}
