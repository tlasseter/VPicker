
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.TimeSeries;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.Toolbars.*;

import java.awt.*;

public class StsTimeSeriesPlot extends StsWizardStep
{
    StsTimeSeriesPlotPanel panel;
    StsHeaderPanel header;
    boolean abort = false;
    Thread threadVCR = null;
    Runnable runVCR = null;
    StsTimeSeriesToolbar tb = null;
    StsMovie movie = null;
    int delay = 0;

    private boolean stop = false;
    private boolean pause = false;
    private boolean forward = true;
    private boolean killIt = false;
    private boolean debug = false;

    public StsTimeSeriesPlot(StsWizard wizard)
    {
        super(wizard);
        panel = new StsTimeSeriesPlotPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, null);
        panel.setPreferredSize(new Dimension(150, 800));
        header.setTitle("Time Series Plot");
        header.setSubtitle("Time Series Graphics");
        header.setInfoText(wizardDialog,"\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#TimeSeries");
        runVCR = new Runnable()
        {
            public void run()
            {
                processVCR();
            }
        };

    }

    public boolean start()
    {
        stop = false;
        pause = false;
        forward = true;
        killIt = false;
        movie = null;
		threadVCR = new Thread(runVCR);
        tb = new StsTimeSeriesToolbar(model.win3d, this);
        wizardDialog.addToolbar(tb);
        wizard.enableFinish();
        disablePrevious();
        disableCancel();
        panel.initialize();
        movie = ((StsTimeSeriesWizard)wizard).getMovie();
        delay = movie.getDelay();
        return panel.initializeMovie();
    }

    public boolean end()
    {
        if(threadVCR.isAlive())
            killIt = true;
        stopAction();
        return true;
    }

    public void destroyChart()
    {
        if(threadVCR.isAlive())
            killIt = true;
        wizardDialog.removeToolbar(tb);
        panel.destroyChart();
    }

    public void startAction()
    {
        delay = movie.getDelay();
		panel.startFrame();
    }
    public void reverseAction()
    {
        killIt = false;
		if(!threadVCR.isAlive())
            threadVCR.start();
        forward = false;
        pause = false;
        stop = false;
    }
    public void stopAction()
    {
        stop = true;
        pause = false;
        delay = movie.getDelay();
        panel.stop();
    }
    public void pauseAction()
    {
        delay = movie.getDelay();
        if(pause)
            pause = false;
        else
            pause = true;
        stop = false;
    }
    public void playAction()
    {
        killIt = false;
        if(!threadVCR.isAlive())
            threadVCR.start();
        stop = false;
        pause = false;
        forward = true;

    }
    public void endAction()
    {
        panel.endFrame();
    }

    public void processVCR()
	{
        long time = 0L;
        panel.startFrame();
		while(true)
		{
			try
			{
//                time = System.currentTimeMillis();
                if (debug) System.out.println("Start time:=" + time);

				if((!pause) && (!stop))
				{
					if(forward)
						panel.nextFrame();
					else
						panel.previousFrame();

                    // Compute any remaining delay.
/*                    long cTime = System.currentTimeMillis();
                    if(time != 0L)
                        time = (long)delay - (cTime - time);
                    else
                        time = delay;

                    if(time < 0)
                        time = delay;
*/
                    // Increase the delay if that supplied with movie is not adaquate
//                    if(time <  0L)
//                        delay = delay + Math.abs( (int) time);

                    // Remove any excess delay beyond that set in movie.
//                    if((time > 100) && (delay > (movie.getDelay()+100)))
//                        delay = delay - 100;

//                    if(time > 0)
                        Thread.sleep(delay);
				}
				if((stop) || (pause))
				{
                    Thread.sleep(2000L);
				}
                if(killIt)
                    break;
			}
			catch (Exception e)
			{
			}
		}
        return;
	}
}


