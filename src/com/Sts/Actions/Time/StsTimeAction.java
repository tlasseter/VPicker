package com.Sts.Actions.Time;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: Movie action class. Performs all the actions associated with the running of a movie.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author Stuat A. Jackson
 * @version 1.0
 */

import com.Sts.Actions.*;
import com.Sts.MVC.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;

public class StsTimeAction extends StsAction
{
    private StsModel model = null;
    StsTimeActionToolbar tb = null;
    StsProject project = null;
    private vcr playThread = null;
//    private int timeIncrement = 250; // How often to send update to model and toolbar clock (250 milliseconds)
    private boolean pause = false;

    /**
     * Movie action constructor
     * @param actionManager the current action manager
     */
    public StsTimeAction(StsActionManager actionManager, StsTimeActionToolbar toolbar, StsProject project)
    {
        super(actionManager);
        this.model = actionManager.getModel();
        tb = toolbar;
        this.project = project;
        return;
    }

    /**
     * Not used - required by StsAction3d
     * @return false since never run
     */
    public boolean start() { return false; }
    /**
     * <p>Title: VCR Controls thread</p>
     * <p>Description: Thread class used to run the movie while still allowing full
     * control of the application. Each frame the movie state is checked and adjustments
     * are made or the movie is stopped.</p>
     * <p>Copyright: Copyright (c) 2002</p>
     * <p>Company: </p>
     * @author Stuart A. Jackson
     * @version 1.0
     */
    class vcr extends Thread
    {
        /**
         * constructor. Automatically instantiates the object and starts the movie.
         */
        vcr() { super(); }
        /**
         * Run the movie. Automatically executed on the instantiation of the VCR object
         */
        public void run()
        {
            play();
            return;
        }
        /**
         * Run the movie
         */
        public void play()
        {
            while(true)
            {
                try { this.sleep(project.getTimeUpdateRate()); }
                catch (Exception e) {}

                StsToolkit.runLaterOnEventThread( new Runnable() { public void run()
                {
                    if(!pause)
                    {
                        if(!project.incrementProjectTime(project.getTimeUpdateRate() * tb.getSpeed()))
                            tb.runRealTime();
                        tb.updateTime();
                    }
                } } );
            }
        }
    }

    public void playAction()
    {
        pause = false;
        chkThread();
        return;
    }

    public boolean isRunning()
    {
        return !pause;
    }

    public void stopAction()
    {
        pause = true;
        return;
    }

    private void chkThread()
    {
        if(playThread == null)
        {
            playThread = new vcr();
            playThread.start();
        }
        return;
    }
    
    public boolean end()
    {
    	killThread();
    	return true;
    }
    
    public void killThread()
    {
    	stopAction();
    	if(playThread != null)
    	{
    		playThread = null;
    	}
    }
}
