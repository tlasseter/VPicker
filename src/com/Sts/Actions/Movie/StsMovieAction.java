package com.Sts.Actions.Movie;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: Movie action class. Performs all the actions associated with the running of a movie.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author Stuat A. Jackson
 * @version 1.0
 */

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;

public class StsMovieAction extends StsAction
{
    private StsModel model = null;
    private StsMovie movie = null;
    private StsSliderBean slider = null;
    // private StsCursor3d cursor = null;
    private vcr playThread = null;
    private boolean pause = false;
    private int maxPosition = -1;
    private int minPosition = -1;
    private int pausePosition = -1;
    private int relativeIncrement = -1;
    private int relativeStart = -1;
    private int relativeEnd = -1;
    private boolean stop = true;
    private boolean paused = false;
    private boolean first = true;
    private boolean save = true;
    private float[] range = null;
    private boolean reverse = false;
    private int orientation = movie.INCREMENT;
    private boolean init = false;

    // StsMovieActionToolbar tb = null;

    /**
     * Movie action constructor
     */
    public StsMovieAction(StsModel model, StsMovie movie)
    {
        this.model = model;
        setMovie(movie);
    }

    public void setMovie(StsMovie movie)
    {
        if(this.movie != null)
        {
            if(this.movie.isPlaying())
                stopAction();
            // movie.setToolbar(this.movie.getToolbar());
            if(slider == null)
                movie.enableAnimation(false);

        }
        this.movie = movie;
        init = false;
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
                try { this.sleep(200); }
                catch (Exception e) {}

                if(pause) continue;
                if(((!reverse) && (orientation == movie.INCREMENT)) || ((reverse) && (orientation == movie.DECREMENT)))
                {
                    for(int i=relativeStart; i<=relativeEnd; i=i+relativeIncrement)
                    {
                        try
                        {
                            if(movie.getDelay() != -1)
                                this.sleep(movie.getDelay());
                        }
                        catch (Exception e) {}
                        if(stop) {
                            movie.setFrame(minPosition);
                            break;
                        }
                        movie.setFrame(i);
                        if(slider != null)
                        {
                            slider.setValue(i);
                        }
                        if(pause)
                        {
                            pausePosition = i;
                            break;
                        }
                        if(movie.getElevationAnimation())
                            movie.setElevationAndAzimuth(movie.INCREMENT, movie.ELEVATION);
                        if(movie.getAzimuthAnimation())
                            movie.setElevationAndAzimuth(movie.INCREMENT, movie.AZIMUTH);
                        if(save)
                        {
                            // Output gif screen capture
                            ;
                        }
                    }
                    relativeStart = minPosition;
                    relativeEnd = maxPosition;
                }
                else
                {
                    for(int i=relativeStart; i>=relativeEnd; i=i+relativeIncrement)
                    {
                        try
                        {
                            if(movie.getDelay() != -1)
                                this.sleep(movie.getDelay());
                        }
                        catch (Exception e) {}
                        if(stop)
                        {
                            movie.setFrame(minPosition);
                            break;
                        }
                        movie.setFrame(i);
                        if(slider != null) slider.setValue(i);
                        if(pause)
                        {
                            pausePosition = i;
                            break;
                        }
                        if(movie.getElevationAnimation()) movie.setElevationAndAzimuth(movie.DECREMENT, movie.ELEVATION);
                        if(movie.getAzimuthAnimation()) movie.setElevationAndAzimuth(movie.DECREMENT, movie.AZIMUTH);
                        if(save)
                        {
                            // Output gif screen capture
                            ;
                        }
                    }
                    relativeStart = minPosition;
                    relativeEnd = maxPosition;
                }
                if(pause)
                    continue;
                if(!movie.getLoop())
                    stopAction();
                else if(movie.getCycleVolumes())
                    movie.nextVolume();
            }
        }
    }

    public void initMovie()
    {
        StsMovieClass movieClass = (StsMovieClass)model.getStsClass(StsMovie.class);
        movie = (StsMovie)movieClass.getCurrentObject();
        if(movie == null) return;
        // if(init) return;
        // tb = (StsMovieActionToolbar) movie.getToolbar();
        range = movie.getRange();
        maxPosition = (int)range[1];
        minPosition = (int)range[0];
        movie.setFrame(minPosition);
        movie.prepareCursor();
        slider = movie.prepareSlider();
        movie.enableAnimation(true);
    }
    /**
     * Play the movie from the current location
     * @see #reverseAction to play the movie backward
     */
    public void playAction()
    {
        initMovie();
        // tb.setState(tb.PLAY);

        pause = true;

        if(first)
        {
            startAction();
            first = false;
        }
        relativeIncrement = (int) movie.getIncrement();
        if(movie.getCurrentFrame() == maxPosition)
            relativeStart = minPosition;
        else
            relativeStart = (int) movie.getCurrentFrame();
        relativeEnd = maxPosition;

        if(relativeIncrement > 0)
            orientation = movie.INCREMENT;
        else
            orientation = movie.DECREMENT;

        chkThread();

        reverse = false;
        pause = false;
        stop = false;
        movie.setPlaying(true);

        // Disable the reverse button
        // tb.setState(tb.PLAY);
        return;
    }
    /**
     * Play the movie in reverse direction from the current location
     * @see #playAction to play the movie forward
     */
    public void reverseAction()
    {
        initMovie();
        // tb.setState(tb.REVERSE);
        pause = true;

        if(first)
        {
            endAction();
            first = false;
        }

        relativeIncrement = (int) -movie.getIncrement();
        if(movie.getCurrentFrame() == minPosition)
            relativeStart = maxPosition;
        else
            relativeStart = (int) movie.getCurrentFrame();
        relativeEnd = minPosition;

        if(relativeIncrement < 0)
            orientation = movie.INCREMENT;
        else
            orientation = movie.DECREMENT;

        chkThread();

        reverse = true;
        pause = false;
        stop = false;
        movie.setPlaying(true);

        // Disable the reverse button
        // tb.setState(tb.REVERSE);

        return;
    }
    /**
     * Move to the first frame in the movie and start playing
     * @see #stopAction()
     */
    public void startAction()
    {
        initMovie();
        pause = true;
        movie.setPlaying(false);
        if(movie.getType() != StsMovie.VOLUME)
        slider.setValue(minPosition);
        movie.setFrame(minPosition);
        // tb.setState(tb.START);
    }
    /**
     * Stop and move to the last frame of the movie
     * @see #startAction()
     */
    public void endAction()
    {
        initMovie();
        pause = true;
        movie.setPlaying(false);
        movie.prepareCursor();
        slider = movie.prepareSlider();
        if(slider != null)
            slider.setValue(maxPosition);
        movie.setFrame(maxPosition);
        // tb.setState(tb.END);
        return;
    }
    /**
     * Stop and return to the start position
     * @see #pauseAction() to stop at current location
     */
    public void stopAction()
    {
        initMovie();
        pause = true;
        stop = true;
        movie.setPlaying(false);
        if(slider != null)
            slider.setValue(minPosition);
        movie.setFrame(minPosition);
        // tb.setState(tb.STOP);
        return;
    }
    /**
     * Pause at the current frame. Play button will restart from current location.
     */
    public void pauseAction()
    {
        initMovie();
        movie.setPlaying(false);
        pause = true;
        // tb.setState(tb.PAUSE);
        return;
    }
    /**
     * Output each frame to a gif file (not enabled yet)
     */
    public void saveAction()
    {
        save = true;
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
}