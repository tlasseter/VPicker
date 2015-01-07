
/**
 * <p>Title: S2S Development</p>
 * <p>Description: Movie toolbar class. Defines the toolbar that is presented
 * once movies have been defined via the movie workflow step.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author Stuat A. Jackson
 * @version 1.0
 */

package com.Sts.UI.Toolbars;

import com.Sts.Actions.Wizards.TimeSeries.*;
import com.Sts.DB.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.*;

import javax.swing.*;
import javax.swing.border.*;

public class StsTimeSeriesToolbar extends StsToolbar implements StsSerializable
{
    /** Toolbar name used to reference the toolbar throughout the code */
    public static final String NAME = "Time Series Toolbar";
    public static final boolean defaultFloatable = true;

    /** Unique identifier button name to start the movie - moviePlay */
    public static final String PLAY_ACTION = "moviePlay";
    public static final int PLAY = 0;
    /** Unique identifier button name to stop the movie and return to start - movieStop */
    public static final String STOP_ACTION = "movieStop";
    public static final int STOP = 1;
    /** Unique identifier button name to move to beginning of movie - movieStart */
    public static final String START_ACTION = "movieStart";
    public static final int START = 2;
    /** Unique identifier button name to move to end of movie - movieEnd */
    public static final String END_ACTION = "movieEnd";
    public static final int END = 3;
    /** Unique identifier button name to pause the move - moviePause */
    public static final String PAUSE_ACTION = "moviePause";
    public static final int PAUSE = 4;
    /** Unique identifier button name to play movie in reverse - movieReverese */
    public static final String REVERSE_ACTION = "movieReverse";
    public static final int REVERSE = 5;

    transient public StsModel model = null;
    transient private StsTimeSeriesPlot wizardStep;

    /** Developer Controlled Buttons */
    transient private StsToggleButton revBtn = null;
    transient private StsToggleButton playBtn = null;
    transient private StsToggleButton pauseBtn = null;
    transient private StsToggleButton stopBtn = null;
    transient ButtonGroup btnGroup = new ButtonGroup();

    /**
     * Movie toolbar constructor
     */
    public StsTimeSeriesToolbar()
     {
         super(NAME);
     }

    public StsTimeSeriesToolbar(StsWin3dBase win3d, StsTimeSeriesPlot wizardStep)
    {
        super(NAME);
        this.wizardStep = wizardStep;
        initialize(win3d);
    }

    public boolean initialize(StsWin3dBase win3d)
    {
        this.model = win3d.getModel();

        Border border = BorderFactory.createEtchedBorder();
        setBorder(border);

        Object[] nullObject = new Object[0];
        add(new StsButton(START_ACTION, "Return to Start of Movie", wizardStep, "startAction"));
        revBtn = new StsToggleButton(REVERSE_ACTION, "Play Movie in Reverse", wizardStep, "reverseAction");
        add(revBtn);
        stopBtn = new StsToggleButton(STOP_ACTION, "Stop and Return to Start of Movie", wizardStep, "stopAction");
        add(stopBtn);
        pauseBtn = new StsToggleButton(PAUSE_ACTION, "Pause on Current Frame", wizardStep, "pauseAction");
        add(pauseBtn);
        playBtn = new StsToggleButton(PLAY_ACTION, "Play from Current Frame", wizardStep , "playAction");
        add(playBtn);
        add(new StsButton(END_ACTION, "Go to the End of the Movie", wizardStep, "endAction"));

        btnGroup.add(playBtn);
        btnGroup.add(revBtn);
//        btnGroup.add(pauseBtn);
        btnGroup.add(stopBtn);
        setMinimumSize();
        return true;
    }

    private JTextField createTextField()
    {
        return null;
    }

/*
    public void setState(int state)
    {
        switch(state)
        {
            case START:
            case END:
                stopBtn.setSelected(true);
                revBtn.setEnabled(true);
                playBtn.setEnabled(true);
                break;
            case STOP:
                stopBtn.setSelected(true);
                revBtn.setEnabled(true);
                playBtn.setEnabled(true);
                break;
            case PAUSE:
                revBtn.setEnabled(true);
                playBtn.setEnabled(true);
                break;
            case PLAY:
                playBtn.setSelected(true);
                revBtn.setEnabled(true);
                break;
            case REVERSE:
                revBtn.setSelected(true);
                playBtn.setEnabled(true);
                break;
        }
        repaint();
    }
    */
}

