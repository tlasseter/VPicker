
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Movie;

import com.Sts.Actions.Movie.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;

public class StsMovieCreate extends StsWizardStep
{
    public StsMovieWizard wizard;
    public StsStatusPanel panel;
    private StsHeaderPanel header;
    // private StsMovieActionToolbar movieToolbar;
    private StsMovie movie;

    public StsMovieCreate(StsWizard wizard)
    {
        super(wizard, new StsStatusPanel(), null, new StsHeaderPanel());

//        movie = StsMovie.getInstance();
        this.wizard = (StsMovieWizard)wizard;
        panel = (StsStatusPanel) getContainer();
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Movie Definition");
        header.setSubtitle("Create Movie");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Movie");        
    }

    public boolean start()
    {
        panel.setTitle("Creating Movie");
        run();
        return true;
    }

    public void run()
    {
        try
        {
            disableFinish();
            // Setup the Movie Object
            panel.setProgress(0.0f);
            panel.setText("Configuring Movie");

            movie = new StsMovie(model);
            movie.setType(StsMovie.LINE);
            movie.setDirection(wizard.getDirection());
            movie.setRange(wizard.getRange());
            movie.setIncrement(wizard.getIncrement());
            movie.setName(wizard.getMovieName());
            movie.setDelay(wizard.getDelay());
            movie.setLoop(wizard.getLoop());
            movie.setCycleVolumes(wizard.getCycleVolumes());
            movie.setElevationAnimation(wizard.getElevation());
            movie.setElevationIncrement(wizard.getElevationIncrement());
            movie.setElevationStart(wizard.getElevationStart());
            movie.setAzimuthAnimation(wizard.getAzimuth());
            movie.setAzimuthIncrement(wizard.getAzimuthIncrement());
            movie.setAzimuthStart(wizard.getAzimuthStart());

            panel.setText("Configuring Movie Toolbar");
            panel.setProgress(50.0f);

            // Enable the Toolbar
            StsMovieActionToolbar tb = (StsMovieActionToolbar)model.win3d.getToolbarNamed(StsMovieActionToolbar.NAME);
            if(tb == null)
            {
                StsMovieAction movieAction = new StsMovieAction(model, movie);
                tb = new StsMovieActionToolbar(model.win3d, movieAction);
                // movie.setToolbar(movieToolbar);
                model.win3d.addToolbar(tb);
            }
            else
            {
                // movie.setToolbar(tb);
                tb.getMovieAction().setMovie(movie);
            }
            wizard.setSelectedMovie(movie);

            panel.setText("Successfully created movie. Use toolbar to play.");
            panel.setProgress(100.0f);
            enableFinish();
            disableCancel();
            wizard.getActionManager().endCurrentAction();
            movie.addToModel();
        }
        catch(Exception e)
        {
            success = false;
        }
    }

    public boolean end()
    {
        return true;
    }

}
