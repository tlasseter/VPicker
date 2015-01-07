
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Movie;

import com.Sts.Actions.Movie.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;

public class StsMovieWizard extends StsWizard
{
    public StsMovieSelect movieSelect = null;
    public StsMovieType movieType = null;
    public StsMovieDefine movieDefine = null;
//    public StsMovieVolumeDefine volumeDefine = null;
//    public StsMovieDefine movieDefine = null;
//    public StsMovieDefineVolume volumeMovieDefine = null;
    public StsMovieCreate movieCreate = null;

    private StsMovie selectedMovie = null;

    private StsWizardStep[] mySteps =
    {
        movieSelect = new StsMovieSelect(this),
        movieType = new StsMovieType(this),
        movieDefine = new StsMovieDefine(this),
        movieCreate = new StsMovieCreate(this)
    };

    public StsMovieWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("Movie Setup");
        dialog.getContentPane().setSize(300, 300);

        if(!super.start()) return false;

        StsMovieClass movieClass = StsMovie.getMovieClass();
        int nMovies = movieClass.getSize();

        StsMovie currentMovie = movieClass.getCurrentMovie();
        if(currentMovie == null && nMovies > 0)
        {
            currentMovie = (StsMovie)movieClass.getLast();
            movieClass.setCurrentObject(currentMovie);
        }
        return true;
    }

    public void checkAddToolbar()
    {
        if(currentStep == null)
            return;
        StsMovieActionToolbar movieToolbar = (StsMovieActionToolbar)model.win3d.getToolbarNamed(StsMovieActionToolbar.NAME);
        if(movieToolbar == null)
        {
            StsMovieAction movieAction = new StsMovieAction(model, selectedMovie);
            movieToolbar = new StsMovieActionToolbar(model.win3d, movieAction);
            model.win3d.addToolbar(movieToolbar);
            // selectedMovie.setToolbar(movieToolbar);
        }
        movieToolbar.comboBoxSetItem(selectedMovie);
    }

    public void createNewMovie()
    {
        gotoStep(movieType);
        enableFinish();
    }

    public boolean end()
    {
        return super.end();
    }

    public void previous()
    {
       gotoPreviousStep();
    }

    public void next()
    {
        if(currentStep == movieSelect)
        {
            selectedMovie = movieSelect.getSelectedMovie();
            if(selectedMovie == null)
            {
                disableFinish();
                new StsMessage(frame, StsMessage.INFO, "No movie selected.");
                return;
            }
            else
            {
                enableFinish();
                model.setCurrentObject(selectedMovie);
                checkAddToolbar();
                super.finish();
            }
        }
        else
        {
            gotoNextStep();
        }

    }

    public void finish()
    {
        checkAddToolbar();
        super.finish();
    }

    public void setSelectedMovie(StsMovie movie) { selectedMovie = movie; }
    public int getDirection() { return movieType.getDirection(); }
    public String getMovieName() {  return movieType.getMovieName(); }
    public float[] getRange()
    {
        return movieDefine.getRange();
    }
    public float getIncrement()
    {
        return movieDefine.getIncrement();
    }
    public int getDelay()
    {
        return movieDefine.getDelay();
    }
    public boolean getLoop()
    {
        return movieDefine.getLoop();
    }
    public boolean getCycleVolumes()
    {
        return movieDefine.getCycleVolumes();
    }
    public boolean getElevation()
    {
        return movieDefine.getElevation();
    }
    public boolean getAzimuth()
    {
        return movieDefine.getAzimuth();
    }
    public int getElevationStart()
    {
        return movieDefine.getElevationStart();
    }
    public int getAzimuthStart()
    {
        return movieDefine.getAzimuthStart();
    }
    public int getAzimuthIncrement()
    {
        return movieDefine.getAzimuthIncrement();
    }
    public int getElevationIncrement()
    {
        return movieDefine.getElevationIncrement();
    }
}
