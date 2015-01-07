package com.Sts.Actions.Wizards.Movie;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;

import javax.swing.*;
import javax.swing.border.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsMovieSelectPanel extends JPanel implements ActionListener
{
    private StsMovieWizard wizard;
    private StsMovieSelect wizardStep;

    private StsModel model = null;
    private StsMovie selectedMovie = null;

    JList movieList = new JList();
    DefaultListModel movieListModel = new DefaultListModel();

    JButton newMovieButton = new JButton();
    Border border1;
    GridBagLayout gridBagLayout1 = new GridBagLayout();

    StsMovie[] movies;

    public StsMovieSelectPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsMovieWizard)wizard;
        this.wizardStep = (StsMovieSelect)wizardStep;
        try
        {
            jbInit();
            initialize();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    private void initialize()
    {
        model = wizard.getModel();
        movies = (StsMovie[])model.getCastObjectList(StsMovie.class);
        int nMovies = movies.length;
        for(int n = 0; n < nMovies; n++)
            movieListModel.addElement(movies[n].getName());
        movieList.setModel(movieListModel);
    }

    public StsMovie getSelectedMovie()
    {
        if(movieList.isSelectionEmpty()) return null;
        int selectedIndex = movieList.getSelectedIndex();
        return movies[selectedIndex];
    }

    void jbInit() throws Exception
    {
        border1 = BorderFactory.createCompoundBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED,Color.white,Color.white,new Color(178, 178, 178),new Color(124, 124, 124)),BorderFactory.createEmptyBorder(10,10,10,10));
        this.setLayout(gridBagLayout1);
        newMovieButton.setText("New Movie");
        movieList.setBorder(BorderFactory.createEtchedBorder());
        movieList.setMaximumSize(new Dimension(200, 200));
        movieList.setMinimumSize(new Dimension(50, 50));
        movieList.setPreferredSize(new Dimension(200, 200));
        this.add(movieList,   new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(4, 4, 4, 4), 0, 0));
        this.add(newMovieButton,   new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(4, 4, 4, 4), 0, 0));
        newMovieButton.addActionListener(this);
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();

        // select new directory
        if(source == newMovieButton)
        {
            wizard.createNewMovie();
        }
    }
}
