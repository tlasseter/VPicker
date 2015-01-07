package com.Sts.Actions.Wizards.TimeSeries;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;

import javax.swing.*;
import javax.swing.border.*;
import java.awt.*;
import java.awt.event.*;
import java.text.*;
import java.util.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsTimeSeriesSelectPanel extends JPanel implements ActionListener
{
    private StsTimeSeriesWizard wizard;
    private StsTimeSeriesSelect wizardStep;

    private StsModel model = null;
    private StsMovie selectedMovie = null;

    DefaultListModel movieListModel = new DefaultListModel();

    Border border1;

    StsMovie[] movies;
    int nMovies = 0;
    ButtonGroup typeRadio = new ButtonGroup();
    JLabel jLabel2 = new JLabel();
    JTextField startTimeTxt = new JTextField();
    JLabel jLabel3 = new JLabel();
    JTextField frameIncTxt = new JTextField();
    JPanel jPanel2 = new JPanel();
    JLabel jLabel4 = new JLabel();
    JRadioButton elapseRadio = new JRadioButton();
    JRadioButton clockRadio = new JRadioButton();
    ButtonGroup timeGroup = new ButtonGroup();
    JLabel jLabel6 = new JLabel();
    JComboBox movieCombo = new JComboBox();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    JLabel jLabel1 = new JLabel();
    JLabel jLabel5 = new JLabel();
    GridBagLayout gridBagLayout1 = new GridBagLayout();

    public StsTimeSeriesSelectPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsTimeSeriesWizard)wizard;
        this.wizardStep = (StsTimeSeriesSelect)wizardStep;
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
        movieCombo.removeAllItems();
        movies = (StsMovie[])wizard.getModel().getCastObjectList(StsMovie.class);
        nMovies = movies.length;
        for(int i=0; i<nMovies; i++)
            movieCombo.addItem(movies[i]);
    }

    public StsMovie getSelectedMovie()
    {
        return movies[movieCombo.getSelectedIndex()];
    }

    void jbInit() throws Exception
    {
        border1 = BorderFactory.createCompoundBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED,Color.white,Color.white,new Color(178, 178, 178),new Color(124, 124, 124)),BorderFactory.createEmptyBorder(10,10,10,10));
        this.setLayout(gridBagLayout2);
        jLabel2.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel2.setText("Start Time:");
        startTimeTxt.setText("0");
        frameIncTxt.setText("5");
        jLabel4.setText("Frame Increment:");
        jLabel4.setHorizontalAlignment(SwingConstants.RIGHT);
        jPanel2.setBorder(BorderFactory.createEtchedBorder());
        jPanel2.setLayout(gridBagLayout1);
        elapseRadio.setSelected(true);
        elapseRadio.setText("Elapsed");
        elapseRadio.addActionListener(this);
        clockRadio.addActionListener(this);
        clockRadio.setText("Clock");
        jLabel6.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel6.setText("Movie:");
        movieCombo.addActionListener(this);
        jLabel1.setText("Seconds");
        jLabel5.setText("Milliseconds");
        jPanel2.add(movieCombo,  new GridBagConstraints(1, 0, 4, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(4, 0, 0, 2), 158, 0));
        jPanel2.add(jLabel6,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(4, 4, 0, 0), 0, 0));
        jPanel2.add(frameIncTxt,  new GridBagConstraints(2, 3, 2, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 6, 0), 63, 3));
        jPanel2.add(startTimeTxt,  new GridBagConstraints(2, 2, 2, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(11, 0, 0, 0), 63, 3));
        jPanel2.add(jLabel2,  new GridBagConstraints(1, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(11, 0, 0, 0), 9, 9));
        jPanel2.add(jLabel4,  new GridBagConstraints(0, 3, 2, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 4, 6, 0), 3, 9));
        jPanel2.add(jLabel1,  new GridBagConstraints(4, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(11, 0, 0, 2), 18, 9));
        jPanel2.add(jLabel5,  new GridBagConstraints(4, 3, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 6, 2), 18, 9));
        jPanel2.add(elapseRadio,  new GridBagConstraints(1, 1, 2, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, -3));
        jPanel2.add(clockRadio,  new GridBagConstraints(3, 1, 2, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(6, 49, 0, 22), 0, -3));
        timeGroup.add(elapseRadio);
        timeGroup.add(clockRadio);
        this.add(jPanel2,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(6, 4, 6, 4), 49, 0));
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
        if((source == elapseRadio) || (source == clockRadio))
        {
            if(elapseRadio.isSelected())
                jLabel1.setText("Seconds");
            else
                jLabel1.setText("Time(hh:mm:ss)");
        }
    }

    public boolean isClockTime()
    {
        return clockRadio.isSelected();
    }

    public long getStartTime()
    {
        Date date = new Date(0L);
        try
        {
        if(isClockTime())
        {
            SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd hh:mm:ss.S");
            String time = startTimeTxt.getText();
            time = wizard.cleanTimeString(time);
            if(time == null)
            {
                new StsMessage(wizard.frame, StsMessage.ERROR,
                    "Invalid time value (##:##:##) in movie selection step.\n" +
                    "\n   Solution: Return to step and re-enter valid start time.");
                return 0L;
            }
            date = format.parse(time);
            return date.getTime();
        }
        else
            return (long)(Long.valueOf(startTimeTxt.getText()).intValue() * 1000);
        }
        catch(Exception e)
        {
            new StsMessage(wizard.frame, StsMessage.ERROR, "Unable to convert supplied moive start to time.\n" +
                           "\n     Syntax: ##:##:## \n");
            return 0L;
        }
    }

    public int getIncrement()
    {
        return Integer.valueOf(frameIncTxt.getText()).intValue();
    }
    public StsMovie getMovie()
    {
        return (StsMovie)movieCombo.getSelectedItem();
    }

}
