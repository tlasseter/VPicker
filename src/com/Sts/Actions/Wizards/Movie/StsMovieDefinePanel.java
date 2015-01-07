package com.Sts.Actions.Wizards.Movie;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.StsMovie;
import com.Sts.MVC.*;
import com.Sts.Types.*;

import javax.swing.*;
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

public class StsMovieDefinePanel extends JPanel implements ActionListener
{
    private StsMovieWizard wizard;
    private StsMovieDefine wizardStep;
    private StsModel model = null;

//    private StsMovie movie;

    JPanel jPanel1 = new JPanel();
    JLabel jLabel1 = new JLabel();
    JLabel jLabel2 = new JLabel();
    JLabel jLabel10 = new JLabel();
    JTextField minFrame = new JTextField();
    JTextField maxFrame = new JTextField();
    JTextField frameInc = new JTextField();
    JTextField cumInc = new JTextField();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    JPanel jPanel2 = new JPanel();
    JCheckBox movieLoopChk = new JCheckBox();
    JCheckBox fullSpeedChk = new JCheckBox();
    JTextField frameSpeed = new JTextField();
    JLabel jLabel3 = new JLabel();
    JLabel jLabel6 = new JLabel();
    JPanel jPanel3 = new JPanel();
    JPanel jPanel4 = new JPanel();
    JCheckBox elevationChk = new JCheckBox();
    JCheckBox azimuthChk = new JCheckBox();
    JTextField elevationStart = new JTextField();
    JTextField elevationInc = new JTextField();
    JTextField azimuthStart = new JTextField();
    JTextField azimuthInc = new JTextField();
    JLabel jLabel4 = new JLabel();
    JLabel jLabel5 = new JLabel();
    GridBagLayout gridBagLayout3 = new GridBagLayout();
    JCheckBox cycleChkBox = new JCheckBox();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    GridBagLayout gridBagLayout4 = new GridBagLayout();

    public StsMovieDefinePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsMovieWizard)wizard;
        this.wizardStep = (StsMovieDefine)wizardStep;
        this.model = wizard.getModel();

        try
        {
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        StsRotatedGridBoundingBox boundingBox = model.getProject().getRotatedBoundingBox();
        float[] range = boundingBox.getBoundingBoxRangeData(wizard.getDirection());        
        minFrame.setText(Float.toString(range[0]));
        maxFrame.setText(Float.toString(range[1]));
        frameInc.setText(Float.toString(range[2]));
        wizard.rebuild();
    }

    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout4);
        this.setMinimumSize(new Dimension(300, 400));
        setPreferredSize(new Dimension(300, 400));
        jPanel1.setBorder(BorderFactory.createEtchedBorder());
        jPanel1.setLayout(gridBagLayout1);
        jLabel1.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel1.setText("Range:");
        jLabel2.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel2.setText("Increment:");
        minFrame.setText("1");
        minFrame.setHorizontalAlignment(SwingConstants.RIGHT);
        maxFrame.setText("100");
        maxFrame.setHorizontalAlignment(SwingConstants.RIGHT);
        frameInc.setText("1");
        frameInc.setHorizontalAlignment(SwingConstants.RIGHT);
        cumInc.setText("2");
        cumInc.setHorizontalAlignment(SwingConstants.RIGHT);
        jPanel2.setBorder(BorderFactory.createEtchedBorder());
        jPanel2.setLayout(gridBagLayout2);
        movieLoopChk.setText("Loop");
        fullSpeedChk.setSelected(true);
        fullSpeedChk.setText("Full Speed");
        fullSpeedChk.addActionListener(this);
        frameSpeed.setText(new Integer(StsMovie.FULL_SPEED).toString());
        frameSpeed.setHorizontalAlignment(SwingConstants.RIGHT);
        frameSpeed.setEnabled(false);
        jLabel3.setText("Delay (ms)");
        jLabel10.setText("# of Frames");
        jPanel3.setBorder(BorderFactory.createEtchedBorder());
        jPanel3.setLayout(gridBagLayout3);
        elevationChk.setText("Elevation");
        azimuthChk.setText("Azimuth");
        elevationStart.setText("0");
        elevationInc.setText("5");
        azimuthStart.setText("0");
        azimuthInc.setText("5");
        jLabel4.setText("Start");
        jLabel5.setText("Increment");
        cycleChkBox.setText("Cycle Volumes");
        cycleChkBox.setSelected(false);
        cycleChkBox.addActionListener(this);

        this.add(jPanel1,           new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.NORTH, GridBagConstraints.HORIZONTAL, new Insets(2, 4, 2, 6), 0, 0));
        jPanel1.add(minFrame,  new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(5, 0, 0, 0), 20, 0));
        jPanel1.add(maxFrame,  new GridBagConstraints(2, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(5, 0, 0, 11), 20, 0));
        jPanel1.add(frameInc,  new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(5, 0, 0, 0), 20, 0));
        jPanel1.add(jLabel1,    new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 0, 5, 0), 0, 0));
        jPanel1.add(jLabel2, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 0, 5, 0), 0, 0));

        this.add(jPanel2,  new GridBagConstraints(0, 3, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 4, 6, 5), 0, 0));
        jPanel2.add(jLabel3,  new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(4, 6, 0, 6), 0, 0));
        jPanel2.add(fullSpeedChk,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(4, 5, 0, 12), 0, 0));
        jPanel2.add(frameSpeed,  new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(4, 0, 0, 0), 0, 0));
        jPanel2.add(movieLoopChk,  new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 5, 0, 5), 0, 0));
        jPanel2.add(cycleChkBox,  new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
                ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 5, 0, 0), 0, 0));

        this.add(jPanel3,  new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 4, 0, 5), 0, 0));
        jPanel3.add(azimuthChk,  new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 5, 0, 0), 0, 0));
        jPanel3.add(elevationChk,  new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 5, 0, 0), 0, 0));
        jPanel3.add(elevationInc,  new GridBagConstraints(2, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 0, 0), 20, 0));
        jPanel3.add(azimuthStart,  new GridBagConstraints(1, 2, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 0, 0), 20, 0));
        jPanel3.add(elevationStart,  new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 0, 0), 20, 0));
        jPanel3.add(azimuthInc,  new GridBagConstraints(2, 2, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 0, 0), 20, 0));
        jPanel3.add(jLabel4,  new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 20, 0, 0), 0, 0));
        jPanel3.add(jLabel5,  new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 20, 0, 0), 1, 1));
    }

    public void actionPerformed(ActionEvent e)
    {
        int i;
        Object source = e.getSource();
        if(source == fullSpeedChk)
        {
            if(fullSpeedChk.isSelected())
                frameSpeed.setEnabled(false);
            else
                frameSpeed.setEnabled(true);
        }
        if(source == cycleChkBox)
        {
            if(cycleChkBox.isSelected())
            {
                movieLoopChk.setSelected(true);
            }
        }
    }

    public float[] getRange()
    {
        float[] range = new float[2];

        Float fp = new Float(minFrame.getText());
        range[0] = fp.floatValue();

        fp = new Float(maxFrame.getText());
        range[1] = fp.floatValue();

        return range;
    }

    public float getIncrement()
    {
        Float fp = new Float(frameInc.getText());
        return fp.floatValue();
    }
 
    public int getDelay()
    {
        Integer ip = new Integer(frameSpeed.getText());
        if(fullSpeedChk.isSelected()) return StsMovie.FULL_SPEED;
        return ip.intValue();
    }
    public boolean getLoop()
    {
        return movieLoopChk.isSelected();
    }
    public boolean getCycleVolumes()
    {
        return cycleChkBox.isSelected();
    }
    public boolean getElevation()
    {
        return elevationChk.isSelected();
    }
    public boolean getAzimuth()
    {
        return azimuthChk.isSelected();
    }
    public int getAzimuthStart()
    {
        Integer as = new Integer(azimuthStart.getText());
        return as.intValue();
    }
    public int getElevationStart()
    {
        Integer es = new Integer(elevationStart.getText());
        return es.intValue();
    }
    public int getAzimuthIncrement()
    {
        Integer ai = new Integer(azimuthInc.getText());
        return ai.intValue();
    }
    public int getElevationIncrement()
    {
        Integer ei = new Integer(elevationInc.getText());
        return ei.intValue();
    }

}