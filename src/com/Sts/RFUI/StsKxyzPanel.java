

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.RFUI;

import javax.swing.*;
import java.awt.*;


public class StsKxyzPanel extends JPanel
{
    JLabel kxLabel = new JLabel();
    JTextField kxField = new JTextField();
    JLabel kyLabel = new JLabel();
    JTextField kyField = new JTextField();
    JLabel kzLabel = new JLabel();
    JTextField kzField = new JTextField();
    GridBagLayout gridBagLayout1 = new GridBagLayout();


    public StsKxyzPanel()
    {
        try
        {
            jbInit();
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }

    private void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout1);
        kxLabel.setHorizontalAlignment(SwingConstants.RIGHT);
        kxLabel.setText("X Permeability factor");
        kyLabel.setHorizontalAlignment(SwingConstants.RIGHT);
        kyLabel.setText("Y Permeability factor");
        kzLabel.setHorizontalAlignment(SwingConstants.RIGHT);
        kzLabel.setText("Z Permeability factor");

        this.add(kxLabel, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 5, 0, 5), 0, 0));
        this.add(kyLabel, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 5, 0, 5), 0, 0));
        this.add(kzLabel, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 5, 0, 5), 0, 0));
        this.add(kxField, new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
        this.add(kyField, new GridBagConstraints(1, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
        this.add(kzField, new GridBagConstraints(1, 2, 1, 1, 1.0, 1.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
    }


    public void setDouble(JTextField field, double d)
    {
        field.setText(Double.toString(d));
    }

    public double getDouble(JTextField f)
    {
        return Double.valueOf(f.getText()).doubleValue();
    }

    public void setKx(double d) { setDouble(kxField, d); }
    public void setKy(double d) { setDouble(kyField, d); }
    public void setKz(double d) { setDouble(kzField, d); }

    public double getKx() { return getDouble(kxField); }
    public double getKy() { return getDouble(kyField); }
    public double getKz() { return getDouble(kzField); }



    public static void main(String[] args)
    {
        StsKxyzPanel kxyPanel = new StsKxyzPanel();
        JFrame frame = new JFrame();
        frame.getContentPane().add(kxyPanel);
        frame.pack();
        frame.setVisible(true);
    }
}
