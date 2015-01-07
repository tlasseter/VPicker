
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.MVC.View3d;

import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Icons.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

public class StsWin3d_AboutBox extends Dialog implements ActionListener
{

    JPanel panel1 = new JPanel();
    JPanel panel2 = new JPanel();
    JPanel insetsPanel1 = new JPanel();
    //JPanel insetsPanel2 = new JPanel();
    JPanel insetsPanel3 = new JPanel();
    JButton button1 = new JButton();
    //ImageControl imageControl1 = new ImageControl();
    JLabel label1 = new JLabel();
    JLabel label2 = new JLabel();
    JLabel label2a = new JLabel();
    JLabel label3 = new JLabel();
    JLabel label4 = new JLabel();
    JLabel label4a = new JLabel();
    JLabel label4b = new JLabel();
    JLabel label4c = new JLabel(); 
    JLabel label4d = new JLabel(); 
    JLabel label4e = new JLabel();
    JLabel label4f = new JLabel(); 
    JLabel label4g = new JLabel();    
    FlowLayout flowLayout1 = new FlowLayout();
    FlowLayout flowLayout2 = new FlowLayout();
    ImageIcon imageIcon = StsIcon.createIcon("S2Slogo.gif");
    String product = "S2S";
    String version = "Version: " + StsModel.version;
    String build = "Build: ";
    String copyright = "Copyright (c) 2001";
    StsJPanel logoPanel = new StsJPanel();
    JLabel logoLabel = new JLabel();
//    XYLayout xYLayout3 = new XYLayout();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    GridBagLayout gridBagLayout3 = new GridBagLayout();

    public StsWin3d_AboutBox(Frame parent)
    {
        super(parent);
        enableEvents(AWTEvent.WINDOW_EVENT_MASK);
        try
        {
            jbInit();
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
        pack();
    }

    private void jbInit() throws Exception
    {
        this.setTitle("About S2S");
        setResizable(false);
        panel1.setLayout(gridBagLayout3);
        panel2.setLayout(gridBagLayout2);
        insetsPanel1.setLayout(flowLayout1);
        //insetsPanel1.setBorder(...FLAT);
        label1.setFont(new java.awt.Font("Dialog", 1, 12));
        label1.setText("S2S Systems");

        label2.setFont(new java.awt.Font("Dialog", 0, 11));
        label2.setText(version + " 1.0"); // + StsModel.VERSION);

        label2a.setFont(new java.awt.Font("Dialog", 0, 11));
        label2a.setText(build + StsModel.version);

        label3.setFont(new java.awt.Font("Dialog", 0, 11));
        label3.setText("Copyright (c) 2001-2011");

        label4.setFont(new java.awt.Font("Dialog", 0, 11));
        label4.setText("S2S Systems International, Inc.");

        label4a.setFont(new java.awt.Font("Dialog", 0, 11));
        label4a.setText("www.s2ssystems.com");
        
        label4b.setFont(new java.awt.Font("Dialog", 0, 9));
        label4b.setText("Some 2D graphics utilize - Chart2D");
        label4c.setFont(new java.awt.Font("Dialog", 0, 9));
        label4c.setText("Copyright (C) Free Software Foundation, Inc.");
        label4d.setFont(new java.awt.Font("Dialog", 0, 9));
        label4d.setText("Available at Chart2D.sourceforge.net");
        
        label4e.setFont(new java.awt.Font("Dialog", 0, 9));
        label4e.setText("Limited DXF Support via Kabeja");        
        label4f.setFont(new java.awt.Font("Dialog", 0, 9));
        label4f.setText("Copyright 2008 Simon Mieth");
        label4g.setFont(new java.awt.Font("Dialog", 0, 9));
        label4g.setText("Available at kabeja.sourceforge.net");        

        insetsPanel3.setLayout(gridBagLayout1);
        
//        insetsPanel3.setMargins(new Insets(10, 60, 10, 10));
//        insetsPanel3.setBevelInner(JPanel.FLAT);
        button1.setText("OK");
        button1.addActionListener(this);
        //imageControl1.setImageName("");
        //insetsPanel2.add(imageControl1, null);
        //panel2.add(insetsPanel2, BorderLayout.WEST);
//    logoPanel.setLayout(xYLayout3);
        logoLabel.setHorizontalAlignment(SwingConstants.RIGHT);
    logoLabel.setIcon(imageIcon);
    panel2.setBorder(BorderFactory.createEtchedBorder());
    logoPanel.setBorder(null);
    insetsPanel3.setFont(new java.awt.Font("Dialog", 0, 10));
    this.add(panel1, null);
        insetsPanel3.add(label1,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 149, 2));
        insetsPanel3.add(label2,  new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 146, 2));
        insetsPanel3.add(label2a,  new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 146, 2));
        insetsPanel3.add(label3,  new GridBagConstraints(0, 3, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 91, 2));
        insetsPanel3.add(label4,  new GridBagConstraints(0, 4, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 122, 2));
        insetsPanel3.add(label4a,  new GridBagConstraints(0, 5, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 1, 0), 0, 0));
        insetsPanel3.add(label4b,  new GridBagConstraints(0, 6, 1, 1, 0.0, 0.0
                ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 0, 1, 0), 0, 0));
        insetsPanel3.add(label4c,  new GridBagConstraints(0, 7, 1, 1, 0.0, 0.0
                ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 1, 0), 0, 0));
        insetsPanel3.add(label4d,  new GridBagConstraints(0, 8, 1, 1, 0.0, 0.0
                ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 1, 0), 0, 0));
        insetsPanel3.add(label4e,  new GridBagConstraints(0, 9, 1, 1, 0.0, 0.0
                ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 0, 1, 0), 0, 0));
        insetsPanel3.add(label4f,  new GridBagConstraints(0, 10, 1, 1, 0.0, 0.0
                ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 1, 0), 0, 0));
        insetsPanel3.add(label4g,  new GridBagConstraints(0, 11, 1, 1, 0.0, 0.0
                ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 1, 0), 0, 0));        
		panel1.add(insetsPanel1,  new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 4, 5), 188, -3));
		insetsPanel1.add(button1, null);
        panel2.add(logoPanel,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(3, 8, 3, 0), 9, 7));
	    logoPanel.add(logoLabel);
        panel1.add(panel2,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(2, 0, 0, 5), 4, 0));
        panel2.add(insetsPanel3,  new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(3, 0, 3, 5), 0, 0));
    }

    protected void processWindowEvent(WindowEvent e)
    {
        if (e.getID() == WindowEvent.WINDOW_CLOSING)
        {
            cancel();
        }
        super.processWindowEvent(e);
    }

    void cancel()
    {
        dispose();
    }

    public void actionPerformed(ActionEvent e)
    {
        if (e.getSource() == button1)
        {
            cancel();
        }
    }
}



