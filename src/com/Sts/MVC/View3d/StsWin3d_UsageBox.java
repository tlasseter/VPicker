
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.MVC.View3d;

import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.s2s.client.usage.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

public class StsWin3d_UsageBox extends Dialog implements ActionListener
{

    JPanel panel1 = new JPanel();
    JPanel insetsPanel1 = new JPanel();
    //JPanel insetsPanel2 = new JPanel();
    JButton button1 = new JButton();
    //ImageControl imageControl1 = new ImageControl();
    FlowLayout flowLayout2 = new FlowLayout();
    String product = "S2S";
    String version = "Version: ";
    String copyright = "Copyright (c) 2001";
    String comments = "S2S Systems LLC";
    String comments2 = "   ";
    JPanel jPanel1 = new JPanel();
    JTextArea jTextArea1 = new JTextArea();

    UsageManager usageManager = null;
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    JButton reticketBtn = new JButton();
    JTextArea jTextArea2 = new JTextArea();

    Frame parent = null;
    JPanel jPanel2 = new JPanel();
    JPanel jPanel3 = new JPanel();
    JLabel jLabel2 = new JLabel();
    JLabel jLabel1 = new JLabel();
    JLabel jLabel3 = new JLabel();
    JTextField ticketValidUntil = new JTextField();
    JTextField ticketUsageHoursRemaining = new JTextField();
    JTextField licenseValidUntil = new JTextField();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    GridBagLayout gridBagLayout4 = new GridBagLayout();
    JTextArea jTextArea3 = new JTextArea();
    JPanel jPanel4 = new JPanel();
    JButton transferTktBtn = new JButton();
    GridBagLayout gridBagLayout5 = new GridBagLayout();
    GridBagLayout gridBagLayout6 = new GridBagLayout();
    JLabel jLabel4 = new JLabel();
    JTextField licensedTo = new JTextField();
    GridBagLayout gridBagLayout3 = new GridBagLayout();

    public StsWin3d_UsageBox(Frame parent)
    {
        super(parent);
        this.parent = parent;
        enableEvents(AWTEvent.WINDOW_EVENT_MASK);
        try
        {
            jbInit();

            setTicketUsageHoursRemaining();
            setTicketDeath();
            setLicenseDeath();
            setLicenseTo();
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
        pack();
    }

    private void jbInit() throws Exception
    {
        this.setTitle("S2S License Status");
        setResizable(false);
        insetsPanel1.setLayout(gridBagLayout2);
        button1.setText("OK");
        button1.addActionListener(this);
        insetsPanel1.setBorder(BorderFactory.createEtchedBorder());
        jPanel1.setLayout(gridBagLayout5);
        jPanel1.setBorder(BorderFactory.createEtchedBorder());
        jPanel1.setPreferredSize(new Dimension(400, 450));
        jTextArea1.setBackground(UIManager.getColor("Menu.background"));
        jTextArea1.setFont(new java.awt.Font("Dialog", 0, 10));
        jTextArea1.setBorder(null);
        jTextArea1.setToolTipText("");
        jTextArea1.setEditable(false);
        jTextArea1.setText("Usage is tracked on the server and on the client. After downloading " +
                           "the application the client is issued a ticket which has a maximum " +
                           "usage and lifespan. These values are used only if the client system " +
                           "is detached from the server. An example would be for travel. The " +
                           "following three values will let you better manage your usage and " +
                           "request different licensing from your software provider if these " +
                           "values are to restrictive. If your system remains connected to the " +
                           "server the License Term is the only value that should be of concern. ");
        jTextArea1.setLineWrap(true);
        jTextArea1.setWrapStyleWord(true);
        panel1.setLayout(gridBagLayout6);
        reticketBtn.setText("ReIssue Ticket");
        reticketBtn.addActionListener(this);
        jTextArea2.setWrapStyleWord(true);
        jTextArea2.setLineWrap(true);
        jTextArea2.setText("Reissue the ticket if you plan on running disconnected from the Internet. " +
                           "This action will ensure that your ticket life and usage are at the " +
                           "maximum allowed.");
        jTextArea2.setEditable(false);
        jTextArea2.setToolTipText("");
        jTextArea2.setBorder(null);
        jTextArea2.setFont(new java.awt.Font("Dialog", 0, 10));
        jTextArea2.setBackground(UIManager.getColor("Menu.background"));
        jPanel2.setBorder(BorderFactory.createEtchedBorder());
        jPanel2.setDebugGraphicsOptions(0);
        jPanel2.setLayout(gridBagLayout1);
        jPanel3.setBorder(BorderFactory.createEtchedBorder());
        jPanel3.setDebugGraphicsOptions(0);
        jPanel3.setLayout(gridBagLayout3);
        jLabel2.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel2.setHorizontalTextPosition(SwingConstants.RIGHT);
        jLabel2.setText("Ticket Valid Until:");
        jLabel1.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel1.setHorizontalTextPosition(SwingConstants.RIGHT);
        jLabel1.setText("License Valid Until:");
        jLabel3.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel3.setHorizontalTextPosition(SwingConstants.RIGHT);
        jLabel3.setText("Ticket Usage Hours Remaining:");
        ticketValidUntil.setBorder(BorderFactory.createEtchedBorder());
        ticketValidUntil.setToolTipText("Reconnect to Server to automatically renew");
        ticketValidUntil.setEditable(false);
        ticketValidUntil.setText("NA");
        ticketValidUntil.setHorizontalAlignment(SwingConstants.RIGHT);
        ticketUsageHoursRemaining.setBorder(BorderFactory.createEtchedBorder());
        ticketUsageHoursRemaining.setToolTipText("Reconnect to Server to automatically renew");
        ticketUsageHoursRemaining.setText("NA");
        ticketUsageHoursRemaining.setHorizontalAlignment(SwingConstants.RIGHT);
        ticketUsageHoursRemaining.setEnabled(false);
        licenseValidUntil.setBorder(BorderFactory.createEtchedBorder());
        licenseValidUntil.setToolTipText("Contact software provider to get license renewed.");
        licenseValidUntil.setEditable(false);
        licenseValidUntil.setText("NA");
        licenseValidUntil.setHorizontalAlignment(SwingConstants.RIGHT);
        jTextArea3.setBackground(UIManager.getColor("Menu.background"));
        jTextArea3.setFont(new java.awt.Font("Dialog", 0, 10));
        jTextArea3.setBorder(null);
        jTextArea3.setToolTipText("");
        jTextArea3.setEditable(false);
        jTextArea3.setText("Transfer the current ticket so that your license can be moved to " +
                           "another system. Tickets are node locked, therefore the existing ticket " +
                           "must be terminated before another ticket can be issued on this or " +
                           "another computer. PRESSING THIS BUTTON WILL AUTOMATICALLY EXIT THE " +
                           "APPLICATION.");
        jTextArea3.setLineWrap(true);
        jTextArea3.setWrapStyleWord(true);
        jPanel4.setLayout(gridBagLayout4);
        jPanel4.setDebugGraphicsOptions(0);
        jPanel4.setBorder(BorderFactory.createEtchedBorder());
        transferTktBtn.setActionCommand("ReIssue Ticket");
        transferTktBtn.setText("Transfer Ticket");
        transferTktBtn.addActionListener(this);
        jLabel4.setText("Licensed To:");
        jLabel4.setHorizontalTextPosition(SwingConstants.RIGHT);
        jLabel4.setHorizontalAlignment(SwingConstants.RIGHT);
        licensedTo.setHorizontalAlignment(SwingConstants.RIGHT);
        licensedTo.setText("Unavailable");
        licensedTo.setEditable(false);
        licensedTo.setSelectionStart(11);
        licensedTo.setToolTipText("Contact software provider to get license renewed.");
        licensedTo.setBorder(BorderFactory.createEtchedBorder());
        this.add(panel1, null);
        panel1.add(jPanel1,   new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 0, 5), 0, 0));
        jPanel2.add(jTextArea2,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 3, 0, 7), 0, 10));
        jPanel2.add(reticketBtn,  new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 104, 8, 99), 68, 2));
        jPanel1.add(jPanel4,     new GridBagConstraints(0, 2, 1, 1, 1.0, 0.5
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 1, 4, 2), 0, 0));
        jPanel4.add(jTextArea3, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 3, 0, 7), 0, 10));
        jPanel4.add(transferTktBtn, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 104, 8, 99), 68, 2));
    panel1.add(insetsPanel1,   new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 6, 5), -9, 0));
    insetsPanel1.add(button1, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 171, 2, 171), 0, 0));
        jPanel1.add(jPanel3,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(2, 1, 0, 2), 0, 10));
    jPanel3.add(jTextArea1,  new GridBagConstraints(0, 0, 2, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(2, 2, 0, 5), 0, 5));
        jPanel3.add(jLabel3,  new GridBagConstraints(0, 4, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 10, 14, 0), 12, 0));
        jPanel3.add(jLabel1,  new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(13, 13, 0, 0), 67, 0));
        jPanel3.add(licenseValidUntil,  new GridBagConstraints(1, 2, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 16), 169, 6));
        jPanel3.add(jLabel2,  new GridBagConstraints(0, 3, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(11, 14, 0, 0), 76, 0));
        jPanel3.add(ticketValidUntil,  new GridBagConstraints(1, 3, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 16), 169, 6));
        jPanel3.add(ticketUsageHoursRemaining,  new GridBagConstraints(1, 4, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 14, 16), 169, 6));
        jPanel3.add(jLabel4,  new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(17, 12, 0, 0), 96, 0));
        jPanel3.add(licensedTo,  new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(9, 0, 0, 16), 129, 6));
    jPanel1.add(jPanel2,   new GridBagConstraints(0, 1, 1, 1, 1.0, 0.2
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 1, 0, 2), -2, 1));
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
        else if(e.getSource() == reticketBtn)
        {
            makeReticketRequest();
        }
        else if(e.getSource() == transferTktBtn)
        {
            transferTicket();
        }
    }

    public void transferTicket()
    {
        try
        {
            if(Main.usageTracking)
            {
                usageManager = UsageManager.getInstance();
                try
                {
                    usageManager.transferTicket();
                    new StsMessage(parent, StsMessage.INFO, "Exiting application, re-start on desired computer to complete transfer.");
                    Main.logUsage(true);
                    System.exit(0);
                }
                catch (Exception ex)
                {
                    new StsMessage(parent, StsMessage.WARNING, "Unable to reissue ticket: Probably disconnected from Internet, reconnect and try again");
                    ex.printStackTrace();
                }
            }
            else
                new StsMessage(parent, StsMessage.WARNING, "Usage management not enabled.");
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void makeReticketRequest()
    {
        try
        {
            if(Main.usageTracking)
            {
                usageManager = UsageManager.getInstance();
                try
                {
                    usageManager.resetTicket();
                    Main.logUsage();
                    setTicketUsageHoursRemaining();
                    setTicketDeath();
                }
                catch (Exception ex)
                {
                    new StsMessage(parent, StsMessage.WARNING, "Unable to reissue ticket: Probably disconnected from Internet, reconnect and try again");
                    ex.printStackTrace();
                }
            }
            else
                new StsMessage(parent, StsMessage.WARNING, "Usage management not enabled.");
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void setLicenseTo()
    {
        try
        {
            if(Main.usageTracking)
            {
                usageManager = UsageManager.getInstance();
                licensedTo.setText(usageManager.getLicensedTo());
            }
            else
                licensedTo.setText("Usage Tracking Not Active");
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void setLicenseDeath()
    {
        try
        {
            if(Main.usageTracking)
            {
                usageManager = UsageManager.getInstance();
                licenseValidUntil.setText(usageManager.getLicenseDeathDate().toString());
            }
            else
                licenseValidUntil.setText("Usage Tracking Not Active");
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void setTicketDeath()
    {
        try
        {
            if(Main.usageTracking)
            {
                usageManager = UsageManager.getInstance();
                ticketValidUntil.setText(usageManager.getTicketDeathDate().toString());
            }
            else
                ticketValidUntil.setText("Usage Tracking Not Active");
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void setTicketUsageHoursRemaining()
    {
        try
        {
            if(Main.usageTracking)
            {
                usageManager = UsageManager.getInstance();
                Float timeLeft = new Float(usageManager.getRemainingUsage()/60.f);
                ticketUsageHoursRemaining.setText(timeLeft.toString());
            }
            else
                ticketUsageHoursRemaining.setText("Usage Tracking Not Active");
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }
}



