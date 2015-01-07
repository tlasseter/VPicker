//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

// dialog box to display a list

package com.Sts.UI;

import com.Sts.UI.Beans.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

public class StsStringBeanDialog extends JDialog
{
    private boolean useCancelBtn;

    private JPanel buttonPanel = new JPanel();
    private JLabel label = new JLabel();
    protected StsStringFieldBean stringBean;
    protected JButton cancelButton = null;
    protected JButton okayButton = new JButton();



    static private final Dimension MIN_SIZE = new Dimension(200, 200);
    GridBagLayout gridBagLayout1 = new GridBagLayout();

    public StsStringBeanDialog(JFrame parent, String title, String text, StsStringFieldBean stringBean)
    {
        this(parent, title, text, stringBean, true, false);
    }

    public StsStringBeanDialog(Frame parent, String title, String text, StsStringFieldBean stringBean, boolean modal, boolean useCancelBtn)
    {
        super(parent, title, modal);
        label.setText(text);
        this.stringBean = stringBean;
        this.useCancelBtn = useCancelBtn;
        jbInit();
        pack();
        if(parent != null) setLocationRelativeTo(parent);
        setVisible(true);
    }

    private void jbInit()
    {
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);

        Container contentPane = getContentPane();
        contentPane.setLayout(gridBagLayout1);
        buttonPanel.setLayout(new FlowLayout());

        okayButton.setText("Okay");
        buttonPanel.add(okayButton);
        if (useCancelBtn)
        {
            cancelButton = new JButton();
            cancelButton.setText("Cancel");
            buttonPanel.add(cancelButton);
        }
        contentPane.add(label,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 10, 10, 10), 0, 0));
        contentPane.add(stringBean,  new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(10, 10, 10, 10), 0, 0));
        contentPane.add(buttonPanel,  new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(10, 10, 10, 10), 0, 0));

        okayButton.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(ActionEvent e)
            {
                okayButtonAction(e);
            }
        });
        if (cancelButton != null)
        {
            cancelButton.addActionListener(new java.awt.event.ActionListener()
            {
                public void actionPerformed(ActionEvent e)
                {
                    cancelButtonAction(e);
                }
            });
        }
    }

    protected void okayButtonAction(ActionEvent e)
    {
         setVisible(false);
         dispose();
    }

    protected void cancelButtonAction(ActionEvent e)
    {
        setVisible(false);
        dispose();
    }

    public StsStringFieldBean getStringFieldBean()
    {
        return stringBean;
    }

    public void setComboBoxFieldBean(StsStringFieldBean stringBean) { this.stringBean = stringBean; }

    public String getTestItem() { return ""; }
    public void setTestItem(String testItem) { System.out.println("setTestItem: " + testItem); }

    public static void main(String[] args)
    {
        StsStringFieldBean stringFieldBean = new StsStringFieldBean(StsStringBeanDialog.class, "testItem", "value", true, "Test Items");
        StsStringBeanDialog d = new StsStringBeanDialog(null, "StsStringBeanDialog test", "text:", stringFieldBean);
        stringFieldBean.setBeanObject(d);
        d.setModal(true);
        d.pack();
        d.setVisible(true);
        System.exit(0);
    }
}
