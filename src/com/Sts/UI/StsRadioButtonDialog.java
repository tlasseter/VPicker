package com.Sts.UI;

import com.Sts.UI.Beans.StsGroupBox;
import com.Sts.UI.Beans.StsJPanel;
import com.Sts.UI.Beans.StsRadioButtonFieldBean;
import com.Sts.Utilities.StsToolkit;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public class StsRadioButtonDialog extends JDialog {
    StsGroupBox buttonBox = new StsGroupBox();
    StsRadioButtonFieldBean[] buttons;
    ButtonGroup buttonGroup = new ButtonGroup();
    String[] items = null;
    String selectedBtn = null;
    protected JButton okayButton = new JButton();
    private BorderLayout paneLayout = new BorderLayout();
       
    public StsRadioButtonDialog(Frame frame, String title, boolean modal, String[] items)
    {
        super(frame, title, modal);
        try
        {
            this.items = items;
            selectedBtn = items[0];
            jbInit();
            pack();
			StsToolkit.centerComponentOnScreen(this);
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    private void jbInit() throws Exception
    {
        StsJPanel panel = new StsJPanel();
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
        buttons = new StsRadioButtonFieldBean[items.length];
        for(int i=0; i<items.length; i++)
        {
            buttons[i] = new StsRadioButtonFieldBean(this, "selectBtn", items[i], buttonGroup);
            buttonBox.addEndRow(buttons[i]);
        }
        panel.addEndRow(buttonBox);
        okayButton.setText("Okay");
        panel.addEndRow(okayButton);

        okayButton.addActionListener(
            new java.awt.event.ActionListener()
            {
                public void actionPerformed(ActionEvent e)
                {
                    okayButton_actionPerformed(e);
                }
            });
        getContentPane().add(panel);
    }

    public String getSelectedButton() { return selectedBtn; }
    public void setSelectBtn(String val)
    {
        selectedBtn = val;
    }
    public String getSelectBtn()
    {
        return selectedBtn;
    }

    protected void okayButton_actionPerformed(ActionEvent e)
    {
        setVisible(false);
    }
}