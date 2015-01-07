package com.Sts.UI;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.1
 */

import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.lang.reflect.*;

public class StsEditColorListPanel extends JPanel
{
    StsColor[] colors = StsColor.colors8;
    String[] names = StsColor.colorNames8;

    int selectedColorIndex = 0;
    String selectedText = null;
    Method selectListenerMethod = null;

    StsColorListComboBox comboBox;
    GridBagConstraints constraints;
    JPanel selectPanel;
    StsColorListComboBox colorComboBox;
    JTextField textField;
    JButton okButton;

    public StsEditColorListPanel(String name)
    {
        try
        {
            setName(name);
            setLayout(new GridBagLayout());

            comboBox = new StsColorListComboBox("new", null, 30, 10, true);
            comboBox.setToolTipText("Select polygon type or create a new type.");

            selectPanel = new JPanel();
            selectPanel.setLayout(new BorderLayout());
//            selectPanel.setLayout(new GridBagLayout());
            colorComboBox = new StsColorListComboBox(null, colors, 15, 15, false);
            colorComboBox.setToolTipText("Select a color to represent the lithologic type.");
            colorComboBox.setMaximumSize(new Dimension(30, 15));
            selectPanel.add(colorComboBox, BorderLayout.WEST);
        /*
            GridBagConstraints constraints1 = new GridBagConstraints(0, 0, 1, 1, 0.0, 1.0,
                GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(4, 4, 4, 4), 0, 0);
            selectPanel.add(colorComboBox, constraints1);
        */
            textField = new JTextField();
            textField.setToolTipText("Specify the lithologic type.");
            textField.setEditable(true);
            textField.setPreferredSize(new Dimension(60, 10));
//            textField.setFont(new java.awt.Font("Dialog", 0, 10));
            selectPanel.add(textField, BorderLayout.CENTER);
        /*
            GridBagConstraints constraints2 = new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0,
                GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(4, 4, 4, 4), 0, 0);
            selectPanel.add(textField, constraints2);
        */
            okButton = new JButton();
            okButton.setToolTipText("Accept the current lithologic definition (Color/Type).");
            okButton.setText("OK");
//            okButton.setFont(new java.awt.Font("Dialog", 0, 10));
            okButton.setMaximumSize(new Dimension(15, 15));
            selectPanel.add(okButton, BorderLayout.EAST);
        /*
            GridBagConstraints constraints3 = new GridBagConstraints(2, 0, 1, 1, 0.0, 1.0,
                GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(4, 4, 4, 4), 0, 0);
            selectPanel.add(okButton, constraints3);
        */
            constraints = new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0,
                GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(4, 4, 4, 4), 0, 0);

            this.add(comboBox, constraints); // selectPanel will be swapped in/out with action
//            this.add(selectPanel, constraints);

            colorComboBox.setSelectedIndex(0);

            okButton.addActionListener(new java.awt.event.ActionListener()
            {
                public void actionPerformed(ActionEvent e)
                {
                    okButtonAction(e);
                }
            });

            comboBox.addActionListener(new java.awt.event.ActionListener()
            {
                public void actionPerformed(ActionEvent e)
                {
                    comboBoxAction(e);
                }
            });

            colorComboBox.addActionListener(new java.awt.event.ActionListener()
            {
                public void actionPerformed(ActionEvent e)
                {
                    colorComboBoxAction(e);
                }
            });

            textField.addActionListener(new java.awt.event.ActionListener()
            {
                public void actionPerformed(ActionEvent e)
                {
                    textFieldAction();
                }
            });

            textField.addFocusListener(new java.awt.event.FocusAdapter()
            {
                public void focusLost(FocusEvent e)
                {
                    textFieldAction();
                }
            });
//        createItem("new", null); // this will toggle display to button, so toggle back
            comboBox.setSelectedIndex(0);
            toggleDisplay(true);
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }

    public StsEditColorListPanel(String name, Class listenerClass, String listenerMethodName, Class[] args)
    {
        this(name);

        try
        {
            selectListenerMethod = listenerClass.getDeclaredMethod(listenerMethodName, args);
            if(selectListenerMethod == null)
            {
                StsException.systemError("StsEditColorListPanel(listener...) failed: couldn't find class.method: " +
                    listenerClass.getName() + "." + listenerMethodName);
                return;
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsEditColorListPanel(listener...) failed.",
                e, StsException.WARNING);
        }
    }

    public StsColorListComboBox getComboBox() { return comboBox; }

    private void createItem(String name, StsColor color)
    {
        comboBox.addItem(name, color);
        int nItems = comboBox.getItemCount();
        comboBox.setSelectedIndex(nItems-1);
    }

    private void notifyColorListener(StsColor color, String name)
    {
        try
        {
            if(selectListenerMethod != null) selectListenerMethod.invoke(null, new Object[] { color, name } );
        }
        catch(Exception e)
        {
            StsException.systemError("StsEditColorListPanel.createItem() failed.\n" +
                "Couldn't invoke method: " + selectListenerMethod.getName());
        }
    }

    public void toggleDisplay(boolean showComboBox)
    {
        if(showComboBox)
        {
            remove(selectPanel);
            add(comboBox, constraints);
        }
        else
        {
            remove(comboBox);
            add(selectPanel, constraints);
        }
        this.validate();
        this.repaint();
    }

    void comboBoxAction(ActionEvent e)
    {
        int index = comboBox.getSelectedIndex();
        if(index == 0)
            toggleDisplay(false);
        else
            notifyColorListener(comboBox.getSelectedColor(), comboBox.getSelectedName());
    }

    void colorComboBoxAction(ActionEvent e)
    {
        selectedColorIndex = colorComboBox.getSelectedIndex();
        textField.setText(names[selectedColorIndex]);
        textField.requestFocus();
        notifyColorListener(null, null);
    }

    void textFieldAction()  // invoked by ActionEvent and FocusEvent(lost)
    {
        selectedText = textField.getText();
        notifyColorListener(null, null);
    }

    void okButtonAction(ActionEvent e)
    {
        createItem(selectedText, colors[selectedColorIndex]);
        notifyColorListener(colors[selectedColorIndex], selectedText);
        toggleDisplay(true);
    }

    public void setEditable(boolean enabled)
    {
        comboBox.setEnabled(enabled);
        colorComboBox.setEnabled(enabled);
        okButton.setEnabled(enabled);
    }

    static void main(String[] args)
    {
        StsEditColorListPanel panel = new StsEditColorListPanel(null);
        StsToolkit.createDialog(panel);
    }
}
