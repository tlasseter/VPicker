//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

// dialog box to display a list

package com.Sts.UI;

import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.lang.reflect.*;

public class StsComboBoxDialog
{
    JDialog dialog;
    private boolean addCancelButton;
    protected Object[] items = null;
    protected Object selectedItem = null;
    protected StsComboBoxFieldBean comboBoxBean;
    private StsJPanel panel = new StsJPanel();
    private StsJPanel buttonPanel = new StsJPanel();
    private Object instance;
    private String fieldName;
    private Field field;
    protected StsButton cancelButton = new StsButton("Cancel", "Cancel operation.", this, "cancel");
    protected StsButton okayButton = new StsButton("OK", "Complete operation.", this, "ok");

    static private final Dimension MIN_SIZE = new Dimension(200, 200);

    public StsComboBoxDialog(JFrame parent, String title, String label, Object[] items, Object instance, String fieldName)
    {
        this(parent, title, label, items, false, instance, fieldName);
    }

    public StsComboBoxDialog(Frame parent, String title, String label, Object[] items, boolean addCancelButton, Object instance, String fieldName)
    {
        this.instance = instance;
        this.fieldName = fieldName;

        StsJPanel comboBoxPanel = new StsJPanel();
        comboBoxBean = new StsComboBoxFieldBean(this, "selectedItem", label, items);
        comboBoxPanel.add(comboBoxBean);
        this.addCancelButton = addCancelButton;
        constructPanel();
        try
        {
            field = instance.getClass().getDeclaredField(fieldName);
            dialog = StsToolkit.createDialog(title, panel, false);
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "constructor", "Failed to construct field setter for " + StsToolkit.getSimpleClassname(instance) + "." + fieldName, e);
        }
        if (parent != null) dialog.setLocationRelativeTo(parent);
    }

    private void constructPanel()
    {
        StsJPanel comboBoxPanel = new StsJPanel();
        comboBoxPanel.add(comboBoxBean);
        buttonPanel.gbc.fill = GridBagConstraints.HORIZONTAL;
        buttonPanel.addToRow(okayButton);
        if (addCancelButton)
            buttonPanel.addEndRow(cancelButton);
        panel.add(comboBoxPanel);
        panel.add(buttonPanel);
    }

    public void ok()
    {
        end();
    }

    public void cancel()
    {
        selectedItem = null;
        end();
    }

    private void end()
    {
        try
        {
            field.set(instance, selectedItem);
        }
        catch (Exception e)
        {
            StsException.systemError(this, "ok", "Failed to set field for " + StsToolkit.getSimpleClassname(instance) + "." + fieldName);
        }

        dialog.setVisible(false);
    }    

    public Object getSelectedItem() { return selectedItem; }
    public void setSelectedItem(Object item) { this.selectedItem = item; }


    public static void main(String[] args)
    {
        TestComboBoxDialog test = new TestComboBoxDialog();
        StsComboBoxDialog d = new StsComboBoxDialog(null, "StsComboBoxDialog test", "list of items:", test.items, true, test, "item");
    }
}

class TestComboBoxDialog
{
    String[] items = new String[]{"One", "Two"};
    String item = items[0];

    public String getItem() { return item; }
    public void setItem(String item)
    {
        this.item = item;
        System.out.println("Selected item: " + item);
    }
}
