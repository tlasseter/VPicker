
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI.Beans;

import com.Sts.DBTypes.StsColor;
import com.Sts.UI.StsTextAreaScrollPane;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;


public class StsListFieldBean extends StsFieldBean implements ListSelectionListener
{
    JList list = new JList();
    DefaultListModel listModel = new DefaultListModel();
    Object[] items;

    public StsListFieldBean()
    {
		initializeList(true);
   }

    public StsListFieldBean(Class c, String fieldName, String fieldLabel, Object[] items)
    {
        super();
        this.initialize(c, fieldName, fieldLabel, items);
    }

    public StsListFieldBean(Object instance, String fieldName, String fieldLabel, Object[] items)
    {
        super();
        this.initialize(instance, fieldName, fieldLabel, items);
    }

    public void initialize(Class c, String fieldName, String fieldLabel, Object[] items)
    {
        classInitialize(c, fieldName, fieldLabel);
        setListItems(items);
		initializeList(true);
        layoutBean();
    }

	public void initialize(Object beanObject, String fieldName, boolean editable)
	{
		this.initialize(beanObject, fieldName, editable, null);
	}

	// not sure what to do with editable.
	public void initialize(Object beanObject, String fieldName, boolean editable, String label)
	{
		super.initialize(beanObject, fieldName, editable, label);
		layoutBean();
		initializeList(editable);
    }

    public void initialize(Object instance, String fieldName, String fieldLabel, Object[] items)
    {
        super.initialize(instance, fieldName, true, fieldLabel);
        setListItems(items);
        layoutBean();
		initializeList(true);
    }

	private void initializeList(boolean editable)
	{
		list.setModel(listModel);
		list.setBorder(BorderFactory.createLoweredBevelBorder());
		if(editable) list.addListSelectionListener(this);
	}

    public void setSingleSelect()
    {
        list.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    }

    protected void layoutBean()
    {
		labelAndComponentLayout(label, list);
    }
    
    // Override to make sure the list extends to edges of provided space. Otherwise one one item is in list it
    // is not obvious that it is a list.
	protected void labelAndComponentLayout(Component label, Component component)
	{
		if (label != null)
		{
			gbc.anchor = GridBagConstraints.WEST;
			addToRow(label);
		}
		gbc.fill = GridBagConstraints.BOTH;
        gbc.anchor = GridBagConstraints.EAST;
		addToRow(component);
		gbc.anchor = GridBagConstraints.CENTER;
		gbc.fill = GridBagConstraints.NONE;
	}
	
    public void valueChanged(ListSelectionEvent e)
    {
        if(e.getValueIsAdjusting()) return;
        setValueInPanelObjects();
    }

	public void setValueInPanelObjects()
	{
		Object[] panelObjects = getPanelObjects();
		if (panelObjects == null)return;
		Object panelObject = null;
		Object value = null;

		try
		{
			if (set != null)
			{
				value = getValueObject();
				if (value == null)return;
                if(value.getClass().isArray())
                {
                    Object[] values = (Object[])value;
                    for (int n = 0; n < panelObjects.length; n++)
                    {
                        panelObject = panelObjects[n];
                        if (panelObject != null)
                        for(Object valueObject : values)
                            set.invoke(panelObject, valueObject);
                    }
                }
                else
                {
                    for (int n = 0; n < panelObjects.length; n++)
                    {
                        panelObject = panelObjects[n];
                        if (panelObject != null) set.invoke(panelObject, value);
                    }
                }
            }
			else if (field != null)
			{
				value = getValueObject();
				if (value == null)return;
                if(value.getClass().isArray())
                {
                    Object[] values = (Object[])value;
                    for (int n = 0; n < panelObjects.length; n++)
                    {
                        panelObject = panelObjects[n];
                        if (panelObject != null)
                        for(Object valueObject : values)
                            field.set(panelObject, valueObject);
                    }
                }
                else
                {
                    for (int n = 0; n < panelObjects.length; n++)
                    {
                        panelObject = panelObjects[n];
                        if (panelObject != null) field.set(panelObject, value);
                    }
                }
			}
		}
		catch (Exception e)
		{
			outputFieldException(e, panelObject, value);
		}
	}

    public Component[] getBeanComponents() { return new Component[] { label, list }; }
	public Component getMainComponent() { return list; }
    public Object getValueObject()
    {
        Object[] selectedValues = list.getSelectedValues();
        if(selectedValues.length == 0) return null;
        return selectedValues; 
    }
//	public Object getDefaultValueObject() { return list.getModel().getElementAt(0); }
//    public Object getValueObject() { return getString(); }
//    public String getString() { return (String) comboBox.getSelectedItem(); }

    public void setListItems(Object[] items)
    {
        if(items == null || items.length == 0) return;
        this.items = items;
        //listModel.removeAllElements();
        listModel.clear();
        {
            for( int i=0; i<items.length; i++ )
                listModel.addElement(items[i]);
        }
//        list.setSelectedIndex(0);
    }

    public void clearSelections()
    {
        list.clearSelection();
    }

    public void doSetValueObject(Object object)
    {
        Object item = getItem(object);
        if(item == null) return;
        list.setSelectedValue(item, true);
    }

    // override in subclass if object is inside an item
    public Object getItem(Object object)
    {
        if(items == null) return null;
        else return object;
    }

	public String toString()
	{
		return list.getSelectedValue().toString();
	}

	public Object fromString(String string)
	{
		if(items == null) return null;
		for(int n = 0; n < items.length; n++)
			if(items[n].toString().equals(string))
				return items[n];
		return null;
	}

    public void setToolTipText(String tip)
    {
        list.setToolTipText(tip);
    }

    public int getSelectedIndex()
    {
        return list.getSelectedIndex();
    }

    public void setSelectedIndex(int index)
    {
        list.setSelectedIndex(index);
    }

    public void setSelectedAll()
    {
        if(items == null) return;
        list.setSelectionInterval(0, items.length-1);
    }

    public void setSelectedValue(Object object)
    {
        list.setSelectedValue(object, true);
    }
    
    public void setSelectedValue(Object[] objects)
    {
        list.setSelectedValue(objects, true);
    }
    
    public Object[] getSelectedObjects()
    {
        return list.getSelectedValues();
    }
/*
    public void setBorder(Border border)
    {
		if(list == null) return;
        list.setBorder(border);
    }
*/
    public void setMaximumSize(Dimension size)
    {
        list.setMaximumSize(size);
    }

    public void setMinimumSize(Dimension size)
    {
        list.setMaximumSize(size);
    }

    public void setPreferredSize(Dimension size)
    {
        list.setMaximumSize(size);
    }

    static public void main(String[] args)
    {
        TestList testList = new TestList();
        StsListFieldBean listBean = new StsListFieldBean(testList, "itemsSelected", "Select items", testList.items);
        JScrollPane itemScrollPane = new JScrollPane();
        itemScrollPane.getViewport().add(listBean, null);
        StsJPanel itemsBox = new StsGroupBox("Box of items");
        itemsBox.setPreferredSize(200, 100);
        itemsBox.gbc.fill = GridBagConstraints.BOTH;
        itemsBox.addEndRow(itemScrollPane);
        StsJPanel panel = new StsJPanel();
        panel.gbc.fill = GridBagConstraints.BOTH;
        panel.add(itemsBox);

        StsTextAreaScrollPane scrollPane = new StsTextAreaScrollPane();
        JTextArea fileTextArea = scrollPane.textArea;
        // fileViewPanel.gbc.weighty = 0.5;
        StsJPanel fileViewPanel = new StsGroupBox("View of item");
        fileViewPanel.gbc.fill = GridBagConstraints.BOTH;
        fileViewPanel.add(scrollPane);
        panel.add(fileViewPanel);
        // listBean.setSingleSelect();
        StsToolkit.createDialog(panel);
    }
}

class TestList
{
    Object[] items = new String[] { "item1", "item2", "item3"};
    Object[] itemsSelected;

    TestList()
    {
        items = new String[100];
        for(int n = 0; n < 100; n++)
            items[n] = new String("item " + n);
    }

    public Object[] getItemsSelected() { return itemsSelected; }
    public void setItemsSelected(Object[] itemsSelected)
    {
        if(itemsSelected == null) return;
        this.itemsSelected = itemsSelected;
        int nItemsSelected = itemsSelected.length;
        if(nItemsSelected == 0) return;
        System.out.print("Selected items: ");
        for(int n = 0; n < nItemsSelected; n++)
            System.out.print((String)itemsSelected[n] + " ");
        System.out.println();
    }
}
